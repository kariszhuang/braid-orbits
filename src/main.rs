use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

// Legacy exact orders fit u128 through 34!. The explicit wide smooth build
// uses u64 point masks through degree63 and does not claim overflowing orders.
// The dense legacy engine still enforces its separate small-group limits.
#[cfg(not(smooth_wide_keys))]
const MAX_D: usize = 34;
#[cfg(smooth_wide_keys)]
const MAX_D: usize = 63;
mod smart;
const FACT_U64: [u64; 17] = [
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5_040,
    40_320,
    362_880,
    3_628_800,
    39_916_800,
    479_001_600,
    6_227_020_800,
    87_178_291_200,
    1_307_674_368_000,
    20_922_789_888_000,
];
const NONE_U32: u32 = u32::MAX;
const EMPTY_KEY: u64 = u64::MAX;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum Ambient {
    Alternating,
    Symmetric,
}

impl Ambient {
    fn parse(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "a" | "alt" | "alternating" => Ok(Self::Alternating),
            "s" | "sym" | "symmetric" => Ok(Self::Symmetric),
            _ => Err(format!("unknown group {s:?}; use A or S")),
        }
    }
    fn short(self) -> &'static str {
        match self { Self::Alternating => "A", Self::Symmetric => "S" }
    }
    fn target_order(self, d: usize) -> u128 {
        let f = factorial_u128(d);
        match self { Self::Alternating => f / 2, Self::Symmetric => f }
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Hash)]
struct Perm {
    img: [u8; MAX_D],
}

impl std::fmt::Debug for Perm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("Perm").field(&self.img).finish()
    }
}

impl Perm {
    fn identity(d: usize) -> Self {
        let mut img = [0u8; MAX_D];
        for (i, x) in img.iter_mut().enumerate().take(d) { *x = i as u8; }
        Self { img }
    }

    // Group convention: (a*b)(x) = a(b(x)).
    fn compose(self, b: Self, d: usize) -> Self {
        let mut out = [0u8; MAX_D];
        for x in 0..d { out[x] = self.img[b.img[x] as usize]; }
        Self { img: out }
    }

    fn inverse(self, d: usize) -> Self {
        let mut out = [0u8; MAX_D];
        for x in 0..d { out[self.img[x] as usize] = x as u8; }
        Self { img: out }
    }

    fn is_identity(self, d: usize) -> bool {
        (0..d).all(|i| self.img[i] as usize == i)
    }

    fn is_even(self, d: usize) -> bool {
        let mut inv = 0usize;
        for i in 0..d {
            for j in (i + 1)..d {
                if self.img[i] > self.img[j] { inv ^= 1; }
            }
        }
        inv == 0
    }

    fn conjugate_by(self, h: Self, d: usize) -> Self {
        h.compose(self, d).compose(h.inverse(d), d)
    }

    fn cycle_string(self, d: usize) -> String {
        let mut seen = [false; MAX_D];
        let mut parts = Vec::new();
        for i in 0..d {
            if seen[i] || self.img[i] as usize == i { continue; }
            let mut cyc = Vec::new();
            let mut x = i;
            while !seen[x] {
                seen[x] = true;
                cyc.push((x + 1).to_string());
                x = self.img[x] as usize;
            }
            parts.push(format!("({})", cyc.join(" ")));
        }
        if parts.is_empty() { "()".to_string() } else { parts.join("") }
    }
}

fn factorial_u128(n: usize) -> u128 {
    (2..=n).fold(1u128, |a, x| a * x as u128)
}

fn factorial_usize(n: usize) -> usize {
    assert!(n < FACT_U64.len(), "legacy factorial rank only supports d <= 16");
    FACT_U64[n] as usize
}

fn lex_rank(p: Perm, d: usize) -> usize {
    assert!(d < FACT_U64.len(), "legacy permutation rank only supports d <= 16");
    let mut rank = 0usize;
    for i in 0..d {
        let mut smaller = 0usize;
        for j in (i + 1)..d {
            if p.img[j] < p.img[i] { smaller += 1; }
        }
        rank += smaller * FACT_U64[d - 1 - i] as usize;
    }
    rank
}

fn next_permutation(a: &mut [u8]) -> bool {
    if a.len() < 2 { return false; }
    let mut i = a.len() - 2;
    loop {
        if a[i] < a[i + 1] { break; }
        if i == 0 { a.reverse(); return false; }
        i -= 1;
    }
    let mut j = a.len() - 1;
    while a[j] <= a[i] { j -= 1; }
    a.swap(i, j);
    a[(i + 1)..].reverse();
    true
}

#[derive(Clone)]
struct SplitMix64 { state: u64 }
impl SplitMix64 {
    fn new(seed: u64) -> Self { Self { state: seed } }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        assert!(n > 0);
        // Rejection removes modulo bias.
        let zone = u64::MAX - (u64::MAX % n as u64);
        loop {
            let x = self.next_u64();
            if x < zone { return (x % n as u64) as usize; }
        }
    }
    fn shuffle<T>(&mut self, a: &mut [T]) {
        for i in (1..a.len()).rev() {
            let j = self.below(i + 1);
            a.swap(i, j);
        }
    }
}

fn random_perm(d: usize, ambient: Ambient, rng: &mut SplitMix64) -> Perm {
    let mut p = Perm::identity(d);
    rng.shuffle(&mut p.img[..d]);
    if ambient == Ambient::Alternating && !p.is_even(d) {
        p.img.swap(d - 2, d - 1); // fixed bijection odd -> even, hence still uniform on A_d
    }
    p
}

fn fixed_split_perm(d: usize, split: usize, ambient: Ambient, rng: &mut SplitMix64) -> Perm {
    assert!(split > 0 && split < d);
    let mut p = Perm::identity(d);
    rng.shuffle(&mut p.img[..split]);
    rng.shuffle(&mut p.img[split..d]);
    if ambient == Ambient::Alternating && !p.is_even(d) {
        if split >= 2 { p.img.swap(0, 1); }
        else { p.img.swap(split, split + 1); }
    }
    p
}

fn wreath_perm(d: usize, block_size: usize, ambient: Ambient, rng: &mut SplitMix64) -> Perm {
    assert!(block_size >= 2 && block_size < d && d % block_size == 0);
    let blocks = d / block_size;
    let mut block_perm: Vec<usize> = (0..blocks).collect();
    rng.shuffle(&mut block_perm);
    let mut out = Perm::identity(d);
    for src_b in 0..blocks {
        let mut local: Vec<usize> = (0..block_size).collect();
        rng.shuffle(&mut local);
        let dst_b = block_perm[src_b];
        for j in 0..block_size {
            out.img[src_b * block_size + j] = (dst_b * block_size + local[j]) as u8;
        }
    }
    if ambient == Ambient::Alternating && !out.is_even(d) {
        // Two positions in one source block map into one target block, so this stays in the wreath subgroup.
        out.img.swap(0, 1);
    }
    out
}

fn nontrivial_divisors(d: usize) -> Vec<usize> {
    (2..d).filter(|s| d % s == 0).collect()
}

// ---------- Exact subgroup order on points: deterministic Schreier recursion ----------

#[cfg(test)]
fn dedup_nonidentity(gens: Vec<Perm>, d: usize) -> Vec<Perm> {
    let mut set = HashSet::with_capacity(gens.len() * 2 + 1);
    let mut out = Vec::new();
    for g in gens {
        if !g.is_identity(d) && set.insert(g) { out.push(g); }
    }
    out
}

fn subgroup_order_schreier(d: usize, gens0: &[Perm]) -> u128 {
    smart::subgroup_order(d, gens0)
}

fn generators_transitive(d: usize, generators: &[Perm]) -> bool {
    let full=(1u64<<d)-1;
    let (mut seen,mut todo)=(1u64,1u64);
    while todo!=0 {
        let x=todo.trailing_zeros() as usize;
        todo &= todo-1;
        for g in generators {
            let bit=1u64<<g.img[x];
            if seen&bit==0 {seen|=bit; todo|=bit;}
        }
        if seen==full {return true;}
    }
    false
}

#[cfg(test)]
fn subgroup_order_reference(d: usize, gens0: &[Perm]) -> u128 {
    let id = Perm::identity(d);
    let mut current = dedup_nonidentity(gens0.to_vec(), d);
    let mut order = 1u128;

    for base in 0..d {
        if current.is_empty() { break; }

        let mut action = current.clone();
        action.extend(current.iter().copied().map(|g| g.inverse(d)));
        action = dedup_nonidentity(action, d);

        let mut trans: [Option<Perm>; MAX_D] = [None; MAX_D];
        let mut q = VecDeque::new();
        trans[base] = Some(id);
        q.push_back(base);
        let mut orbit = Vec::new();

        while let Some(x) = q.pop_front() {
            orbit.push(x);
            let ux = trans[x].unwrap();
            for &s in &action {
                let y = s.img[x] as usize;
                if trans[y].is_none() {
                    trans[y] = Some(s.compose(ux, d));
                    q.push_back(y);
                }
            }
        }

        order *= orbit.len() as u128;

        let mut next = Vec::new();
        let mut seen = HashSet::new();
        for &x in &orbit {
            let ux = trans[x].unwrap();
            for &s in &action {
                let y = s.img[x] as usize;
                let uy = trans[y].unwrap();
                // uy^{-1} * s * ux fixes base.
                let h = uy.inverse(d).compose(s.compose(ux, d), d);
                if !h.is_identity(d) && seen.insert(h) { next.push(h); }
            }
        }
        current = next;
    }
    order
}

fn tuple_generates_dynamic(d: usize, ambient: Ambient, t: &[Perm; 4]) -> bool {
    generators_transitive(d,&t[..3])
        && subgroup_order_schreier(d, &t[..3]) == ambient.target_order(d)
}

fn pair_contains_ad_dynamic(d: usize, ambient: Ambient, a: Perm, b: Perm) -> bool {
    if !generators_transitive(d,&[a,b]) {return false;}
    let ord = subgroup_order_schreier(d, &[a, b]);
    match ambient {
        Ambient::Alternating => ord == factorial_u128(d) / 2,
        Ambient::Symmetric => {
            if a.is_even(d) && b.is_even(d) {
                ord == factorial_u128(d) / 2
            } else {
                ord == factorial_u128(d)
            }
        }
    }
}

fn point_orbit_count_perms(d: usize, a: Perm, b: Perm) -> usize {
    let mut parent = [0u8; MAX_D];
    for i in 0..d { parent[i] = i as u8; }
    fn find(parent: &mut [u8; MAX_D], mut x: usize) -> usize {
        while parent[x] as usize != x {
            let p = parent[x] as usize;
            parent[x] = parent[p];
            x = parent[x] as usize;
        }
        x
    }
    fn union(parent: &mut [u8; MAX_D], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb { parent[rb] = ra as u8; }
    }
    for x in 0..d {
        union(&mut parent, x, a.img[x] as usize);
        union(&mut parent, x, b.img[x] as usize);
    }
    let mut roots = [false; MAX_D];
    for x in 0..d { let r = find(&mut parent, x); roots[r] = true; }
    roots[..d].iter().filter(|&&z| z).count()
}

fn product_one_fourth(d: usize, a: Perm, b: Perm, c: Perm) -> Perm {
    a.compose(b, d).compose(c, d).inverse(d)
}

fn hurwitz_move_dynamic(d: usize, t: [Perm; 4], mv: usize) -> [Perm; 4] {
    let mut u = t;
    match mv {
        0 | 1 | 2 => {
            let i = mv;
            let a = t[i];
            let b = t[i + 1];
            u[i] = b.conjugate_by(a, d); // a*b*a^{-1}
            u[i + 1] = a;
        }
        3 | 4 | 5 => {
            let i = mv - 3;
            let a = t[i];
            let b = t[i + 1];
            u[i] = b;
            u[i + 1] = a.conjugate_by(b.inverse(d), d); // b^{-1}*a*b
        }
        _ => unreachable!(),
    }
    u
}

fn tuple_product_is_identity(d: usize, t: &[Perm; 4]) -> bool {
    t[0].compose(t[1], d).compose(t[2], d).compose(t[3], d).is_identity(d)
}

// ---------- Dense finite group engine for exact quotient certification ----------

struct FiniteGroup {
    d: usize,
    ambient: Ambient,
    elems: Vec<Perm>,
    rank_to_local: Vec<u32>,
    inv: Vec<u32>,
    mul: Vec<u32>,
    conj: Vec<u32>,
    identity: u32,
    rep_of: Vec<u32>,
    to_rep: Vec<u32>,
    reps: Vec<u32>,
    rep_slot: Vec<u32>,
    centralizers: Vec<Vec<u32>>,
}

impl FiniteGroup {
    fn build(d: usize, ambient: Ambient, threads: usize) -> Self {
        assert!(d <= MAX_D);
        let full_n = factorial_usize(d);
        eprintln!("enumerating {}{} inside S_{} ({} full permutations)...", ambient.short(), d, d, full_n);
        let mut elems = Vec::with_capacity(match ambient { Ambient::Alternating => full_n / 2, Ambient::Symmetric => full_n });
        let mut rank_to_local = vec![NONE_U32; full_n];
        let mut a: Vec<u8> = (0..d as u8).collect();
        let mut full_rank = 0usize;
        loop {
            let mut img = [0u8; MAX_D];
            img[..d].copy_from_slice(&a);
            let p = Perm { img };
            if ambient == Ambient::Symmetric || p.is_even(d) {
                rank_to_local[full_rank] = elems.len() as u32;
                elems.push(p);
            }
            full_rank += 1;
            if !next_permutation(&mut a) { break; }
        }
        let m = elems.len();
        eprintln!("|G| = {m}");
        // S8 needs two 6.05-GiB tables plus the sharded canonical-state set.
        if m > 40_320 {
            panic!("exact quotient mode currently supports groups of order <= 40320 (through S_8)");
        }
        if m > 6000 {
            eprintln!("HIGH-MEMORY exact mode: |G|={m}; dense mul+conj tables alone use about {:.2} GiB",
                (2.0 * (m as f64) * (m as f64) * 4.0) / (1024.0_f64.powi(3)));
        }

        let identity = rank_to_local[0];
        assert_ne!(identity, NONE_U32);
        let mut inv = vec![0u32; m];
        for (i, &p) in elems.iter().enumerate() {
            let r = lex_rank(p.inverse(d), d);
            inv[i] = rank_to_local[r];
            assert_ne!(inv[i], NONE_U32);
        }

        let threads = threads.max(1).min(m.max(1));
        eprintln!("building {m} x {m} multiplication table with {threads} worker(s)...");
        let t0 = Instant::now();
        let mul = build_mul_table_parallel(d, &elems, &rank_to_local, threads);
        eprintln!("multiplication table: {:.2?}", t0.elapsed());

        eprintln!("building conjugation table...");
        let t1 = Instant::now();
        let mut conj = vec![0u32; m * m];
        for h in 0..m {
            let hi = inv[h] as usize;
            for x in 0..m {
                let hx = mul[h * m + x] as usize;
                conj[h * m + x] = mul[hx * m + hi];
            }
        }
        eprintln!("conjugation table: {:.2?}", t1.elapsed());

        eprintln!("computing conjugacy classes and centralizers...");
        let mut rep_of = vec![NONE_U32; m];
        let mut to_rep = vec![NONE_U32; m];
        let mut reps = Vec::new();
        for x in 0..m {
            if rep_of[x] != NONE_U32 { continue; }
            let r = x as u32;
            reps.push(r);
            for h in 0..m {
                let y = conj[h * m + x] as usize;
                if rep_of[y] == NONE_U32 {
                    rep_of[y] = r;
                    to_rep[y] = inv[h]; // if y=h*x*h^-1, then h^-1 maps y back to x=r
                }
            }
        }
        let mut rep_slot = vec![NONE_U32; m];
        let mut centralizers = Vec::with_capacity(reps.len());
        for (slot, &r) in reps.iter().enumerate() {
            rep_slot[r as usize] = slot as u32;
            let mut c = Vec::new();
            for h in 0..m {
                if conj[h * m + r as usize] == r { c.push(h as u32); }
            }
            centralizers.push(c);
        }
        eprintln!("{} conjugacy classes", reps.len());

        Self { d, ambient, elems, rank_to_local, inv, mul, conj, identity, rep_of, to_rep, reps, rep_slot, centralizers }
    }

    #[inline] fn m(&self) -> usize { self.elems.len() }
    #[inline] fn mul_id(&self, a: u32, b: u32) -> u32 { self.mul[a as usize * self.m() + b as usize] }
    #[inline] fn conj_id(&self, h: u32, x: u32) -> u32 { self.conj[h as usize * self.m() + x as usize] }

    fn fourth_id(&self, a: u32, b: u32, c: u32) -> u32 {
        let ab = self.mul_id(a, b);
        let abc = self.mul_id(ab, c);
        self.inv[abc as usize]
    }

    fn encode3(&self, a: u32, b: u32, c: u32) -> u64 {
        let m = self.m() as u64;
        (a as u64 * m + b as u64) * m + c as u64
    }

    fn decode3(&self, key: u64) -> (u32, u32, u32) {
        let m = self.m() as u64;
        let c = (key % m) as u32;
        let q = key / m;
        let b = (q % m) as u32;
        let a = (q / m) as u32;
        (a, b, c)
    }

    // Canonical representative for diagonal conjugation by G.
    fn canonical_key(&self, a: u32, b: u32, c: u32) -> u64 {
        let r = self.rep_of[a as usize];
        let h0 = self.to_rep[a as usize];
        let b1 = self.conj_id(h0, b);
        let c1 = self.conj_id(h0, c);
        let slot = self.rep_slot[r as usize] as usize;
        let cent = &self.centralizers[slot];
        let mut bb_best = b1;
        let mut cc_best = c1;
        for &z in cent {
            let bb = self.conj_id(z, b1);
            let cc = self.conj_id(z, c1);
            if (bb, cc) < (bb_best, cc_best) {
                bb_best = bb;
                cc_best = cc;
            }
        }
        self.encode3(r, bb_best, cc_best)
    }

    fn braid_neighbor_canonical(&self, key: u64, mv: usize) -> u64 {
        let (a, b, c) = self.decode3(key);
        let (x, y, z) = match mv {
            0 => (self.conj_id(a, b), a, c),
            1 => (a, self.conj_id(b, c), b),
            2 => {
                let d = self.fourth_id(a, b, c);
                (a, b, self.conj_id(c, d))
            }
            3 => (b, self.conj_id(self.inv[b as usize], a), c),
            4 => (a, c, self.conj_id(self.inv[c as usize], b)),
            5 => {
                let d = self.fourth_id(a, b, c);
                (a, b, d)
            }
            _ => unreachable!(),
        };
        self.canonical_key(x, y, z)
    }

    fn tuple_string_ids(&self, a: u32, b: u32, c: u32) -> String {
        let d4 = self.fourth_id(a, b, c);
        format!("[{}, {}, {}, {}]",
            self.elems[a as usize].cycle_string(self.d),
            self.elems[b as usize].cycle_string(self.d),
            self.elems[c as usize].cycle_string(self.d),
            self.elems[d4 as usize].cycle_string(self.d))
    }

    fn class_signature(&self, a: u32, b: u32, c: u32) -> Vec<u32> {
        let d4 = self.fourth_id(a, b, c);
        let mut v = vec![self.rep_of[a as usize], self.rep_of[b as usize], self.rep_of[c as usize], self.rep_of[d4 as usize]];
        v.sort_unstable();
        v
    }

    fn signature_bucket(&self, a: u32, b: u32, c: u32, shard_count: usize) -> usize {
        if shard_count <= 1 { return 0; }
        let sig = self.class_signature(a,b,c);
        let mut h = 0x243F6A8885A308D3u64;
        for x in sig {
            h ^= (x as u64).wrapping_add(0x9E3779B97F4A7C15).wrapping_add(h << 6).wrapping_add(h >> 2);
            h = U64Set::hash(h);
        }
        (h % shard_count as u64) as usize
    }
}

fn build_mul_table_parallel(d: usize, elems: &[Perm], rank_to_local: &[u32], threads: usize) -> Vec<u32> {
    let m = elems.len();
    let rows_per = (m + threads - 1) / threads;
    let mut pieces: Vec<(usize, Vec<u32>)> = thread::scope(|scope| {
        let mut handles = Vec::new();
        for tid in 0..threads {
            let start = tid * rows_per;
            let end = (start + rows_per).min(m);
            if start >= end { continue; }
            handles.push(scope.spawn(move || {
                let mut part = Vec::with_capacity((end - start) * m);
                for a in start..end {
                    for b in 0..m {
                        let p = elems[a].compose(elems[b], d);
                        let r = lex_rank(p, d);
                        let id = rank_to_local[r];
                        assert_ne!(id, NONE_U32);
                        part.push(id);
                    }
                }
                (start, part)
            }));
        }
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    pieces.sort_by_key(|x| x.0);
    let mut out = Vec::with_capacity(m * m);
    for (_, mut p) in pieces { out.append(&mut p); }
    assert_eq!(out.len(), m * m);
    out
}

struct GroupClosureWorkspace {
    mark: Vec<u32>,
    epoch: u32,
    queue: Vec<u32>,
}
impl GroupClosureWorkspace {
    fn new(m: usize) -> Self { Self { mark: vec![0; m], epoch: 1, queue: Vec::with_capacity(m) } }
    fn reset(&mut self) {
        self.queue.clear();
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 { self.mark.fill(0); self.epoch = 1; }
    }
    fn reaches(&mut self, g: &FiniteGroup, gens: &[u32], target: usize) -> bool {
        self.reset();
        let mut action = Vec::with_capacity(gens.len() * 2);
        for &x in gens {
            action.push(x);
            let xi = g.inv[x as usize];
            if xi != x { action.push(xi); }
        }
        self.mark[g.identity as usize] = self.epoch;
        self.queue.push(g.identity);
        let mut head = 0usize;
        while head < self.queue.len() {
            let x = self.queue[head];
            head += 1;
            for &s in &action {
                let y = g.mul_id(x, s);
                if self.mark[y as usize] != self.epoch {
                    self.mark[y as usize] = self.epoch;
                    self.queue.push(y);
                    if self.queue.len() >= target { return true; }
                }
            }
        }
        self.queue.len() >= target
    }
}

fn point_orbit_count_ids(g: &FiniteGroup, a: u32, b: u32) -> usize {
    point_orbit_count_perms(g.d, g.elems[a as usize], g.elems[b as usize])
}

fn tuple_generates_ids(g: &FiniteGroup, a: u32, b: u32, c: u32, ws: &mut GroupClosureWorkspace) -> bool {
    ws.reaches(g, &[a, b, c], g.m())
}

fn pair_contains_ad_ids(g: &FiniteGroup, a: u32, b: u32, ws: &mut GroupClosureWorkspace) -> bool {
    match g.ambient {
        Ambient::Alternating => ws.reaches(g, &[a, b], g.m()),
        Ambient::Symmetric => {
            let ae = g.elems[a as usize].is_even(g.d);
            let be = g.elems[b as usize].is_even(g.d);
            if ae && be { ws.reaches(g, &[a, b], g.m() / 2) }
            else { ws.reaches(g, &[a, b], g.m()) }
        }
    }
}

struct BitSet { words: Vec<u64>, nbits: usize }
impl BitSet {
    fn new(nbits: usize) -> Self { Self { words: vec![0; (nbits + 63) / 64], nbits } }
    #[inline] fn get(&self, i: usize) -> bool { debug_assert!(i < self.nbits); (self.words[i >> 6] >> (i & 63)) & 1 != 0 }
    #[inline] fn set(&mut self, i: usize) { debug_assert!(i < self.nbits); self.words[i >> 6] |= 1u64 << (i & 63); }
}

// Fast integer set; avoids SipHash overhead for millions of canonical state IDs.
struct U64Set {
    slots: Vec<u64>,
    len: usize,
}
impl U64Set {
    fn with_expected(expected: usize) -> Self {
        let need = ((expected as f64) / 0.70).ceil() as usize;
        let cap = need.max(16).next_power_of_two();
        Self { slots: vec![EMPTY_KEY; cap], len: 0 }
    }
    fn len(&self) -> usize { self.len }
    #[inline] fn hash(k: u64) -> u64 {
        let mut z = k.wrapping_add(0x9E3779B97F4A7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn contains(&self, key: u64) -> bool {
        let mask = self.slots.len() - 1;
        let mut i = Self::hash(key) as usize & mask;
        loop {
            let x = self.slots[i];
            if x == EMPTY_KEY { return false; }
            if x == key { return true; }
            i = (i + 1) & mask;
        }
    }
    fn insert(&mut self, key: u64) -> bool {
        assert_ne!(key, EMPTY_KEY);
        if (self.len + 1) * 10 >= self.slots.len() * 7 { self.grow(); }
        self.insert_no_grow(key)
    }
    fn insert_no_grow(&mut self, key: u64) -> bool {
        let mask = self.slots.len() - 1;
        let mut i = Self::hash(key) as usize & mask;
        loop {
            let x = self.slots[i];
            if x == EMPTY_KEY {
                self.slots[i] = key;
                self.len += 1;
                return true;
            }
            if x == key { return false; }
            i = (i + 1) & mask;
        }
    }
    fn grow(&mut self) {
        let new_len = self.slots.len() * 2;
        let old = std::mem::replace(&mut self.slots, vec![EMPTY_KEY; new_len]);
        self.len = 0;
        for x in old { if x != EMPTY_KEY { self.insert_no_grow(x); } }
    }
    fn keys(&self) -> impl Iterator<Item=u64> + '_ { self.slots.iter().copied().filter(|&x| x != EMPTY_KEY) }
}

#[derive(Default)]
struct ExactStats {
    joint_components: u64,
    generating_components: u64,
    strong_pass: u64,
    weak_pass: u64,
    strong_fail: u64,
    weak_fail: u64,
    generating_states_quotient: u64,
}

fn save_checkpoint(path: &Path, g: &FiniteGroup, seen: &U64Set, stats: &ExactStats, shard_count: usize, shard_index: usize) -> io::Result<()> {
    let tmp = path.with_extension("tmp");
    let mut w = BufWriter::new(File::create(&tmp)?);
    w.write_all(b"HBQ2")?;
    w.write_all(&[g.d as u8, match g.ambient { Ambient::Alternating => 0, Ambient::Symmetric => 1 }])?;
    w.write_all(&(shard_count as u32).to_le_bytes())?;
    w.write_all(&(shard_index as u32).to_le_bytes())?;
    for x in [stats.joint_components, stats.generating_components, stats.strong_pass, stats.weak_pass, stats.strong_fail, stats.weak_fail, stats.generating_states_quotient, seen.len() as u64] {
        w.write_all(&x.to_le_bytes())?;
    }
    for k in seen.keys() { w.write_all(&k.to_le_bytes())?; }
    w.flush()?;
    fs::rename(tmp, path)?;
    Ok(())
}

fn load_checkpoint(path: &Path, g: &FiniteGroup, shard_count: usize, shard_index: usize) -> io::Result<Option<(U64Set, ExactStats)>> {
    if !path.exists() { return Ok(None); }
    let mut r = File::open(path)?;
    let mut magic = [0u8; 4]; r.read_exact(&mut magic)?;
    if &magic != b"HBQ2" { return Err(io::Error::new(io::ErrorKind::InvalidData, "bad/old checkpoint magic")); }
    let mut hdr = [0u8; 2]; r.read_exact(&mut hdr)?;
    let amb = if hdr[1] == 0 { Ambient::Alternating } else { Ambient::Symmetric };
    if hdr[0] as usize != g.d || amb != g.ambient { return Err(io::Error::new(io::ErrorKind::InvalidData, "checkpoint group mismatch")); }
    let mut b4=[0u8;4];
    r.read_exact(&mut b4)?; let ck_sc=u32::from_le_bytes(b4) as usize;
    r.read_exact(&mut b4)?; let ck_si=u32::from_le_bytes(b4) as usize;
    if ck_sc != shard_count || ck_si != shard_index { return Err(io::Error::new(io::ErrorKind::InvalidData, "checkpoint shard mismatch")); }
    fn read_u64(r: &mut File) -> io::Result<u64> { let mut b=[0u8;8]; r.read_exact(&mut b)?; Ok(u64::from_le_bytes(b)) }
    let mut stats = ExactStats::default();
    stats.joint_components = read_u64(&mut r)?;
    stats.generating_components = read_u64(&mut r)?;
    stats.strong_pass = read_u64(&mut r)?;
    stats.weak_pass = read_u64(&mut r)?;
    stats.strong_fail = read_u64(&mut r)?;
    stats.weak_fail = read_u64(&mut r)?;
    stats.generating_states_quotient = read_u64(&mut r)?;
    let count = read_u64(&mut r)? as usize;
    let mut set = U64Set::with_expected(count.max(16));
    for _ in 0..count { set.insert(read_u64(&mut r)?); }
    Ok(Some((set, stats)))
}

fn exact_quotient(
    d: usize,
    ambient: Ambient,
    threads: usize,
    out_dir: &Path,
    checkpoint_every: u64,
    resume: bool,
    stop_on_counterexample: bool,
    shard_count: usize,
    shard_index: usize,
    max_components: u64,
) -> io::Result<()> {
    fs::create_dir_all(out_dir)?;
    if d < 4 || d > 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "exact dense quotient mode is enabled for A_4..A_8 and S_4..S_8",
        ));
    }
    if shard_count == 0 || shard_index >= shard_count {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "require 0 <= shard-index < shard-count"));
    }
    if d == 8 && ambient == Ambient::Symmetric && shard_count < 8 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "S8 requires at least 8 shards on this 64-GiB machine"));
    }
    eprintln!("exact signature shard {}/{}", shard_index, shard_count);
    let g = FiniteGroup::build(d, ambient, threads);
    let m = g.m();
    let ckpt = out_dir.join("checkpoint.bin");
    let (mut global_seen, mut stats) = if resume {
        if let Some((s, st)) = load_checkpoint(&ckpt, &g, shard_count, shard_index)? {
            eprintln!("resumed checkpoint: {} canonical states already complete", s.len());
            (s, st)
        } else {
            (U64Set::with_expected(((m * m + shard_count - 1) / shard_count).min(20_000_000)), ExactStats::default())
        }
    } else {
        (U64Set::with_expected(((m * m + shard_count - 1) / shard_count).min(20_000_000)), ExactStats::default())
    };

    let csv_path = out_dir.join("generating_components.csv");
    let new_csv = !csv_path.exists() || !resume;
    let csv_file = OpenOptions::new().create(true).write(true).append(resume).truncate(!resume).open(&csv_path)?;
    let mut csv = BufWriter::new(csv_file);
    if new_csv {
        writeln!(csv, "component,quotient_states,strong,weak,max_pair_rank,signature,seed,strong_witness,weak_witness")?;
    }

    // Certified failures are never overwritten.  --keep-going may discover more
    // than one bad component, so keep both a machine-readable append-only ledger
    // and one human-readable file per failing component.  On a fresh (non-resume)
    // run we truncate the ledger; on resume we append after the checkpoint.
    let ce_path = out_dir.join("CERTIFIED_COUNTEREXAMPLES.tsv");
    let ce_new = !ce_path.exists() || !resume;
    let ce_file = OpenOptions::new().create(true).write(true).append(resume).truncate(!resume).open(&ce_path)?;
    let mut ce_log = BufWriter::new(ce_file);
    if ce_new {
        writeln!(ce_log, "component\tgroup\tweak\tquotient_states\tmax_pair_rank\tsignature\tseed")?;
        ce_log.flush()?;
    }

    let mut closure = GroupClosureWorkspace::new(m);
    let mut queue = Vec::<u64>::new();
    let t_all = Instant::now();
    let pair_space = m.checked_mul(m).expect("m^2 overflow");

    for (rep_idx, &a_rep) in g.reps.iter().enumerate() {
        let slot = g.rep_slot[a_rep as usize] as usize;
        let cent = &g.centralizers[slot];
        eprintln!("class {}/{}: first entry rep {}, |C_G(rep)|={}", rep_idx + 1, g.reps.len(), g.elems[a_rep as usize].cycle_string(d), cent.len());
        let mut pair_seen = BitSet::new(pair_space);

        for pair in 0..pair_space {
            if pair_seen.get(pair) { continue; }
            let b0 = (pair / m) as u32;
            let c0 = (pair % m) as u32;

            // Enumerate the centralizer orbit of (b,c), mark it, and take lexicographic minimum.
            let mut best = (b0, c0);
            for &z in cent {
                let bb = g.conj_id(z, b0);
                let cc = g.conj_id(z, c0);
                pair_seen.set(bb as usize * m + cc as usize);
                if (bb, cc) < best { best = (bb, cc); }
            }
            let seed_key = g.encode3(a_rep, best.0, best.1);
            if g.signature_bucket(a_rep, best.0, best.1, shard_count) != shard_index { continue; }
            if global_seen.contains(seed_key) { continue; }

            let (sa, sb, sc) = g.decode3(seed_key);
            let generating = tuple_generates_ids(&g, sa, sb, sc, &mut closure);
            stats.joint_components += 1;

            queue.clear();
            queue.push(seed_key);
            global_seen.insert(seed_key);
            let mut head = 0usize;
            let mut comp_size = 0u64;
            let mut strong = false;
            let mut weak = false;
            let mut max_rank = 0usize;
            let mut strong_witness: Option<String> = None;
            let mut weak_witness: Option<String> = None;

            while head < queue.len() {
                let key = queue[head]; head += 1; comp_size += 1;
                if comp_size % 1_000_000 == 0 {
                    eprintln!("  current quotient component: {comp_size} states; total seen={}", global_seen.len());
                }
                let (a,b,c) = g.decode3(key);

                if generating {
                    let oc = point_orbit_count_ids(&g, a, b);
                    let rank = d - oc;
                    if rank > max_rank { max_rank = rank; }
                    if oc == 1 {
                        if !weak {
                            weak = true;
                            weak_witness = Some(g.tuple_string_ids(a,b,c));
                        }
                        if !strong && pair_contains_ad_ids(&g, a, b, &mut closure) {
                            strong = true;
                            strong_witness = Some(g.tuple_string_ids(a,b,c));
                        }
                    }
                }

                for mv in 0..6 {
                    let nk = g.braid_neighbor_canonical(key, mv);
                    if global_seen.insert(nk) { queue.push(nk); }
                }
            }

            if generating {
                stats.generating_components += 1;
                stats.generating_states_quotient += comp_size;
                if strong { stats.strong_pass += 1; } else { stats.strong_fail += 1; }
                if weak { stats.weak_pass += 1; } else { stats.weak_fail += 1; }
                let signature = g.class_signature(sa,sb,sc).iter().map(|x| x.to_string()).collect::<Vec<_>>().join(";");
                let seed_str = g.tuple_string_ids(sa,sb,sc).replace('"', "''");
                let sw = strong_witness.clone().unwrap_or_default().replace('"', "''");
                let ww = weak_witness.clone().unwrap_or_default().replace('"', "''");
                writeln!(csv, "{},{},{},{},{},\"{}\",\"{}\",\"{}\",\"{}\"",
                    stats.generating_components, comp_size, strong, weak, max_rank, signature, seed_str, sw, ww)?;
                csv.flush()?;

                if !strong {
                    let ce_no = stats.strong_fail;
                    let ce_name = format!("CERTIFIED_COUNTEREXAMPLE_{ce_no:06}_component_{:06}.txt", stats.generating_components);
                    let mut w = BufWriter::new(File::create(out_dir.join(&ce_name))?);
                    writeln!(w, "CERTIFIED failure of the strong statement in quotient B4 x diagonal-conjugation component")?;
                    writeln!(w, "group: {}_{}", ambient.short(), d)?;
                    writeln!(w, "generating component number: {}", stats.generating_components)?;
                    writeln!(w, "weak property found: {weak}")?;
                    writeln!(w, "quotient component size: {comp_size}")?;
                    writeln!(w, "signature: {signature}")?;
                    writeln!(w, "seed: {}", g.tuple_string_ids(sa,sb,sc))?;
                    writeln!(w, "max pair connectivity rank d-#orbits: {max_rank} / {}", d-1)?;
                    if !weak {
                        writeln!(w, "This also certifies failure of the weaker transitivity statement.")?;
                    }
                    w.flush()?;

                    // Append to the ledger as well; this survives any later discoveries.
                    let seed_tsv = g.tuple_string_ids(sa,sb,sc).replace('\t', " ").replace('\n', " ");
                    writeln!(ce_log, "{}\t{}_{}\t{}\t{}\t{}\t{}\t{}",
                        stats.generating_components, ambient.short(), d, weak, comp_size, max_rank, signature, seed_tsv)?;
                    ce_log.flush()?;

                    eprintln!("*** CERTIFIED COUNTEREXAMPLE LOGGED (continuing={}): strong=false, weak={weak}: {}",
                        !stop_on_counterexample, g.tuple_string_ids(sa,sb,sc));
                    if stop_on_counterexample {
                        save_checkpoint(&ckpt, &g, &global_seen, &stats, shard_count, shard_index)?;
                        return Ok(());
                    }
                }
            }

            if checkpoint_every > 0 && stats.joint_components % checkpoint_every == 0 {
                save_checkpoint(&ckpt, &g, &global_seen, &stats, shard_count, shard_index)?;
                eprintln!("checkpoint: components={}, generating={}, seen={} ({:.2?})",
                    stats.joint_components, stats.generating_components, global_seen.len(), t_all.elapsed());
            }
            if max_components > 0 && stats.joint_components >= max_components {
                save_checkpoint(&ckpt, &g, &global_seen, &stats, shard_count, shard_index)?;
                eprintln!("PILOT STOP: components={}, generating={}, seen={}, elapsed={:.2?}; shard incomplete, resume with --resume",
                    stats.joint_components, stats.generating_components, global_seen.len(), t_all.elapsed());
                return Ok(());
            }
        }
    }

    save_checkpoint(&ckpt, &g, &global_seen, &stats, shard_count, shard_index)?;
    let summary = out_dir.join("SUMMARY.txt");
    let mut w = BufWriter::new(File::create(&summary)?);
    writeln!(w, "group = {}_{}", ambient.short(), d)?;
    writeln!(w, "signature shard = {}/{}", shard_index, shard_count)?;
    writeln!(w, "joint B4 + diagonal-conjugation components = {}", stats.joint_components)?;
    writeln!(w, "generating components = {}", stats.generating_components)?;
    writeln!(w, "strong pass = {}, strong fail = {}", stats.strong_pass, stats.strong_fail)?;
    writeln!(w, "weak pass = {}, weak fail = {}", stats.weak_pass, stats.weak_fail)?;
    writeln!(w, "generating quotient states visited = {}", stats.generating_states_quotient)?;
    writeln!(w, "all canonical states visited = {}", global_seen.len())?;
    writeln!(w, "elapsed = {:.2?}", t_all.elapsed())?;
    if stats.strong_fail == 0 {
        writeln!(w, "RESULT: every generating component in this signature shard contains a pair generating/containing A_d.")?;
    }
    if stats.weak_fail == 0 {
        writeln!(w, "RESULT: every generating component in this signature shard contains a transitive first pair.")?;
    }
    w.flush()?;
    eprintln!("DONE. Summary written to {}", summary.display());
    Ok(())
}

// ---------- Massive randomized/adversarial search for d >= 7 ----------

#[derive(Clone, Copy, Debug)]
enum SeedFamily { Uniform, Intransitive, Imprimitive, Cyclic, Affine, Mixed }
impl SeedFamily {
    fn parse(s: &str) -> Result<Self,String> {
        match s.to_ascii_lowercase().as_str() {
            "uniform" => Ok(Self::Uniform),
            "intransitive" | "intrans" => Ok(Self::Intransitive),
            "imprimitive" | "imprim" => Ok(Self::Imprimitive),
            "cyclic" | "transitive-proper" | "proper-transitive" => Ok(Self::Cyclic),
            "affine" | "agl" => Ok(Self::Affine),
            "mixed" => Ok(Self::Mixed),
            _ => Err("family must be uniform|intransitive|imprimitive|cyclic|affine|mixed".into()),
        }
    }
}

fn random_d_cycle(d: usize, ambient: Ambient, rng: &mut SplitMix64) -> Option<Perm> {
    // A single d-cycle is transitive.  Its parity is (-1)^(d-1), so it belongs
    // to A_d exactly when d is odd.  For S_d it is always available.
    if ambient == Ambient::Alternating && d % 2 == 0 { return None; }
    let mut pts: Vec<usize> = (0..d).collect();
    rng.shuffle(&mut pts);
    let mut p = Perm::identity(d);
    for i in 0..d { p.img[pts[i]] = pts[(i + 1) % d] as u8; }
    debug_assert!(ambient == Ambient::Symmetric || p.is_even(d));
    Some(p)
}

fn perm_power(mut a: Perm, mut k: usize, d: usize) -> Perm {
    let mut out = Perm::identity(d);
    while k > 0 {
        if k & 1 == 1 { out = out.compose(a, d); }
        k >>= 1;
        if k > 0 { a = a.compose(a, d); }
    }
    out
}

fn is_prime_usize(n: usize) -> bool {
    if n < 2 { return false; }
    if n % 2 == 0 { return n == 2; }
    let mut q = 3usize;
    while q * q <= n {
        if n % q == 0 { return false; }
        q += 2;
    }
    true
}

fn random_affine_pair_prime(d: usize, ambient: Ambient, rng: &mut SplitMix64) -> Option<(Perm, Perm)> {
    // Prime degree p: translation x->x+1 plus scaling x->a*x lies in AGL(1,p).
    // The pair is transitive but generates a subgroup of order <= p(p-1), far
    // smaller than A_p for p>=5.  This targets prime degrees where there is no
    // nontrivial imprimitive block system.
    if d < 5 || !is_prime_usize(d) { return None; }

    // Randomize coordinates so A_p experiments sample both split p-cycle classes.
    let mut coord: Vec<usize> = (0..d).collect();
    rng.shuffle(&mut coord);

    let mut translation = Perm::identity(d);
    for i in 0..d {
        translation.img[coord[i]] = coord[(i + 1) % d] as u8;
    }
    debug_assert!(translation.is_even(d));

    for _ in 0..64 {
        let a = 2 + rng.below(d - 2); // 2,...,p-1; avoid identity scaling
        let mut scaling = Perm::identity(d);
        for i in 0..d {
            scaling.img[coord[i]] = coord[(a * i) % d] as u8;
        }
        if ambient == Ambient::Symmetric || scaling.is_even(d) {
            debug_assert_eq!(point_orbit_count_perms(d, translation, scaling), 1);
            return Some((translation, scaling));
        }
    }
    None
}

fn make_seed(d: usize, ambient: Ambient, fam: SeedFamily, rng: &mut SplitMix64) -> Option<[Perm;4]> {
    let actual = match fam {
        SeedFamily::Mixed => {
            let divs = nontrivial_divisors(d);
            let cyclic_ok = ambient == Ambient::Symmetric || d % 2 == 1;
            let affine_ok = is_prime_usize(d) && d >= 5;
            let r = rng.below(100);
            if affine_ok && r < 30 { SeedFamily::Affine }
            else if !divs.is_empty() && r < 30 { SeedFamily::Imprimitive }
            else if r < 65 { SeedFamily::Intransitive }
            else if cyclic_ok && r < 85 { SeedFamily::Cyclic }
            else { SeedFamily::Uniform }
        }
        x => x,
    };

    for _ in 0..200 {
        let (a,b) = match actual {
            SeedFamily::Uniform => (random_perm(d,ambient,rng), random_perm(d,ambient,rng)),
            SeedFamily::Intransitive => {
                let split = 1 + rng.below(d-1);
                (fixed_split_perm(d,split,ambient,rng), fixed_split_perm(d,split,ambient,rng))
            }
            SeedFamily::Imprimitive => {
                let divs = nontrivial_divisors(d);
                if divs.is_empty() { continue; }
                let s = divs[rng.below(divs.len())];
                let mut pair = None;
                for _ in 0..50 {
                    let x = wreath_perm(d,s,ambient,rng);
                    let y = wreath_perm(d,s,ambient,rng);
                    if point_orbit_count_perms(d,x,y) == 1 { pair = Some((x,y)); break; }
                }
                match pair { Some(p) => p, None => continue }
            }
            SeedFamily::Cyclic => {
                let a = match random_d_cycle(d, ambient, rng) { Some(x) => x, None => continue };
                // <a,b> = <a> is transitive but tiny, hence definitely does not
                // contain A_d for d>=4.  This directly targets the strong claim.
                let k = 1 + rng.below(d - 1);
                let b = perm_power(a, k, d);
                (a,b)
            }
            SeedFamily::Affine => {
                match random_affine_pair_prime(d, ambient, rng) {
                    Some(pair) => pair,
                    None => continue,
                }
            }
            SeedFamily::Mixed => unreachable!(),
        };
        let c = random_perm(d,ambient,rng);
        let d4 = product_one_fourth(d,a,b,c);
        let t = [a,b,c,d4];
        if tuple_generates_dynamic(d,ambient,&t) { return Some(t); }
    }
    None
}

#[derive(Default)]
struct WalkResult {
    weak_step: Option<u64>,
    strong_step: Option<u64>,
    best_rank: usize,
    steps: u64,
    timed_out: bool,
}

fn search_orbit_walk(
    d: usize,
    ambient: Ambient,
    seed: [Perm;4],
    rng: &mut SplitMix64,
    steps: u64,
    greedy_every: u64,
    deadline: Option<Instant>,
) -> WalkResult {
    let mut t = seed;
    let mut out = WalkResult::default();
    for step in 0..=steps {
        // Keep time-budgeted overnight campaigns responsive even when one seed gets
        // promoted to a very deep walk. Checking every 4096 iterations makes the
        // timing overhead negligible.
        if (step & 4095) == 0 {
            if let Some(dl) = deadline {
                if Instant::now() >= dl {
                    out.steps = step;
                    out.timed_out = true;
                    return out;
                }
            }
        }
        let oc = point_orbit_count_perms(d,t[0],t[1]);
        out.best_rank = out.best_rank.max(d - oc);
        if oc == 1 {
            if out.weak_step.is_none() { out.weak_step = Some(step); }
            if pair_contains_ad_dynamic(d,ambient,t[0],t[1]) {
                out.strong_step = Some(step);
                out.steps = step;
                return out;
            }
        }
        if step == steps { break; }

        let mv = if greedy_every > 0 && step % greedy_every == 0 {
            // Lexicographic objective:
            //   (1) maximize point connectivity of <g1,g2>;
            //   (2) once transitive, maximize the exact subgroup order.
            // This deliberately pushes a transitive pair out of small primitive/
            // imprimitive subgroups toward A_d/S_d instead of wandering blindly.
            let mut best_rank = 0usize;
            let mut best_order = 0u128;
            let mut choices = [0usize; 6];
            let mut nchoices = 0usize;
            for cand in 0..6 {
                let u = hurwitz_move_dynamic(d,t,cand);
                let rank = d - point_orbit_count_perms(d,u[0],u[1]);
                let ord = if rank == d - 1 {
                    subgroup_order_schreier(d, &[u[0], u[1]])
                } else {
                    0
                };
                if rank > best_rank || (rank == best_rank && ord > best_order) {
                    best_rank = rank;
                    best_order = ord;
                    choices[0] = cand;
                    nchoices = 1;
                } else if rank == best_rank && ord == best_order {
                    choices[nchoices] = cand;
                    nchoices += 1;
                }
            }
            choices[rng.below(nchoices)]
        } else {
            rng.below(6)
        };
        t = hurwitz_move_dynamic(d,t,mv);
    }
    out.steps = steps;
    out
}

fn tuple_string_dynamic(d: usize, t: &[Perm;4]) -> String {
    format!("[{}, {}, {}, {}]", t[0].cycle_string(d),t[1].cycle_string(d),t[2].cycle_string(d),t[3].cycle_string(d))
}

#[derive(Default, Debug)]
struct RawCertResult {
    exhausted: bool,
    strong_found: bool,
    weak_seen: bool,
    states: usize,
    max_rank: usize,
    timed_out: bool,
}

fn certify_raw_orbit_bounded(
    d: usize,
    ambient: Ambient,
    seed: [Perm;4],
    max_states: usize,
    deadline: Option<Instant>,
) -> RawCertResult {
    let mut out = RawCertResult::default();
    if max_states == 0 { return out; }

    // This is deliberately the *raw* B_4 orbit, with no diagonal-conjugacy
    // quotient.  Therefore, if the queue empties, exhaustion is an immediate
    // exact certificate for this seed.  The cap only controls resource use.
    let mut seen: HashSet<[Perm;4]> = HashSet::with_capacity(max_states.min(1_000_000));
    let mut queue = VecDeque::<[Perm;4]>::new();
    seen.insert(seed);
    queue.push_back(seed);
    let mut popped = 0usize;

    while let Some(t) = queue.pop_front() {
        popped += 1;
        if (popped & 4095) == 0 {
            if let Some(dl) = deadline {
                if Instant::now() >= dl {
                    out.states = seen.len();
                    out.timed_out = true;
                    return out;
                }
            }
        }

        let oc = point_orbit_count_perms(d, t[0], t[1]);
        out.max_rank = out.max_rank.max(d - oc);
        if oc == 1 {
            out.weak_seen = true;
            if pair_contains_ad_dynamic(d, ambient, t[0], t[1]) {
                out.strong_found = true;
                out.states = seen.len();
                return out;
            }
        }

        for mv in 0..6 {
            let u = hurwitz_move_dynamic(d, t, mv);
            if seen.contains(&u) { continue; }
            if seen.len() >= max_states {
                out.states = seen.len();
                return out; // capped: inconclusive
            }
            seen.insert(u);
            queue.push_back(u);
        }
    }

    out.states = seen.len();
    out.exhausted = true;
    out
}

fn hunt(
    d: usize,
    ambient: Ambient,
    threads_n: usize,
    seeds_per_thread: u64,
    walk_steps: u64,
    restarts: u64,
    deep_multiplier: u64,
    family: SeedFamily,
    out_dir: &Path,
    master_seed: u64,
    max_seconds: u64,
    certify_cap: usize,
) -> io::Result<()> {
    if d > MAX_D { return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("d>{MAX_D} not compiled; raise MAX_D"))); }
    fs::create_dir_all(out_dir)?;
    let candidate_file = Arc::new(Mutex::new(BufWriter::new(File::create(out_dir.join("hard_candidates.txt"))?)));
    let certified_file = Arc::new(Mutex::new(BufWriter::new(File::create(out_dir.join("CERTIFIED_RAW_COUNTEREXAMPLES.tsv"))?)));
    {
        let mut w = certified_file.lock().unwrap();
        writeln!(w, "number\tgroup\tweak\traw_orbit_states\tmax_pair_rank\tseed")?;
        w.flush()?;
    }
    let total = Arc::new(AtomicU64::new(0));
    let generated = Arc::new(AtomicU64::new(0));
    let weak_pass = Arc::new(AtomicU64::new(0));
    let strong_pass = Arc::new(AtomicU64::new(0));
    let certified_strong_fail = Arc::new(AtomicU64::new(0));
    let certified_weak_fail = Arc::new(AtomicU64::new(0));
    let hard = Arc::new(AtomicU64::new(0));
    let t0 = Instant::now();
    let deadline = if max_seconds == 0 {
        None
    } else {
        Some(t0 + Duration::from_secs(max_seconds))
    };

    thread::scope(|scope| {
        for tid in 0..threads_n.max(1) {
            let candidate_file = Arc::clone(&candidate_file);
            let certified_file = Arc::clone(&certified_file);
            let total = Arc::clone(&total);
            let generated = Arc::clone(&generated);
            let weak_pass = Arc::clone(&weak_pass);
            let strong_pass = Arc::clone(&strong_pass);
            let certified_strong_fail = Arc::clone(&certified_strong_fail);
            let certified_weak_fail = Arc::clone(&certified_weak_fail);
            let hard = Arc::clone(&hard);
            scope.spawn(move || {
                let mut rng = SplitMix64::new(master_seed ^ ((tid as u64 + 1).wrapping_mul(0xD1B54A32D192ED03)));
                for _ in 0..seeds_per_thread {
                    if let Some(dl) = deadline {
                        if Instant::now() >= dl { break; }
                    }
                    let attempt_no = total.fetch_add(1, Ordering::Relaxed) + 1;
                    if tid == 0 && attempt_no % 100 == 0 {
                        eprintln!("hunt progress: attempted={} generating={} strong={} hard={}",
                            attempt_no, generated.load(Ordering::Relaxed), strong_pass.load(Ordering::Relaxed), hard.load(Ordering::Relaxed));
                    }
                    let seed = match make_seed(d,ambient,family,&mut rng) { Some(t) => t, None => continue };
                    generated.fetch_add(1, Ordering::Relaxed);
                    let mut best = WalkResult::default();
                    let mut solved = false;
                    for _ in 0..restarts.max(1) {
                        let r = search_orbit_walk(d,ambient,seed,&mut rng,walk_steps,16,deadline);
                        best.best_rank = best.best_rank.max(r.best_rank);
                        best.weak_step = match (best.weak_step,r.weak_step) { (None,x)=>x,(x,None)=>x,(Some(a),Some(b))=>Some(a.min(b)) };
                        best.strong_step = match (best.strong_step,r.strong_step) { (None,x)=>x,(x,None)=>x,(Some(a),Some(b))=>Some(a.min(b)) };
                        if r.strong_step.is_some() { solved = true; break; }
                        if r.timed_out { break; }
                    }
                    if best.weak_step.is_some() { weak_pass.fetch_add(1,Ordering::Relaxed); }
                    if solved {
                        strong_pass.fetch_add(1,Ordering::Relaxed);
                        continue;
                    }

                    // Automatically spend much more effort on anything that survives the normal budget.
                    let deep_steps = walk_steps.saturating_mul(deep_multiplier.max(1));
                    let had_weak_before_deep = best.weak_step.is_some();
                    let deep = search_orbit_walk(d,ambient,seed,&mut rng,deep_steps,8,deadline);
                    best.best_rank = best.best_rank.max(deep.best_rank);
                    if best.weak_step.is_none() { best.weak_step = deep.weak_step; }
                    if !had_weak_before_deep && deep.weak_step.is_some() {
                        weak_pass.fetch_add(1,Ordering::Relaxed);
                    }
                    if deep.strong_step.is_some() {
                        strong_pass.fetch_add(1,Ordering::Relaxed);
                        continue;
                    }

                    // Before calling this merely "hard", try an honest bounded
                    // exhaustive BFS of the raw Hurwitz orbit.  If the queue closes
                    // before certify_cap, this is a genuine certificate even for d>7.
                    let mut certification_timed_out = false;
                    if certify_cap > 0 {
                        let cert = certify_raw_orbit_bounded(d, ambient, seed, certify_cap, deadline);
                        certification_timed_out = cert.timed_out;
                        best.best_rank = best.best_rank.max(cert.max_rank);
                        if cert.weak_seen && best.weak_step.is_none() {
                            best.weak_step = Some(0); // only used as a boolean in reporting below
                            weak_pass.fetch_add(1,Ordering::Relaxed);
                        }
                        if cert.strong_found {
                            strong_pass.fetch_add(1,Ordering::Relaxed);
                            continue;
                        }
                        if cert.exhausted {
                            let no = certified_strong_fail.fetch_add(1,Ordering::Relaxed) + 1;
                            if !cert.weak_seen { certified_weak_fail.fetch_add(1,Ordering::Relaxed); }
                            let seed_s = tuple_string_dynamic(d,&seed).replace('\t', " ").replace('\n', " ");
                            {
                                let mut w = certified_file.lock().unwrap();
                                writeln!(w, "{}\t{}_{}\t{}\t{}\t{}\t{}",
                                    no,ambient.short(),d,cert.weak_seen,cert.states,cert.max_rank,seed_s).unwrap();
                                w.flush().unwrap();
                            }
                            eprintln!("*** CERTIFIED RAW-ORBIT COUNTEREXAMPLE #{} on {}_{}: weak={} states={} seed={}",
                                no,ambient.short(),d,cert.weak_seen,cert.states,tuple_string_dynamic(d,&seed));
                            // Deliberately CONTINUE.  An overnight campaign should
                            // collect every discovery it can instead of terminating.
                            continue;
                        }
                    }

                    hard.fetch_add(1,Ordering::Relaxed);
                    let mut w = candidate_file.lock().unwrap();
                    writeln!(w, "group={}_{} thread={} best_rank={}/{} weak_seen={} normal_restarts={} normal_steps={} deep_steps={} certify_cap={} budget_cutoff={} seed={}",
                        ambient.short(),d,tid,best.best_rank,d-1,best.weak_step.is_some(),restarts,walk_steps,deep_steps,certify_cap,(deep.timed_out || certification_timed_out),tuple_string_dynamic(d,&seed)).unwrap();
                    w.flush().unwrap();
                }
            });
        }
    });

    let summary = out_dir.join("HUNT_SUMMARY.txt");
    let mut w = BufWriter::new(File::create(&summary)?);
    writeln!(w,"group = {}_{}",ambient.short(),d)?;
    writeln!(w,"family = {:?}",family)?;
    writeln!(w,"threads = {threads_n}")?;
    writeln!(w,"time budget seconds = {max_seconds} (0 means unlimited)")?;
    writeln!(w,"bounded raw-orbit certification cap = {certify_cap} states per deep survivor")?;
    writeln!(w,"attempted seeds = {}",total.load(Ordering::Relaxed))?;
    writeln!(w,"accepted generating seeds = {}",generated.load(Ordering::Relaxed))?;
    writeln!(w,"weak success seen = {}",weak_pass.load(Ordering::Relaxed))?;
    writeln!(w,"strong success seen = {}",strong_pass.load(Ordering::Relaxed))?;
    writeln!(w,"CERTIFIED raw-orbit strong failures = {}",certified_strong_fail.load(Ordering::Relaxed))?;
    writeln!(w,"CERTIFIED raw-orbit weak failures = {}",certified_weak_fail.load(Ordering::Relaxed))?;
    writeln!(w,"hard unsolved candidates = {}",hard.load(Ordering::Relaxed))?;
    writeln!(w,"elapsed = {:.2?}",t0.elapsed())?;
    writeln!(w,"IMPORTANT: an unsolved hunt candidate is NOT a counterexample; only exhaustive orbit certification can prove failure.")?;
    w.flush()?;
    eprintln!("hunt complete: {}", summary.display());
    Ok(())
}

// ---------- Self tests ----------

fn cycle_perm(d: usize, cyc: &[usize]) -> Perm {
    let mut p = Perm::identity(d);
    if cyc.len() >= 2 {
        for i in 0..cyc.len() { p.img[cyc[i]] = cyc[(i+1)%cyc.len()] as u8; }
    }
    p
}

fn selftest() {
    for d in 4..=7 {
        let id = Perm::identity(d);
        let a = cycle_perm(d,&[0,1,2]);
        let b = cycle_perm(d,&[0,1]);
        assert_eq!(a.compose(a.inverse(d),d),id);
        assert_eq!(b.compose(b.inverse(d),d),id);
        let c = cycle_perm(d,&[0,d-1]);
        let d4 = product_one_fourth(d,a,b,c);
        let t=[a,b,c,d4];
        assert!(tuple_product_is_identity(d,&t));
        for mv in 0..6 {
            let u=hurwitz_move_dynamic(d,t,mv);
            assert!(tuple_product_is_identity(d,&u));
            let back = if mv<3 { hurwitz_move_dynamic(d,u,mv+3) } else { hurwitz_move_dynamic(d,u,mv-3) };
            assert_eq!(back,t);
        }
    }
    let s5_a = cycle_perm(5,&[0,1]);
    let s5_b = cycle_perm(5,&[0,1,2,3,4]);
    assert_eq!(subgroup_order_schreier(5,&[s5_a,s5_b]),120);
    let a5_a = cycle_perm(5,&[0,1,2]);
    let a5_b = cycle_perm(5,&[0,3,4]);
    assert_eq!(subgroup_order_schreier(5,&[a5_a,a5_b]),60);

    // New adversarial seed-family checks.  A d-cycle gives a transitive but
    // proper cyclic subgroup (for A_d only when d is odd).
    let mut rng = SplitMix64::new(0x1234_5678_9ABC_DEF0);
    let cyc_a5 = random_d_cycle(5, Ambient::Alternating, &mut rng).unwrap();
    let cyc_a5_b = perm_power(cyc_a5, 2, 5);
    assert_eq!(point_orbit_count_perms(5,cyc_a5,cyc_a5_b),1);
    assert!(!pair_contains_ad_dynamic(5,Ambient::Alternating,cyc_a5,cyc_a5_b));
    let cyc_s6 = random_d_cycle(6, Ambient::Symmetric, &mut rng).unwrap();
    let cyc_s6_b = perm_power(cyc_s6, 5, 6);
    assert_eq!(point_orbit_count_perms(6,cyc_s6,cyc_s6_b),1);
    assert!(!pair_contains_ad_dynamic(6,Ambient::Symmetric,cyc_s6,cyc_s6_b));
    assert!(random_d_cycle(6, Ambient::Alternating, &mut rng).is_none());

    // Prime-degree affine family: transitive from the start, but trapped in a
    // tiny proper affine subgroup before braid moves are applied.
    for ambient in [Ambient::Alternating, Ambient::Symmetric] {
        let (aa,bb) = random_affine_pair_prime(5, ambient, &mut rng).expect("affine prime seed");
        assert_eq!(point_orbit_count_perms(5,aa,bb),1);
        assert!(!pair_contains_ad_dynamic(5,ambient,aa,bb));
    }

    // Bounded raw-orbit certifier: a known strong A5 pair must succeed
    // immediately; the all-identity tuple has a one-state exhausted orbit.
    let id5 = Perm::identity(5);
    let strong_t = [a5_a,a5_b,id5,product_one_fourth(5,a5_a,a5_b,id5)];
    let strong_cert = certify_raw_orbit_bounded(5,Ambient::Alternating,strong_t,100,None);
    assert!(strong_cert.strong_found);
    let trivial_t = [id5,id5,id5,id5];
    let trivial_cert = certify_raw_orbit_bounded(5,Ambient::Alternating,trivial_t,10,None);
    assert!(trivial_cert.exhausted && !trivial_cert.strong_found && !trivial_cert.weak_seen);
    assert_eq!(trivial_cert.states,1);

    // Exhaustively validate the diagonal-conjugation canonicalizer on A4.
    // Python reference enumeration gives 178 conjugacy classes of triples in A4^3.
    let g = FiniteGroup::build(4, Ambient::Alternating, 1);
    let m = g.m();
    let mut canonical = HashSet::new();
    for a in 0..m as u32 {
        for b in 0..m as u32 {
            for c in 0..m as u32 {
                let k = g.canonical_key(a,b,c);
                canonical.insert(k);
                for h in 0..m as u32 {
                    let aa=g.conj_id(h,a); let bb=g.conj_id(h,b); let cc=g.conj_id(h,c);
                    assert_eq!(g.canonical_key(aa,bb,cc),k);
                }
                for mv in 0..3 {
                    let k2=g.braid_neighbor_canonical(k,mv);
                    let k3=g.braid_neighbor_canonical(k2,mv+3);
                    assert_eq!(k3,k);
                }
            }
        }
    }
    assert_eq!(canonical.len(),178);
    eprintln!("selftest passed");
}

// ---------- CLI ----------

fn get_arg(args: &[String], name: &str, default: Option<&str>) -> Result<String,String> {
    if let Some(i)=args.iter().position(|x|x==name) {
        args.get(i+1).cloned().ok_or_else(||format!("missing value after {name}"))
    } else if let Some(d)=default { Ok(d.to_string()) }
    else { Err(format!("missing required {name}")) }
}
fn has_flag(args:&[String],name:&str)->bool { args.iter().any(|x|x==name) }
fn parse_num<T:std::str::FromStr>(s:String,name:&str)->Result<T,String> where T::Err:std::fmt::Display { s.parse::<T>().map_err(|e|format!("bad {name}: {e}")) }

fn usage() {
    eprintln!(r#"hurwitz_breaker -- std-only search and certification for the four-tuple question

Commands:
  selftest

  smart-hunt --group A|S --d N --out NEW_DIRECTORY
             [--threads 8 --seconds 300 --family low-index|mixed|intransitive|passport|deform|braid-fixed]
             [--passport 3.3/3.3/3.3/3.3 --parents FILE]
             [--involution-pairs M --genus G] (braid-fixed family only)
             [--cap 200000 --candidate-seconds 10 --steps 64 --weak-only]
             [--replay FILE --replay-only --seeds 1000000000 --seed INTEGER]
    Sparse exact orbit exploration modulo S_d conjugation after short probes.
    Closed components are certified; capped/time-limited cases are candidates.

  braid-restrict --certificate FILE --out NEW_DIRECTORY
                 [--seconds 300 --cap 200000 --max-involutions 1000000]
    Restrict a known three-involution sigma_1-fixed seed to invariant subsets
    of <a,b,h>, then exhaust compatible third involutions commuting with h.
    Counts describe this construction only, never an entire degree/passport.

  rigid-glue --certificate FILE --blocks 9,9,16 --out NEW_DIRECTORY
             [--entry-order 2 --seconds 300 --max-nodes 1000000000 --cap 64]
    Glue braid-symmetry blocks and solve a^order=1, ak=ka, (a h^-1 k^-1)^3=1.
    Entry order may be 2..6; degree may be 4..34.
    Blocks may include t20:0 (four-fiber central lift, degree 20, phase 0).
  hexagon-search (--k 2|3|4|5 | --pair-blocks 16,2) --out NEW_DIRECTORY
                 [--all-h --all-types --seconds 600 --cap 200000]
    Enumerate sigma_1 relabelings and commuting perfect-matching thirds for a
    fixed pair graph. --pair-blocks searches all fourth types and relabelings.
    smart-hunt --family braid-fixed also accepts --exceptional-type 4.4.4
    to require the fourth entry's cycle type (unlisted points are fixed).
    Full generation and every claimed weak failure are checked exactly.

  orbit-map --group A|S --d N --passport TYPES --out NEW_DIRECTORY
            [--seconds 600 --seeds 1000000000 --cap 1024 --max-states 250000]
            [--seed INTEGER --replay FILE]
    Map bounded S_d-quotient neighborhoods even after finding weak witnesses.
    State overlap proves same orbit; disjoint partial maps prove nothing.
    Closed weak failures produce independently verifiable certificates.

  certify --group A|S --d N --tuple '[(), (), (), ()]' --out NEW_FILE
          [--cap 200000 --seconds 60]
    The tuple must generate the specified A_d/S_d. Writes a certificate only
    if the quotient closes without a strong pair.

  passport-exact --group A|S --d N --passport 3.3/3.3/3.3/3.3 --out NEW_DIRECTORY
                 [--start 0 --count 500000 --seconds 300 --threads 8]
    d<=11 only. Scan every g2/g3 in the specified classes with fixed g1.
    SUMMARY records scan next for resume; capped orbits remain unresolved.

  verify --certificate FILE
    Rechecks generation, pair orders, and closure of the stored certificate.

  spectrum --components A5.csv,S5.csv [--certificates DIR_OR_FILE[,DIR_OR_FILE...]]
           [--out DIR --threads 10 --full-max 256 --lanczos-steps 128
            --degree5-full-max 1024 --baseline-reps 8
            --baseline-kinds random_6_regular,random_schreier,cycle_matched_schreier,loop_matched_schreier
            --trace-order 8 --trace-probes 16 --maxcut-restarts 4 --canonical-max 24
            --switch-attempts-per-edge 20 --lanczos-checks 96,192,384
            --select-case-ids FILE --summary-only]
    Enumerates simultaneous-S_d quotient graphs, writes six braid transitions,
    genus, exact small spectra, iterative large spectra, and matched baselines.
    Summary-only omits the potentially huge per-state and per-edge tables.

  research --components A5.csv,S5.csv [--certificates DIR_OR_FILE[,DIR_OR_FILE...]]
           --out NEW_DIRECTORY [--only CASE_ID,CASE_ID --max-lift 30000]
           [--null-max 128 --null-reps 64 --full-max 32 --standard-max 256
            --subset-max-k 0 --ordered-max-k 0 --subset-cap 1000000
            --lanczos-steps 128 --trace-probes 8]
    Computes exact permutation holonomy and labeled-orbit sizes, then labeled
    spectra when feasible.  Distinct move-labeled quotient graphs each get
    cycle-matched, loop-position-matched, Schreier, and simple-regular nulls.
    Writes progress tables and explicit edge-voltage permutations.

  analyze-nulls --inputs RUN_DIR,RUN_DIR --out NEW_DIRECTORY
    Class-balanced exploratory comparison of completed research null samples.

  rayleigh --certificate FILE --ordered-k 3 --iterations 6000 --out NEW_DIRECTORY
           [--seed INTEGER | --vector EXISTING_VECTOR_FILE]
    Writes an integer test vector and exact rational Rayleigh lower bound on
    labeled lambda2 (hence an upper bound on the labeled relaxation gap).

  mechanism-audit --certificates FILE[,FILE...] --out NEW_DIRECTORY
    Compares verified weak components by pair-orbit partitions, sigma1/sigma3
    cusp classes, fixed braid moves, and seed relabeling symmetries.

  s18-fiber-census --certificate FILE --out NEW_DIRECTORY
    Exhausts degree-ten braid pairs of type 7+2+1 compatible with the
    S18 position fiber; groups them by the fixed generator's centralizer.

  For a resumable campaign: ./run_campaign.py (builds/tests, detaches, resumes jobs).

  exact --group A|S --d 7 [--threads 16] [--out runs/exact/A7]
        [--checkpoint-every 100] [--max-components 0] [--resume] [--keep-going]
        [--shard-count N --shard-index I]

    Exact modulo simultaneous conjugation.  This is enough to certify the
    yes/no property because diagonal conjugation commutes with Hurwitz moves
    and preserves transitivity / containing A_d.  Dense exact mode supports
    A4..A8 and S4..S8.  S8 requires at least 8 shards and high memory.

  hunt --group A|S --d 8 [--threads 16] [--seeds 1000]
       [--seconds 0] [--steps 50000] [--restarts 4] [--deep-multiplier 20]
       [--certify-cap 500000]
       [--family mixed|uniform|intransitive|imprimitive|cyclic|affine]
       [--seed 123456789] [--out runs/experiments/hunt-A8]

    Massive adversarial/random search. --seconds gives a wall-clock budget for
    the whole hunt (0 = unlimited).  Deep survivors get a bounded exhaustive
    raw-orbit BFS of up to --certify-cap states.  If that BFS exhausts without
    success, the failure is CERTIFIED and logged; capped survivors remain only
    candidates.

Recommended:
  RUSTFLAGS="-C target-cpu=native" cargo run --release -- selftest
  RUSTFLAGS="-C target-cpu=native" cargo run --release -- exact --group A --d 7 --threads 16 --out runs/exact/A7
  RUSTFLAGS="-C target-cpu=native" cargo run --release -- exact --group S --d 7 --threads 16 --out runs/exact/S7
  RUSTFLAGS="-C target-cpu=native" cargo run --release -- hunt --group A --d 8 --threads 16 --seeds 5000 --steps 100000 --restarts 8 --deep-multiplier 50 --family mixed --out runs/experiments/hunt-A8
"#);
}

fn real_main() -> Result<(),String> {
    let args:Vec<String>=env::args().collect();
    if args.len()<2 { usage(); return Ok(()); }
    match args[1].as_str() {
        "smart-hunt" | "orbit-map" | "passport-exact" | "braid-restrict" | "rigid-glue" | "braid-library" | "hexagon-search" | "certify" | "verify" | "spectrum" | "research" | "rayleigh" | "position-spectrum" | "braid-symmetry" | "mechanism-audit" | "s18-fiber-census" | "analyze-nulls" | "defect-census" | "certificate-profile" | "defect-mission" => smart::run(&args),
        "selftest" => { selftest(); Ok(()) }
        "exact" => {
            let ambient=Ambient::parse(&get_arg(&args,"--group",None)?)?;
            let d:usize=parse_num(get_arg(&args,"--d",None)?,"--d")?;
            let threads:usize=parse_num(get_arg(&args,"--threads",Some("8"))?,"--threads")?;
            let out=PathBuf::from(get_arg(&args,"--out",Some("runs/exact"))?);
            let ce:u64=parse_num(get_arg(&args,"--checkpoint-every",Some("100"))?,"--checkpoint-every")?;
            let resume=has_flag(&args,"--resume");
            let stop=!has_flag(&args,"--keep-going");
            let shard_count:usize=parse_num(get_arg(&args,"--shard-count",Some("1"))?,"--shard-count")?;
            let shard_index:usize=parse_num(get_arg(&args,"--shard-index",Some("0"))?,"--shard-index")?;
            let max_components:u64=parse_num(get_arg(&args,"--max-components",Some("0"))?,"--max-components")?;
            exact_quotient(d,ambient,threads,&out,ce,resume,stop,shard_count,shard_index,max_components).map_err(|e|e.to_string())
        }
        "hunt" => {
            let ambient=Ambient::parse(&get_arg(&args,"--group",None)?)?;
            let d:usize=parse_num(get_arg(&args,"--d",None)?,"--d")?;
            let threads_n:usize=parse_num(get_arg(&args,"--threads",Some("8"))?,"--threads")?;
            let seeds:u64=parse_num(get_arg(&args,"--seeds",Some("1000"))?,"--seeds")?;
            let seconds:u64=parse_num(get_arg(&args,"--seconds",Some("0"))?,"--seconds")?;
            let steps:u64=parse_num(get_arg(&args,"--steps",Some("50000"))?,"--steps")?;
            let restarts:u64=parse_num(get_arg(&args,"--restarts",Some("4"))?,"--restarts")?;
            let deep:u64=parse_num(get_arg(&args,"--deep-multiplier",Some("20"))?,"--deep-multiplier")?;
            let family=SeedFamily::parse(&get_arg(&args,"--family",Some("mixed"))?)?;
            let seed:u64=parse_num(get_arg(&args,"--seed",Some("7640891576956012809"))?,"--seed")?;
            let certify_cap:usize=parse_num(get_arg(&args,"--certify-cap",Some("500000"))?,"--certify-cap")?;
            let out=PathBuf::from(get_arg(&args,"--out",Some("runs/experiments/hunt"))?);
            hunt(d,ambient,threads_n,seeds,steps,restarts,deep,family,&out,seed,seconds,certify_cap).map_err(|e|e.to_string())
        }
        _ => { usage(); Err(format!("unknown command {:?}",args[1])) }
    }
}

fn main() {
    if let Err(e)=real_main() { eprintln!("ERROR: {e}"); std::process::exit(2); }
}
