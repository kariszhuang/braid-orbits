//! All-zero disproof attempt (Nov 8). Conjecture: a positive-genus braid orbit
//! whose group is primitive and whose defect is 0 at EVERY state is a C_p orbit.
//!
//! For each group G in groups.txt (generators as 0-based image lists): enumerate all
//! identity-free product-one 4-tuples in G with positive genus, up to G-conjugacy, with g1 a
//! class representative of minimal index among the four entries (any entry can be braided to
//! position 1) and g2 up to conjugation by C(g1). Defect at {12|34} is
//! cycles(g1g2) - orbits<g1,g2> - orbits<g3,g4> + 1. A tuple that is transitive with defect 0
//! at the three cusps of its state is closed under the six moves (states modulo simultaneous
//! S_d conjugation), stopping at the first state of positive defect. Every all-zero orbit
//! found is reported with: genus, number of states, whether its tuple generates a primitive
//! group, whether that group is abelian (hence C_p if primitive), and its order.
//! Conventions as src/main.rs: rightmost-first composition, g1 g2 g3 g4 = 1.
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

const NM: usize = 64;
type P = [u8; NM];
type T = [P; 4];

#[derive(Clone, Copy)]
struct Cx { n: usize }
impl Cx {
    fn id(&self) -> P { let mut p = [0u8; NM]; for i in 0..NM { p[i] = i as u8; } p }
    #[inline] fn mul(&self, a: &P, b: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[i] = a[b[i] as usize]; } c }
    #[inline] fn inv(&self, a: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[a[i] as usize] = i as u8; } c }
    #[inline] fn is_id(&self, a: &P) -> bool { (0..self.n).all(|i| a[i] as usize == i) }
    #[inline] fn ncyc(&self, a: &P) -> usize { let mut seen: u64 = 0; let mut k = 0;
        for s in 0..self.n { if seen >> s & 1 == 1 { continue; } k += 1; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u64 << x; x = a[x] as usize; } } k }
    #[inline] fn ind(&self, a: &P) -> usize { self.n - self.ncyc(a) }
    /// number of orbits of <gens>
    #[inline] fn norb(&self, gens: &[&P]) -> usize {
        let mut seen: u64 = 0; let mut k = 0; let mut st = [0u8; NM];
        for s in 0..self.n { if seen >> s & 1 == 1 { continue; } k += 1; seen |= 1u64 << s; let mut sp = 1; st[0] = s as u8;
            while sp > 0 { sp -= 1; let x = st[sp] as usize; for g in gens { let y = g[x] as usize; if seen >> y & 1 == 0 { seen |= 1u64 << y; st[sp] = y as u8; sp += 1; } } } }
        k
    }
    fn fmt(&self, p: &P) -> String {
        let mut seen = [false; NM]; let mut s = String::new();
        for i in 0..self.n { if seen[i] { continue; } let mut c = vec![]; let mut x = i;
            while !seen[x] { seen[x] = true; c.push(x + 1); x = p[x] as usize; }
            if c.len() > 1 { s.push('('); s.push_str(&c.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ")); s.push(')'); } }
        if s.is_empty() { "()".into() } else { s }
    }
}
/// NEON kernels for permutations of up to 64 points (aarch64): a permutation is four 16-byte
/// registers and composition is four `tbl4` lookups. Entries beyond n are fixed points.
#[cfg(target_arch = "aarch64")]
mod n64 {
    use std::arch::aarch64::*;
    use super::P;
    #[derive(Clone, Copy)] pub struct V(pub [uint8x16_t; 4]);
    #[inline(always)] pub unsafe fn ld(p: &P) -> V { let q = p.as_ptr(); V([vld1q_u8(q), vld1q_u8(q.add(16)), vld1q_u8(q.add(32)), vld1q_u8(q.add(48))]) }
    static ID: [u8; 64] = { let mut a = [0u8; 64]; let mut i = 0; while i < 64 { a[i] = i as u8; i += 1; } a };
    #[inline(always)] pub unsafe fn idv() -> V { ld(&ID) }
    /// (a o b)[i] = a[b[i]]
    #[inline(always)] pub unsafe fn comp(a: V, b: V) -> V {
        let t = uint8x16x4_t(a.0[0], a.0[1], a.0[2], a.0[3]);
        V([vqtbl4q_u8(t, b.0[0]), vqtbl4q_u8(t, b.0[1]), vqtbl4q_u8(t, b.0[2]), vqtbl4q_u8(t, b.0[3])])
    }
    #[inline(always)] unsafe fn vmin(a: V, b: V) -> V { V([vminq_u8(a.0[0], b.0[0]), vminq_u8(a.0[1], b.0[1]), vminq_u8(a.0[2], b.0[2]), vminq_u8(a.0[3], b.0[3])]) }
    #[inline(always)] unsafe fn same(a: V, b: V) -> bool {
        let x = vorrq_u8(vorrq_u8(veorq_u8(a.0[0], b.0[0]), veorq_u8(a.0[1], b.0[1])), vorrq_u8(veorq_u8(a.0[2], b.0[2]), veorq_u8(a.0[3], b.0[3]))); vmaxvq_u8(x) == 0 }
    /// number of i < n with l[i] == i (l = least-element labels)
    #[inline(always)] unsafe fn leaders(l: V, n: usize) -> usize {
        let id = idv(); let one = vdupq_n_u8(1); let mut c = 0usize;
        for k in 0..4 { c += vaddvq_u8(vandq_u8(vceqq_u8(l.0[k], id.0[k]), one)) as usize; }
        c - (64 - n)
    }
    /// number of cycles of q
    #[inline(always)] pub unsafe fn ncyc(q: V, n: usize) -> usize {
        let mut m = idv(); let mut r = q;
        for _ in 0..6 { m = vmin(m, comp(m, r)); r = comp(r, r); }
        leaders(m, n)
    }
    /// number of orbits of <a, b>
    #[inline(always)] pub unsafe fn norb(a: V, b: V, n: usize) -> usize {
        let mut l = idv();
        loop { let n1 = vmin(l, vmin(comp(l, a), comp(l, b))); let n2 = comp(n1, n1); if same(n2, l) { return leaders(l, n); } l = n2; }
    }
}

#[inline] fn defect(cx: &Cx, t: &T) -> i64 {
    let p = cx.mul(&t[0], &t[1]);
    cx.ncyc(&p) as i64 - cx.norb(&[&t[0], &t[1]]) as i64 - cx.norb(&[&t[2], &t[3]]) as i64 + 1
}
fn mv(cx: &Cx, t: &T, m: usize) -> T {
    let j = m % 3; let a = t[j]; let b = t[j + 1]; let mut u = *t;
    if m < 3 { u[j] = cx.mul(&a, &cx.mul(&b, &cx.inv(&a))); u[j + 1] = a; } else { u[j] = b; u[j + 1] = cx.mul(&cx.inv(&b), &cx.mul(&a, &b)); }
    u
}
fn rot(t: &T) -> T { [t[1], t[2], t[3], t[0]] }
fn canon(cx: &Cx, t: &T) -> Vec<u8> {
    let n = cx.n;
    let mut ls = [[0u8; NM]; 4];
    for j in 0..4 { let a = &t[j]; let mut seen: u64 = 0;
        for s in 0..n { if seen >> s & 1 == 1 { continue; } let mut l = 0u8; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u64 << x; x = a[x] as usize; l += 1; }
            let mut x = s; for _ in 0..l { ls[j][x] = l; x = a[x] as usize; } } }
    let mut sig = [0u32; NM]; for r in 0..n { sig[r] = ((ls[0][r] as u32) << 24) | ((ls[1][r] as u32) << 16) | ((ls[2][r] as u32) << 8) | ls[3][r] as u32; }
    let mut srt = sig; srt[..n].sort_unstable();
    let (mut best_sig, mut best_cnt) = (0u32, usize::MAX); let mut i = 0;
    while i < n { let mut j = i; while j < n && srt[j] == srt[i] { j += 1; } if j - i < best_cnt { best_cnt = j - i; best_sig = srt[i]; } i = j; }
    let mut best = vec![u8::MAX; 3 * n]; let mut have = false;
    for r in 0..n {
        if sig[r] != best_sig { continue; }
        let mut lab = [u8::MAX; NM]; let mut order = [0u8; NM]; let mut len = 1; lab[r] = 0; order[0] = r as u8; let mut q = 0;
        while q < len { let x = order[q] as usize; q += 1; for g in &t[..3] { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y as u8; len += 1; } } }
        let mut state = if have { 0i8 } else { -1 }; let mut key = vec![0u8; 3 * n];
        'o: for j in 0..3 { for k in 0..n { let v = lab[t[j][order[k] as usize] as usize]; key[j * n + k] = v;
            if state == 0 { let b = best[j * n + k]; if v < b { state = -1; } else if v > b { state = 1; break 'o; } } } }
        if state == -1 { best = key; have = true; }
    }
    best
}
fn unpack(cx: &Cx, k: &[u8]) -> T {
    let n = cx.n; let mut t = [cx.id(); 4];
    for j in 0..3 { for i in 0..n { t[j][i] = k[j * n + i]; } }
    t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2])));
    t
}
fn h128(k: &[u8]) -> u128 {
    let (mut a, mut b) = (0xcbf2_9ce4_8422_2325u64, 0x9e37_79b9_7f4a_7c15u64);
    for &x in k { a = (a ^ x as u64).wrapping_mul(0x0000_0100_0000_01b3); b = (b.rotate_left(7) ^ x as u64).wrapping_mul(0x517c_c1b7_2722_0a95); }
    ((a as u128) << 64) | (b as u128)
}
/// minimal nontrivial blocks: true if <gens> (transitive) is primitive
fn primitive(cx: &Cx, gens: &[P]) -> bool {
    let n = cx.n;
    for j in 1..n {
        let mut par: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, x: usize) -> usize { let mut r = x; while p[r] != r { r = p[r]; } let mut y = x; while p[y] != r { let z = p[y]; p[y] = r; y = z; } r }
        let mut q = vec![(0usize, j)];
        while let Some((x, y)) = q.pop() { let (a, b) = (find(&mut par, x), find(&mut par, y)); if a == b { continue; } par[a] = b; for g in gens { q.push((g[x] as usize, g[y] as usize)); } }
        let r0 = find(&mut par, 0); if (0..n).filter(|&x| find(&mut par, x) == r0).count() < n { return false; }
    }
    true
}
fn group_order(cx: &Cx, gens: &[P], cap: usize) -> usize {
    let mut seen: HashSet<P> = HashSet::new(); let id = cx.id(); seen.insert(id); let mut q = vec![id];
    while let Some(x) = q.pop() { for g in gens { let y = cx.mul(g, &x); if seen.insert(y) { if seen.len() > cap { return 0; } q.push(y); } } }
    seen.len()
}

struct Found { states: usize }
fn close(cx: &Cx, t0: &T, cap: usize) -> (bool, Vec<u128>, usize) {
    let root = canon(cx, t0); let mut idx: HashSet<u128> = HashSet::new(); idx.insert(h128(&root));
    let mut st = vec![root]; let mut i = 0;
    while i < st.len() {
        let u = unpack(cx, &st[i]); i += 1;
        if defect(cx, &u) != 0 { return (false, st.iter().map(|k| h128(k)).collect(), st.len()); }
        for m in 0..6 { let v = mv(cx, &u, m); let k = canon(cx, &v); if idx.insert(h128(&k)) { st.push(k); if st.len() > cap { return (false, vec![], st.len()); } } }
    }
    let hs = st.iter().map(|k| h128(k)).collect(); let n = st.len(); (true, hs, n)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let file = get("--groups").expect("--groups FILE");
    let only = get("--only");
    let maxo: usize = get("--max-order").map(|x| x.parse().unwrap()).unwrap_or(120_000);
    let threads: usize = get("--threads").map(|x| x.parse().unwrap()).unwrap_or(10);
    let min_genus: usize = get("--min-genus").map(|x| x.parse().unwrap()).unwrap_or(1);
    let cap: usize = get("--cap").map(|x| x.parse().unwrap()).unwrap_or(2_000_000);
    let out = Mutex::new(std::fs::OpenOptions::new().create(true).append(true).open(get("--out").expect("--out FILE")).unwrap());
    for line in std::fs::read_to_string(&file).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect(); if f.len() < 3 { continue; }
        let name = f[0]; if let Some(o) = &only { if !name.contains(o.as_str()) { continue; } }
        let n: usize = f[1].parse().unwrap(); assert!(n <= NM); let cx = Cx { n };
        let gens: Vec<P> = f[2].split(';').map(|g| { let mut p = cx.id(); for (i, x) in g.split(' ').enumerate() { p[i] = x.parse().unwrap(); } p }).collect();
        let t0 = std::time::Instant::now();
        if cx.norb(&gens.iter().collect::<Vec<_>>()) != 1 || !primitive(&cx, &gens) {
            writeln!(out.lock().unwrap(), "{{\"group\":{:?},\"d\":{},\"skipped\":\"not primitive\"}}", name, n).unwrap(); eprintln!("{} d{}: not primitive, skipped", name, n); continue; }
        // elements
        let mut elems: Vec<P> = vec![cx.id()]; let mut ix: HashMap<P, u32> = HashMap::new(); ix.insert(cx.id(), 0); let mut qi = 0; let mut too_big = false;
        while qi < elems.len() { let x = elems[qi]; qi += 1; for g in &gens { let y = cx.mul(g, &x); if !ix.contains_key(&y) { ix.insert(y, elems.len() as u32); elems.push(y); if elems.len() > maxo { too_big = true; break; } } } if too_big { break; } }
        if too_big { writeln!(out.lock().unwrap(), "{{\"group\":{:?},\"d\":{},\"skipped\":\"order > {}\"}}", name, n, maxo).unwrap(); eprintln!("{} d{}: order > {}, skipped", name, n, maxo); continue; }
        let ord = elems.len();
        let ind: Vec<u8> = elems.iter().map(|e| cx.ind(e) as u8).collect();
        let ginv: Vec<P> = gens.iter().map(|g| cx.inv(g)).collect();
        // conjugacy classes (union-find under conjugation by generators)
        let mut par: Vec<u32> = (0..ord as u32).collect();
        fn find(p: &mut Vec<u32>, x: u32) -> u32 { let mut r = x; while p[r as usize] != r { r = p[r as usize]; } let mut y = x; while p[y as usize] != r { let z = p[y as usize]; p[y as usize] = r; y = z; } r }
        for e in 0..ord { for (g, gi) in gens.iter().zip(&ginv) { let c = cx.mul(g, &cx.mul(&elems[e], gi)); let j = ix[&c]; let (a, b) = (find(&mut par, e as u32), find(&mut par, j)); if a != b { par[a as usize] = b; } } }
        let mut reps: Vec<usize> = (0..ord).filter(|&e| find(&mut par, e as u32) == e as u32 && ind[e] > 0).collect();
        reps.sort_by_key(|&e| ind[e]);
        let nclasses = reps.len() + 1;
        // jobs: (class rep r, list of g2 representatives modulo C(r))
        let mut rng: u64 = 0x9e3779b97f4a7c15;
        let mut jobs: Vec<(usize, Vec<u32>)> = vec![];
        for &r in &reps {
            let cent: Vec<usize> = (0..ord).filter(|&c| cx.mul(&elems[c], &elems[r]) == cx.mul(&elems[r], &elems[c])).collect();
            let mut cg: Vec<usize> = vec![]; for _ in 0..8 { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; cg.push(cent[(rng % cent.len() as u64) as usize]); }
            let cgi: Vec<P> = cg.iter().map(|&c| cx.inv(&elems[c])).collect();
            let mut p2: Vec<u32> = (0..ord as u32).collect();
            for e in 0..ord { if ind[e] < ind[r] { continue; } for (k, &c) in cg.iter().enumerate() { let y = cx.mul(&elems[c], &cx.mul(&elems[e], &cgi[k])); let j = ix[&y]; let (a, b) = (find(&mut p2, e as u32), find(&mut p2, j)); if a != b { p2[a as usize] = b; } } }
            let g2s: Vec<u32> = (0..ord).filter(|&e| ind[e] >= ind[r] && find(&mut p2, e as u32) == e as u32).map(|e| e as u32).collect();
            for ch in g2s.chunks(64) { jobs.push((r, ch.to_vec())); }
        }
        let (tested, d0, pass3, closures, cyc, prim, imprim) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0));
        let done: Mutex<HashSet<u128>> = Mutex::new(HashSet::new()); let probe_rej = AtomicU64::new(0);
        let next = AtomicUsize::new(0);
        let elems_r = &elems; let ind_r = &ind;
        std::thread::scope(|s| { for _ in 0..threads { s.spawn(|| { loop {
            let j = next.fetch_add(1, Ordering::Relaxed); if j >= jobs.len() { break; }
            let (r, ref g2s) = jobs[j]; let g1 = elems_r[r]; let i1 = ind_r[r] as usize;
            let (mut lt, mut ld0, mut lp3) = (0u64, 0u64, 0u64);
            for &k2 in g2s {
                let g2 = elems_r[k2 as usize]; let p12 = cx.mul(&g1, &g2);
                let c12 = cx.ncyc(&p12) as i64; let v12 = cx.norb(&[&g1, &g2]) as i64;
                #[cfg(target_arch = "aarch64")] let p12v = unsafe { n64::ld(&p12) };
                let i2 = ind_r[k2 as usize] as usize;
                for k3 in 0..elems_r.len() {
                    let i3 = ind_r[k3] as usize; if i3 < i1 { continue; }
                    let g3 = &elems_r[k3];
                    // i4 = ind(g4) = ind(g1 g2 g3); orbits <g3,g4> = orbits <g3, g1 g2> (g4 = (g1g2g3)^-1)
                    #[cfg(target_arch = "aarch64")]
                    let (i4, v34) = unsafe { let g3v = n64::ld(g3); let q = n64::comp(p12v, g3v); let i4 = n - n64::ncyc(q, n);
                        if i4 < i1 || i1 + i2 + i3 + i4 < 2 * n - 2 + 2 * min_genus { continue; }
                        (i4, n64::norb(g3v, p12v, n)) };
                    #[cfg(not(target_arch = "aarch64"))]
                    let (i4, v34) = { let q = cx.mul(&p12, g3); let i4 = cx.ind(&q); if i4 < i1 || i1 + i2 + i3 + i4 < 2 * n - 2 + 2 * min_genus { continue; } (i4, cx.norb(&[g3, &p12])) };
                    let _ = i4;
                    lt += 1;
                    if c12 - v12 - v34 as i64 + 1 != 0 { continue; }
                    let g4 = cx.inv(&cx.mul(&p12, g3));
                    if cx.norb(&[&g1, &g2, g3]) != 1 { continue; }
                    ld0 += 1;
                    let t: T = [g1, g2, *g3, g4];
                    if defect(&cx, &rot(&t)) != 0 || defect(&cx, &mv(&cx, &t, 1)) != 0 || defect(&cx, &mv(&cx, &t, 4)) != 0 { continue; }
                    lp3 += 1;
                    // probe BEFORE the canonical form: 64 random moves; any state with positive
                    // defect rejects the orbit (sound). Seeded from the tuple bytes.
                    { let mut r = { let mut b = Vec::with_capacity(3 * n); for q in &t[..3] { b.extend_from_slice(&q[..n]); } h128(&b) as u64 | 1 }; let mut u = t; let mut bad = false;
                      for _ in 0..64 { r ^= r << 13; r ^= r >> 7; r ^= r << 17; u = mv(&cx, &u, (r % 6) as usize); if defect(&cx, &u) != 0 { bad = true; break; } }
                      if bad { probe_rej.fetch_add(1, Ordering::Relaxed); continue; } }
                    let key = h128(&canon(&cx, &t));
                    if done.lock().unwrap().contains(&key) { continue; }
                    closures.fetch_add(1, Ordering::Relaxed);
                    let (ok, hs, ns) = close(&cx, &t, cap);
                    { let mut dn = done.lock().unwrap(); for h in &hs { dn.insert(*h); } }
                    if ns > cap { writeln!(out.lock().unwrap(), "{{\"group\":{:?},\"d\":{},\"capped\":true,\"tuple\":{:?}}}", name, n, t.iter().map(|p| cx.fmt(p)).collect::<Vec<_>>()).unwrap(); continue; }
                    if !ok { continue; }
                    let tv: Vec<P> = t.to_vec();
                    let pr = primitive(&cx, &tv);
                    let ab = (0..4).all(|a| (0..4).all(|b| cx.mul(&t[a], &t[b]) == cx.mul(&t[b], &t[a])));
                    let go = group_order(&cx, &tv, 2_000_000);
                    let genus = (t.iter().map(|p| cx.ind(p)).sum::<usize>() as i64 - 2 * n as i64 + 2) / 2;
                    if ab && pr { cyc.fetch_add(1, Ordering::Relaxed); } else if pr { prim.fetch_add(1, Ordering::Relaxed); } else { imprim.fetch_add(1, Ordering::Relaxed); }
                    let line = format!("{{\"group\":{:?},\"d\":{},\"genus\":{},\"states\":{},\"primitive\":{},\"abelian\":{},\"gen_order\":{},\"tuple\":{:?}}}", name, n, genus, ns, pr, ab, go, t.iter().map(|p| cx.fmt(p)).collect::<Vec<_>>());
                    if pr && !ab { eprintln!("!!! PRIMITIVE NON-ABELIAN ALL-ZERO ORBIT: {}", line); }
                    writeln!(out.lock().unwrap(), "{}", line).unwrap();
                }
            }
            tested.fetch_add(lt, Ordering::Relaxed); d0.fetch_add(ld0, Ordering::Relaxed); pass3.fetch_add(lp3, Ordering::Relaxed);
        } }); } });
        let ld = |a: &AtomicU64| a.load(Ordering::Relaxed);
        let summ = format!("{{\"group\":{:?},\"d\":{},\"order\":{},\"classes\":{},\"tuples_tested\":{},\"defect0_first_cusp_transitive\":{},\"defect0_three_cusps\":{},\"probe_rejects\":{},\"closures\":{},\"allzero_cyclic\":{},\"allzero_primitive_nonabelian\":{},\"allzero_imprimitive_subgroup\":{},\"seconds\":{:.1},\"complete\":true}}",
            name, n, ord, nclasses, ld(&tested), ld(&d0), ld(&pass3), ld(&probe_rej), ld(&closures), ld(&cyc), ld(&prim), ld(&imprim), t0.elapsed().as_secs_f64());
        eprintln!("{}", summ); writeln!(out.lock().unwrap(), "{}", summ).unwrap();
    }
}
