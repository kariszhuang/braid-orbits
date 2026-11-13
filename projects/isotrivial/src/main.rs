//! isotrivial: decide whether the family of curves of a braid orbit is isotrivial.
//!
//! Input: a product-one tuple (g1, g2, g3, g4) of permutations (or g1, g2, g3).
//! The program
//!   1. checks the tuple and computes the genus g (Riemann-Hurwitz);
//!   2. closes the braid orbit (states up to simultaneous conjugation);
//!   3. checks that every state is one-piece (defect 0, at most one piece of positive genus);
//!   4. computes the Hodge degree per sheet exactly:
//!        12 L = d - sum_i rho(g_i) + 3 * mean over states of rho(g1 g2),
//!      where rho(p) = sum over the cycles c of p of 1/|c|;
//!   5. decides: not one-piece => some limit curve is singular => not isotrivial;
//!      one-piece and g >= 2 => isotrivial iff L = 0.
//!
//! Conventions: points are 1..d, products are composed right to left, (x y)(i) = x(y(i)).
//! No dependencies.

use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};
use std::time::Instant;

type Perm = Vec<u8>;

// ---------- permutations ----------

fn compose(a: &Perm, b: &Perm) -> Perm {
    b.iter().map(|&x| a[x as usize]).collect()
}

fn inverse(a: &Perm) -> Perm {
    let mut q = vec![0u8; a.len()];
    for (i, &x) in a.iter().enumerate() {
        q[x as usize] = i as u8;
    }
    q
}

fn cycle_lengths(p: &Perm) -> Vec<usize> {
    let mut seen = vec![false; p.len()];
    let mut out = Vec::new();
    for s in 0..p.len() {
        if seen[s] {
            continue;
        }
        let (mut x, mut l) = (s, 0);
        while !seen[x] {
            seen[x] = true;
            x = p[x] as usize;
            l += 1;
        }
        out.push(l);
    }
    out
}

fn index(p: &Perm) -> usize {
    p.len() - cycle_lengths(p).len()
}

fn is_identity(p: &Perm) -> bool {
    p.iter().enumerate().all(|(i, &x)| i == x as usize)
}

fn transitive(gens: &[&Perm]) -> bool {
    let d = gens[0].len();
    let mut seen = vec![false; d];
    let mut stack = vec![0usize];
    seen[0] = true;
    let mut count = 1;
    while let Some(x) = stack.pop() {
        for g in gens {
            let y = g[x] as usize;
            if !seen[y] {
                seen[y] = true;
                count += 1;
                stack.push(y);
            }
        }
    }
    count == d
}

// ---------- exact rationals ----------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Q {
    n: i128,
    d: i128,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

impl Q {
    fn new(n: i128, d: i128) -> Q {
        let g = gcd(n, d).max(1);
        let s = if d < 0 { -1 } else { 1 };
        Q { n: s * n / g, d: s * d / g }
    }
    fn int(n: i128) -> Q {
        Q { n, d: 1 }
    }
    fn add(self, o: Q) -> Q {
        let g = gcd(self.d, o.d);
        let l = self.d / g * o.d;
        Q::new(self.n * (l / self.d) + o.n * (l / o.d), l)
    }
    fn neg(self) -> Q {
        Q { n: -self.n, d: self.d }
    }
    fn mul_int(self, k: i128) -> Q {
        Q::new(self.n * k, self.d)
    }
    fn div_int(self, k: i128) -> Q {
        Q::new(self.n, self.d * k)
    }
    fn sign(self) -> i32 {
        self.n.signum() as i32
    }
}

impl std::fmt::Display for Q {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.d == 1 {
            write!(f, "{}", self.n)
        } else {
            write!(f, "{}/{}", self.n, self.d)
        }
    }
}

/// rho summed from a histogram: counts[l] = number of cycles of length l.
fn rho_from_counts(counts: &[u64]) -> Q {
    let mut r = Q::int(0);
    for (l, &c) in counts.iter().enumerate() {
        if c > 0 {
            r = r.add(Q::new(c as i128, l as i128));
        }
    }
    r
}

fn rho(p: &Perm) -> Q {
    let mut counts = vec![0u64; p.len() + 1];
    for l in cycle_lengths(p) {
        counts[l] += 1;
    }
    rho_from_counts(&counts)
}

// ---------- canonical form up to simultaneous conjugation ----------
//
// For a root r, number the points in the order a breadth-first search from r reaches them
// through g1, g2, g3, and record (label of g1 x, label of g2 x, label of g3 x) for the points x in
// that order. The smallest such record over all roots is the key. Only roots with the smallest
// conjugation-invariant signature (cycle lengths through the point in g1..g4) are tried; this set
// of roots is itself invariant, so the key is still a canonical form.

const NONE: u8 = u8::MAX;

struct Canon {
    lab: Vec<u8>,
    order: Vec<u8>,
    cur: Vec<u8>,
    best: Vec<u8>,
}

impl Canon {
    fn new(d: usize) -> Canon {
        Canon { lab: vec![NONE; d], order: Vec::with_capacity(d), cur: Vec::with_capacity(3 * d), best: Vec::with_capacity(3 * d) }
    }

    fn key(&mut self, t: &[Perm; 4]) -> Vec<u8> {
        let d = t[0].len();
        // signature per point
        let mut sig = vec![[0u8; 4]; d];
        for (k, p) in t.iter().enumerate() {
            let mut seen = vec![false; d];
            for s in 0..d {
                if seen[s] {
                    continue;
                }
                let mut cyc = Vec::new();
                let mut x = s;
                while !seen[x] {
                    seen[x] = true;
                    cyc.push(x);
                    x = p[x] as usize;
                }
                for &y in &cyc {
                    sig[y][k] = cyc.len() as u8;
                }
            }
        }
        let min_sig = *sig.iter().min().unwrap();
        self.best.clear();
        let mut have_best = false;
        for r in 0..d {
            if sig[r] != min_sig {
                continue;
            }
            for v in self.lab.iter_mut() {
                *v = NONE;
            }
            self.order.clear();
            self.cur.clear();
            self.lab[r] = 0;
            self.order.push(r as u8);
            // 0: equal so far, -1: already smaller than best, (abort if larger)
            let mut state = if have_best { 0i8 } else { -1 };
            let mut aborted = false;
            let mut i = 0;
            'outer: while i < self.order.len() {
                let x = self.order[i] as usize;
                i += 1;
                for p in &t[..3] {
                    let y = p[x] as usize;
                    if self.lab[y] == NONE {
                        self.lab[y] = self.order.len() as u8;
                        self.order.push(y as u8);
                    }
                    let v = self.lab[y];
                    if state == 0 {
                        let b = self.best[self.cur.len()];
                        if v > b {
                            aborted = true;
                            break 'outer;
                        }
                        if v < b {
                            state = -1;
                        }
                    }
                    self.cur.push(v);
                }
            }
            if !aborted && state == -1 {
                std::mem::swap(&mut self.best, &mut self.cur);
                have_best = true;
            }
        }
        self.best.clone()
    }
}

fn unpack(key: &[u8]) -> [Perm; 4] {
    let d = key.len() / 3;
    let mut g = [vec![0u8; d], vec![0u8; d], vec![0u8; d]];
    for i in 0..d {
        for j in 0..3 {
            g[j][i] = key[3 * i + j];
        }
    }
    let g4 = inverse(&compose(&g[0], &compose(&g[1], &g[2])));
    let [a, b, c] = g;
    [a, b, c, g4]
}

fn moves(t: &[Perm; 4]) -> Vec<[Perm; 4]> {
    let mut out = Vec::with_capacity(6);
    for j in 0..3 {
        let (a, b) = (&t[j], &t[j + 1]);
        let ai = inverse(a);
        let bi = inverse(b);
        let mut u = t.clone();
        u[j] = compose(a, &compose(b, &ai));
        u[j + 1] = a.clone();
        out.push(u);
        let mut v = t.clone();
        v[j] = b.clone();
        v[j + 1] = compose(&bi, &compose(a, b));
        out.push(v);
    }
    out
}

// ---------- one-piece test at the collision of positions 1, 2 ----------

fn find(par: &mut [usize], mut x: usize) -> usize {
    while par[x] != x {
        par[x] = par[par[x]];
        x = par[x];
    }
    x
}

/// Returns (defect, number of pieces of positive genus).
fn collision(t: &[Perm; 4]) -> (i64, usize) {
    let d = t[0].len();
    let p = compose(&t[0], &t[1]);
    let nodes = cycle_lengths(&p).len() as i64;
    let mut pos = 0;
    let mut pieces = 0i64;
    for (x, y) in [(&t[0], &t[1]), (&t[2], &t[3])] {
        let mut par: Vec<usize> = (0..d).collect();
        for i in 0..d {
            for g in [x, y] {
                let (a, b) = (find(&mut par, i), find(&mut par, g[i] as usize));
                if a != b {
                    par[a] = b;
                }
            }
        }
        let xy = compose(x, y);
        let mut size = vec![0i64; d];
        let mut cyc = vec![0i64; d];
        for i in 0..d {
            let r = find(&mut par, i);
            size[r] += 1;
        }
        for q in [x, y, &xy] {
            let mut seen = vec![false; d];
            for s in 0..d {
                if seen[s] {
                    continue;
                }
                let mut z = s;
                while !seen[z] {
                    seen[z] = true;
                    z = q[z] as usize;
                }
                let r = find(&mut par, s);
                cyc[r] += 1;
            }
        }
        for r in 0..d {
            if size[r] > 0 {
                pieces += 1;
                // 2 g(O) = |O| + 2 - c_O(x) - c_O(y) - c_O(xy)
                if size[r] + 2 - cyc[r] > 0 {
                    pos += 1;
                }
            }
        }
    }
    (nodes - pieces + 1, pos)
}

fn additive(t: &[Perm; 4]) -> bool {
    let z = index(&compose(&t[0], &t[1]));
    z == index(&t[0]) + index(&t[1]) || z == index(&t[2]) + index(&t[3])
}

// ---------- input ----------

fn parse_tuple(s: &str) -> Result<[Perm; 4], String> {
    let parts: Vec<&str> = s.split(';').map(|x| x.trim()).filter(|x| !x.is_empty()).collect();
    if parts.len() != 3 && parts.len() != 4 {
        return Err(format!("expected 3 or 4 permutations separated by ';', got {}", parts.len()));
    }
    let mut cycles_all: Vec<Vec<Vec<usize>>> = Vec::new();
    let mut d = 0;
    for part in &parts {
        let mut cycles = Vec::new();
        let mut rest = *part;
        while let Some(open) = rest.find('(') {
            let close = rest[open..].find(')').ok_or("missing ')'")? + open;
            let pts: Result<Vec<usize>, _> =
                rest[open + 1..close].split(|c: char| c == ' ' || c == ',').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>()).collect();
            let pts = pts.map_err(|e| format!("bad point in {}: {}", part, e))?;
            for &x in &pts {
                if x == 0 {
                    return Err("points are numbered from 1".into());
                }
                d = d.max(x);
            }
            cycles.push(pts);
            rest = &rest[close + 1..];
        }
        if !rest.trim().is_empty() {
            return Err(format!("unexpected text '{}'", rest.trim()));
        }
        cycles_all.push(cycles);
    }
    if d > 250 {
        return Err("degree above 250 is not supported".into());
    }
    let mut perms = Vec::new();
    for cycles in &cycles_all {
        let mut p: Perm = (0..d as u8).collect();
        let mut used = vec![false; d];
        for c in cycles {
            for (k, &x) in c.iter().enumerate() {
                if used[x - 1] {
                    return Err(format!("point {} appears twice in one permutation", x));
                }
                used[x - 1] = true;
                p[x - 1] = (c[(k + 1) % c.len()] - 1) as u8;
            }
        }
        perms.push(p);
    }
    if perms.len() == 3 {
        let g4 = inverse(&compose(&perms[0], &compose(&perms[1], &perms[2])));
        perms.push(g4);
    }
    let t: [Perm; 4] = [perms[0].clone(), perms[1].clone(), perms[2].clone(), perms[3].clone()];
    let prod = compose(&t[0], &compose(&t[1], &compose(&t[2], &t[3])));
    if is_identity(&prod) {
        return Ok(t);
    }
    let rev = compose(&t[3], &compose(&t[2], &compose(&t[1], &t[0])));
    if is_identity(&rev) {
        // product 1 in the other order: use the reversed tuple (the same family, mirrored)
        return Ok([t[3].clone(), t[2].clone(), t[1].clone(), t[0].clone()]);
    }
    Err("the product g1 g2 g3 g4 is not 1 (in either order)".into())
}

// ---------- the test ----------

struct Report {
    degree: usize,
    genus: i64,
    states: usize,
    one_piece: bool,
    all_smooth: bool,
    first_bad_defect: i64,
    first_bad_positive: usize,
    avg: [Q; 3],
    twelve_l: Option<Q>,
    verdict: &'static str,
    reason: &'static str,
    seconds: f64,
}

fn run(t: [Perm; 4], max_states: usize, quiet: bool) -> Result<Report, String> {
    let start = Instant::now();
    let d = t[0].len();
    if t.iter().any(is_identity) {
        return Err("an entry is the identity".into());
    }
    if !transitive(&[&t[0], &t[1], &t[2]]) {
        return Err("the tuple is not transitive".into());
    }
    let sum_ind: i64 = t.iter().map(|p| index(p) as i64).sum();
    let genus = (2 - 2 * d as i64 + sum_ind) / 2;
    let s_total = t.iter().fold(Q::int(0), |acc, p| acc.add(rho(p)));

    // close the orbit
    let mut canon = Canon::new(d);
    let mut seen: HashSet<Vec<u8>, BuildHasherDefault<Fnv>> = HashSet::default();
    let mut queue: Vec<Vec<u8>> = Vec::new();
    let root = canon.key(&t);
    seen.insert(root.clone());
    queue.push(root);
    let mut counts = [vec![0u64; d + 1], vec![0u64; d + 1], vec![0u64; d + 1]];
    let (mut one_piece, mut all_smooth) = (true, true);
    let (mut bad_def, mut bad_pos) = (0i64, 0usize);
    let mut i = 0;
    while i < queue.len() {
        let u = unpack(&queue[i]);
        i += 1;
        // per-state data
        let p1 = compose(&u[0], &u[1]);
        let p2 = compose(&u[1], &u[2]);
        let p3 = compose(&u[0], &compose(&u[1], &compose(&u[2], &inverse(&u[1]))));
        for (k, p) in [p1, p2, p3].iter().enumerate() {
            for l in cycle_lengths(p) {
                counts[k][l] += 1;
            }
        }
        let (def, pos) = collision(&u);
        if one_piece && (def != 0 || pos > 1) {
            one_piece = false;
            bad_def = def;
            bad_pos = pos;
        }
        if !additive(&u) {
            all_smooth = false;
        }
        for v in moves(&u) {
            let k = canon.key(&v);
            if !seen.contains(&k) {
                if seen.len() >= max_states {
                    return Err(format!("orbit has more than {} states (raise --max-states)", max_states));
                }
                seen.insert(k.clone());
                queue.push(k);
            }
        }
        if !quiet && i % 100_000 == 0 {
            eprintln!("  ... {} states processed, {} found", i, queue.len());
        }
    }
    let n = queue.len();
    let avg = [
        rho_from_counts(&counts[0]).div_int(n as i128),
        rho_from_counts(&counts[1]).div_int(n as i128),
        rho_from_counts(&counts[2]).div_int(n as i128),
    ];
    let twelve_l = Q::int(d as i128).add(s_total.neg()).add(avg[0].mul_int(3));

    let (verdict, reason, tl) = if genus <= 0 {
        ("isotrivial", "genus 0: every curve is the projective line", None)
    } else if !one_piece {
        ("not isotrivial", "some state is not one-piece, so some limit curve is singular", None)
    } else if genus == 1 {
        ("isotrivial", "genus 1 and one-piece: the j-invariant of a complete smooth family is constant", Some(twelve_l))
    } else if twelve_l.sign() > 0 {
        ("not isotrivial", "one-piece and L > 0", Some(twelve_l))
    } else if twelve_l.sign() == 0 {
        ("isotrivial", "one-piece and L = 0", Some(twelve_l))
    } else {
        ("error", "L < 0, which is impossible for a one-piece orbit", Some(twelve_l))
    };
    Ok(Report {
        degree: d,
        genus,
        states: n,
        one_piece,
        all_smooth,
        first_bad_defect: bad_def,
        first_bad_positive: bad_pos,
        avg,
        twelve_l: tl,
        verdict,
        reason,
        seconds: start.elapsed().as_secs_f64(),
    })
}

fn print_report(r: &Report, json: bool) {
    let l_str = r.twelve_l.map(|q| q.to_string());
    if json {
        println!(
            "{{\"degree\":{},\"genus\":{},\"states\":{},\"one_piece\":{},\"all_smooth\":{},\"avg_rho_p1\":\"{}\",\"avg_rho_p2\":\"{}\",\"avg_rho_p3\":\"{}\",\"twelve_L\":{},\"L\":{},\"verdict\":\"{}\",\"reason\":\"{}\",\"seconds\":{:.3}}}",
            r.degree,
            r.genus,
            r.states,
            r.one_piece,
            r.all_smooth,
            r.avg[0],
            r.avg[1],
            r.avg[2],
            l_str.clone().map(|s| format!("\"{}\"", s)).unwrap_or("null".into()),
            r.twelve_l.map(|q| format!("\"{}\"", q.div_int(12))).unwrap_or("null".into()),
            r.verdict,
            r.reason,
            r.seconds
        );
        return;
    }
    println!("degree {}, genus {}", r.degree, r.genus);
    println!("orbit: {} states", r.states);
    if r.one_piece {
        println!("one-piece at every state: yes (all-smooth: {})", if r.all_smooth { "yes" } else { "no" });
    } else {
        println!(
            "one-piece at every state: no (a state has defect {} and {} pieces of positive genus)",
            r.first_bad_defect, r.first_bad_positive
        );
    }
    let agree = r.avg[0] == r.avg[1] && r.avg[1] == r.avg[2];
    println!(
        "average of rho over the orbit: p1 {}, p2 {}, p3 {} ({})",
        r.avg[0],
        r.avg[1],
        r.avg[2],
        if agree { "agree" } else { "DO NOT AGREE: report this" }
    );
    if let Some(q) = r.twelve_l {
        println!("12L = {}   L = {}", q, q.div_int(12));
    }
    println!("verdict: {} ({})", r.verdict.to_uppercase(), r.reason);
    println!("time: {:.3} s", r.seconds);
}

// ---------- self-test ----------

fn selftest() -> bool {
    let cases: [(&str, usize, &str, &str); 3] = [
        ("(1 2 4)(3 5 7)(6 8 9);(1 3 6);(1 2 7)(3 8 4)(5 9 6);(1 8 4)(2 9 3)(5 7 6)", 12, "0", "isotrivial"),
        (
            "(1 2 3);(1 4 7 10 13)(2 5 8 11 14)(3 6 9 12 15);(1 5 9 12 15)(2 6 8 11 14)(3 4 7 10 13);(1 12 5 14 8)(2 10 4 15 9)(3 11 6 13 7)",
            540,
            "16/15",
            "not isotrivial",
        ),
        (
            "(1 2 4)(3 5 6);(1 7 13)(2 8 14)(3 9 15)(4 10 16)(5 11 17)(6 12 18);(1 7 18)(2 9 16 6 10 15)(3 12 17)(4 8 13)(5 11 14);(1 9 17 5 10 13)(2 8 15)(3 7 18)(4 11 14)(6 12 16)",
            1296,
            "2/3",
            "not isotrivial",
        ),
    ];
    let mut ok = true;
    for (tuple, states, l12, verdict) in cases {
        let t = parse_tuple(tuple).expect("parse");
        let r = run(t, 10_000_000, true).expect("run");
        let got = r.twelve_l.map(|q| q.to_string()).unwrap_or_default();
        let pass = r.states == states && got == l12 && r.verdict == verdict && r.avg[0] == r.avg[1] && r.avg[1] == r.avg[2];
        println!("{} degree {}: states {} (want {}), 12L {} (want {}), {}", if pass { "PASS" } else { "FAIL" }, r.degree, r.states, states, got, l12, r.verdict);
        ok &= pass;
    }
    ok
}

// ---------- FNV hasher ----------

#[derive(Default)]
struct Fnv(u64);

impl Hasher for Fnv {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        let mut h = if self.0 == 0 { 0xcbf29ce484222325 } else { self.0 };
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        self.0 = h;
    }
}

// ---------- main ----------

const USAGE: &str = "usage:
  isotrivial \"(1 2 4)(3 5 6);(...);(...);(...)\"   four permutations with product 1
  isotrivial \"(...);(...);(...)\"                  three permutations (g4 is computed)
  isotrivial --file FILE                            one tuple per line
options:
  --json             one JSON line per tuple
  --max-states N     give up above N states (default 20000000)
  --selftest         run the three built-in examples";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!("{}", USAGE);
        return;
    }
    if args.iter().any(|a| a == "--selftest") {
        std::process::exit(if selftest() { 0 } else { 1 });
    }
    let json = args.iter().any(|a| a == "--json");
    let mut max_states = 20_000_000usize;
    let mut tuples: Vec<String> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => {}
            "--max-states" => {
                i += 1;
                max_states = args.get(i).and_then(|x| x.parse().ok()).unwrap_or_else(|| {
                    eprintln!("--max-states needs a number");
                    std::process::exit(2)
                });
            }
            "--file" => {
                i += 1;
                let path = args.get(i).unwrap_or_else(|| {
                    eprintln!("--file needs a path");
                    std::process::exit(2)
                });
                let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
                    eprintln!("cannot read {}: {}", path, e);
                    std::process::exit(2)
                });
                tuples.extend(text.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty() && !l.starts_with('#')));
            }
            other => tuples.push(other.to_string()),
        }
        i += 1;
    }
    let mut failed = false;
    for (k, s) in tuples.iter().enumerate() {
        if k > 0 && !json {
            println!();
        }
        match parse_tuple(s).and_then(|t| run(t, max_states, json)) {
            Ok(r) => print_report(&r, json),
            Err(e) => {
                failed = true;
                if json {
                    println!("{{\"error\":\"{}\",\"input\":\"{}\"}}", e, s);
                } else {
                    println!("error: {} (input: {})", e, s);
                }
            }
        }
    }
    if failed {
        std::process::exit(1);
    }
}
