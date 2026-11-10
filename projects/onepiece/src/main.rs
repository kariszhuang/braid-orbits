//! One-piece (geometrically smooth) cusps, clusters and censuses for
//! four-point covers.
//!
//! Conventions (same as ../geosmooth):
//! composition is rightmost-first, (x*y)[i] = x[y[i]], g1 g2 g3 g4 = 1. The six
//! moves are sigma_j: (a,b) -> (a b a^-1, a), m = j-1 in 0..3, and the inverses
//! (a,b) -> (b, b^-1 a b), m = j+2. States are taken modulo simultaneous
//! conjugation in S_n.
//!
//! At the cusp {12|34} of a state, the limit source curve has one piece for
//! each orbit of <g1,g2> and of <g3,g4>, and one node for each cycle of
//! p = g1 g2. A piece on orbit O has genus 2g_O = |O| + 2 - c_O(x) - c_O(y) - c_O(p).
//! ONE-PIECE: the dual graph is a tree (defect 0) with exactly one positive-genus
//! piece; equivalently the stable limit is smooth (Lemma G1 of C22).
//! COMBINATORIALLY SMOOTH (the index equality): ind(g1g2) = ind g1 + ind g2
//! or ind(g3g4) = ind g3 + ind g4. Smooth implies one-piece.
//!
//! Clusters: orbits of <sigma1, sigma3> on the states. Every invariant
//! of the {12|34} cusp is constant on a cluster.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

#[derive(Default)]
struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 { self.0 }
    fn write(&mut self, b: &[u8]) {
        for c in b.chunks(8) {
            let mut v = 0u64;
            for (i, &x) in c.iter().enumerate() { v |= (x as u64) << (8 * i); }
            self.0 = (self.0.rotate_left(5) ^ v).wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }
}
type FxB = BuildHasherDefault<Fx>;

#[cfg(not(feature = "wide"))]
const M: usize = 16;
#[cfg(feature = "wide")]
const M: usize = 24;
type P = [u8; M];
type T = [P; 4];
type Key = [u8; 3 * M];
/// lcm(1..M): node thicknesses 1/e are stored as LL/e.
#[cfg(not(feature = "wide"))]
const LL: i64 = 720720;
#[cfg(feature = "wide")]
const LL: i64 = 5354228880;

#[derive(Clone, Copy)]
struct Cx { n: usize }

impl Cx {
    fn id(&self) -> P { let mut p = [0u8; M]; for i in 0..M { p[i] = i as u8; } p }
    #[inline] fn mul(&self, a: &P, b: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[i] = a[b[i] as usize]; } c }
    #[inline] fn inv(&self, a: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[a[i] as usize] = i as u8; } c }
    fn is_id(&self, a: &P) -> bool { (0..self.n).all(|i| a[i] as usize == i) }
    #[inline] fn ncyc(&self, a: &P) -> usize {
        let mut seen = 0u32; let mut c = 0;
        for s in 0..self.n { if seen >> s & 1 == 1 { continue; } c += 1; let mut x = s;
            while seen >> x & 1 == 0 { seen |= 1 << x; x = a[x] as usize; } }
        c
    }
    #[inline] fn ind(&self, a: &P) -> usize { self.n - self.ncyc(a) }
    fn cyc_lens(&self, a: &P) -> Vec<usize> {
        let mut seen = [false; M]; let mut out = vec![];
        for s in 0..self.n { if seen[s] { continue; } let mut l = 0; let mut x = s;
            while !seen[x] { seen[x] = true; x = a[x] as usize; l += 1; } out.push(l); }
        out
    }
    /// LL * sum over cycles of 1/length
    fn h(&self, a: &P) -> i64 { self.cyc_lens(a).iter().map(|&l| LL / l as i64).sum() }
    /// per-point length of its cycle
    fn len_at(&self, a: &P) -> [u8; M] {
        let mut out = [0u8; M]; let mut seen = [false; M];
        for s in 0..self.n { if seen[s] { continue; } let mut c = vec![]; let mut x = s;
            while !seen[x] { seen[x] = true; c.push(x); x = a[x] as usize; }
            for &y in &c { out[y] = c.len() as u8; } }
        out
    }
    fn orbits(&self, gens: &[&P]) -> ([u8; M], usize) {
        let mut lab = [u8::MAX; M]; let mut k = 0u8; let mut stack = [0u8; M];
        for s in 0..self.n {
            if lab[s] != u8::MAX { continue; }
            let mut sp = 1; stack[0] = s as u8; lab[s] = k;
            while sp > 0 { sp -= 1; let x = stack[sp] as usize;
                for g in gens { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = k; stack[sp] = y as u8; sp += 1; } } }
            k += 1;
        }
        (lab, k as usize)
    }
    fn cyc_per_orbit(&self, a: &P, lab: &[u8; M], no: usize) -> Vec<i64> {
        let mut out = vec![0i64; no]; let mut seen = [false; M];
        for s in 0..self.n { if seen[s] { continue; } out[lab[s] as usize] += 1; let mut x = s;
            while !seen[x] { seen[x] = true; x = a[x] as usize; } }
        out
    }
    fn fmt(&self, p: &P) -> String {
        let mut seen = [false; M]; let mut s = String::new();
        for i in 0..self.n { if seen[i] { continue; } let mut c = vec![]; let mut x = i;
            while !seen[x] { seen[x] = true; c.push(x + 1); x = p[x] as usize; }
            if c.len() > 1 { s.push('('); s.push_str(&c.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ")); s.push(')'); } }
        if s.is_empty() { "()".into() } else { s }
    }
    fn ctype(&self, p: &P) -> String { let mut c = self.cyc_lens(p); c.retain(|&l| l > 1); c.sort_unstable_by(|a, b| b.cmp(a)); format!("{:?}", c) }
}

// ---------------------------------------------------------------- cusps

/// Adds, for each orbit label, the number of cycles of `a` in that orbit.
#[inline]
fn add_cycles(cx: &Cx, a: &P, lab: &[u8; M], out: &mut [i32; M]) {
    let mut seen = 0u32;
    for s in 0..cx.n { if seen >> s & 1 == 1 { continue; } out[lab[s] as usize] += 1; let mut x = s;
        while seen >> x & 1 == 0 { seen |= 1 << x; x = a[x] as usize; } }
}
/// For the pair (x, y) with product p and orbit labels lab (no orbits):
/// (number of positive-genus orbits, total genus). Allocation-free.
#[inline]
fn side(cx: &Cx, x: &P, y: &P, p: &P, lab: &[u8; M], no: usize) -> (u32, i32) {
    let mut c = [0i32; M];
    for i in 0..cx.n { c[lab[i] as usize] += 1; }
    for v in 0..no { c[v] += 2; }
    let mut d = [0i32; M];
    add_cycles(cx, x, lab, &mut d); add_cycles(cx, y, lab, &mut d); add_cycles(cx, p, lab, &mut d);
    let (mut pos, mut tot) = (0u32, 0i32);
    for v in 0..no { let tg = c[v] - d[v]; if tg > 0 { pos += 1; tot += tg / 2; } }
    (pos, tot)
}

/// NEON version of one_piece for M = 16 (aarch64): composition is one `tbl`, orbit and cycle
/// labels (least element) are computed by min-and-jump iterations. Checked against
/// one_piece_scalar by `--selftest`.
#[cfg(all(target_arch = "aarch64", not(feature = "wide")))]
mod nfast {
    use std::arch::aarch64::*;
    use super::{Cx, P, T};
    #[inline(always)] unsafe fn ld(p: &P) -> uint8x16_t { vld1q_u8(p.as_ptr()) }
    #[inline(always)] unsafe fn idv() -> uint8x16_t { let a: [u8; 16] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; vld1q_u8(a.as_ptr()) }
    /// least element of each cycle of q (doubling: 4 rounds cover 16 points)
    #[inline(always)] unsafe fn cyc_min(q: uint8x16_t) -> uint8x16_t {
        let mut m = idv(); let mut r = q;
        m = vminq_u8(m, vqtbl1q_u8(m, r)); r = vqtbl1q_u8(r, r);
        m = vminq_u8(m, vqtbl1q_u8(m, r)); r = vqtbl1q_u8(r, r);
        m = vminq_u8(m, vqtbl1q_u8(m, r)); r = vqtbl1q_u8(r, r);
        vminq_u8(m, vqtbl1q_u8(m, r))
    }
    /// least element of each orbit of <a, b>
    #[inline(always)] unsafe fn orb_min(a: uint8x16_t, b: uint8x16_t) -> uint8x16_t {
        let mut l = idv();
        loop {
            let n1 = vminq_u8(l, vminq_u8(vqtbl1q_u8(l, a), vqtbl1q_u8(l, b)));
            let n2 = vqtbl1q_u8(n1, n1);
            if vmaxvq_u8(veorq_u8(n2, l)) == 0 { return n2; }
            l = n2;
        }
    }
    /// defect 0 at the cusp: cycles(xy) = orbits<x,y> + orbits<z,w> - 1
    #[inline(always)] pub unsafe fn defect0(cx: &Cx, t: &T) -> bool {
        let n = cx.n;
        let (x, y, z, w) = (ld(&t[0]), ld(&t[1]), ld(&t[2]), ld(&t[3]));
        let p = vqtbl1q_u8(x, y); let id = idv(); let ones = vdupq_n_u8(1);
        let cnt = |l: uint8x16_t| -> usize { vaddvq_u8(vandq_u8(vceqq_u8(l, id), ones)) as usize - (16 - n) };
        cnt(cyc_min(p)) + 1 == cnt(orb_min(x, y)) + cnt(orb_min(z, w))
    }
    /// is the cusp of the pair sides (x,y | z,w) one-piece? p = x y.
    #[inline(always)] pub unsafe fn one_piece(cx: &Cx, t: &T) -> bool {
        let n = cx.n;
        let (x, y, z, w) = (ld(&t[0]), ld(&t[1]), ld(&t[2]), ld(&t[3]));
        let p = vqtbl1q_u8(x, y);
        let id = idv();
        let la = orb_min(x, y); let lb = orb_min(z, w); let lp = cyc_min(p);
        let ones = vdupq_n_u8(1);
        let cnt = |l: uint8x16_t| -> usize { vaddvq_u8(vandq_u8(vceqq_u8(l, id), ones)) as usize - (16 - n) };
        let (na, nb, e) = (cnt(la), cnt(lb), cnt(lp));
        if e + 1 != na + nb { return false; }
        // per-orbit 2g - 2 = sum over the orbit of (1 - [leader of x] - [leader of y] - [leader of p])
        let lpm = vandq_u8(vceqq_u8(lp, id), ones);
        let side = |l: uint8x16_t, a: uint8x16_t, b: uint8x16_t| -> u32 {
            let la_ = vandq_u8(vceqq_u8(cyc_min(a), id), ones); let lb_ = vandq_u8(vceqq_u8(cyc_min(b), id), ones);
            let lead = vaddq_u8(vaddq_u8(la_, lb_), lpm);
            let mut lv = [0u8; 16]; let mut cv = [0u8; 16]; vst1q_u8(lv.as_mut_ptr(), l); vst1q_u8(cv.as_mut_ptr(), lead);
            let mut acc = [0i8; 16];
            for i in 0..n { acc[lv[i] as usize] += 1 - cv[i] as i8; }
            let mut pos = 0; for i in 0..n { if lv[i] as usize == i && acc[i] + 2 > 0 { pos += 1; } } pos
        };
        let pa = side(la, x, y); if pa > 1 { return false; }
        pa + side(lb, z, w) <= 1
    }
}

/// --all-zero mode: the per-state test is "defect 0" instead of "one-piece".
static ALLZERO: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
/// --no-s3norm turns the sigma3-orbit normalization of the half-cusp tables off (for checks)
static S3OFF: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
// Set when the census prints its results to the terminal: progress lines are then hidden.
static QUIET: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
macro_rules! progress {
    ($($arg:tt)*) => { if !QUIET.load(Ordering::Relaxed) { eprintln!($($arg)*); } };
}
#[inline] fn az() -> bool { ALLZERO.load(Ordering::Relaxed) }
#[inline] fn defect0_scalar(cx: &Cx, t: &T) -> bool {
    let p = cx.mul(&t[0], &t[1]);
    cx.ncyc(&p) + 1 == cx.orbits(&[&t[0], &t[1]]).1 + cx.orbits(&[&t[2], &t[3]]).1
}
#[inline] fn defect0(cx: &Cx, t: &T) -> bool {
    #[cfg(all(target_arch = "aarch64", not(feature = "wide")))]
    { return unsafe { nfast::defect0(cx, t) }; }
    #[allow(unreachable_code)]
    defect0_scalar(cx, t)
}
/// the per-state property the census looks for
#[inline] fn good(cx: &Cx, t: &T) -> bool { if az() { defect0(cx, t) } else { one_piece(cx, t) } }

#[inline]
fn one_piece(cx: &Cx, t: &T) -> bool {
    #[cfg(all(target_arch = "aarch64", not(feature = "wide")))]
    { return unsafe { nfast::one_piece(cx, t) }; }
    #[allow(unreachable_code)]
    one_piece_scalar(cx, t)
}

/// Fast test: is the {12|34} cusp of t one-piece?
#[inline]
fn one_piece_scalar(cx: &Cx, t: &T) -> bool {
    let p = cx.mul(&t[0], &t[1]);
    let e = cx.ncyc(&p);
    let (la, na) = cx.orbits(&[&t[0], &t[1]]);
    let (lb, nb) = cx.orbits(&[&t[2], &t[3]]);
    if e + 1 != na + nb { return false; }
    let (pa, _) = side(cx, &t[0], &t[1], &p, &la, na);
    if pa > 1 { return false; }
    let (pb, _) = side(cx, &t[2], &t[3], &p, &lb, nb);
    pa + pb <= 1
}

#[derive(Clone, Debug, Default)]
struct Cusp {
    defect: i64,
    /// genera of the pieces on side {12} and side {34}
    gen_a: Vec<i64>, gen_b: Vec<i64>,
    one_piece: bool,
    comb: bool,
    /// LL * (delta_0, delta_1, delta_2, ...) of the stable limit, thickness 1/e per node
    delta: [i64; 9],
    /// edges of the stable dual graph that survive (nodes of the stable limit)
    stable_nodes: usize,
}

fn cusp_full(cx: &Cx, t: &T) -> Cusp {
    let n = cx.n;
    let p = cx.mul(&t[0], &t[1]);
    let ip = cx.ind(&p);
    let comb = ip == cx.ind(&t[0]) + cx.ind(&t[1]) || ip == cx.ind(&t[2]) + cx.ind(&t[3]);
    let (la, na) = cx.orbits(&[&t[0], &t[1]]);
    let (lb, nb) = cx.orbits(&[&t[2], &t[3]]);
    let e = cx.ncyc(&p) as i64;
    let defect = e - (na + nb) as i64 + 1;
    let mut gens = vec![0i64; na + nb];
    let mut size = vec![0i64; na + nb];
    for x in 0..n { size[la[x] as usize] += 1; size[na + lb[x] as usize] += 1; }
    let ca = [cx.cyc_per_orbit(&t[0], &la, na), cx.cyc_per_orbit(&t[1], &la, na), cx.cyc_per_orbit(&p, &la, na)];
    let cb = [cx.cyc_per_orbit(&t[2], &lb, nb), cx.cyc_per_orbit(&t[3], &lb, nb), cx.cyc_per_orbit(&p, &lb, nb)];
    for v in 0..na { let tg = size[v] + 2 - ca[0][v] - ca[1][v] - ca[2][v]; assert!(tg >= 0 && tg % 2 == 0); gens[v] = tg / 2; }
    for v in 0..nb { let tg = size[na + v] + 2 - cb[0][v] - cb[1][v] - cb[2][v]; assert!(tg >= 0 && tg % 2 == 0); gens[na + v] = tg / 2; }
    let npos = gens.iter().filter(|&&g| g > 0).count();
    let one = defect == 0 && npos <= 1;
    // dual graph: one edge per cycle of p
    let mut edges: Vec<(usize, usize, i64, bool)> = vec![]; // (u, v, thickness*LL, alive)
    let mut seen = [false; M];
    for s in 0..n { if seen[s] { continue; } let mut l = 0; let mut x = s;
        while !seen[x] { seen[x] = true; x = p[x] as usize; l += 1; }
        edges.push((la[s] as usize, na + lb[s] as usize, LL / l as i64, true)); }
    let nv = na + nb;
    let mut alive = vec![true; nv];
    // stabilize: drop genus-0 leaves, smooth genus-0 bivalent vertices
    loop {
        let mut deg = vec![0usize; nv];
        for &(u, v, _, a) in &edges { if a { deg[u] += 1; deg[v] += 1; } }
        let mut changed = false;
        for w in 0..nv {
            if !alive[w] || gens[w] > 0 { continue; }
            if deg[w] == 0 && alive.iter().filter(|&&a| a).count() > 1 { alive[w] = false; changed = true; break; }
            if deg[w] == 1 {
                for ed in edges.iter_mut() { if ed.3 && (ed.0 == w || ed.1 == w) { ed.3 = false; } }
                alive[w] = false; changed = true; break;
            }
            if deg[w] == 2 {
                let inc: Vec<usize> = (0..edges.len()).filter(|&i| edges[i].3 && (edges[i].0 == w || edges[i].1 == w)).collect();
                if inc.len() == 2 {
                    let other = |i: usize| if edges[i].0 == w { edges[i].1 } else { edges[i].0 };
                    let (a, b) = (other(inc[0]), other(inc[1]));
                    let th = edges[inc[0]].2 + edges[inc[1]].2;
                    edges[inc[0]].3 = false; edges[inc[1]].3 = false;
                    edges.push((a, b, th, true)); alive[w] = false; changed = true; break;
                }
            }
        }
        if !changed { break; }
    }
    let total_g: i64 = gens.iter().sum::<i64>() + defect;
    let mut delta = [0i64; 9];
    let live: Vec<usize> = (0..edges.len()).filter(|&i| edges[i].3).collect();
    for &i in &live {
        // remove edge i and find the component of edges[i].0
        let (u, v) = (edges[i].0, edges[i].1);
        let mut comp = vec![false; nv]; let mut st = vec![u]; comp[u] = true;
        while let Some(x) = st.pop() {
            for &j in &live { if j == i { continue; } let (a, b) = (edges[j].0, edges[j].1);
                if a == x && !comp[b] { comp[b] = true; st.push(b); }
                if b == x && !comp[a] { comp[a] = true; st.push(a); } }
        }
        if comp[v] { delta[0] += edges[i].2; } else {
            let gv: i64 = (0..nv).filter(|&w| comp[w] && alive[w]).map(|w| gens[w]).sum();
            let ev = live.iter().filter(|&&j| j != i && comp[edges[j].0]).count() as i64;
            let vv = (0..nv).filter(|&w| comp[w] && alive[w]).count() as i64;
            let h = gv + ev - vv + 1;
            let k = h.min(total_g - h) as usize; delta[k.min(8)] += edges[i].2;
        }
    }
    let gen_a = gens[..na].to_vec(); let gen_b = gens[na..].to_vec();
    Cusp { defect, gen_a, gen_b, one_piece: one, comb, delta, stable_nodes: live.len() }
}

// ---------------------------------------------------------------- moves and keys

#[inline]
fn mv(cx: &Cx, t: &T, m: usize) -> T {
    let j = m % 3; let a = t[j]; let b = t[j + 1]; let mut u = *t;
    if m < 3 { u[j] = cx.mul(&a, &cx.mul(&b, &cx.inv(&a))); u[j + 1] = a; }
    else { u[j] = b; u[j + 1] = cx.mul(&cx.inv(&b), &cx.mul(&a, &b)); }
    u
}
fn rot(t: &T) -> T { [t[1], t[2], t[3], t[0]] }

/// Canonical form modulo simultaneous S_n conjugation. Only base points with
/// the lexicographically least (g1,g2,g3)-cycle-length signature are tried.
fn canon(cx: &Cx, t: &T) -> Key {
    let n = cx.n;
    // Base points: the signature class (cycle lengths in g1..g4) with the FEWEST
    // points, ties broken by the signature itself. Conjugation-invariant, so canonical.
    let l0 = cx.len_at(&t[0]); let l1 = cx.len_at(&t[1]); let l2 = cx.len_at(&t[2]); let l3 = cx.len_at(&t[3]);
    let sig = |r: usize| ((l0[r] as u32) << 24) | ((l1[r] as u32) << 16) | ((l2[r] as u32) << 8) | l3[r] as u32;
    let mut sigs = [0u32; M]; for r in 0..n { sigs[r] = sig(r); }
    let mut best_sig = u32::MAX; let mut best_cnt = usize::MAX;
    for r in 0..n { let c = (0..n).filter(|&q| sigs[q] == sigs[r]).count(); if c < best_cnt || (c == best_cnt && sigs[r] < best_sig) { best_cnt = c; best_sig = sigs[r]; } }
    let sig = |r: usize| sigs[r];
    let mut best = [u8::MAX; 3 * M];
    for r in 0..n {
        if sig(r) != best_sig { continue; }
        let mut lab = [u8::MAX; M]; let mut order = [0u8; M]; let mut len = 1; lab[r] = 0; order[0] = r as u8;
        let mut i = 0;
        while i < len { let x = order[i] as usize; i += 1;
            for g in &t[..3] { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y as u8; len += 1; } } }
        debug_assert_eq!(len, n);
        // build key in lex order with early abort
        let mut key = [u8::MAX; 3 * M]; let mut state = 0i8; // 0 equal so far, -1 less, 1 greater
        'outer: for j in 0..3 { for k in 0..n {
            let v = lab[t[j][order[k] as usize] as usize]; key[j * M + k] = v;
            if state == 0 { let b = best[j * M + k]; if v < b { state = -1; } else if v > b { state = 1; break 'outer; } }
        } }
        if state == -1 { best = key; }
    }
    best
}
fn unpack(cx: &Cx, k: &Key) -> T {
    let mut t = [cx.id(); 4];
    for j in 0..3 { for i in 0..cx.n { t[j][i] = k[j * M + i]; } }
    t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2])));
    t
}

// ---------------------------------------------------------------- orbits

/// A state can be a seed of census_fast only if g1 has the minimal index of the four
/// entries and its {12|34} cusp is one-piece. Only such states need to be remembered
/// for de-duplication (memory: the shared table holds a small fraction of all states).
fn seed_like(cx: &Cx, k: &Key) -> bool {
    let t = unpack(cx, k); let i1 = cx.ind(&t[0]);
    if t[1..].iter().any(|q| cx.ind(q) < i1) { return false; }
    good(cx, &t)
}

struct Orbit { states: Vec<Key>, nbr: Vec<[u32; 6]>, all_one: bool, capped: bool, dup: bool }

/// Six-move closure. If `early`, stop at the first state whose {12|34} cusp is not
/// one-piece (all_one = false; states then holds the visited part). The one-piece test
/// runs before canonicalization, so a bad neighbour costs no canonical form. With
/// `sh = Some((table, id))`, visited states are published every 256 new states and the
/// closure aborts (dup = true) on meeting a state owned by a smaller id.
fn close(cx: &Cx, t0: &T, cap: usize, early: bool, sh: Option<(&Shared, u32)>) -> Orbit {
    let root = canon(cx, t0);
    let mut idx: HashMap<Key, u32, FxB> = HashMap::default(); idx.insert(root, 0);
    let mut st = vec![root]; let mut nbr: Vec<[u32; 6]> = vec![]; let mut i = 0; let mut all_one = true;
    let mut flushed = 0usize;
    let ret = |st: Vec<Key>, nbr, all_one, capped, dup| Orbit { states: st, nbr, all_one, capped, dup };
    if early && !good(cx, &unpack(cx, &root)) { return ret(st, nbr, false, false, false); }
    while i < st.len() {
        let u = unpack(cx, &st[i]); let mut row = [0u32; 6];
        for m in 0..6 {
            let v = mv(cx, &u, m);
            if early && !good(cx, &v) { return ret(st, nbr, false, false, false); }
            let k = canon(cx, &v);
            let id = match idx.get(&k) { Some(&j) => j, None => {
                if let Some((tab, me)) = sh { if st.len() % 16 == 0 { if let Some(o) = tab.get(&k) { if o < me { return ret(st, nbr, false, false, true); } } } }
                let j = st.len() as u32; idx.insert(k, j); st.push(k);
                if st.len() > cap { return ret(st, nbr, false, true, false); }
                if let Some((tab, me)) = sh { if st.len() - flushed >= 256 { let keep: Vec<Key> = st[flushed..].iter().filter(|k| seed_like(cx, k)).cloned().collect(); tab.add(&keep, me); flushed = st.len(); } }
                j } };
            row[m] = id;
        }
        nbr.push(row); i += 1;
    }
    if !early { all_one = (0..st.len()).all(|j| good(cx, &unpack(cx, &st[j]))); }
    ret(st, nbr, all_one, false, false)
}

fn components(nbr: &[[u32; 6]], moves: &[usize]) -> (Vec<u32>, usize) {
    let n = nbr.len(); let mut lab = vec![u32::MAX; n]; let mut k = 0u32;
    for s in 0..n { if lab[s] != u32::MAX { continue; } lab[s] = k; let mut stack = vec![s];
        while let Some(x) = stack.pop() { for &m in moves { let y = nbr[x][m] as usize; if lab[y] == u32::MAX { lab[y] = k; stack.push(y); } } }
        k += 1; }
    (lab, k as usize)
}

/// Sizes of the nontrivial block systems of <t>, as a sorted list of block sizes.
fn block_sizes(cx: &Cx, t: &T) -> Vec<usize> {
    let n = cx.n; let mut out = vec![];
    for j in 1..n {
        let mut par: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, x: usize) -> usize { let mut r = x; while p[r] != r { r = p[r]; } let mut y = x; while p[y] != r { let z = p[y]; p[y] = r; y = z; } r }
        let mut queue = vec![(0usize, j)];
        while let Some((x, y)) = queue.pop() {
            let (rx, ry) = (find(&mut par, x), find(&mut par, y)); if rx == ry { continue; }
            par[rx] = ry;
            for g in t.iter() { queue.push((g[x] as usize, g[y] as usize)); }
        }
        let r0 = find(&mut par, 0); let sz = (0..n).filter(|&x| find(&mut par, x) == r0).count();
        if sz < n { out.push(sz); }
    }
    out.sort(); out.dedup(); out
}

fn genus_of(cx: &Cx, t: &T) -> i64 { let s: usize = t.iter().map(|p| cx.ind(p)).sum(); (s as i64 - 2 * cx.n as i64 + 2) / 2 }

/// Full report of a closed orbit as one JSON object.
fn report(cx: &Cx, o: &Orbit) -> String { report2(cx, o).0 }
fn report2(cx: &Cx, o: &Orbit) -> (String, usize) {
    let n = o.states.len();
    let ts: Vec<T> = o.states.iter().map(|k| unpack(cx, k)).collect();
    let t0 = ts[0]; let g = genus_of(cx, &t0);
    // clusters: <sigma1, sigma3>. The {12|34} cusp is constant on a cluster and
    // the {23|41} cusp on a <sigma2>-orbit, so each is computed once per class.
    let (clab, nclu) = components(&o.nbr, &[0, 2, 3, 5]);
    let mut per_clu: Vec<Option<Cusp>> = vec![None; nclu];
    for s in 0..n { let c = clab[s] as usize; if per_clu[c].is_none() { per_clu[c] = Some(cusp_full(cx, &ts[s])); } }
    let cusps: Vec<Cusp> = (0..n).map(|s| per_clu[clab[s] as usize].clone().unwrap()).collect();
    let (s2lab, ns2) = components(&o.nbr, &[1, 4]);
    let mut per_s2: Vec<Option<[i64; 9]>> = vec![None; ns2];
    for s in 0..n { let c = s2lab[s] as usize; if per_s2[c].is_none() { per_s2[c] = Some(cusp_full(cx, &rot(&ts[s])).delta); } }
    let (s1lab, ns1) = components(&o.nbr, &[0, 3]);
    let mut clu_bad = vec![false; nclu]; let mut clu_defect = vec![false; nclu]; let mut clu_sep = vec![false; nclu];
    let mut clu_size = vec![0usize; nclu]; let mut clu_comb = vec![false; nclu];
    for s in 0..n { let c = clab[s] as usize; clu_size[c] += 1;
        if !cusps[s].one_piece { clu_bad[c] = true; if cusps[s].defect > 0 { clu_defect[c] = true; } else { clu_sep[c] = true; } }
        if cusps[s].comb { clu_comb[c] = true; } }
    // consistency: one-piece status constant on clusters
    for s in 0..n { assert_eq!(one_piece(cx, &ts[s]), !clu_bad[clab[s] as usize], "one-piece status not constant on a cluster"); }
    let mut s1_bad = vec![false; ns1]; for s in 0..n { if !cusps[s].one_piece { s1_bad[s1lab[s] as usize] = true; } }
    let bad_states = cusps.iter().filter(|c| !c.one_piece).count();
    let comb_states = cusps.iter().filter(|c| c.comb).count();
    let defect_states = cusps.iter().filter(|c| c.defect > 0).count();
    // lambda (C18 formula, as in geosmooth): 12 * LL * deg lambda on the N-state curve
    let base: i64 = LL * cx.n as i64 - t0.iter().map(|p| cx.h(p)).sum::<i64>();
    let mut lam = base * n as i64; let mut delta = [0i64; 9];
    for (s, u) in ts.iter().enumerate() {
        let ab = cx.mul(&u[0], &u[1]); let bc = cx.mul(&u[1], &u[2]);
        let third = cx.mul(&u[0], &cx.mul(&u[1], &cx.mul(&u[2], &cx.inv(&u[1]))));
        lam += cx.h(&ab) + cx.h(&bc) + cx.h(&third);
        // the same three cusps: {12|34} at s, {23|41} at s (= {12|34} of the rotation), {12|34} of sigma2(s)
        let c2 = per_s2[s2lab[s] as usize].unwrap(); let c3 = &cusps[o.nbr[s][1] as usize];
        for k in 0..9 { delta[k] += cusps[s].delta[k] + c2[k] + c3.delta[k]; }
    }
    let blocks = block_sizes(cx, &t0);
    let nbad = clu_bad.iter().filter(|&&b| b).count();
    let ntwo = (0..nclu).filter(|&c| clu_bad[c] && !clu_defect[c]).count();
    let ndef = (0..nclu).filter(|&c| clu_defect[c]).count();
    let mut bad_sizes: Vec<usize> = (0..nclu).filter(|&c| clu_bad[c]).map(|c| clu_size[c]).collect(); bad_sizes.sort();
    let types: Vec<String> = t0.iter().map(|p| cx.ctype(p)).collect();
    let tuple: Vec<String> = t0.iter().map(|p| cx.fmt(p)).collect();
    // Mumford check in genus 2: 10 lambda = delta0 + 2 delta1  <=>  10*lam = 12*(d0 + 2 d1)
    let mumford = if g == 2 { format!("{}", 10 * lam == 12 * (delta[0] + 2 * delta[1])) } else { "null".into() };
    (format!("{{\"d\":{},\"genus\":{},\"states\":{},\"clusters\":{},\"bad_clusters\":{},\"bad_defect_clusters\":{},\"bad_separating_clusters\":{},\"bad_cluster_sizes\":{:?},\"sigma1_orbits\":{},\"bad_sigma1_orbits\":{},\"bad_states\":{},\"defect_states\":{},\"comb_smooth_states\":{},\"all_one_piece\":{},\"all_zero\":{},\"all_smooth\":{},\"primitive\":{},\"block_sizes\":{:?},\"lambda12_LL\":{},\"delta_LL\":{:?},\"mumford_g2\":{},\"LL\":{},\"types\":{:?},\"tuple\":{:?}}}",
        cx.n, g, n, nclu, nbad, ndef, ntwo, bad_sizes, ns1, s1_bad.iter().filter(|&&b| b).count(), bad_states, defect_states, comb_states,
        bad_states == 0, defect_states == 0, comb_states == n, blocks.is_empty(), blocks, lam, delta, mumford, LL, types, tuple), nbad)
}

// ---------------------------------------------------------------- parsing

fn parse_perm(n: usize, s: &str) -> P {
    let mut p = [0u8; M]; for i in 0..M { p[i] = i as u8; }
    for cyc in s.split(')') { let cyc = cyc.trim().trim_start_matches('('); if cyc.is_empty() { continue; }
        let v: Vec<usize> = cyc.split(|c: char| c == ' ' || c == ',').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap() - 1).collect();
        for i in 0..v.len() { p[v[i]] = v[(i + 1) % v.len()] as u8; } }
    let _ = n; p
}
fn parse_tuple(s: &str) -> (Cx, T, bool) {
    let parts: Vec<&str> = s.split(|c| c == ';' || c == '|').collect(); assert_eq!(parts.len(), 4, "need four permutations separated by ;");
    let n = s.split(|c: char| !c.is_ascii_digit()).filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap()).max().unwrap();
    if n > M {
        eprintln!("error: this tuple moves {} points, but this build handles at most {}. Build with `cargo build --release --features wide` for up to 24 points, or use projects/geosmooth/geo_check.py, which has no limit.", n, M);
        std::process::exit(2);
    }
    let cx = Cx { n };
    let mut t = [cx.id(); 4]; for j in 0..4 { t[j] = parse_perm(n, parts[j]); }
    let prod = cx.mul(&t[0], &cx.mul(&t[1], &cx.mul(&t[2], &t[3])));
    if cx.is_id(&prod) { return (cx, t, false); }
    // left-to-right convention: reverse the tuple
    let r = [t[3], t[2], t[1], t[0]];
    let prod2 = cx.mul(&r[0], &cx.mul(&r[1], &cx.mul(&r[2], &r[3])));
    assert!(cx.is_id(&prod2), "product is not 1 in either convention");
    (cx, r, true)
}

// ---------------------------------------------------------------- census

fn all_perms(n: usize) -> Vec<P> {
    let mut out = vec![]; let mut a: Vec<u8> = (0..n as u8).collect();
    fn rec(k: usize, a: &mut Vec<u8>, out: &mut Vec<P>) {
        if k == a.len() { let mut p = [0u8; M]; for i in 0..M { p[i] = i as u8; } for i in 0..a.len() { p[i] = a[i]; } out.push(p); return; }
        for i in k..a.len() { a.swap(k, i); rec(k + 1, a, out); a.swap(k, i); }
    }
    rec(0, &mut a, &mut out); out
}
fn partitions(n: usize, maxp: usize) -> Vec<Vec<usize>> {
    if n == 0 { return vec![vec![]]; }
    let mut out = vec![];
    for p in (1..=maxp.min(n)).rev() { for mut rest in partitions(n - p, p) { rest.insert(0, p); out.push(rest); } }
    out
}
fn class_rep(n: usize, parts: &[usize]) -> P {
    let mut p = [0u8; M]; for i in 0..M { p[i] = i as u8; }
    let mut s = 0; for &l in parts { for i in 0..l { p[s + i] = (s + (i + 1) % l) as u8; } s += l; } assert_eq!(s, n); p
}
fn rank(n: usize, p: &P) -> usize {
    // Lehmer code
    let mut r = 0usize; let mut used = 0u32;
    for i in 0..n { let v = p[i] as u32; let smaller = (used ^ ((1u32 << n) - 1)) & ((1u32 << v) - 1); r = r * (n - i) + smaller.count_ones() as usize; used |= 1 << v; }
    r
}

/// Shared table of visited states, each with the id of the closure that owns it.
/// Concurrent closures of the same orbit: the one with the larger id aborts on meeting
/// a state owned by a smaller id; the smallest id never aborts, so every orbit is
/// completed by exactly one closure.
struct Shared { known: Vec<Mutex<HashMap<Key, u32, FxB>>> }
impl Shared {
    fn new() -> Self { Shared { known: (0..64).map(|_| Mutex::new(HashMap::default())).collect() } }
    fn shard(k: &Key) -> usize { let mut h = Fx::default(); h.write(k); (h.finish() >> 58) as usize }
    fn get(&self, k: &Key) -> Option<u32> { self.known[Self::shard(k)].lock().unwrap().get(k).copied() }
    fn contains(&self, k: &Key) -> bool { self.get(k).is_some() }
    fn add(&self, ks: &[Key], owner: u32) {
        for k in ks { let mut m = self.known[Self::shard(k)].lock().unwrap(); let e = m.entry(*k).or_insert(owner); if owner < *e { *e = owner; } }
    }
}

fn census(n: usize, genus: usize, all_orbits: bool, need_good: bool, max_bad: usize, threads: usize, cap: usize, out: &str) {
    let cx = Cx { n };
    let t_start = std::time::Instant::now();
    let perms = all_perms(n);
    let pind: Vec<u8> = perms.iter().map(|p| cx.ind(p) as u8).collect();
    let target = 2 * n - 2 + 2 * genus;
    // work items: (g1 class rep of minimal index, g2 representative under conjugation by C(g1))
    let mut items: Vec<(P, P)> = vec![];
    for parts in partitions(n, n) {
        let g1 = class_rep(n, &parts); let i1 = cx.ind(&g1);
        if i1 == 0 || 4 * i1 > target { continue; }
        let cent: Vec<P> = perms.iter().filter(|c| cx.mul(c, &g1) == cx.mul(&g1, c)).cloned().collect();
        let mut seen = vec![false; perms.len()];
        let mut cnt = 0;
        for (k, g2) in perms.iter().enumerate() {
            if seen[k] { continue; }
            for c in &cent { let h = cx.mul(c, &cx.mul(g2, &cx.inv(c))); seen[rank(n, &h)] = true; }
            let i2 = pind[k] as usize;
            if i2 < i1 || i1 + i2 + 2 * i1 > target || i2 + 3 * i1 > target { continue; }
            if !all_orbits {
                // <g1,g2> may carry at most one positive-genus piece
                let p = cx.mul(&g1, g2); let (lab, no) = cx.orbits(&[&g1, g2]);
                let mut size = [0i64; M]; for i in 0..n { size[lab[i] as usize] += 1; }
                let a = cx.cyc_per_orbit(&g1, &lab, no); let b = cx.cyc_per_orbit(g2, &lab, no); let c = cx.cyc_per_orbit(&p, &lab, no);
                if (0..no).filter(|&v| size[v] + 2 - a[v] - b[v] - c[v] > 0).count() > 1 { continue; }
            }
            items.push((g1, *g2)); cnt += 1;
        }
        progress!("class {:?}: |C|={} g2 reps kept {}", parts, cent.len(), cnt);
    }
    progress!("items {} setup {:.1}s", items.len(), t_start.elapsed().as_secs_f64());
    let shared = Shared::new();
    let next = AtomicUsize::new(0);
    let seeds = AtomicU64::new(0); let passed = AtomicU64::new(0); let closed = AtomicU64::new(0); let capped = AtomicU64::new(0);
    let found = AtomicU64::new(0);
    let file = Mutex::new(std::fs::File::create(out).unwrap());
    let ids: Mutex<HashSet<Key, FxB>> = Mutex::new(HashSet::default());
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| {
                loop {
                    let it = next.fetch_add(1, Ordering::Relaxed); if it >= items.len() { break; }
                    if it % 2000 == 0 { progress!("item {}/{} {:.0}s seeds {} passed {} closed {} found {} capped {}", it, items.len(), t_start.elapsed().as_secs_f64(),
                        seeds.load(Ordering::Relaxed), passed.load(Ordering::Relaxed), closed.load(Ordering::Relaxed), found.load(Ordering::Relaxed), capped.load(Ordering::Relaxed)); }
                    let (g1, g2) = items[it]; let i1 = cx.ind(&g1); let i2 = cx.ind(&g2);
                    let p12 = cx.mul(&g1, &g2);
                    let (la, na) = cx.orbits(&[&g1, &g2]);
                    let e = cx.ncyc(&p12);
                    // p-cycle id per point, for the tree test via a DSU on cycles
                    let mut cid = [0u8; M]; { let mut seen = [false; M]; let mut k = 0u8;
                        for s in 0..n { if seen[s] { continue; } let mut x = s; while !seen[x] { seen[x] = true; cid[x] = k; x = p12[x] as usize; } k += 1; } }
                    let need = (e + 1) as isize - na as isize;
                    if (!all_orbits || need_good) && need < 1 { continue; }
                    let _ = la;
                    let (mut s_seeds, mut s_pass) = (0u64, 0u64);
                    for (k3, g3) in perms.iter().enumerate() {
                        let i3 = pind[k3] as usize; if i3 < i1 || i1 + i2 + i3 + i1 > target { continue; }
                        if !all_orbits || need_good {
                            let mut par = [0u8; M]; for i in 0..e { par[i] = i as u8; }
                            let mut comps = e;
                            for i in 0..n { let (mut a, mut b) = (cid[i] as usize, cid[g3[i] as usize] as usize);
                                while par[a] as usize != a { a = par[a] as usize; } while par[b] as usize != b { b = par[b] as usize; }
                                if a != b { par[a] = b as u8; comps -= 1; } }
                            if comps as isize != need { continue; }
                        }
                        let q = cx.mul(&p12, g3); let i4 = cx.ind(&q);
                        if i1 + i2 + i3 + i4 != target || i4 < i1 { continue; }
                        let t = [g1, g2, *g3, cx.inv(&q)];
                        if cx.orbits(&[&t[0], &t[1], &t[2]]).1 != 1 { continue; }
                        s_seeds += 1;
                        if !all_orbits {
                            if !one_piece(&cx, &t) || !one_piece(&cx, &rot(&t)) || !one_piece(&cx, &mv(&cx, &t, 1)) || !one_piece(&cx, &mv(&cx, &t, 4)) { continue; }
                        } else if need_good && !one_piece(&cx, &t) { continue; }
                        s_pass += 1;
                        let key = canon(&cx, &t);
                        if shared.contains(&key) { continue; }
                        let o = close(&cx, &t, cap, !all_orbits, None);
                        closed.fetch_add(1, Ordering::Relaxed);
                        shared.add(&o.states, 0);
                        if o.capped { capped.fetch_add(1, Ordering::Relaxed); let line = format!("{{\"capped\":true,\"d\":{},\"tuple\":{:?}}}\n", n, t.iter().map(|p| cx.fmt(p)).collect::<Vec<_>>()); file.lock().unwrap().write_all(line.as_bytes()).unwrap(); continue; }
                        if all_orbits || o.all_one {
                            let id = *o.states.iter().min().unwrap();
                            if !ids.lock().unwrap().insert(id) { continue; }
                            let o2 = if o.nbr.len() == o.states.len() { o } else { close(&cx, &t, cap, false, None) };
                            let (line, nbad) = report2(&cx, &o2);
                            if nbad > max_bad { continue; }
                            found.fetch_add(1, Ordering::Relaxed);
                            let line = line + "\n"; file.lock().unwrap().write_all(line.as_bytes()).unwrap();
                        }
                    }
                    seeds.fetch_add(s_seeds, Ordering::Relaxed); passed.fetch_add(s_pass, Ordering::Relaxed);
                }
            });
        }
    });
    let summary = format!("{{\"summary\":true,\"d\":{},\"genus\":{},\"mode\":\"{}\",\"need_good\":{},\"max_bad\":{},\"items\":{},\"seeds\":{},\"passed_prefilter\":{},\"closures\":{},\"orbits_reported\":{},\"capped\":{},\"seconds\":{:.1},\"complete\":true}}\n",
        n, genus, if all_orbits { "all_orbits" } else { "all_one_piece" }, need_good, max_bad, items.len(), seeds.load(Ordering::Relaxed), passed.load(Ordering::Relaxed),
        closed.load(Ordering::Relaxed), found.load(Ordering::Relaxed), capped.load(Ordering::Relaxed), t_start.elapsed().as_secs_f64());
    file.lock().unwrap().write_all(summary.as_bytes()).unwrap();
    if !QUIET.load(Ordering::Relaxed) { eprint!("{}", summary); }
    if std::env::var_os("ONEPIECE_PROFILE").is_some() { eprintln!("PROFILE pass-pattern (rot,s2,s2inv) counts: {:?}", PROF.iter().map(|a| a.load(Ordering::Relaxed)).collect::<Vec<_>>()); }
}

/// v3 census (one-piece and need-good modes). Same population and normal form as
/// `census`: g1 a class representative of minimal index, g2 modulo C(g1), every g3.
/// New: the {12|34} cusp condition on the (g3,g4) half depends only on p = g1 g2.
/// For each cycle type lambda of p we precompute, once, every g3 such that, with
/// p_lambda the class representative and g4 = (p_lambda g3)^-1, the side <g3,g4> has
/// at most one positive-genus orbit and total genus 0 or G. Buckets are keyed by
/// (lambda, number of <g3,g4> orbits, genus>0) and sorted by min(ind g3, ind g4).
/// Each item is conjugated so that p = p_lambda; the tree condition fixes the bucket,
/// and the min-index condition is a prefix of it. Items whose <g1,g2> side has genus
/// strictly between 0 and G are impossible (two positive pieces) and are skipped.
static PROF: [AtomicU64; 8] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
fn census_fast(n: usize, genus: usize, need_good: bool, max_bad: usize, threads: usize, cap: usize, out: &str, start_item: usize) {
    let cx = Cx { n };
    let t_start = std::time::Instant::now();
    let perms = all_perms(n);
    let pind: Vec<u8> = perms.iter().map(|p| cx.ind(p) as u8).collect();
    let target = 2 * n - 2 + 2 * genus;
    let g = genus as i32;
    let parts_all = partitions(n, n);
    let lam_index: HashMap<Vec<usize>, usize> = parts_all.iter().enumerate().map(|(i, q)| (q.clone(), i)).collect();
    let reps: Vec<P> = parts_all.iter().map(|q| class_rep(n, q)).collect();
    // ---- half-cusp tables
    let s3norm = !S3OFF.load(Ordering::Relaxed);
    let idx_of: Vec<u32> = if s3norm { let mut v = vec![0u32; perms.len()]; for (k, q) in perms.iter().enumerate() { v[rank(n, q)] = k as u32; } v } else { vec![] };
    type Bucket = Vec<(u32, u8)>; // (rank of g3, min(ind g3, ind g4))
    let next_l = AtomicUsize::new(0);
    let tables: Vec<Mutex<HashMap<(u8, u8), Bucket>>> = (0..parts_all.len()).map(|_| Mutex::new(HashMap::new())).collect();
    let entries = AtomicU64::new(0);
    std::thread::scope(|sc| {
        for _ in 0..threads { sc.spawn(|| loop {
            let li = next_l.fetch_add(1, Ordering::Relaxed); if li >= parts_all.len() { break; }
            let p = reps[li]; let pinv = cx.inv(&p);
            let mut map: HashMap<(u8, u8), Bucket> = HashMap::new();
            // sigma3-normalization: sigma3 (g3,g4) -> (g3 g4 g3^-1, g3) keeps g1, g2, p, <g3,g4> and the
            // cluster, so only one g3 per sigma3-orbit is needed (the one with the least index).
            let mut vis = if s3norm { vec![false; perms.len()] } else { vec![] };
            for (k3, g3) in perms.iter().enumerate() {
                if s3norm && vis[k3] { continue; }
                let i3 = pind[k3]; if i3 == 0 { continue; }
                let g4 = cx.mul(&cx.inv(g3), &pinv); // g3 g4 = p^-1
                let i4 = cx.ind(&g4) as u8; if i4 == 0 { continue; }
                let (lab, no) = cx.orbits(&[g3, &p]);
                let (pos, tot) = side(&cx, g3, &g4, &p, &lab, no);
                if (!az() && (pos > 1 || (tot != 0 && tot != g))) || tot > g { continue; }
                let mut kmin = k3;
                if s3norm {
                    let (mut a, mut b) = (*g3, g4);
                    loop { let na = cx.mul(&a, &cx.mul(&b, &cx.inv(&a))); b = a; a = na; if a == *g3 { break; }
                        let k = idx_of[rank(n, &a)] as usize; vis[k] = true; if k < kmin { kmin = k; } }
                    vis[k3] = true;
                }
                map.entry((no as u8, tot as u8)).or_default().push((kmin as u32, i3.min(i4)));
            }
            let mut cnt = 0u64;
            for b in map.values_mut() { b.sort_unstable_by(|a, b| b.1.cmp(&a.1)); cnt += b.len() as u64; }
            entries.fetch_add(cnt, Ordering::Relaxed);
            *tables[li].lock().unwrap() = map;
        }); }
    });
    let tables: Vec<HashMap<(u8, u8), Bucket>> = tables.into_iter().map(|m| m.into_inner().unwrap()).collect();
    progress!("half-cusp tables: {} entries, {:.1}s", entries.load(Ordering::Relaxed), t_start.elapsed().as_secs_f64());
    // ---- items: (g1 class rep of minimal index, g2 modulo C(g1))
    let mut items: Vec<(P, P)> = vec![];
    for parts in &parts_all {
        let g1 = class_rep(n, parts); let i1 = cx.ind(&g1);
        if i1 == 0 || 4 * i1 > target { continue; }
        let cent: Vec<P> = perms.iter().filter(|c| cx.mul(c, &g1) == cx.mul(&g1, c)).cloned().collect();
        let mut seen = vec![false; perms.len()]; let mut cnt = 0;
        for (k, g2) in perms.iter().enumerate() {
            if seen[k] { continue; }
            for c in &cent { let h = cx.mul(c, &cx.mul(g2, &cx.inv(c))); seen[rank(n, &h)] = true; }
            let i2 = pind[k] as usize;
            if i2 < i1 || i2 + 3 * i1 > target { continue; }
            let p = cx.mul(&g1, g2); let (lab, no) = cx.orbits(&[&g1, g2]);
            let (pos, tot) = side(&cx, &g1, g2, &p, &lab, no);
            if (!az() && (pos > 1 || (tot != 0 && tot != g))) || tot > g { continue; }
            items.push((g1, *g2)); cnt += 1;
        }
        progress!("class {:?}: |C|={} g2 reps kept {}", parts, cent.len(), cnt);
    }
    progress!("items {} setup {:.1}s", items.len(), t_start.elapsed().as_secs_f64());
    let shared = Shared::new();
    // --start-item S resumes a run: the item list is deterministic, items < S were completed earlier
    let profile = std::env::var_os("ONEPIECE_PROFILE").is_some();
    let next = AtomicUsize::new(start_item);
    if start_item > 0 { eprintln!("resuming at item {} of {}", start_item, items.len()); }
    let scanned = AtomicU64::new(0); let seeds = AtomicU64::new(0); let passed = AtomicU64::new(0);
    let closed = AtomicU64::new(0); let capped = AtomicU64::new(0); let found = AtomicU64::new(0); let probe_rej = AtomicU64::new(0);
    let next_id = std::sync::atomic::AtomicU32::new(1); let dups = AtomicU64::new(0); let stored = AtomicU64::new(0);
    let ns_close = AtomicU64::new(0); let ns_report = AtomicU64::new(0); let ns_canon = AtomicU64::new(0); let closed_states = AtomicU64::new(0);
    let file = Mutex::new(std::fs::File::create(out).unwrap());
    let ids: Mutex<HashSet<Key, FxB>> = Mutex::new(HashSet::default());
    std::thread::scope(|sc| {
        for _ in 0..threads { sc.spawn(|| loop {
            let it = next.fetch_add(1, Ordering::Relaxed); if it >= items.len() { break; }
            if it % 20000 == 0 { progress!("item {}/{} {:.0}s scanned {} seeds {} passed {} closed {} found {} capped {}", it, items.len(), t_start.elapsed().as_secs_f64(),
                scanned.load(Ordering::Relaxed), seeds.load(Ordering::Relaxed), passed.load(Ordering::Relaxed), closed.load(Ordering::Relaxed), found.load(Ordering::Relaxed), capped.load(Ordering::Relaxed)); }
            let (g1, g2) = items[it]; let i1 = cx.ind(&g1) as u8;
            let p = cx.mul(&g1, &g2);
            // conjugate so that p becomes its class representative
            let mut cyc: Vec<Vec<u8>> = vec![]; { let mut seen = 0u32;
                for s0 in 0..n { if seen >> s0 & 1 == 1 { continue; } let mut c = vec![]; let mut x = s0;
                    while seen >> x & 1 == 0 { seen |= 1 << x; c.push(x as u8); x = p[x] as usize; } cyc.push(c); } }
            cyc.sort_by(|a, b| b.len().cmp(&a.len()));
            let lam: Vec<usize> = cyc.iter().map(|c| c.len()).collect();
            let li = lam_index[&lam];
            let mut sig = [0u8; M]; for i in 0..M { sig[i] = i as u8; }
            { let mut base = 0u8; for c in &cyc { for (i, &x) in c.iter().enumerate() { sig[x as usize] = base + i as u8; } base += c.len() as u8; } }
            let sinv = cx.inv(&sig);
            let h1 = cx.mul(&sig, &cx.mul(&g1, &sinv)); let h2 = cx.mul(&sig, &cx.mul(&g2, &sinv));
            let rep = reps[li]; debug_assert!(cx.mul(&h1, &h2) == rep);
            let (la, na) = cx.orbits(&[&h1, &h2]);
            let (_, tot12) = side(&cx, &h1, &h2, &rep, &la, na);
            let e = lam.len();
            let need = e + 1 - na; // e >= na always (each orbit holds a p-cycle)
            if tot12 > g { continue; } let Some(bucket) = tables[li].get(&(need as u8, (g - tot12) as u8)) else { continue; };
            let repinv = cx.inv(&rep);
            let (mut s_scan, mut s_seeds, mut s_pass) = (0u64, 0u64, 0u64);
            for &(k3, m) in bucket.iter() {
                if m < i1 { break; }
                s_scan += 1;
                let g3 = &perms[k3 as usize];
                // transitivity: <g1,g2>-orbits joined by g3
                let mut par = [0u8; M]; for i in 0..na { par[i] = i as u8; } let mut comps = na;
                for x in 0..n { let (mut a, mut b) = (la[x] as usize, la[g3[x] as usize] as usize);
                    while par[a] as usize != a { a = par[a] as usize; } while par[b] as usize != b { b = par[b] as usize; }
                    if a != b { par[a] = b as u8; comps -= 1; if comps == 1 { break; } } }
                if comps != 1 { continue; }
                s_seeds += 1;
                let t = [h1, h2, *g3, cx.mul(&cx.inv(g3), &repinv)];
                debug_assert!(good(&cx, &t));
                if !need_good {
                    if profile {
                        let a = good(&cx, &rot(&t)); let b = good(&cx, &mv(&cx, &t, 1)); let c = good(&cx, &mv(&cx, &t, 4));
                        PROF[(a as usize) | (b as usize) << 1 | (c as usize) << 2].fetch_add(1, Ordering::Relaxed);
                        if !(a && b && c) { continue; }
                    } else if !good(&cx, &rot(&t)) || !good(&cx, &mv(&cx, &t, 1)) || !good(&cx, &mv(&cx, &t, 4)) { continue; }
                }
                s_pass += 1;
                // probe BEFORE the canonical form (not in need-good mode, which closes every orbit):
                // a random walk of 64 moves; a state without the property rejects the orbit.
                // Sound: it never rejects a good orbit. Seeded from the tuple bytes (deterministic).
                if !need_good {
                    let mut r: u64 = { let mut h = Fx::default(); for q in &t[..3] { h.write(&q[..n]); } h.finish() | 1 }; let mut u = t; let mut bad = false;
                    for _ in 0..64 { r ^= r << 13; r ^= r >> 7; r ^= r << 17; u = mv(&cx, &u, (r % 6) as usize); if !good(&cx, &u) { bad = true; break; } }
                    if bad { probe_rej.fetch_add(1, Ordering::Relaxed); continue; }
                }
                let t0 = std::time::Instant::now();
                let key = canon(&cx, &t);
                let hit = shared.contains(&key);
                ns_canon.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
                if hit { continue; }
                let t0 = std::time::Instant::now();
                let me = next_id.fetch_add(1, Ordering::Relaxed);
                let o = close(&cx, &t, cap, !need_good, Some((&shared, me)));
                ns_close.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
                closed_states.fetch_add(o.states.len() as u64, Ordering::Relaxed);
                closed.fetch_add(1, Ordering::Relaxed);
                { let keep: Vec<Key> = o.states.iter().filter(|k| seed_like(&cx, k)).cloned().collect(); shared.add(&keep, me); stored.fetch_add(keep.len() as u64, Ordering::Relaxed); }
                if o.dup { dups.fetch_add(1, Ordering::Relaxed); continue; }
                if o.capped { capped.fetch_add(1, Ordering::Relaxed); let line = format!("{{\"capped\":true,\"d\":{},\"tuple\":{:?}}}\n", n, t.iter().map(|q| cx.fmt(q)).collect::<Vec<_>>()); file.lock().unwrap().write_all(line.as_bytes()).unwrap(); continue; }
                if need_good || o.all_one {
                    let id = *o.states.iter().min().unwrap();
                    if !ids.lock().unwrap().insert(id) { continue; }
                    let t0 = std::time::Instant::now();
                    let (line, nbad) = report2(&cx, &o);
                    ns_report.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    if nbad > max_bad { continue; }
                    found.fetch_add(1, Ordering::Relaxed);
                    let line = line + "\n"; file.lock().unwrap().write_all(line.as_bytes()).unwrap();
                }
            }
            scanned.fetch_add(s_scan, Ordering::Relaxed); seeds.fetch_add(s_seeds, Ordering::Relaxed); passed.fetch_add(s_pass, Ordering::Relaxed);
        }); }
    });
    let summary = format!("{{\"summary\":true,\"version\":\"v9\",\"d\":{},\"genus\":{},\"mode\":\"{}\",\"need_good\":{},\"max_bad\":{},\"table_entries\":{},\"items\":{},\"start_item\":{},\"scanned\":{},\"seeds\":{},\"passed_prefilter\":{},\"probe_rejects\":{},\"closures\":{},\"orbits_reported\":{},\"capped\":{},\"dup_aborts\":{},\"stored_states\":{},\"closed_states\":{},\"thread_s_canon\":{:.1},\"thread_s_close\":{:.1},\"thread_s_report\":{:.1},\"seconds\":{:.1},\"complete\":true}}\n",
        n, genus, if need_good { "all_orbits" } else if az() { "all_zero" } else { "all_one_piece" }, need_good, max_bad, entries.load(Ordering::Relaxed), items.len(), start_item, scanned.load(Ordering::Relaxed), seeds.load(Ordering::Relaxed),
        passed.load(Ordering::Relaxed), probe_rej.load(Ordering::Relaxed), closed.load(Ordering::Relaxed), found.load(Ordering::Relaxed), capped.load(Ordering::Relaxed),
        dups.load(Ordering::Relaxed), stored.load(Ordering::Relaxed), closed_states.load(Ordering::Relaxed), ns_canon.load(Ordering::Relaxed) as f64 * 1e-9, ns_close.load(Ordering::Relaxed) as f64 * 1e-9, ns_report.load(Ordering::Relaxed) as f64 * 1e-9, t_start.elapsed().as_secs_f64());
    file.lock().unwrap().write_all(summary.as_bytes()).unwrap();
    if !QUIET.load(Ordering::Relaxed) { eprint!("{}", summary); }
    if std::env::var_os("ONEPIECE_PROFILE").is_some() { eprintln!("PROFILE pass-pattern (rot,s2,s2inv) counts: {:?}", PROF.iter().map(|a| a.load(Ordering::Relaxed)).collect::<Vec<_>>()); }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let has = |k: &str| args.iter().any(|a| a == k);
    if let Some(s) = get("--orbit") {
        let (cx, t, rev) = parse_tuple(&s);
        if rev { eprintln!("note: product was 1 in left-to-right order; using the reversed tuple"); }
        let cap: usize = get("--cap").map(|x| x.parse().unwrap()).unwrap_or(50_000_000);
        let o = close(&cx, &t, cap, false, None);
        if o.capped { println!("{{\"capped\":true,\"states_seen\":{}}}", o.states.len()); return; }
        println!("{}", report(&cx, &o));
        if has("--clusters") {
            let ts: Vec<T> = o.states.iter().map(|k| unpack(&cx, k)).collect();
            let (clab, nclu) = components(&o.nbr, &[0, 2, 3, 5]);
            for c in 0..nclu {
                let s = (0..ts.len()).find(|&s| clab[s] as usize == c).unwrap();
                let cu = cusp_full(&cx, &ts[s]); let size = clab.iter().filter(|&&x| x as usize == c).count();
                println!("cluster {} size {} one_piece {} comb {} defect {} genera12 {:?} genera34 {:?} delta_LL {:?} rep [{}]", c, size, cu.one_piece, cu.comb, cu.defect, cu.gen_a, cu.gen_b, &cu.delta[..3],
                    ts[s].iter().map(|p| cx.fmt(p)).collect::<Vec<_>>().join("; "));
            }
        }
    } else if has("--selftest") {
        // NEON one_piece vs scalar on random tuples of every degree 2..M (and their moves)
        let mut r: u64 = 0x2545F4914F6CDD1D; let mut rnd = |k: usize| { r ^= r << 13; r ^= r >> 7; r ^= r << 17; (r % k as u64) as usize };
        let (mut cnt, mut yes) = (0u64, 0u64);
        for _ in 0..3_000_000 {
            let n = 2 + rnd(M - 1); let cx = Cx { n }; let mut t = [cx.id(); 4];
            // bias towards one-piece: sometimes use few-cycle elements
            for j in 0..3 { for i in (1..n).rev() { if rnd(3) == 0 { continue; } let k = rnd(i + 1); t[j].swap(i, k); } }
            t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2])));
            for u in [t, mv(&cx, &t, 1), rot(&t)] { let a = one_piece(&cx, &u); assert_eq!(a, one_piece_scalar(&cx, &u), "mismatch at {:?}", u.iter().map(|q| cx.fmt(q)).collect::<Vec<_>>()); cnt += 1; yes += a as u64;
                assert_eq!(defect0(&cx, &u), defect0_scalar(&cx, &u), "defect0 mismatch"); }
        }
        println!("selftest OK: {} states, {} one-piece, fast one_piece == scalar", cnt, yes);
    } else if has("--census") {
        let n: usize = get("--d").unwrap().parse().unwrap();
        let g: usize = get("--genus").unwrap().parse().unwrap();
        let th: usize = get("--threads").map(|x| x.parse().unwrap()).unwrap_or(8);
        let cap: usize = get("--cap").map(|x| x.parse().unwrap()).unwrap_or(20_000_000);
        // Without --out the results go to the terminal and the progress lines are hidden.
        let out = match get("--out") { Some(f) => f, None => { QUIET.store(true, Ordering::Relaxed); "/dev/stdout".to_string() } };
        if has("--all-zero") { ALLZERO.store(true, Ordering::Relaxed); }
        if has("--no-s3norm") { S3OFF.store(true, Ordering::Relaxed); }
        let mb: usize = get("--max-bad").map(|x| x.parse().unwrap()).unwrap_or(usize::MAX);
        if has("--old") || (has("--all-orbits") && !has("--need-good")) || g == 0 {
            census(n, g, has("--all-orbits"), has("--need-good"), mb, th, cap, &out);
        } else {
            census_fast(n, g, has("--need-good"), mb, th, cap, &out, get("--start-item").map(|x| x.parse().unwrap()).unwrap_or(0));
        }
    } else {
        eprintln!("usage:\n  onepiece --orbit \"(..);(..);(..);(..)\" [--clusters]\n  onepiece --census --d D --genus G [--out FILE] [--all-orbits [--need-good] [--max-bad B]] [--threads K] [--cap N]\n  --need-good: only orbits with at least one one-piece cluster (seeds must be one-piece at {{12|34}})\n  without --out, results print to the terminal");
    }
}
