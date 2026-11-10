//! quotient_exhaust: EXHAUSTIVE enumeration of source-genus-3 four-point tuples whose
//! cover factors C -> Y -> P^1 through a FIXED quotient Y of genus h (h = 1 or 0), i.e.
//! an imprimitive tuple with a block system on which the first entry s acts trivially.
//!
//! Every such tuple, up to simultaneous S_d conjugacy and up to moving s to position 1,
//! is (s, t2, t3, t4) with s = (id; sigma), t2 = (x; f), t3 = (y; g), t4 = (s t2 t3)^-1,
//! where (x, y, z) is a degree-m dessin of genus h and (sigma, f, g) a flat S_k-connection
//! on its 2-complex (vertices = blocks, faces = cycles of x, y, xy and the m "sigma faces").
//! Gauge: a primal spanning tree of edges is set to the identity (all vertex gauge used
//! except a global S_k conjugation, broken by fixing the first nontrivial element to a
//! class representative). The remaining edges are the 2h handle edges (free in S_k) and the
//! dual-tree edges, which are SOLVED from prescribed face holonomies. Enumerating all
//! handles, all sigma's and all non-root face holonomies with total relative index <= R
//! (R = 4 - k(2h-2), Riemann-Hurwitz for source genus 3) and checking the root face gives
//! every connection exactly once up to that global conjugation. Dessins are taken one per
//! B3-orbit (moves not involving s act on (x,y,z) only and reach the whole quotient orbit).
//!
//! Each tuple: random-walk smoothness filter (a non-smooth state is a rigorous NEGATIVE),
//! then exact six-move closure modulo S_d with exact Hodge degree
//!   12 deg(lambda)/N = d - sum h(g_i) + mean[h(g1g2)+h(g2g3)+h(g1g2g3g2^-1)],
//! lambda > 0 <=> non-isotrivial for an all-smooth orbit. Caps are UNRESOLVED.
//! Conventions: compose(p,q)[x] = p[q[x]] (rightmost first); g1 g2 g3 g4 = 1.

use std::collections::{HashMap, HashSet};
use std::fs::{self, OpenOptions};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, Mutex};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

// ---------------- fast hash ----------------
#[derive(Default)]
struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 { self.0 }
    fn write(&mut self, bytes: &[u8]) {
        for c in bytes.chunks(8) { let mut b = [0u8; 8]; b[..c.len()].copy_from_slice(c); self.0 = (self.0.rotate_left(5) ^ u64::from_le_bytes(b)).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95); }
    }
}
type FxB = BuildHasherDefault<Fx>;
type FxSet<K> = HashSet<K, FxB>;

// ---------------- permutations (Vec, small) ----------------
type P = Vec<u8>;
fn compose(p: &[u8], q: &[u8]) -> P { q.iter().map(|&x| p[x as usize]).collect() }
fn inverse(p: &[u8]) -> P { let mut q = vec![0u8; p.len()]; for (i, &x) in p.iter().enumerate() { q[x as usize] = i as u8; } q }
fn ident(n: usize) -> P { (0..n as u8).collect() }
fn cycles(p: &[u8]) -> Vec<usize> {
    let n = p.len(); let mut seen = vec![false; n]; let mut out = Vec::new();
    for s in 0..n { if seen[s] { continue; } let mut l = 0; let mut x = s; while !seen[x] { seen[x] = true; x = p[x] as usize; l += 1; } out.push(l); }
    out
}
fn ctype(p: &[u8]) -> Vec<usize> { let mut c: Vec<usize> = cycles(p).into_iter().filter(|&l| l > 1).collect(); c.sort_unstable_by(|a, b| b.cmp(a)); c }
fn index(p: &[u8]) -> usize { p.len() - cycles(p).len() }

// ---------------- fixed-array tuples for the hot path ----------------
const MX: usize = 64;
type A = [u8; MX];
#[inline] fn acomp(p: &A, q: &A, n: usize) -> A { let mut r = [0u8; MX]; for i in 0..n { r[i] = p[q[i] as usize]; } r }
#[inline] fn ainv(p: &A, n: usize) -> A { let mut r = [0u8; MX]; for i in 0..n { r[p[i] as usize] = i as u8; } r }
#[inline] fn aindex(p: &A, n: usize) -> usize {
    let mut seen: u64 = 0; let mut c = 0;
    for s in 0..n { if seen >> s & 1 == 1 { continue; } c += 1; let mut x = s; while seen >> x & 1 == 0 { seen |= 1 << x; x = p[x] as usize; } }
    n - c
}
fn to_a(p: &[u8]) -> A { let mut a = [0u8; MX]; a[..p.len()].copy_from_slice(p); a }
#[derive(Clone, Copy)]
struct Tp { g: [A; 4], idx: [usize; 4], n: usize }
impl Tp {
    fn new(t: &[P; 4]) -> Tp { let n = t[0].len(); let g = [to_a(&t[0]), to_a(&t[1]), to_a(&t[2]), to_a(&t[3])]; Tp { g, idx: [index(&t[0]), index(&t[1]), index(&t[2]), index(&t[3])], n } }
    #[inline] fn mv(&mut self, m: usize) {
        let n = self.n; let j = m % 3; let a = self.g[j]; let b = self.g[j + 1];
        if m < 3 { self.g[j] = acomp(&a, &acomp(&b, &ainv(&a, n), n), n); self.g[j + 1] = a; }
        else { self.g[j] = b; self.g[j + 1] = acomp(&ainv(&b, n), &acomp(&a, &b, n), n); }
        self.idx.swap(j, j + 1);
    }
    #[inline] fn smooth(&self) -> bool {
        // g3 g4 = (g1 g2)^-1 has the same index as g1 g2
        let z = aindex(&acomp(&self.g[0], &self.g[1], self.n), self.n);
        z == self.idx[0] + self.idx[1] || z == self.idx[2] + self.idx[3]
    }
    fn vecs(&self) -> [P; 4] { let n = self.n; [self.g[0][..n].to_vec(), self.g[1][..n].to_vec(), self.g[2][..n].to_vec(), self.g[3][..n].to_vec()] }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0; z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
}

/// rigorous negative if false: some state reached by literal moves is non-smooth
fn walk_ok(t: &Tp, steps: usize, r: &mut Rng) -> bool {
    let mut u = *t; if !u.smooth() { return false; }
    for _ in 0..steps { u.mv(r.below(6)); if !u.smooth() { return false; } }
    true
}

// canonical form modulo simultaneous conjugation (BFS relabelling from every root)
fn canonical(g: &[A; 4], n: usize) -> Vec<u8> {
    let inv = [ainv(&g[0], n), ainv(&g[1], n), ainv(&g[2], n), ainv(&g[3], n)];
    let mut best: Vec<u8> = Vec::new(); let mut lab = [u8::MAX; MX]; let mut order = [0usize; MX];
    let mut key = Vec::with_capacity(3 * n);
    for r in 0..n {
        for v in lab[..n].iter_mut() { *v = u8::MAX; } lab[r] = 0; order[0] = r; let mut len = 1; let mut h = 0;
        while h < len { let x = order[h]; h += 1;
            for p in g.iter().chain(inv.iter()) { let y = p[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y; len += 1; } } }
        key.clear(); let mut state = if best.is_empty() { 1 } else { 0 }; // 0 undecided, 1 better, 2 worse
        'o: for p in g[..3].iter() { for &x in order[..n].iter() {
            let v = lab[p[x] as usize];
            if state == 0 { let b = best[key.len()]; if v < b { state = 1; } else if v > b { state = 2; break 'o; } }
            key.push(v);
        } }
        if state == 1 { best = key.clone(); }
    }
    best
}
fn unpack(key: &[u8], n: usize) -> [A; 4] {
    let a = to_a(&key[..n]); let b = to_a(&key[n..2 * n]); let c = to_a(&key[2 * n..3 * n]);
    let d = ainv(&acomp(&a, &acomp(&b, &c, n), n), n); [a, b, c, d]
}
fn lcm_upto(n: usize) -> u128 { fn g(a: u128, b: u128) -> u128 { if b == 0 { a } else { g(b, a % b) } } let mut l = 1u128; for i in 1..=n as u128 { l = l / g(l, i) * i; } l }
fn hl(p: &A, n: usize, l: u128) -> u128 {
    let mut seen: u64 = 0; let mut s = 0u128;
    for st in 0..n { if seen >> st & 1 == 1 { continue; } let mut c = 0u128; let mut x = st; while seen >> x & 1 == 0 { seen |= 1 << x; x = p[x] as usize; c += 1; } s += l / c; }
    s
}
enum Closure { Closed { states: usize, lam_num: i128, lam_den: u128, keys: Vec<Vec<u8>> }, Nonsmooth(Vec<Vec<u8>>), Capped(usize) }
fn close(t: &Tp, cap: usize, seconds: f64, keep_bad: bool) -> Closure {
    let n = t.n; let l = lcm_upto(n); let start = Instant::now();
    let root = canonical(&t.g, n); let mut seen: FxSet<Vec<u8>> = FxSet::default(); seen.insert(root.clone());
    let mut q = vec![root]; let mut h = 0; let mut s: u128 = 0;
    while h < q.len() {
        let g = unpack(&q[h], n); h += 1;
        let u = Tp { g, idx: [aindex(&g[0], n), aindex(&g[1], n), aindex(&g[2], n), aindex(&g[3], n)], n };
        if !u.smooth() { return Closure::Nonsmooth(if keep_bad { q } else { Vec::new() }); }
        let p12 = acomp(&g[0], &g[1], n); let p23 = acomp(&g[1], &g[2], n);
        let p13 = acomp(&g[0], &acomp(&g[1], &acomp(&g[2], &ainv(&g[1], n), n), n), n);
        s += hl(&p12, n, l) + hl(&p23, n, l) + hl(&p13, n, l);
        for m in 0..6 { let mut w = u; w.mv(m); let k = canonical(&w.g, n); if !seen.contains(&k) { seen.insert(k.clone()); q.push(k); } }
        if q.len() > cap || (h & 255 == 0 && start.elapsed().as_secs_f64() > seconds) { return Closure::Capped(q.len()); }
    }
    let nst = q.len() as i128; let base = (n as i128) * l as i128 - t.g.iter().map(|p| hl(p, n, l) as i128).sum::<i128>();
    Closure::Closed { states: q.len(), lam_num: nst * base + s as i128, lam_den: l, keys: if keep_bad { q } else { Vec::new() } }
}

// ---------------- dessins ----------------
fn all_perms(m: usize) -> Vec<P> {
    let mut out = Vec::new(); let mut a = ident(m);
    fn heap(k: usize, a: &mut P, out: &mut Vec<P>) { if k == 1 { out.push(a.clone()); return; } heap(k - 1, a, out); for i in 0..k - 1 { if k % 2 == 0 { a.swap(i, k - 1) } else { a.swap(0, k - 1) } heap(k - 1, a, out); } }
    heap(m, &mut a, &mut out); out
}
fn transitive_gen(gs: &[&P]) -> bool {
    let n = gs[0].len(); let mut seen = vec![false; n]; seen[0] = true; let mut st = vec![0usize]; let mut c = 1;
    while let Some(x) = st.pop() { for p in gs.iter() { let y = p[x] as usize; if !seen[y] { seen[y] = true; c += 1; st.push(y); } } }
    c == n
}
fn canon3(x: &P, y: &P) -> Vec<u8> { let n = x.len(); let z = inverse(&compose(x, y)); canonical(&[to_a(x), to_a(y), to_a(&z), to_a(&ident(n))], n) }
/// triple moves (x,y,z) with x y z = 1
fn mv3(t: &(P, P, P), m: usize) -> (P, P, P) {
    let v = [t.0.clone(), t.1.clone(), t.2.clone()]; let j = m % 2; let (a, b) = (&v[j], &v[j + 1]); let mut u = v.clone();
    if m < 2 { u[j] = compose(a, &compose(b, &inverse(a))); u[j + 1] = a.clone(); } else { u[j] = b.clone(); u[j + 1] = compose(&inverse(b), &compose(a, b)); }
    (u[0].clone(), u[1].clone(), u[2].clone())
}
/// one representative (x, y, z) per B3-orbit of transitive genus-h triples with all entries nontrivial
fn dessins(m: usize, h: usize) -> Vec<(P, P, P)> {
    let perms = all_perms(m); let mut seen: FxSet<Vec<u8>> = FxSet::default(); let mut out = Vec::new();
    let mut reps: Vec<P> = Vec::new(); let mut tseen = HashSet::new();
    for x in perms.iter() { let c = ctype(x); if tseen.insert(c) { reps.push(x.clone()); } }
    for x in reps.iter() { if index(x) == 0 { continue; } for y in perms.iter() {
        let z = inverse(&compose(x, y));
        if index(x) + index(y) + index(&z) != 2 * m + 2 * h - 2 { continue; }
        if index(y) == 0 || index(&z) == 0 { continue; }
        if !transitive_gen(&[x, y]) { continue; }
        let key = canon3(x, y); if seen.contains(&key) { continue; }
        // whole B3 orbit
        let mut orb = vec![(x.clone(), y.clone(), z.clone())]; let mut ks = vec![key.clone()]; seen.insert(key); let mut i = 0;
        while i < orb.len() { let t = orb[i].clone(); i += 1; for mm in 0..4 { let u = mv3(&t, mm); let k = canon3(&u.0, &u.1); if seen.insert(k.clone()) { ks.push(k); orb.push(u); } } }
        // representative: minimal key
        let (bi, _) = ks.iter().enumerate().min_by(|a, b| a.1.cmp(b.1)).unwrap();
        out.push(orb[bi].clone());
    } }
    out
}

// ---------------- 2-complex of a dessin ----------------
#[derive(Clone, Copy, PartialEq, Debug)]
enum Face { X(usize), Y(usize), Z(usize) }
enum Fac { E(usize), S(usize) }
struct Dessin { m: usize, x: P, y: P, xy: P, xcyc: Vec<Vec<usize>>, ycyc: Vec<Vec<usize>>, zcyc: Vec<Vec<usize>>, xface_of: Vec<usize>, yface_of: Vec<usize>, zface_of: Vec<usize> }
fn cyc_list(p: &[u8]) -> (Vec<Vec<usize>>, Vec<usize>) {
    let n = p.len(); let mut of = vec![usize::MAX; n]; let mut out = Vec::new();
    for s in 0..n { if of[s] != usize::MAX { continue; } let mut c = Vec::new(); let mut x = s; while of[x] == usize::MAX { of[x] = out.len(); c.push(x); x = p[x] as usize; } out.push(c); }
    (out, of)
}
impl Dessin {
    fn new(x: &P, y: &P) -> Dessin {
        let xy = compose(x, y); let (xcyc, xface_of) = cyc_list(x); let (ycyc, yface_of) = cyc_list(y); let (zcyc, zface_of) = cyc_list(&xy);
        Dessin { m: x.len(), x: x.clone(), y: y.clone(), xy, xcyc, ycyc, zcyc, xface_of, yface_of, zface_of }
    }
    fn faces(&self) -> Vec<Face> { (0..self.xcyc.len()).map(Face::X).chain((0..self.ycyc.len()).map(Face::Y)).chain((0..self.zcyc.len()).map(Face::Z)).collect() }
    fn face_id(&self, f: Face) -> usize { match f { Face::X(c) => c, Face::Y(c) => self.xcyc.len() + c, Face::Z(c) => self.xcyc.len() + self.ycyc.len() + c } }
    // x-edge j (f_j: block j -> x j), y-edge j (g_j: block j -> y j)
    fn edge_faces(&self, e: usize) -> (usize, usize) {
        let m = self.m;
        if e < m { let j = e; let jz = inverse(&self.y)[j] as usize; (self.face_id(Face::X(self.xface_of[j])), self.face_id(Face::Z(self.zface_of[jz]))) }
        else { let j = e - m; (self.face_id(Face::Y(self.yface_of[j])), self.face_id(Face::Z(self.zface_of[j]))) }
    }
    // holonomy words, rightmost applied first, transported back to the starting block
    fn word(&self, f: Face) -> Vec<Fac> {
        match f {
            Face::X(c) => self.xcyc[c].iter().map(|&j| Fac::E(j)).collect(),
            Face::Y(c) => self.ycyc[c].iter().map(|&j| Fac::E(self.m + j)).collect(),
            Face::Z(c) => { let mut w = Vec::new(); for &j in self.zcyc[c].iter() { w.push(Fac::E(self.m + j)); w.push(Fac::E(self.y[j] as usize)); w.push(Fac::S(self.xy[j] as usize)); } w }
        }
    }
}
fn eval(w: &[Fac], edges: &[P], sig: &[P], k: usize, skip: usize) -> (P, P) {
    let mut right = ident(k); let mut left = ident(k); let mut after = false;
    for fac in w.iter() {
        let p: &P = match fac { Fac::E(e) => { if *e == skip { after = true; continue; } &edges[*e] } Fac::S(j) => &sig[*j] };
        if after { left = compose(p, &left); } else { right = compose(p, &right); }
    }
    (left, right)
}

/// Gauge data: primal tree edges (identity), dual tree solve order, handle edges, root face
struct Frame { faces: Vec<Face>, handles: Vec<usize>, order: Vec<(usize, usize)>, root: usize, words: Vec<Vec<Fac>>, h2form: Option<(Vec<(usize, i8)>, Vec<(usize, i8)>, Vec<(usize, i8)>, i8)> }
fn frame(d: &Dessin) -> Option<Frame> {
    let m = d.m; let mut tree = vec![false; 2 * m]; let mut seen = vec![false; m]; seen[0] = true; let mut grow = true;
    while grow { grow = false; for e in 0..2 * m { if tree[e] { continue; }
        let (a, b) = if e < m { (e, d.x[e] as usize) } else { (e - m, d.y[e - m] as usize) };
        if seen[a] != seen[b] { seen[a] = true; seen[b] = true; tree[e] = true; grow = true; } } }
    let faces = d.faces(); let nf = faces.len();
    let mut uf: Vec<usize> = (0..nf).collect();
    fn rt(u: &mut Vec<usize>, mut x: usize) -> usize { while u[x] != x { u[x] = u[u[x]]; x = u[x]; } x }
    let mut dual = Vec::new(); let mut handles = Vec::new();
    for e in 0..2 * m { if tree[e] { continue; } let (a, b) = d.edge_faces(e); let (ra, rb) = (rt(&mut uf, a), rt(&mut uf, b)); if ra != rb { uf[ra] = rb; dual.push(e); } else { handles.push(e); } }
    if dual.len() != nf - 1 { return None; }
    let root = d.face_id(Face::Z(d.zcyc.len() - 1));
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); nf];
    for &e in dual.iter() { let (a, b) = d.edge_faces(e); adj[a].push((b, e)); adj[b].push((a, e)); }
    let mut parent = vec![(usize::MAX, usize::MAX); nf]; let mut ord = vec![root]; let mut vis = vec![false; nf]; vis[root] = true; let mut h = 0;
    while h < ord.len() { let f = ord[h]; h += 1; for &(g, e) in adj[f].iter() { if !vis[g] { vis[g] = true; parent[g] = (f, e); ord.push(g); } } }
    if ord.len() != nf { return None; }
    let order: Vec<(usize, usize)> = ord.iter().rev().filter(|&&f| f != root).map(|&f| (f, parent[f].1)).collect();
    let words: Vec<Vec<Fac>> = faces.iter().map(|&f| d.word(f)).collect();
    // symbolic root holonomy in the free group on H_i (0,1), T_f (2+f), S_j (2+nf+j); product order = compose order
    type W = Vec<(usize, i8)>;
    fn red(w: W) -> W { let mut o: W = Vec::new(); for x in w { if let Some(l) = o.last() { if l.0 == x.0 && l.1 == -x.1 { o.pop(); continue; } } o.push(x); } o }
    fn winv(w: &W) -> W { w.iter().rev().map(|&(a, e)| (a, -e)).collect() }
    fn wmul(a: &W, b: &W) -> W { let mut c = a.clone(); c.extend_from_slice(b); red(c) }
    let mut val: Vec<W> = vec![Vec::new(); 2 * m];
    for (i, &e) in handles.iter().enumerate() { val[e] = vec![(i, 1)]; }
    let sv = |j: usize| -> W { vec![(2 + nf + j, 1)] };
    let ev = |w: &Vec<Fac>, val: &Vec<W>, skip: usize| -> (W, W) { let mut l: W = Vec::new(); let mut r: W = Vec::new(); let mut after = false;
        for f in w.iter() { let p: W = match f { Fac::E(e) => { if *e == skip { after = true; continue; } val[*e].clone() } Fac::S(j) => sv(*j) };
            if after { l = wmul(&p, &l); } else { r = wmul(&p, &r); } } (l, r) };
    let words: Vec<Vec<Fac>> = words;
    for &(f, e) in order.iter() { let (l, r) = ev(&words[f], &val, e); val[e] = wmul(&winv(&l), &wmul(&vec![(2 + f, 1)], &winv(&r))); }
    let (_, rw) = ev(&words[root], &val, usize::MAX);
    let pos: Vec<usize> = rw.iter().enumerate().filter(|(_, x)| x.0 == 1).map(|(i, _)| i).collect();
    let h2form = if handles.len() == 2 && pos.len() == 2 && rw[pos[0]].1 == -rw[pos[1]].1 { Some((rw[..pos[0]].to_vec(), rw[pos[0] + 1..pos[1]].to_vec(), rw[pos[1] + 1..].to_vec(), rw[pos[0]].1)) } else { None };
    Some(Frame { faces, handles, order, root, words, h2form })
}

// ---------------- S_k class tables ----------------
struct Classes { k: usize, all: Vec<P>, by_class: Vec<(Vec<usize>, usize, Vec<usize>)>, mul: Vec<u16>, inv: Vec<u16>, ind: Vec<u8>, id: u16, cls_of: Vec<u16>, conj_rep: Vec<u16>, cent: Vec<Vec<u16>> } // (type, index, element ids)
fn classes(k: usize) -> Classes {
    let mut all = all_perms(k); all.sort();
    let mut map: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    for (i, p) in all.iter().enumerate() { map.entry(ctype(p)).or_default().push(i); }
    let mut by_class: Vec<(Vec<usize>, usize, Vec<usize>)> = map.into_iter().map(|(t, ids)| { let ind = t.iter().map(|l| l - 1).sum(); (t, ind, ids) }).collect();
    by_class.sort_by(|a, b| (a.1, &a.0).cmp(&(b.1, &b.0)));
    let ne = all.len(); let pos: HashMap<P, u16> = all.iter().enumerate().map(|(i, p)| (p.clone(), i as u16)).collect();
    // all[] is lexicographically sorted, so the index of a permutation is its Lehmer rank
    let rank = |p: &[u8]| -> u16 { let mut r = 0usize; for i in 0..k { let mut c = 0; for j in i + 1..k { if p[j] < p[i] { c += 1; } } r = r * (k - i) + c; } r as u16 };
    for (i, p) in all.iter().enumerate().take(50) { assert_eq!(rank(p) as usize, i); }
    let mut mul = vec![0u16; ne * ne];
    { let nt = 10usize.min(ne); let chunk = (ne + nt - 1) / nt; let allr = &all;
      std::thread::scope(|sc| { for (ti, part) in mul.chunks_mut(chunk * ne).enumerate() { sc.spawn(move || {
          let mut c = vec![0u8; k];
          for (ai, row) in part.chunks_mut(ne).enumerate() { let a = &allr[ti * chunk + ai]; for (b, slot) in row.iter_mut().enumerate() { let bb = &allr[b]; for x in 0..k { c[x] = a[bb[x] as usize]; } *slot = rank(&c); } } }); } }); }
    for a in (0..ne).step_by(97) { for b in (0..ne).step_by(89) { assert_eq!(mul[a * ne + b], pos[&compose(&all[a], &all[b])]); } }
    let inv: Vec<u16> = all.iter().map(|p| pos[&inverse(p)]).collect(); let ind: Vec<u8> = all.iter().map(|p| index(p) as u8).collect(); let id = pos[&ident(k)];
    // class of each element, a conjugator r_g with r_g rep r_g^-1 = g, centralizer of each class rep
    let mut cls_of = vec![0u16; ne]; let mut conj_rep = vec![0u16; ne]; let mut cent: Vec<Vec<u16>> = Vec::new();
    let sorted_cycles = |p: &P| -> Vec<Vec<usize>> { let n = p.len(); let mut seen = vec![false; n]; let mut cs = Vec::new();
        for st in 0..n { if seen[st] { continue; } let mut c = Vec::new(); let mut x = st; while !seen[x] { seen[x] = true; c.push(x); x = p[x] as usize; } cs.push(c); }
        cs.sort_by(|a, b| b.len().cmp(&a.len())); cs };
    for (ci, (_, _, ids)) in by_class.iter().enumerate() {
        let rep = &all[ids[0]]; let rc = sorted_cycles(rep);
        for &g in ids.iter() { cls_of[g] = ci as u16; let gc = sorted_cycles(&all[g]); let mut r = vec![0u8; k];
            for (a, b) in rc.iter().zip(gc.iter()) { for (x, y) in a.iter().zip(b.iter()) { r[*x] = *y as u8; } }
            let ri = pos[&r]; debug_assert!(compose(&r, &compose(rep, &inverse(&r))) == all[g]); conj_rep[g] = ri; }
        cent.push((0..ne).filter(|&g| compose(&all[g], rep) == compose(rep, &all[g])).map(|g| g as u16).collect());
    }
    for g in 0..ne { let r = &all[conj_rep[g] as usize]; let rep = &all[by_class[cls_of[g] as usize].2[0]]; assert!(compose(r, &compose(rep, &inverse(r))) == all[g]); }
    Classes { k, all, by_class, mul, inv, ind, id, cls_of, conj_rep, cent }
}

// ---------------- the exhaustive job ----------------
#[derive(Default, Clone)]
struct Cnt { built: u64, walkpass: u64, closed: u64, closed_pos: u64, nonsmooth: u64, capped: u64, dup: u64 }
struct JobOut { counts: HashMap<String, Cnt>, closures: Vec<String>, capped: Vec<String>, cand: u64, keys: FxSet<Vec<u8>>, solve_bad: u64 }

fn notation(p: &[u8]) -> String {
    let mut seen = vec![false; p.len()]; let mut s = String::new();
    for a in 0..p.len() { if seen[a] { continue; } let mut c = Vec::new(); let mut x = a; while !seen[x] { seen[x] = true; c.push(x + 1); x = p[x] as usize; }
        if c.len() > 1 { s.push('('); s.push_str(&c.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ")); s.push(')'); } }
    if s.is_empty() { "()".into() } else { s }
}
fn gcd(a: u128, b: u128) -> u128 { if b == 0 { a } else { gcd(b, a % b) } }
fn tstr(t: &[P; 4]) -> String { t.iter().map(|p| format!("\"{}\"", notation(p))).collect::<Vec<_>>().join(",") }

struct Params { walk: usize, cap: usize, secs: f64, stop: Arc<AtomicBool>, collect: bool, solve_h2: bool }

/// all class placements (blocks, then non-root faces) with total index <= R and root remainder <= k-1
fn assignments(d: &Dessin, h: usize, cl: &Classes) -> Vec<(Vec<usize>, i64)> {
    let k = cl.k; let m = d.m; let r_total: i64 = 4 - (k as i64) * (2 * h as i64 - 2);
    let Some(fr) = frame(d) else { return Vec::new() };
    let npos = m + fr.faces.len() - 1; let ncl = cl.by_class.len();
    let usable: Vec<usize> = (1..ncl).filter(|&c| (cl.by_class[c].1 as i64) <= r_total).collect();
    let mut out = Vec::new(); let mut a = vec![0usize; npos];
    fn rec(pos: usize, left: i64, a: &mut Vec<usize>, usable: &[usize], cl: &Classes, k: usize, out: &mut Vec<(Vec<usize>, i64)>) {
        if pos == a.len() { if left <= k as i64 - 1 { out.push((a.clone(), left)); } return; }
        a[pos] = 0; rec(pos + 1, left, a, usable, cl, k, out);
        for &c in usable { let ind = cl.by_class[c].1 as i64; if ind <= left { a[pos] = c; rec(pos + 1, left - ind, a, usable, cl, k, out); } }
        a[pos] = 0;
    }
    rec(0, r_total, &mut a, &usable, cl, k, &mut out); out
}
fn run_job(d: &Dessin, h: usize, cl: &Classes, prm: &Params, rng_seed: u64, which: &[(Vec<usize>, i64)]) -> Option<JobOut> {
    let k = cl.k; let m = d.m; let n = m * k; if n > MX { return None; }
    let fr = frame(d)?;
    let nf = fr.faces.len();
    let nonroot: Vec<usize> = (0..nf).filter(|&f| f != fr.root).collect();
    let npos = m + nonroot.len();
    let ne = cl.all.len(); let mul = &cl.mul; let inv = &cl.inv; let ind = &cl.ind; let id = cl.id;
    // compiled words: variable ids  e (edge) -> e ; sigma j -> 2m + j
    let words: Vec<Vec<usize>> = fr.words.iter().map(|w| w.iter().map(|f| match f { Fac::E(e) => *e, Fac::S(j) => 2 * m + *j }).collect()).collect();
    // for each solved (face, edge): word split around the edge
    let solve: Vec<(usize, usize, Vec<usize>, Vec<usize>)> = fr.order.iter().map(|&(f, e)| {
        let w = &words[f]; let p = w.iter().position(|&v| v == e).unwrap(); (f, e, w[..p].to_vec(), w[p + 1..].to_vec()) }).collect();
    let root_w = words[fr.root].clone();
    let mut out = JobOut { counts: HashMap::new(), closures: Vec::new(), capped: Vec::new(), cand: 0, keys: FxSet::default(), solve_bad: 0 };
    let mut seen_keys: FxSet<Vec<u8>> = FxSet::default();
    let mut bad: FxSet<Vec<u8>> = FxSet::default();
    let mut rng = Rng(rng_seed);
    let mut vars = vec![id; 3 * m]; // edges 0..2m, sigma 2m..3m
    let mut target = vec![id; nf];
    let solved = prm.solve_h2 && fr.h2form.is_some(); let solved_root_check = solved; let mut h2s: Vec<u16> = Vec::new();
    let allids: Vec<u16> = (0..ne as u16).collect(); let mut vis = vec![false; ne];
    for (asg, left) in which.iter() {
        let (asg, left) = (asg, *left);
        if prm.stop.load(Relaxed) { return None; }
        if left > (k as i64 - 1) { continue; }
        let nh = fr.handles.len();
        let brute2 = nh == 2 && !solved;
        // free variables in order: nontrivial positions, handle 0, (handle 1 when not solved)
        let mut fvars: Vec<(u8, usize)> = Vec::new(); let mut fcands: Vec<Vec<u16>> = Vec::new();
        for p in 0..npos { if asg[p] != 0 { fvars.push((0, p)); fcands.push(cl.by_class[asg[p]].2.iter().map(|&x| x as u16).collect()); } }
        if nh >= 1 { fvars.push((1, 0)); fcands.push(allids.clone()); }
        if brute2 { fvars.push((1, 1)); fcands.push(allids.clone()); }
        let roots: Vec<u16> = (0..ne).filter(|&g| ind[g] as i64 == left).map(|g| g as u16).collect();
        let mut mov: Vec<String> = Vec::new(); let mut fix: Vec<String> = Vec::new();
        for p in 0..npos { if asg[p] != 0 { let t = format!("{:?}", cl.by_class[asg[p]].0).replace(' ', ""); if p < m { mov.push(t) } else { fix.push(t) } } }
        mov.sort(); fix.sort();
        for v in vars.iter_mut() { *v = id; } for t in target.iter_mut() { *t = id; }
        // stabilizer-chain DFS: each free variable only up to conjugation by the stabilizer of earlier choices
        let mk_level = |cand: &Vec<u16>, stab: &Vec<u16>, vis: &mut Vec<bool>| -> Lv {
            if stab.len() <= 1 { return Lv { reps: cand.clone(), stabs: Vec::new(), trivial: true, i: 0 }; }
            let mut reps = Vec::new(); let mut stabs = Vec::new();
            for &e in cand.iter() { vis[e as usize] = false; }
            for &e in cand.iter() { if vis[e as usize] { continue; }
                for &g in stab.iter() { let c = mul[mul[g as usize * ne + e as usize] as usize * ne + inv[g as usize] as usize]; vis[c as usize] = true; }
                reps.push(e); stabs.push(stab.iter().copied().filter(|&g| mul[g as usize * ne + e as usize] == mul[e as usize * ne + g as usize]).collect::<Vec<u16>>()); }
            Lv { reps, stabs, trivial: false, i: 0 }
        };
        let mut stack: Vec<Lv> = Vec::new();
        if !fvars.is_empty() { let lv = mk_level(&fcands[0], &allids, &mut vis); stack.push(lv); }
        let mut once = fvars.is_empty();
        loop {
            if once { once = false; } else {
                let mut found = false;
                loop {
                    let depth = match stack.len() { 0 => break, l => l - 1 };
                    let (e, child) = { let top = stack.last_mut().unwrap();
                        if top.i >= top.reps.len() { stack.pop(); continue; }
                        let e = top.reps[top.i]; let ch = if top.trivial { None } else { Some(top.stabs[top.i].clone()) }; top.i += 1; (e, ch) };
                    match fvars[depth] { (0, p) => { if p < m { vars[2 * m + p] = e; } else { target[nonroot[p - m]] = e; } } (_, i) => { vars[fr.handles[i]] = e; } }
                    if depth + 1 == fvars.len() { found = true; break; }
                    let st = child.unwrap_or_else(|| vec![id]);
                    let lv = mk_level(&fcands[depth + 1], &st, &mut vis); stack.push(lv);
                }
                if !found { break; }
            }
            h2s.clear();
            if nh == 2 {
                if let (true, Some((uw, ww, yw, e2))) = (solved, &fr.h2form) {
                    let symv = |sy: usize, vars: &Vec<u16>, target: &Vec<u16>| -> u16 { if sy == 0 { vars[fr.handles[0]] } else if sy < 2 + nf { target[sy - 2] } else { vars[2 * m + sy - 2 - nf] } };
                    let evw = |w: &Vec<(usize, i8)>, vars: &Vec<u16>, target: &Vec<u16>| -> u16 { let mut v = id; for &(sy, e) in w.iter() { let mut x = symv(sy, vars, target); if e < 0 { x = inv[x as usize]; } v = mul[v as usize * ne + x as usize]; } v };
                    let (u, w, y) = (evw(uw, &vars, &target), evw(ww, &vars, &target), evw(yw, &vars, &target));
                    let (ui, yi) = (inv[u as usize] as usize, inv[y as usize] as usize);
                    for &r in roots.iter() {
                        let v = mul[mul[ui * ne + r as usize] as usize * ne + yi];
                        if cl.cls_of[v as usize] != cl.cls_of[w as usize] { continue; }
                        // e2=+1: H2 w H2^-1 = v ;  e2=-1: H2 v H2^-1 = w
                        let (a, b) = if *e2 > 0 { (v, w) } else { (w, v) }; // want H2 b H2^-1 = a
                        let (ra, rbi) = (cl.conj_rep[a as usize] as usize, inv[cl.conj_rep[b as usize] as usize] as usize);
                        for &c in cl.cent[cl.cls_of[b as usize] as usize].iter() { h2s.push(mul[mul[ra * ne + c as usize] as usize * ne + rbi]); }
                    }
                } else { h2s.push(vars[fr.handles[1]]); }
            } else { h2s.push(u16::MAX); }
            for hi in 0..h2s.len() { if nh == 2 { vars[fr.handles[1]] = h2s[hi]; }
            for (f, ed, pre, post) in solve.iter() {
                // holonomy = L * e * R with R = product of pre (applied first), L = product of post
                let mut r = id; for &v in pre.iter() { r = mul[vars[v] as usize * ne + r as usize]; }
                let mut l = id; for &v in post.iter() { l = mul[vars[v] as usize * ne + l as usize]; }
                let e = mul[mul[inv[l as usize] as usize * ne + target[*f] as usize] as usize * ne + inv[r as usize] as usize];
                vars[*ed] = e;
            }
            let mut rh = id; for &v in root_w.iter() { rh = mul[vars[v] as usize * ne + rh as usize]; }
            if ind[rh as usize] as i64 == left {
                out.cand += 1;
                let mut tp = Tp { g: [[0u8; MX]; 4], idx: [0; 4], n };
                for j in 0..m { let sg = &cl.all[vars[2 * m + j] as usize]; let f = &cl.all[vars[j] as usize]; let g = &cl.all[vars[m + j] as usize];
                    for i in 0..k { tp.g[0][j * k + i] = (j * k + sg[i] as usize) as u8; tp.g[1][j * k + i] = (d.x[j] as usize * k + f[i] as usize) as u8; tp.g[2][j * k + i] = (d.y[j] as usize * k + g[i] as usize) as u8; } }
                tp.g[3] = ainv(&acomp(&tp.g[0], &acomp(&tp.g[1], &tp.g[2], n), n), n);
                for i in 0..4 { tp.idx[i] = aindex(&tp.g[i], n); }
                let ok = tp.idx.iter().all(|&x| x > 0) && tp.idx.iter().sum::<usize>() == 2 * n + 4 && atransitive(&tp);
                if ok && prm.collect {
                    out.keys.insert(canonical(&tp.g, n));
                    for &(f, _) in fr.order.iter() { let mut hh = id; for &v in words[f].iter() { hh = mul[vars[v] as usize * ne + hh as usize]; } if hh != target[f] { out.solve_bad += 1; } }
                } else if ok {
                    let label = if left == 0 { format!("k{} mov[{}] fix[{}]", k, mov.join(","), fix.join(",")) } else {
                        let mut f2 = fix.clone(); f2.push(format!("{:?}", ctype(&cl.all[rh as usize])).replace(' ', "")); f2.sort(); format!("k{} mov[{}] fix[{}]", k, mov.join(","), f2.join(",")) };
                    let c = out.counts.entry(label.clone()).or_default(); c.built += 1;
                    if walk_ok(&tp, prm.walk, &mut rng) {
                        c.walkpass += 1;
                        let key = canonical(&tp.g, n);
                        if seen_keys.contains(&key) || bad.contains(&key) { c.dup += 1; }
                        else {
                            seen_keys.insert(key.clone());
                            let tup = tp.vecs();
                            match close(&tp, prm.cap, prm.secs, true) {
                                Closure::Nonsmooth(states) => { c.nonsmooth += 1; if bad.len() < 3_000_000 { for st in states { bad.insert(st); } } }
                                Closure::Capped(sz) => { c.capped += 1; out.capped.push(format!("{{\"states_seen\":{},\"class\":\"{}\",\"degree\":{},\"tuple\":[{}]}}", sz, label, n, tstr(&tup))); }
                                Closure::Closed { states, lam_num, lam_den, keys } => {
                                    for kk in keys { seen_keys.insert(kk); }
                                    c.closed += 1; if lam_num > 0 { c.closed_pos += 1; }
                                    let g = gcd(lam_num.unsigned_abs(), lam_den * states as u128); let g = if g == 0 { 1 } else { g };
                                    out.closures.push(format!("{{\"states\":{},\"lam12_per_state\":\"{}/{}\",\"nonisotrivial\":{},\"class\":\"{}\",\"degree\":{},\"tuple\":[{}]}}",
                                        states, lam_num / g as i128, lam_den * states as u128 / g, lam_num > 0, label, n, tstr(&tup)));
                                }
                            }
                        }
                    }
                }
            }
            if solved_root_check && ind[rh as usize] as i64 != left { out.solve_bad += 1; }
            }
        }
    }
    Some(out)
}
struct Lv { reps: Vec<u16>, stabs: Vec<Vec<u16>>, trivial: bool, i: usize }
fn atransitive(t: &Tp) -> bool {
    let n = t.n; let mut seen: u64 = 1; let mut st = [0u8; MX]; let mut sp = 1; st[0] = 0;
    while sp > 0 { sp -= 1; let x = st[sp] as usize; for p in t.g[..3].iter() { let y = p[x] as usize; if seen >> y & 1 == 0 { seen |= 1 << y; st[sp] = y as u8; sp += 1; } } }
    seen.count_ones() as usize == n
}

// ---------------- face-word self check against actual lifted cycle structure ----------------
fn check_faces(d: &Dessin, k: usize, r: &mut Rng) -> bool {
    let m = d.m; let n = m * k;
    let rp = |r: &mut Rng| -> P { let mut v = ident(k); for i in (1..k).rev() { let j = r.below(i + 1); v.swap(i, j); } v };
    let sig: Vec<P> = (0..m).map(|_| rp(r)).collect(); let edges: Vec<P> = (0..2 * m).map(|_| rp(r)).collect();
    let mut s = vec![0u8; n]; let mut t2 = vec![0u8; n]; let mut t3 = vec![0u8; n];
    for j in 0..m { for i in 0..k {
        s[j * k + i] = (j * k + sig[j][i] as usize) as u8; t2[j * k + i] = (d.x[j] as usize * k + edges[j][i] as usize) as u8; t3[j * k + i] = (d.y[j] as usize * k + edges[m + j][i] as usize) as u8; } }
    let st2t3 = compose(&s, &compose(&t2, &t3));
    // lifted cycle type over a quotient cycle c of p equals {l*|c| : l in cycles(hol)}
    let lift_type = |p: &P, blocks: &Vec<usize>| -> Vec<usize> {
        let mut seen = vec![false; n]; let mut out = Vec::new();
        for &b in blocks { for i in 0..k { let st = b * k + i; if seen[st] { continue; } let mut l = 0; let mut x = st; while !seen[x] { seen[x] = true; x = p[x] as usize; l += 1; } out.push(l); } }
        out.sort(); out
    };
    for f in d.faces() {
        let w = d.word(f); let (_, hol) = eval(&w, &edges, &sig, k, usize::MAX);
        let (p, blocks) = match f { Face::X(c) => (&t2, &d.xcyc[c]), Face::Y(c) => (&t3, &d.ycyc[c]), Face::Z(c) => (&st2t3, &d.zcyc[c]) };
        let mut want: Vec<usize> = cycles(&hol).iter().map(|l| l * blocks.len()).collect(); want.sort();
        if lift_type(p, blocks) != want { eprintln!("face {:?} m{} k{} lift {:?} want {:?} hol {:?}", f, m, k, lift_type(p, blocks), want, hol); return false; }
    }
    for j in 0..m { let mut want = cycles(&sig[j]); want.sort(); if lift_type(&s, &vec![j]) != want { eprintln!("sigma mismatch"); return false; } }
    true
}

fn now() -> f64 { SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64() }
fn append(path: &Path, line: &str) { let mut f = OpenOptions::new().create(true).append(true).open(path).unwrap(); f.write_all(line.as_bytes()).unwrap(); f.write_all(b"\n").unwrap(); f.sync_all().unwrap(); }

fn parse(parts: &[&str], n: usize) -> [P; 4] {
    let mut out: Vec<P> = Vec::new();
    for s in parts { let mut p = ident(n); for row in s.split('(').skip(1) { let row = row.split(')').next().unwrap();
        let c: Vec<usize> = row.split_whitespace().map(|v| v.parse::<usize>().unwrap() - 1).collect(); for i in 0..c.len() { p[c[i]] = c[(i + 1) % c.len()] as u8; } } out.push(p); }
    [out[0].clone(), out[1].clone(), out[2].clone(), out[3].clone()]
}
fn selftest() {
    let show = |t: &[P; 4]| match close(&Tp::new(t), 100000, 30.0, false) { Closure::Closed { states, lam_num, lam_den, .. } => { let g = gcd(lam_num.unsigned_abs(), lam_den * states as u128).max(1); format!("closed N={} 12lam/N={}/{}", states, lam_num / g as i128, lam_den * states as u128 / g) } Closure::Nonsmooth(_) => "nonsmooth".into(), Closure::Capped(_) => "capped".into() };
    let t: [P; 4] = [vec![1, 0], vec![1, 0], vec![1, 0], vec![1, 0]];
    println!("legendre: {} (expect 12lam/N = 6 is for nodal family; legendre is not all-smooth)", show(&t));
    let t = parse(&["(5 6 7)", "(1 2 5)(3 7 9)(4 6 8)", "(1 3 7)(2 6 9)(4 5 8)", "(1 4 5)(2 7 8)(3 6 9)"], 9);
    println!("deg9 g2 iso: {} (expect N=12, 0)", show(&t));
    let t = parse(&["(1 4)(2 3)", "(1 5 9 13 17)(2 6 10 14 18)(3 7 11 15 19)(4 8 12 16 20)", "(1 7 12 16 18)(2 5 11 15 17)(3 8 10 14 20)(4 6 9 13 19)", "(1 14 8 19 11)(2 13 6 17 9)(3 16 7 20 12)(4 15 5 18 10)"], 20);
    println!("g6 noniso: {} (expect N=240, 6/5)", show(&t));
    let t = parse(&["(5 6)(15 16)", "(1 13 5 9 21 17)(2 14 6 10 22 18)(3 15 7 11 23 19)(4 16 8 12 24 20)", "(1 12)(2 11)(3 9)(4 10)(5 19)(6 20)(7 18)(8 17)(13 23)(14 24)(15 22)(16 21)", "(1 8 21 3 5 24)(2 7 22 4 6 23)(9 19 13 12 17 16)(10 20 14 11 18 15)"], 24);
    println!("g4 unknown->noniso: {} (expect N=192, 1/2)", show(&t));
    let mut r = Rng(11); let mut bad = 0; let mut tot = 0;
    for h in [0usize, 1] { for m in 2..=5 { for (x, y, _) in dessins(m, h) { let d = Dessin::new(&x, &y); for k in 2..=4 { for _ in 0..3 { tot += 1; if !check_faces(&d, k, &mut r) { bad += 1; } } } } } }
    println!("face-word check: {} random connections, {} mismatches", tot, bad);
    for h in [1usize, 0] { let c: Vec<usize> = (1..=6).map(|m| dessins(m, h).len()).collect(); println!("genus-{} dessins up to B3, m=1..6: {:?}", h, c); }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str, d: &str| -> String { args.iter().position(|a| a == k).map(|i| args[i + 1].clone()).unwrap_or(d.to_string()) };
    if args.iter().any(|a| a == "--selftest") { selftest(); return; }
    if args.iter().any(|a| a == "--crosscheck2") { crosscheck2(get("--genus", "1").parse().unwrap(), get("--m", "3").parse().unwrap(), get("--k", "3").parse().unwrap(), get("--maxd", "1000").parse().unwrap()); return; }
    if args.iter().any(|a| a == "--crosscheck") { crosscheck(get("--genus", "1").parse().unwrap(), get("--m", "3").parse().unwrap(), get("--k", "3").parse().unwrap()); return; }
    let out = PathBuf::from(get("--out", "qx_out")); fs::create_dir_all(&out).unwrap();
    let threads: usize = get("--threads", "8").parse().unwrap();
    let h: usize = get("--genus", "1").parse().unwrap();
    let ms: Vec<usize> = get("--m", "3,4,5").split(',').map(|v| v.parse().unwrap()).collect();
    let ks: Vec<usize> = get("--k", "3").split(',').map(|v| v.parse().unwrap()).collect();
    let walk: usize = get("--walk", "200").parse().unwrap();
    let cap: usize = get("--cap", "400000").parse().unwrap();
    let secs: f64 = get("--closure-seconds", "120").parse().unwrap();
    let stop = Arc::new(AtomicBool::new(false));
    let solve_h2 = !args.iter().any(|a| a == "--brute-h2");
    // jobs: (m, k, dessin index, assignment range); resumable via DONE keys
    let filter = get("--filter", "all");
    let chunk_w: f64 = get("--chunk-weight", "3e6").parse().unwrap();
    let done: HashSet<String> = fs::read_to_string(out.join("done.jsonl")).unwrap_or_default().lines().filter_map(|l| l.split("\"job\":\"").nth(1).map(|s| s.split('"').next().unwrap().to_string())).collect();
    let classes_k: Arc<HashMap<usize, Classes>> = Arc::new(ks.iter().map(|&k| (k, classes(k))).collect());
    let mut jobs: Vec<(usize, usize, usize, usize, usize)> = Vec::new(); let mut ds: HashMap<usize, Arc<Vec<(P, P, P)>>> = HashMap::new();
    let mut total_w = 0f64;
    for &m in ms.iter() { let t0 = Instant::now(); let v = dessins(m, h); eprintln!("genus-{} dessins m={}: {} B3-orbits ({:.1}s)", h, m, v.len(), t0.elapsed().as_secs_f64());
        for &k in ks.iter() { if m * k > MX { continue; } let cl = &classes_k[&k];
            for i in 0..v.len() { let d = Dessin::new(&v[i].0, &v[i].1); let al = assignments(&d, h, cl);
                let hw = (cl.all.len() as f64).powi(2 * h as i32);
                let mut c0 = 0; let mut w = 0f64;
                for (ai, (asg, left)) in al.iter().enumerate() {
                    let keep = match filter.as_str() { "iso" => asg[..m].iter().filter(|&&c| c != 0).count() == 1 && asg[m..].iter().all(|&c| c == 0) && *left == 0, "noniso" => !(asg[..m].iter().filter(|&&c| c != 0).count() == 1 && asg[m..].iter().all(|&c| c == 0) && *left == 0), _ => true };
                    let mut aw = hw; let mut first = true; for &c in asg.iter() { if c != 0 { if first { first = false; } else { aw *= cl.by_class[c].2.len() as f64; } } }
                    if !keep { aw = 0.0; }
                    w += aw; total_w += aw;
                    if w >= chunk_w || ai + 1 == al.len() { if w > 0.0 { jobs.push((m, k, i, c0, ai + 1)); } c0 = ai + 1; w = 0.0; }
                } } } ds.insert(m, Arc::new(v)); }
    jobs.retain(|&(m, k, i, a0, a1)| !done.contains(&format!("h{}m{}k{}d{}a{}-{}", h, m, k, i, a0, a1)));
    jobs.sort_by_key(|&(m, k, i, a0, _)| (m * k, k, m, i, a0));
    eprintln!("jobs pending {} total candidate weight {:.3e}", jobs.len(), total_w);
    if args.iter().any(|a| a == "--plan") { return; }
    append(&out.join("CONFIG.jsonl"), &format!("{{\"started\":{:.0},\"genus\":{},\"m\":{:?},\"k\":{:?},\"walk\":{},\"cap\":{},\"closure_seconds\":{},\"threads\":{},\"filter\":\"{}\",\"jobs_pending\":{},\"weight\":{:.3e}}}", now(), h, ms, ks, walk, cap, secs, threads, filter, jobs.len(), total_w));
    let jobs = Arc::new(jobs); let next = Arc::new(AtomicUsize::new(0)); let lock = Arc::new(Mutex::new(()));
    let mut hs = Vec::new();
    for _ in 0..threads {
        let (jobs, next, ds, out, stop, lock, ck, filter) = (jobs.clone(), next.clone(), ds.clone(), out.clone(), stop.clone(), lock.clone(), classes_k.clone(), filter.clone());
        hs.push(std::thread::spawn(move || {
            let prm = Params { walk, cap, secs, stop: stop.clone(), collect: false, solve_h2 };
            loop {
                let ji = next.fetch_add(1, Relaxed); if ji >= jobs.len() || stop.load(Relaxed) { break; }
                let (m, k, i, a0, a1) = jobs[ji]; let (x, y, z) = &ds[&m][i]; let d = Dessin::new(x, y);
                let t0 = Instant::now();
                let al = assignments(&d, h, &ck[&k]);
                let sel: Vec<(Vec<usize>, i64)> = al[a0..a1].iter().filter(|(asg, left)| match filter.as_str() { "iso" => asg[..m].iter().filter(|&&c| c != 0).count() == 1 && asg[m..].iter().all(|&c| c == 0) && *left == 0, "noniso" => !(asg[..m].iter().filter(|&&c| c != 0).count() == 1 && asg[m..].iter().all(|&c| c == 0) && *left == 0), _ => true }).cloned().collect();
                let res = run_job(&d, h, &ck[&k], &prm, 0x1234_5678 ^ ((i * 1000003 + a0) as u64).wrapping_mul(0x9E3779B97F4A7C15), &sel);
                let Some(o) = res else { if !stop.load(Relaxed) { let _g = lock.lock().unwrap(); append(&out.join("done.jsonl"), &format!("{{\"job\":\"h{}m{}k{}d{}a{}-{}\",\"skipped\":true}}", h, m, k, i, a0, a1)); } continue; };
                let _g = lock.lock().unwrap();
                for c in o.closures.iter() { append(&out.join("closed.jsonl"), &format!("{{\"job\":\"h{}m{}k{}d{}a{}-{}\",{}", h, m, k, i, a0, a1, &c[1..])); if c.contains("\"nonisotrivial\":true") { append(&out.join("POSITIVE_CANDIDATES.jsonl"), &format!("{{\"job\":\"h{}m{}k{}d{}a{}-{}\",{}", h, m, k, i, a0, a1, &c[1..])); } }
                for c in o.capped.iter() { append(&out.join("capped.jsonl"), &format!("{{\"job\":\"h{}m{}k{}d{}a{}-{}\",{}", h, m, k, i, a0, a1, &c[1..])); }
                let mut cs: Vec<_> = o.counts.iter().collect(); cs.sort_by(|a, b| a.0.cmp(b.0));
                let cstr = cs.iter().map(|(l, c)| format!("\"{}\":[{},{},{},{},{},{},{}]", l, c.built, c.walkpass, c.dup, c.nonsmooth, c.capped, c.closed, c.closed_pos)).collect::<Vec<_>>().join(",");
                append(&out.join("done.jsonl"), &format!("{{\"job\":\"h{}m{}k{}d{}a{}-{}\",\"t\":{:.0},\"secs\":{:.2},\"dessin\":[\"{}\",\"{}\",\"{}\"],\"candidates\":{},\"counts_built_walkpass_dup_nonsmooth_capped_closed_closedpos\":{{{}}}}}",
                    h, m, k, i, a0, a1, now(), t0.elapsed().as_secs_f64(), notation(x), notation(y), notation(z), o.cand, cstr));
            }
        }));
    }
    let t0 = Instant::now();
    loop {
        std::thread::sleep(std::time::Duration::from_secs(5));
        if out.join("STOP").exists() || out.parent().map(|p| p.join("STOP").exists()).unwrap_or(false) { stop.store(true, Relaxed); }
        let fin = hs.iter().all(|h| h.is_finished());
        let hb = format!("{{\"t\":{:.0},\"elapsed\":{:.0},\"jobs_total\":{},\"jobs_issued\":{}}}", now(), t0.elapsed().as_secs_f64(), jobs.len(), next.load(Relaxed).min(jobs.len()));
        fs::write(out.join("HEARTBEAT.json"), hb).unwrap();
        if fin { break; }
    }
    for h in hs { h.join().unwrap(); }
    append(&out.join("RECEIPT.jsonl"), &format!("{{\"finished\":{:.0},\"stopped\":{},\"elapsed\":{:.0}}}", now(), stop.load(Relaxed), t0.elapsed().as_secs_f64()));
}

/// brute force over ALL (sigma, f, g) in S_k^{3m} versus the gauge-fixed enumerator
fn crosscheck(h: usize, m: usize, k: usize) {
    let cl = classes(k); let perms = &cl.all; let np = perms.len();
    let r_total: i64 = 4 - (k as i64) * (2 * h as i64 - 2);
    for (di, (x, y, z)) in dessins(m, h).iter().enumerate() {
        let d = Dessin::new(x, y); let n = m * k;
        let prm = Params { walk: 0, cap: 0, secs: 0.0, stop: Arc::new(AtomicBool::new(false)), collect: true, solve_h2: false };
        let al = assignments(&d, h, &cl); let o = run_job(&d, h, &cl, &prm, 1, &al).unwrap();
        let mut brute: FxSet<Vec<u8>> = FxSet::default(); let total = np.pow(3 * m as u32); let mut ctr = vec![0usize; 3 * m];
        for _ in 0..total {
            let mut s = vec![0u8; n]; let mut t2 = vec![0u8; n]; let mut t3 = vec![0u8; n];
            for j in 0..m { let (sg, f, g) = (&perms[ctr[j]], &perms[ctr[m + j]], &perms[ctr[2 * m + j]]); for i in 0..k {
                s[j * k + i] = (j * k + sg[i] as usize) as u8; t2[j * k + i] = (d.x[j] as usize * k + f[i] as usize) as u8; t3[j * k + i] = (d.y[j] as usize * k + g[i] as usize) as u8; } }
            let t4 = inverse(&compose(&s, &compose(&t2, &t3))); let tup: [P; 4] = [s, t2, t3, t4];
            let ok = tup.iter().all(|p| index(p) > 0) && tup.iter().map(|p| index(p)).sum::<usize>() == 2 * n + 4 && transitive_gen(&[&tup[0], &tup[1], &tup[2]]);
            if ok { let tp = Tp::new(&tup); brute.insert(canonical(&tp.g, n)); }
            let mut i = 0; while i < ctr.len() { ctr[i] += 1; if ctr[i] < np { break; } ctr[i] = 0; i += 1; }
        }
        let _ = r_total;
        println!("h{} m{} k{} dessin{} [{} {} {}]: enumerator {} classes (cand {}, solve mismatches {}), brute {} classes, equal={}",
            h, m, k, di, notation(x), notation(y), notation(z), o.keys.len(), o.cand, o.solve_bad, brute.len(), o.keys == brute);
    }
}

/// solved-H2 enumerator versus brute-H2 enumerator: identical canonical tuple sets
fn crosscheck2(h: usize, m: usize, k: usize, maxd: usize) {
    let cl = classes(k); let mut all_eq = true; let mut nd = 0;
    for (di, (x, y, _)) in dessins(m, h).iter().enumerate().take(maxd) {
        let d = Dessin::new(x, y); let al = assignments(&d, h, &cl);
        let mk = |s: bool| Params { walk: 0, cap: 0, secs: 0.0, stop: Arc::new(AtomicBool::new(false)), collect: true, solve_h2: s };
        let t0 = Instant::now(); let a = run_job(&d, h, &cl, &mk(false), 1, &al).unwrap(); let ta = t0.elapsed().as_secs_f64();
        let t1 = Instant::now(); let b = run_job(&d, h, &cl, &mk(true), 1, &al).unwrap(); let tb = t1.elapsed().as_secs_f64();
        let eq = a.keys == b.keys; all_eq &= eq; nd += 1;
        println!("h{} m{} k{} d{}: brute {} classes ({:.2}s, cand {})  solved {} classes ({:.2}s, cand {}, root/solve mismatches {})  equal={}", h, m, k, di, a.keys.len(), ta, a.cand, b.keys.len(), tb, b.cand, b.solve_bad, eq);
    }
    println!("ALL EQUAL over {} dessins: {}", nd, all_eq);
}
