//! Lift search (Nov 7). Lifts of a base four-point orbit (degree b, genus h) to degree
//! m*b and genus g, preserving a system of b blocks of size m whose quotient is the base.
//! Goal: genus-3 orbits that are one-piece at every state, over the genus-1 bases that are
//! one-piece everywhere (one-piece descends to block quotients; see README section 8).
//!
//! Same conventions as src/main.rs: composition rightmost-first, (x*y)[i] = x[y[i]],
//! g1 g2 g3 g4 = 1, six moves, states modulo simultaneous conjugation in S_d.
//!
//! A lift over the base state (a1..a4) is g_i(x, j) = (a_i x, V_i[x] j) for i = 1,2,3
//! (letter x*m + j), and g4 = (g1 g2 g3)^-1. Gauge: V_i[x] -> H[a_i x] V_i[x] H[x]^-1.
//! Voltages on a spanning tree are set to 1; the other 2b+1 (genus-1 base) are unknowns.
//! The return of a base cycle c of a_i is the product of voltages around c; its index e_c
//! satisfies sum_c e_c = 2g - 2 - m(2h - 2) (relative Riemann-Hurwitz). A "pattern" fixes
//! which cycles have a nontrivial return and its cycle type. The trivial-return relations
//! are solved by Tietze elimination (a variable that occurs once in a relator is expressed
//! by the others); designated relators introduce a variable ranging over the prescribed
//! class. The remaining variables are enumerated, the first free one up to conjugacy.
//! `--brute` instead enumerates every non-tree voltage assignment (independent check).
//!
//! Every lift is checked directly on the permutations: genus, transitivity, the cycle
//! types of the returns, one-piece at three cusps, then a six-move closure with early
//! exit at the first state that is not one-piece.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

const NM: usize = 80;
type P = [u8; NM];
type T = [P; 4];
type Key = Vec<u8>;

#[derive(Clone, Copy)]
struct Cx { n: usize }
impl Cx {
    fn id(&self) -> P { let mut p = [0u8; NM]; for i in 0..NM { p[i] = i as u8; } p }
    #[inline] fn mul(&self, a: &P, b: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[i] = a[b[i] as usize]; } c }
    #[inline] fn inv(&self, a: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[a[i] as usize] = i as u8; } c }
    fn is_id(&self, a: &P) -> bool { (0..self.n).all(|i| a[i] as usize == i) }
    #[inline] fn ncyc(&self, a: &P) -> usize {
        let mut seen = [false; NM]; let mut c = 0;
        for s in 0..self.n { if seen[s] { continue; } c += 1; let mut x = s; while !seen[x] { seen[x] = true; x = a[x] as usize; } }
        c
    }
    #[inline] fn ind(&self, a: &P) -> usize { self.n - self.ncyc(a) }
    fn cyc_lens(&self, a: &P) -> Vec<usize> {
        let mut seen = [false; NM]; let mut out = vec![];
        for s in 0..self.n { if seen[s] { continue; } let mut l = 0; let mut x = s; while !seen[x] { seen[x] = true; x = a[x] as usize; l += 1; } out.push(l); }
        out
    }
    fn orbits(&self, gens: &[&P]) -> ([u8; NM], usize) {
        let mut lab = [u8::MAX; NM]; let mut k = 0u8; let mut stack = [0u8; NM];
        for s in 0..self.n {
            if lab[s] != u8::MAX { continue; }
            let mut sp = 1; stack[0] = s as u8; lab[s] = k;
            while sp > 0 { sp -= 1; let x = stack[sp] as usize;
                for g in gens { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = k; stack[sp] = y as u8; sp += 1; } } }
            k += 1;
        }
        (lab, k as usize)
    }
    fn fmt(&self, p: &P) -> String {
        let mut seen = [false; NM]; let mut s = String::new();
        for i in 0..self.n { if seen[i] { continue; } let mut c = vec![]; let mut x = i;
            while !seen[x] { seen[x] = true; c.push(x + 1); x = p[x] as usize; }
            if c.len() > 1 { s.push('('); s.push_str(&c.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ")); s.push(')'); } }
        if s.is_empty() { "()".into() } else { s }
    }
    fn ctype(&self, p: &P) -> Vec<usize> { let mut c = self.cyc_lens(p); c.retain(|&l| l > 1); c.sort_unstable_by(|a, b| b.cmp(a)); c }
}

fn genus_of(cx: &Cx, t: &T) -> i64 { let s: usize = t.iter().map(|p| cx.ind(p)).sum(); (s as i64 - 2 * cx.n as i64 + 2) / 2 }

/// (number of positive-genus orbits, total genus) on one side of the {12|34} cusp.
fn side(cx: &Cx, x: &P, y: &P, p: &P, lab: &[u8; NM], no: usize) -> (u32, i32) {
    let mut c = vec![0i32; no];
    for i in 0..cx.n { c[lab[i] as usize] += 1; }
    for v in 0..no { c[v] += 2; }
    for a in [x, y, p] {
        let mut seen = [false; NM];
        for s in 0..cx.n { if seen[s] { continue; } c[lab[s] as usize] -= 1; let mut z = s; while !seen[z] { seen[z] = true; z = a[z] as usize; } }
    }
    let (mut pos, mut tot) = (0u32, 0i32);
    for v in 0..no { if c[v] > 0 { pos += 1; tot += c[v] / 2; } }
    (pos, tot)
}
fn one_piece_ref(cx: &Cx, t: &T) -> bool {
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
fn comb_smooth(cx: &Cx, t: &T) -> bool {
    let ip = cx.ind(&cx.mul(&t[0], &t[1]));
    ip == cx.ind(&t[0]) + cx.ind(&t[1]) || ip == cx.ind(&t[2]) + cx.ind(&t[3])
}
fn mv_ref(cx: &Cx, t: &T, m: usize) -> T {
    let j = m % 3; let a = t[j]; let b = t[j + 1]; let mut u = *t;
    if m < 3 { u[j] = cx.mul(&a, &cx.mul(&b, &cx.inv(&a))); u[j + 1] = a; } else { u[j] = b; u[j + 1] = cx.mul(&cx.inv(&b), &cx.mul(&a, &b)); }
    u
}
/// Fast move: a b a^-1 via c[a[i]] = a[b[i]] (no inverse); b^-1 a b via c[binv[k]] = binv[a[k]].
#[inline] fn mv(cx: &Cx, t: &T, m: usize) -> T {
    let n = cx.n; let j = m % 3; let mut u = *t;
    if m < 3 { let (a, b) = (&t[j], &t[j + 1]); let mut c = cx.id(); for i in 0..n { c[a[i] as usize] = a[b[i] as usize]; } u[j] = c; u[j + 1] = *a; }
    else { let (a, b) = (&t[j], &t[j + 1]); let mut bi = cx.id(); for i in 0..n { bi[b[i] as usize] = i as u8; }
        let mut c = cx.id(); for k in 0..n { c[bi[k] as usize] = bi[a[k] as usize]; } u[j] = *b; u[j + 1] = c; }
    u
}

// ---- allocation-free fast paths (checked against the _ref versions by --selftest)
#[inline] fn orbits2(n: usize, a: &P, b: &P, lab: &mut [u8; NM]) -> usize {
    for i in 0..n { lab[i] = u8::MAX; } let mut k = 0u8; let mut stack = [0u8; NM];
    for s in 0..n { if lab[s] != u8::MAX { continue; } let mut sp = 1; stack[0] = s as u8; lab[s] = k;
        while sp > 0 { sp -= 1; let x = stack[sp] as usize;
            let y = a[x] as usize; if lab[y] == u8::MAX { lab[y] = k; stack[sp] = y as u8; sp += 1; }
            let y = b[x] as usize; if lab[y] == u8::MAX { lab[y] = k; stack[sp] = y as u8; sp += 1; } }
        k += 1; }
    k as usize
}
#[inline] fn sub_cycles(n: usize, a: &P, lab: &[u8; NM], c: &mut [i32; NM]) -> usize {
    let mut seen: u128 = 0; let mut k = 0;
    for s in 0..n { if seen >> s & 1 == 1 { continue; } k += 1; c[lab[s] as usize] -= 1; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u128 << x; x = a[x] as usize; } }
    k
}
#[inline] fn ncyc_fast(n: usize, a: &P) -> usize {
    let mut seen: u128 = 0; let mut k = 0;
    for s in 0..n { if seen >> s & 1 == 1 { continue; } k += 1; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u128 << x; x = a[x] as usize; } }
    k
}
#[inline] fn side_pos(n: usize, x: &P, y: &P, p: &P, lab: &[u8; NM], no: usize) -> u32 {
    let mut c = [0i32; NM]; for i in 0..n { c[lab[i] as usize] += 1; } for v in 0..no { c[v] += 2; }
    sub_cycles(n, x, lab, &mut c); sub_cycles(n, y, lab, &mut c); sub_cycles(n, p, lab, &mut c);
    (0..no).filter(|&v| c[v] > 0).count() as u32
}
#[inline] fn one_piece(cx: &Cx, t: &T) -> bool {
    let n = cx.n; let mut p = [0u8; NM]; for i in 0..n { p[i] = t[0][t[1][i] as usize]; }
    let mut la = [0u8; NM]; let mut lb = [0u8; NM];
    let na = orbits2(n, &t[0], &t[1], &mut la); let nb = orbits2(n, &t[2], &t[3], &mut lb);
    if ncyc_fast(n, &p) + 1 != na + nb { return false; }
    let pa = side_pos(n, &t[0], &t[1], &p, &la, na); if pa > 1 { return false; }
    pa + side_pos(n, &t[2], &t[3], &p, &lb, nb) <= 1
}
/// Same output as canon_ref: least key over base points of the rarest (len g1..g4) signature.
fn canon(cx: &Cx, t: &T) -> Key {
    let n = cx.n;
    let mut ls = [[0u8; NM]; 4];
    for j in 0..4 { let a = &t[j]; let mut seen: u128 = 0;
        for s in 0..n { if seen >> s & 1 == 1 { continue; } let mut l = 0u8; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u128 << x; x = a[x] as usize; l += 1; }
            let mut x = s; for _ in 0..l { ls[j][x] = l; x = a[x] as usize; } } }
    let mut sig = [0u32; NM]; for r in 0..n { sig[r] = ((ls[0][r] as u32) << 24) | ((ls[1][r] as u32) << 16) | ((ls[2][r] as u32) << 8) | ls[3][r] as u32; }
    let mut srt = sig; srt[..n].sort_unstable();
    let (mut best_sig, mut best_cnt) = (0u32, usize::MAX); let mut i = 0;
    while i < n { let mut j = i; while j < n && srt[j] == srt[i] { j += 1; } if j - i < best_cnt { best_cnt = j - i; best_sig = srt[i]; } i = j; }
    let mut best = [u8::MAX; 3 * NM]; let mut have = false;
    for r in 0..n {
        if sig[r] != best_sig { continue; }
        let mut lab = [u8::MAX; NM]; let mut order = [0u8; NM]; let mut len = 1; lab[r] = 0; order[0] = r as u8; let mut q = 0;
        while q < len { let x = order[q] as usize; q += 1; for g in &t[..3] { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y as u8; len += 1; } } }
        // compare incrementally with early abort
        let mut state = if have { 0i8 } else { -1 }; let mut key = [0u8; 3 * NM];
        'o: for j in 0..3 { for k in 0..n { let v = lab[t[j][order[k] as usize] as usize]; key[j * n + k] = v;
            if state == 0 { let b = best[j * n + k]; if v < b { state = -1; } else if v > b { state = 1; break 'o; } } } }
        if state == -1 { best = key; have = true; }
    }
    best[..3 * n].to_vec()
}
fn rot(t: &T) -> T { [t[1], t[2], t[3], t[0]] }
/// one-piece at {12|34}, {23|41} and at {12|34} of sigma2 (the three cusps of a state)
fn three_cusps(cx: &Cx, t: &T) -> bool { one_piece(cx, t) && one_piece(cx, &rot(t)) && one_piece(cx, &mv(cx, t, 1)) }

fn canon_ref(cx: &Cx, t: &T) -> Key {
    let n = cx.n;
    let ls: Vec<Vec<usize>> = t.iter().map(|p| { let mut o = vec![0; n]; let mut seen = vec![false; n];
        for s in 0..n { if seen[s] { continue; } let mut c = vec![]; let mut x = s; while !seen[x] { seen[x] = true; c.push(x); x = p[x] as usize; } for &y in &c { o[y] = c.len(); } } o }).collect();
    let sig: Vec<u64> = (0..n).map(|r| ((ls[0][r] as u64) << 48) | ((ls[1][r] as u64) << 32) | ((ls[2][r] as u64) << 16) | ls[3][r] as u64).collect();
    let mut cnt: HashMap<u64, usize> = HashMap::new(); for &s in &sig { *cnt.entry(s).or_insert(0) += 1; }
    let best_sig = *cnt.iter().min_by_key(|(s, c)| (**c, **s)).unwrap().0;
    let mut best: Option<Key> = None;
    for r in 0..n {
        if sig[r] != best_sig { continue; }
        let mut lab = vec![u8::MAX; n]; let mut order = vec![r]; lab[r] = 0; let mut i = 0;
        while i < order.len() { let x = order[i]; i += 1; for g in &t[..3] { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = order.len() as u8; order.push(y); } } }
        let mut key = Vec::with_capacity(3 * n);
        for j in 0..3 { for k in 0..n { key.push(lab[t[j][order[k]] as usize]); } }
        if best.as_ref().map_or(true, |b| key < *b) { best = Some(key); }
    }
    best.unwrap()
}
fn unpack(cx: &Cx, k: &Key) -> T {
    let n = cx.n; let mut t = [cx.id(); 4];
    for j in 0..3 { for i in 0..n { t[j][i] = k[j * n + i]; } }
    t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2])));
    t
}

fn h128(k: &Key) -> u128 {
    let (mut a, mut b) = (0xcbf2_9ce4_8422_2325u64, 0x9e37_79b9_7f4a_7c15u64);
    for &x in k { a = (a ^ x as u64).wrapping_mul(0x0000_0100_0000_01b3); b = (b.rotate_left(7) ^ x as u64).wrapping_mul(0x517c_c1b7_2722_0a95); }
    ((a as u128) << 64) | (b as u128)
}
/// States visited by failed closures (their orbits are not one-piece everywhere).
struct Bad { sh: Vec<Mutex<HashSet<u128>>>, n: AtomicU64, max: u64 }
impl Bad {
    fn has(&self, h: u128) -> bool { self.sh[(h as usize) & 63].lock().unwrap().contains(&h) }
    fn add(&self, hs: &[u128]) { if self.n.load(Ordering::Relaxed) >= self.max { return; }
        for &h in hs { if self.sh[(h as usize) & 63].lock().unwrap().insert(h) { self.n.fetch_add(1, Ordering::Relaxed); } } }
}

struct Orbit { states: Vec<Key>, nbr: Vec<[u32; 6]>, all_one: bool, capped: bool }
fn close(cx: &Cx, t0: &T, cap: usize) -> Orbit {
    let root = canon(cx, t0);
    let mut idx: HashMap<Key, u32> = HashMap::new(); idx.insert(root.clone(), 0);
    let mut st = vec![root]; let mut nbr = vec![]; let mut i = 0;
    while i < st.len() {
        let u = unpack(cx, &st[i]); let mut row = [0u32; 6];
        for m in 0..6 {
            let v = mv(cx, &u, m);
            if !one_piece(cx, &v) { return Orbit { states: st, nbr, all_one: false, capped: false }; }
            let k = canon(cx, &v);
            let id = match idx.get(&k) { Some(&j) => j, None => { let j = st.len() as u32; idx.insert(k.clone(), j); st.push(k);
                if st.len() > cap { return Orbit { states: st, nbr, all_one: false, capped: true }; } j } };
            row[m] = id;
        }
        nbr.push(row); i += 1;
    }
    Orbit { states: st, nbr, all_one: true, capped: false }
}
fn components(nbr: &[[u32; 6]], moves: &[usize]) -> usize {
    let n = nbr.len(); let mut lab = vec![u32::MAX; n]; let mut k = 0;
    for s in 0..n { if lab[s] != u32::MAX { continue; } lab[s] = k; let mut stack = vec![s];
        while let Some(x) = stack.pop() { for &m in moves { let y = nbr[x][m] as usize; if lab[y] == u32::MAX { lab[y] = k; stack.push(y); } } } k += 1; }
    k as usize
}
fn gcd(a: i128, b: i128) -> i128 { if b == 0 { a.abs() } else { gcd(b, a % b) } }
fn block_sizes(cx: &Cx, t: &T) -> Vec<usize> {
    let n = cx.n; let mut out = vec![];
    for j in 1..n {
        let mut par: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, x: usize) -> usize { let mut r = x; while p[r] != r { r = p[r]; } let mut y = x; while p[y] != r { let z = p[y]; p[y] = r; y = z; } r }
        let mut queue = vec![(0usize, j)];
        while let Some((x, y)) = queue.pop() { let (rx, ry) = (find(&mut par, x), find(&mut par, y)); if rx == ry { continue; } par[rx] = ry; for g in t.iter() { queue.push((g[x] as usize, g[y] as usize)); } }
        let r0 = find(&mut par, 0); let sz = (0..n).filter(|&x| find(&mut par, x) == r0).count();
        if sz < n { out.push(sz); }
    }
    out.sort(); out.dedup(); out
}
/// JSON report: states, clusters, lambda (C18 formula, exact rational 12*deg(lambda) on the
/// N-state curve, as numerator/denominator), comb-smooth states, block sizes.
fn report(cx: &Cx, o: &Orbit) -> String {
    let ts: Vec<T> = o.states.iter().map(|k| unpack(cx, k)).collect();
    let t0 = ts[0];
    // lcm of all cycle lengths that enter the formula
    let mut ll: i128 = 1;
    let mut take = |p: &P, ll: &mut i128| { for l in cx.cyc_lens(p) { let l = l as i128; *ll = *ll / gcd(*ll, l) * l; } };
    for p in &t0 { take(p, &mut ll); }
    let prods: Vec<[P; 3]> = ts.iter().map(|u| [cx.mul(&u[0], &u[1]), cx.mul(&u[1], &u[2]), cx.mul(&u[0], &cx.mul(&u[1], &cx.mul(&u[2], &cx.inv(&u[1]))))]).collect();
    for q in &prods { for p in q { take(p, &mut ll); } }
    let h = |p: &P| -> i128 { cx.cyc_lens(p).iter().map(|&l| ll / l as i128).sum() };
    let base: i128 = ll * cx.n as i128 - t0.iter().map(|p| h(p)).sum::<i128>();
    let mut lam = base * ts.len() as i128;
    for q in &prods { lam += h(&q[0]) + h(&q[1]) + h(&q[2]); }
    let g = gcd(lam, ll); let (num, den) = (lam / g.max(1), ll / g.max(1));
    let comb = ts.iter().filter(|u| comb_smooth(cx, u)).count();
    let types: Vec<Vec<usize>> = t0.iter().map(|p| cx.ctype(p)).collect();
    let tuple: Vec<String> = t0.iter().map(|p| cx.fmt(p)).collect();
    format!("{{\"d\":{},\"genus\":{},\"states\":{},\"clusters\":{},\"all_one_piece\":true,\"comb_smooth_states\":{},\"all_smooth\":{},\"lambda12_num\":{},\"lambda12_den\":{},\"block_sizes\":{:?},\"types\":{:?},\"tuple\":{:?}}}",
        cx.n, genus_of(cx, &t0), ts.len(), components(&o.nbr, &[0, 2, 3, 5]), comb, comb == ts.len(), num, den, block_sizes(cx, &t0), types, tuple)
}

// ---------------------------------------------------------------- fibre group S_m

type F = [u8; 8];
fn fid() -> F { [0, 1, 2, 3, 4, 5, 6, 7] }
#[inline] fn fmul(m: usize, a: &F, b: &F) -> F { let mut c = fid(); for i in 0..m { c[i] = a[b[i] as usize]; } c }
#[inline] fn finv(m: usize, a: &F) -> F { let mut c = fid(); for i in 0..m { c[a[i] as usize] = i as u8; } c }
fn fctype(m: usize, a: &F) -> Vec<usize> {
    let mut seen = [false; 8]; let mut out = vec![];
    for s in 0..m { if seen[s] { continue; } let mut l = 0; let mut x = s; while !seen[x] { seen[x] = true; x = a[x] as usize; l += 1; } out.push(l); }
    out.sort_unstable_by(|a, b| b.cmp(a)); out
}
fn all_f(m: usize) -> Vec<F> {
    let mut out = vec![]; let mut a: Vec<u8> = (0..m as u8).collect();
    fn rec(k: usize, a: &mut Vec<u8>, out: &mut Vec<F>) { if k == a.len() { let mut p = fid(); for i in 0..a.len() { p[i] = a[i]; } out.push(p); return; }
        for i in k..a.len() { a.swap(k, i); rec(k + 1, a, out); a.swap(k, i); } }
    rec(0, &mut a, &mut out); out
}
fn partitions(n: usize, maxp: usize) -> Vec<Vec<usize>> {
    if n == 0 { return vec![vec![]]; }
    let mut out = vec![]; for p in (1..=maxp.min(n)).rev() { for mut r in partitions(n - p, p) { r.insert(0, p); out.push(r); } } out
}

// ---------------------------------------------------------------- words and Tietze

/// A word in variables, in APPLICATION order: [w1, w2, ...] evaluates to ... w2 w1.
type Word = Vec<(u16, bool)>; // (variable, inverted)
fn winv(w: &Word) -> Word { w.iter().rev().map(|&(v, i)| (v, !i)).collect() }
fn freduce(w: &Word) -> Word {
    let mut out: Word = vec![];
    for &x in w { if let Some(&l) = out.last() { if l.0 == x.0 && l.1 != x.1 { out.pop(); continue; } } out.push(x); }
    out
}
fn creduce(w: &Word) -> Word {
    let mut w = freduce(w);
    while w.len() >= 2 { let (a, b) = (w[0], w[w.len() - 1]); if a.0 == b.0 && a.1 != b.1 { w.remove(0); w.pop(); } else { break; } }
    w
}
fn subst(w: &Word, var: u16, e: &Word) -> Word {
    let ei = winv(e); let mut out = vec![];
    for &(v, i) in w { if v == var { out.extend_from_slice(if i { &ei } else { e }); } else { out.push((v, i)); } }
    freduce(&out)
}
fn eval(m: usize, w: &Word, val: &[F]) -> F {
    let mut r = fid();
    for &(v, i) in w { let x = if i { finv(m, &val[v as usize]) } else { val[v as usize] }; r = fmul(m, &x, &r); }
    r
}

struct Base { cx: Cx, a: T, b: usize, h: i64 }

/// Voltage variable layout: var(i, x) = i*b + x for entries i = 0,1,2.
struct Setup {
    nvar: usize,          // 3b voltage variables (tree ones fixed to 1)
    tree: Vec<bool>,
    cycles: Vec<(usize, Vec<usize>)>, // (entry 0..3, points in order)
    words: Vec<Word>,     // return word of each cycle, tree variables removed
}
fn setup(base: &Base) -> Setup {
    let b = base.b; let a = &base.a;
    let mut tree = vec![false; 3 * b];
    let mut seen = vec![false; b]; seen[0] = true; let mut q = vec![0usize]; let mut qi = 0;
    while qi < q.len() { let y = q[qi]; qi += 1;
        for i in 0..3 { let z = a[i][y] as usize; if !seen[z] { seen[z] = true; tree[i * b + y] = true; q.push(z); }
            let w = (0..b).find(|&w| a[i][w] as usize == y).unwrap(); if !seen[w] { seen[w] = true; tree[i * b + w] = true; q.push(w); } } }
    let mut cycles = vec![]; let mut words = vec![];
    for i in 0..4 {
        let mut s = vec![false; b];
        for x0 in 0..b { if s[x0] { continue; } let mut pts = vec![]; let mut x = x0; while !s[x] { s[x] = true; pts.push(x); x = a[i][x] as usize; }
            let mut w: Word = vec![];
            for &x in &pts {
                if i < 3 { w.push(((i * b + x) as u16, false)); }
                else { // V4[x] = V3[y]^-1 V2[a3 y]^-1 V1[a2 a3 y]^-1, y = a4 x; application order: V1^-1 first
                    let y = a[3][x] as usize; let y3 = a[2][y] as usize; let y23 = a[1][y3] as usize;
                    w.push(((0 * b + y23) as u16, true)); w.push(((1 * b + y3) as u16, true)); w.push(((2 * b + y) as u16, true));
                }
            }
            let w: Word = w.into_iter().filter(|&(v, _)| !tree[v as usize]).collect();
            cycles.push((i, pts)); words.push(freduce(&w));
        }
    }
    Setup { nvar: 3 * b, tree, cycles, words }
}

/// A designated pattern: (cycle index, cycle type in S_m) for each nontrivial return.
type Pattern = Vec<(usize, Vec<usize>)>;
fn patterns(ncyc: usize, m: usize, excess: usize) -> Vec<Pattern> {
    let mut by_e: Vec<Vec<Vec<usize>>> = vec![vec![]; excess + 1];
    for p in partitions(m, m) { let e = m - p.len(); if e >= 1 && e <= excess { by_e[e].push(p); } }
    let mut out = vec![];
    fn rec(c: usize, left: usize, ncyc: usize, by_e: &Vec<Vec<Vec<usize>>>, cur: &mut Pattern, out: &mut Vec<Pattern>) {
        if left == 0 { out.push(cur.clone()); return; }
        if c == ncyc { return; }
        rec(c + 1, left, ncyc, by_e, cur, out);
        for e in 1..=left { for p in &by_e[e] { cur.push((c, p.clone())); rec(c + 1, left - e, ncyc, by_e, cur, out); cur.pop(); } }
    }
    rec(0, excess, ncyc, &by_e, &mut vec![], &mut out); out
}

/// Result of Tietze elimination for one pattern.
struct Plan {
    free: Vec<u16>,                   // free S_m variables (in enumeration order)
    classes: Vec<(u16, Vec<usize>)>,  // class variables and their cycle types
    exprs: Vec<(u16, Word)>,          // eliminated voltage variables, in terms of free + class vars
    checks: Vec<(Word, Option<Vec<usize>>)>, // remaining relators: trivial (None) or class
}
fn tietze(st: &Setup, pat: &Pattern) -> Option<Plan> {
    let nv = st.nvar as u16; let mut next = nv;
    let des: HashMap<usize, Vec<usize>> = pat.iter().cloned().collect();
    let mut rels: Vec<(Word, Option<Vec<usize>>)> = st.words.iter().enumerate().map(|(c, w)| (creduce(w), des.get(&c).cloned())).collect();
    // empty designated relator => identity return, impossible
    if rels.iter().any(|(w, d)| w.is_empty() && d.is_some()) { return None; }
    rels.retain(|(w, _)| !w.is_empty());
    let mut exprs: Vec<(u16, Word)> = vec![]; let mut classes = vec![];
    loop {
        // pick a relator with a voltage variable that occurs exactly once; trivial first, shortest first
        let mut pick: Option<(usize, usize)> = None; let mut best = (2u8, usize::MAX);
        for (ri, (w, d)) in rels.iter().enumerate() {
            let mut cnt: HashMap<u16, usize> = HashMap::new(); for &(v, _) in w { *cnt.entry(v).or_insert(0) += 1; }
            if let Some(pos) = w.iter().position(|&(v, _)| v < nv && cnt[&v] == 1) {
                let key = (if d.is_some() { 1 } else { 0 }, w.len());
                if key < best { best = key; pick = Some((ri, pos)); }
            }
        }
        let Some((ri, pos)) = pick else { break };
        let (w, d) = rels.remove(ri);
        // rotate so that the variable is applied first: w = [u^eps, rest...]
        let mut r = w[pos..].to_vec(); r.extend_from_slice(&w[..pos]);
        let (u, eps_inv) = r[0]; let rest: Word = r[1..].to_vec();
        // rest o u^eps = target  =>  u^eps = rest^-1 o target  (application order: target, then rest^-1)
        let mut ue: Word = vec![];
        if let Some(ct) = d { let s = next; next += 1; classes.push((s, ct)); ue.push((s, false)); }
        ue.extend(winv(&rest));
        let e = if eps_inv { winv(&ue) } else { ue };
        let e = freduce(&e);
        for x in exprs.iter_mut() { x.1 = subst(&x.1, u, &e); }
        for x in rels.iter_mut() { x.0 = creduce(&subst(&x.0, u, &e)); }
        exprs.push((u, e));
        if rels.iter().any(|(w, d)| w.is_empty() && d.is_some()) { return None; }
        rels.retain(|(w, _)| !w.is_empty());
    }
    let elim: HashSet<u16> = exprs.iter().map(|x| x.0).collect();
    let free: Vec<u16> = (0..nv).filter(|&v| !st.tree[v as usize] && !elim.contains(&v)).collect();
    Some(Plan { free, classes, exprs, checks: rels })
}

// ---------------------------------------------------------------- lifting

fn lift(base: &Base, m: usize, volt: &[F]) -> (Cx, T) {
    let b = base.b; let cx = Cx { n: m * b }; let mut t = [cx.id(); 4];
    for i in 0..3 { for x in 0..b { let v = &volt[i * b + x]; let y = base.a[i][x] as usize; for j in 0..m { t[i][x * m + j] = (y * m + v[j] as usize) as u8; } } }
    t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2])));
    (cx, t)
}

struct Stats { cand: AtomicU64, checks_ok: AtomicU64, genus_ok: AtomicU64, trans: AtomicU64, pre: AtomicU64, closures: AtomicU64, found: AtomicU64, capped: AtomicU64, skipped_known: AtomicU64, skipped_bad: AtomicU64, failed_states: AtomicU64, fail_hist: [AtomicU64; 8], probe_rejects: AtomicU64 }
struct Ctx<'a> { base: &'a Base, m: usize, g: i64, cap: usize, no_solve: bool, probe: usize, stats: Stats, known: Mutex<HashSet<Key>>, out: Mutex<std::fs::File>, seeds: Option<Mutex<HashSet<Key>>>, bad: Bad }

fn process(ctx: &Ctx, volt: &[F], pat: Option<&Pattern>, st: &Setup) {
    ctx.stats.checks_ok.fetch_add(1, Ordering::Relaxed);
    let (cx, t) = lift(ctx.base, ctx.m, volt);
    let nth = ctx.stats.checks_ok.load(Ordering::Relaxed);
    // In Tietze mode the genus is forced by the pattern; it and the return types are asserted on
    // every 64th candidate (and on every candidate in brute mode, where genus is a filter).
    let sample = pat.is_none() || nth % 64 == 0;
    if sample && genus_of(&cx, &t) != ctx.g { assert!(pat.is_none(), "pattern lift has wrong genus"); return; }
    ctx.stats.genus_ok.fetch_add(1, Ordering::Relaxed);
    if cx.orbits(&[&t[0], &t[1], &t[2]]).1 != 1 { return; }
    ctx.stats.trans.fetch_add(1, Ordering::Relaxed);
    // independent check of the pattern: cycle types of the returns, read off the permutations
    if let (Some(pat), true) = (pat, sample) {
        let des: HashMap<usize, &Vec<usize>> = pat.iter().map(|(c, p)| (*c, p)).collect();
        for (ci, (i, pts)) in st.cycles.iter().enumerate() {
            let l = pts.len(); let x = pts[0];
            let mut r = fid(); for j in 0..ctx.m { let mut z = x * ctx.m + j; for _ in 0..l { z = t[*i][z] as usize; } assert_eq!(z / ctx.m, x); r[j] = (z % ctx.m) as u8; }
            let ct = fctype(ctx.m, &r);
            match des.get(&ci) { Some(p) => assert_eq!(&&ct, p, "designated return has wrong type"), None => assert!(ct.iter().all(|&c| c == 1), "return should be trivial") }
        }
    }
    if !three_cusps(&cx, &t) { return; }
    ctx.stats.pre.fetch_add(1, Ordering::Relaxed);
    let k = canon(&cx, &t);
    if let Some(s) = &ctx.seeds { s.lock().unwrap().insert(k.clone()); }
    if ctx.known.lock().unwrap().contains(&k) { ctx.stats.skipped_known.fetch_add(1, Ordering::Relaxed); return; }
    let hk = h128(&k);
    if ctx.bad.has(hk) { ctx.stats.skipped_bad.fetch_add(1, Ordering::Relaxed); return; }
    // probe: a random walk of `probe` moves; any state that is not one-piece rejects the orbit
    if ctx.probe > 0 {
        let mut r = (hk as u64) | 1; let mut u = t;
        for _ in 0..ctx.probe { r ^= r << 13; r ^= r >> 7; r ^= r << 17; u = mv(&cx, &u, (r % 6) as usize);
            if !one_piece(&cx, &u) { ctx.stats.probe_rejects.fetch_add(1, Ordering::Relaxed); ctx.bad.add(&[hk]); return; } }
    }
    ctx.stats.closures.fetch_add(1, Ordering::Relaxed);
    let o = close(&cx, &t, ctx.cap);
    if !o.all_one && !o.capped { let hs: Vec<u128> = o.states.iter().map(h128).collect(); ctx.bad.add(&hs);
        ctx.stats.failed_states.fetch_add(o.states.len() as u64, Ordering::Relaxed);
        let b = match o.states.len() { 1 => 0, 2..=3 => 1, 4..=10 => 2, 11..=30 => 3, 31..=100 => 4, 101..=1000 => 5, 1001..=10000 => 6, _ => 7 }; ctx.stats.fail_hist[b].fetch_add(1, Ordering::Relaxed); }
    if o.capped { ctx.stats.capped.fetch_add(1, Ordering::Relaxed);
        let tu: Vec<String> = t.iter().map(|p| cx.fmt(p)).collect();
        writeln!(ctx.out.lock().unwrap(), "{{\"capped\":true,\"d\":{},\"tuple\":{:?}}}", cx.n, tu).unwrap(); return; }
    if !o.all_one { return; }
    {
        let mut kn = ctx.known.lock().unwrap();
        if kn.contains(&k) { return; }
        for s in &o.states { kn.insert(s.clone()); }
    }
    ctx.stats.found.fetch_add(1, Ordering::Relaxed);
    let line = report(&cx, &o);
    writeln!(ctx.out.lock().unwrap(), "{}", line).unwrap();
    eprintln!("FOUND {}", line);
}

/// Conjugator y with y p y^-1 = q (same cycle type assumed): match cycles of equal length.
fn conjugator(m: usize, p: &F, q: &F) -> Option<F> {
    let cyc = |a: &F| { let mut seen = [false; 8]; let mut out: Vec<Vec<u8>> = vec![];
        for s in 0..m { if seen[s] { continue; } let mut c = vec![]; let mut x = s; while !seen[x] { seen[x] = true; c.push(x as u8); x = a[x] as usize; } out.push(c); }
        out.sort_by_key(|c| c.len()); out };
    let (cp, cq) = (cyc(p), cyc(q));
    if cp.len() != cq.len() || cp.iter().zip(&cq).any(|(a, b)| a.len() != b.len()) { return None; }
    let mut y = fid(); for (a, b) in cp.iter().zip(&cq) { for i in 0..a.len() { y[a[i] as usize] = b[i]; } }
    Some(y)
}

/// If the single remaining relator contains a free variable X exactly twice with opposite
/// signs, W = [X^e, U.., X^-e, V..] (application order) means Y^-1 U Y = V^-1 with Y = X^e,
/// so Y runs over y0 C(V^-1) where y0 V^-1 y0^-1 = U. Returns (X, e_inverted, U, V).
fn solvable(plan: &Plan) -> Option<(u16, bool, Word, Word)> {
    if plan.checks.len() != 1 { return None; }
    let w = &plan.checks[0].0;
    for &x in &plan.free {
        let pos: Vec<usize> = (0..w.len()).filter(|&i| w[i].0 == x).collect();
        if pos.len() != 2 || w[pos[0]].1 == w[pos[1]].1 { continue; }
        let mut r = w[pos[0]..].to_vec(); r.extend_from_slice(&w[..pos[0]]);
        let j = pos[1] - pos[0];
        return Some((x, r[0].1, r[1..j].to_vec(), r[j + 1..].to_vec()));
    }
    None
}

fn run_pattern(ctx: &Ctx, st: &Setup, pat: &Pattern, plan: &Plan, all: &[F], reps: &[F], by_type: &HashMap<Vec<usize>, Vec<F>>, cent: &HashMap<Vec<usize>, (F, Vec<F>)>) {
    let m = ctx.m;
    if !ctx.no_solve { if let Some((x, xinv, uw, vw)) = solvable(plan) {
        let others: Vec<u16> = plan.free.iter().cloned().filter(|&v| v != x).collect();
        let doms: Vec<&[F]> = (0..others.len()).map(|i| if i == 0 { reps } else { all }).chain(plan.classes.iter().map(|(_, ct)| by_type[ct].as_slice())).collect();
        let vars: Vec<u16> = others.iter().cloned().chain(plan.classes.iter().map(|c| c.0)).collect();
        if doms.iter().any(|d| d.is_empty()) { return; }
        let mut val = vec![fid(); 1 + st.nvar + plan.classes.len() + 8];
        let mut ix = vec![0usize; vars.len()];
        loop {
            for (k, &v) in vars.iter().enumerate() { val[v as usize] = doms[k][ix[k]]; }
            let u = eval(m, &uw, &val); let vinv = finv(m, &eval(m, &vw, &val));
            // eval(check) = V Y^-1 U Y = t (t = 1, or t in the prescribed class) => Y^-1 U Y = V^-1 t
            let one = [fid()];
            let targets: &[F] = match &plan.checks[0].1 { None => &one, Some(ct) => by_type[ct].as_slice() };
            for t in targets { let vi = fmul(m, &vinv, t);
            if let Some(y0) = conjugator(m, &vi, &u) {
                // C(vi) = z C(r) z^-1 where z r z^-1 = vi
                let (r, cr) = &cent[&fctype(m, &vi)]; let z = conjugator(m, r, &vi).unwrap(); let zi = finv(m, &z);
                for c in cr {
                    let cc = fmul(m, &z, &fmul(m, c, &zi)); let y = fmul(m, &y0, &cc);
                    val[x as usize] = if xinv { finv(m, &y) } else { y };
                    ctx.stats.cand.fetch_add(1, Ordering::Relaxed);
                    debug_assert!({ let r = eval(m, &plan.checks[0].0, &val); match &plan.checks[0].1 { None => (0..m).all(|j| r[j] as usize == j), Some(ct) => &fctype(m, &r) == ct } });
                    for (uu, e) in &plan.exprs { val[*uu as usize] = eval(m, e, &val); }
                    process(ctx, &val[..st.nvar], Some(pat), st);
                }
            } }
            let mut k = 0; loop { if k == vars.len() { return; } ix[k] += 1; if ix[k] < doms[k].len() { break; } ix[k] = 0; k += 1; }
        }
    } }
    let nfree = plan.free.len(); let ncl = plan.classes.len();
    let total = 1 + st.nvar + ncl + 8;
    let mut val = vec![fid(); total];
    let doms: Vec<&[F]> = (0..nfree).map(|i| if i == 0 { reps } else { all }).chain(plan.classes.iter().map(|(_, ct)| by_type[ct].as_slice())).collect();
    let vars: Vec<u16> = plan.free.iter().cloned().chain(plan.classes.iter().map(|c| c.0)).collect();
    let mut ix = vec![0usize; vars.len()];
    if doms.iter().any(|d| d.is_empty()) { return; }
    loop {
        for (k, &v) in vars.iter().enumerate() { val[v as usize] = doms[k][ix[k]]; }
        ctx.stats.cand.fetch_add(1, Ordering::Relaxed);
        let ok = plan.checks.iter().all(|(w, d)| { let r = eval(m, w, &val); match d { None => (0..m).all(|j| r[j] as usize == j), Some(ct) => &fctype(m, &r) == ct } });
        if ok {
            for (u, e) in &plan.exprs { val[*u as usize] = eval(m, e, &val); }
            process(ctx, &val[..st.nvar], Some(pat), st);
        }
        // odometer
        let mut k = 0; loop { if k == vars.len() { return; } ix[k] += 1; if ix[k] < doms[k].len() { break; } ix[k] = 0; k += 1; }
        if vars.is_empty() { return; }
    }
}

/// --dedupe FILE...: orbit id = least canonical state over the whole (re)closed orbit.
/// Prints per-file and global distinct counts; within-file duplicates would be a bug.
fn dedupe(files: &[String]) {
    let mut global: HashMap<u128, String> = HashMap::new(); let mut total = 0usize; let mut cross = 0usize;
    for f in files {
        let txt = std::fs::read_to_string(f).unwrap(); let mut local: HashSet<u128> = HashSet::new(); let mut n = 0usize; let mut dup_in = 0usize; let mut dup_x = 0usize;
        for line in txt.lines() {
            let Some(i) = line.find("\"tuple\":[") else { continue };
            if !line.contains("\"states\"") { continue; }
            let body = &line[i + 9..line.rfind(']').unwrap()];
            let perms: Vec<String> = body.split("\", \"").map(|x| x.trim_matches(|c| c == '"' || c == ' ').to_string()).collect();
            assert_eq!(perms.len(), 4, "{}", body);
            let nn = perms.iter().flat_map(|p| p.split(|c: char| !c.is_ascii_digit())).filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap()).max().unwrap();
            let cx = Cx { n: nn }; let mut t = [cx.id(); 4];
            for j in 0..4 { for cyc in perms[j].split(')') { let cyc = cyc.trim().trim_start_matches('('); if cyc.is_empty() { continue; }
                let v: Vec<usize> = cyc.split(' ').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap() - 1).collect();
                for k in 0..v.len() { t[j][v[k]] = v[(k + 1) % v.len()] as u8; } } }
            // the true degree may exceed the largest moved letter only if some letter is fixed by all four: impossible (transitive)
            assert!(cx.is_id(&cx.mul(&t[0], &cx.mul(&t[1], &cx.mul(&t[2], &t[3])))));
            let o = close(&cx, &t, 50_000_000); assert!(o.all_one && !o.capped, "orbit not one-piece everywhere on reclosure: {}", line);
            let mn = o.states.iter().min().unwrap().clone(); let mut key = vec![nn as u8]; key.extend(mn); let h = h128(&key);
            n += 1;
            if !local.insert(h) { dup_in += 1; }
            match global.get(&h) { Some(g) if g != f => { dup_x += 1; } None => { global.insert(h, f.clone()); } _ => {} }
        }
        total += n; cross += dup_x;
        println!("{}\torbits {}\tdistinct_in_file {}\tdup_within_file {}\tseen_in_earlier_file {}", f, n, local.len(), dup_in, dup_x);
    }
    println!("TOTAL records {} distinct orbits {} cross-file duplicates {}", total, global.len(), cross);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 && args[1] == "--dedupe" { dedupe(&args[2..]); return; }
    if args.len() > 1 && args[1] == "--selftest" {
        // fast one_piece / canon / mv vs the _ref versions: random S_n tuples, random wreath
        // tuples, and every state of the first orbits recorded in the given files
        let mut rng: u64 = 0x9e3779b97f4a7c15; let mut rnd = |k: usize| { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; (rng % k as u64) as usize };
        let mut checked = [0u64; 4]; // tuples, one_piece true, canon, mv
        let mut test = |cx: &Cx, t: &T, checked: &mut [u64; 4]| {
            if cx.orbits(&[&t[0], &t[1], &t[2]]).1 != 1 { return; }
            let a = one_piece(cx, t); assert_eq!(a, one_piece_ref(cx, t), "one_piece mismatch"); checked[1] += a as u64;
            assert_eq!(canon(cx, t), canon_ref(cx, t), "canon mismatch"); checked[2] += 1;
            for m in 0..6 { assert_eq!(mv(cx, t, m), mv_ref(cx, t, m), "mv mismatch"); } checked[3] += 6; checked[0] += 1;
        };
        for _ in 0..20000 { // random S_n
            let n = 4 + rnd(60); let cx = Cx { n }; let mut t = [cx.id(); 4];
            for j in 0..3 { for i in (1..n).rev() { let k = rnd(i + 1); t[j].swap(i, k); } }
            t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2]))); test(&cx, &t, &mut checked);
        }
        for _ in 0..200000 { // random wreath S_m wr S_b (many one-piece cases)
            let m = 2 + rnd(4); let b = 3 + rnd(6); let n = m * b; let cx = Cx { n }; let mut t = [cx.id(); 4];
            for j in 0..3 { let mut top: Vec<usize> = (0..b).collect(); for i in (1..b).rev() { let k = rnd(i + 1); top.swap(i, k); }
                for x in 0..b { let mut f: Vec<usize> = (0..m).collect(); for i in (1..m).rev() { let k = rnd(i + 1); f.swap(i, k); } for q in 0..m { t[j][x * m + q] = (top[x] * m + f[q]) as u8; } } }
            t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2]))); test(&cx, &t, &mut checked);
        }
        for f in &args[2..] { // all states of the first 200 orbits of each file
            let txt = std::fs::read_to_string(f).unwrap(); let mut k = 0;
            for line in txt.lines() { if !line.contains("\"states\"") { continue; } k += 1; if k > 200 { break; }
                let i = line.find("\"tuple\":[").unwrap(); let body = &line[i + 9..line.rfind(']').unwrap()];
                let perms: Vec<String> = body.split("\", \"").map(|x| x.trim_matches(|c| c == '"' || c == ' ').to_string()).collect();
                let nn = perms.iter().flat_map(|p| p.split(|c: char| !c.is_ascii_digit())).filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap()).max().unwrap();
                let cx = Cx { n: nn }; let mut t = [cx.id(); 4];
                for j in 0..4 { for cyc in perms[j].split(')') { let cyc = cyc.trim().trim_start_matches('('); if cyc.is_empty() { continue; }
                    let v: Vec<usize> = cyc.split(' ').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap() - 1).collect();
                    for q in 0..v.len() { t[j][v[q]] = v[(q + 1) % v.len()] as u8; } } }
                let o = close(&cx, &t, 10_000_000); assert!(o.all_one);
                for st in &o.states { let u = unpack(&cx, st); test(&cx, &u, &mut checked); }
            }
        }
        println!("selftest OK: {} tuples ({} one-piece), {} canon, {} moves compared with reference code", checked[0], checked[1], checked[2], checked[3]);
        return;
    }
    if args.len() > 1 && args[1] == "--bench" {
        // --bench FILE: time canon, three_cusps, one_piece, mv on the first recorded tuple
        let txt = std::fs::read_to_string(&args[2]).unwrap();
        let line = txt.lines().find(|l| l.contains("\"states\"")).unwrap();
        let i = line.find("\"tuple\":[").unwrap(); let body = &line[i + 9..line.rfind(']').unwrap()];
        let perms: Vec<String> = body.split("\", \"").map(|x| x.trim_matches(|c| c == '"' || c == ' ').to_string()).collect();
        let nn = perms.iter().flat_map(|p| p.split(|c: char| !c.is_ascii_digit())).filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap()).max().unwrap();
        let cx = Cx { n: nn }; let mut t = [cx.id(); 4];
        for j in 0..4 { for cyc in perms[j].split(')') { let cyc = cyc.trim().trim_start_matches('('); if cyc.is_empty() { continue; }
            let v: Vec<usize> = cyc.split(' ').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap() - 1).collect();
            for k in 0..v.len() { t[j][v[k]] = v[(k + 1) % v.len()] as u8; } } }
        let reps = 200000; let mut acc = 0usize;
        let st = std::time::Instant::now(); for k in 0..reps { let u = mv(&cx, &t, k % 6); acc += canon(&cx, &u)[0] as usize; } let tc = st.elapsed().as_nanos() as f64 / reps as f64;
        let st = std::time::Instant::now(); for k in 0..reps { let u = mv(&cx, &t, k % 6); acc += three_cusps(&cx, &u) as usize; } let t3 = st.elapsed().as_nanos() as f64 / reps as f64;
        let st = std::time::Instant::now(); for k in 0..reps { let u = mv(&cx, &t, k % 6); acc += one_piece(&cx, &u) as usize; } let t1 = st.elapsed().as_nanos() as f64 / reps as f64;
        let st = std::time::Instant::now(); for k in 0..reps { let u = mv(&cx, &t, k % 6); acc += u[0][0] as usize; } let tm = st.elapsed().as_nanos() as f64 / reps as f64;
        println!("n={} ns per call: mv {:.0}, mv+one_piece {:.0}, mv+three_cusps {:.0}, mv+canon {:.0} (acc {})", nn, tm, t1, t3, tc, acc);
        return;
    }
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let has = |k: &str| args.iter().any(|a| a == k);
    let bs = get("--base").expect("--base \"a;b;c;d\"");
    let m: usize = get("--m").expect("--m").parse().unwrap();
    let g: i64 = get("--genus").expect("--genus").parse().unwrap();
    let cap: usize = get("--cap").map(|x| x.parse().unwrap()).unwrap_or(2_000_000);
    let threads: usize = get("--threads").map(|x| x.parse().unwrap()).unwrap_or(1);
    let out = get("--out").unwrap_or("/dev/stdout".into());
    assert!(m >= 2 && m <= 8);
    // base
    let parts: Vec<&str> = bs.split(';').collect(); assert_eq!(parts.len(), 4);
    let b = bs.split(|c: char| !c.is_ascii_digit()).filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap()).max().unwrap();
    assert!(m * b <= NM);
    let bcx = Cx { n: b }; let mut a = [bcx.id(); 4];
    for j in 0..4 { for cyc in parts[j].split(')') { let cyc = cyc.trim().trim_start_matches('('); if cyc.is_empty() { continue; }
        let v: Vec<usize> = cyc.split(|c: char| c == ' ' || c == ',').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().unwrap() - 1).collect();
        for i in 0..v.len() { a[j][v[i]] = v[(i + 1) % v.len()] as u8; } } }
    assert!(bcx.is_id(&bcx.mul(&a[0], &bcx.mul(&a[1], &bcx.mul(&a[2], &a[3])))), "base product is not 1 (rightmost-first)");
    let h = genus_of(&bcx, &a);
    let base = Base { cx: bcx, a, b, h };
    let excess = 2 * g - 2 - m as i64 * (2 * h - 2);
    assert!(excess >= 0, "target genus too small");
    let st = setup(&base);
    let all = all_f(m);
    let mut reps: Vec<F> = vec![]; let mut by_type: HashMap<Vec<usize>, Vec<F>> = HashMap::new();
    for p in &all { by_type.entry(fctype(m, p)).or_default().push(*p); }
    for (_, v) in by_type.iter() { reps.push(v[0]); }
    reps.sort();
    let mut cent: HashMap<Vec<usize>, (F, Vec<F>)> = HashMap::new();
    for r in &reps { let c: Vec<F> = all.iter().filter(|c| fmul(m, c, r) == fmul(m, r, c)).cloned().collect(); cent.insert(fctype(m, r), (*r, c)); }
    let file = std::fs::File::create(&out).unwrap();
    let ctx = Ctx { base: &base, m, g, cap, no_solve: has("--no-solve"), probe: get("--probe").map(|x| x.parse().unwrap()).unwrap_or(64), stats: Stats { cand: AtomicU64::new(0), checks_ok: AtomicU64::new(0), genus_ok: AtomicU64::new(0), trans: AtomicU64::new(0), pre: AtomicU64::new(0), closures: AtomicU64::new(0), found: AtomicU64::new(0), capped: AtomicU64::new(0), skipped_known: AtomicU64::new(0), skipped_bad: AtomicU64::new(0), failed_states: AtomicU64::new(0), fail_hist: Default::default(), probe_rejects: AtomicU64::new(0) },
        known: Mutex::new(HashSet::new()), out: Mutex::new(file), seeds: if has("--seeds") { Some(Mutex::new(HashSet::new())) } else { None },
        bad: Bad { sh: (0..64).map(|_| Mutex::new(HashSet::new())).collect(), n: AtomicU64::new(0), max: get("--bad-max").map(|x| x.parse().unwrap()).unwrap_or(100_000_000) } };
    let t0 = std::time::Instant::now();
    let mode;
    let mut work_est: f64 = 0.0; let mut npat = 0usize; let mut nsolv = 0usize; let mut maxfree = 0usize; let mut nchecks_max = 0usize;
    if has("--brute") {
        mode = "brute";
        // every non-tree voltage assignment; free variable 0 up to conjugacy
        let free: Vec<u16> = (0..st.nvar as u16).filter(|&v| !st.tree[v as usize]).collect();
        work_est = reps.len() as f64 * (all.len() as f64).powi(free.len() as i32 - 1);
        eprintln!("brute: {} free voltages, {:.3e} assignments", free.len(), work_est);
        let next = AtomicUsize::new(0);
        std::thread::scope(|s| { for _ in 0..threads { s.spawn(|| {
            let mut val = vec![fid(); st.nvar];
            loop {
                let job = next.fetch_add(1, Ordering::Relaxed); // job = (rep of var0, value of var1)
                let n1 = if free.len() > 1 { all.len() } else { 1 };
                if job >= reps.len() * n1 { break; }
                val[free[0] as usize] = reps[job / n1]; if free.len() > 1 { val[free[1] as usize] = all[job % n1]; }
                let rest = &free[2.min(free.len())..]; let mut ix = vec![0usize; rest.len()];
                loop {
                    for (k, &v) in rest.iter().enumerate() { val[v as usize] = all[ix[k]]; }
                    ctx.stats.cand.fetch_add(1, Ordering::Relaxed);
                    process(&ctx, &val, None, &st);
                    let mut k = 0; let mut done = false; loop { if k == rest.len() { done = true; break; } ix[k] += 1; if ix[k] < all.len() { break; } ix[k] = 0; k += 1; }
                    if done { break; }
                }
            } }); } });
    } else {
        mode = "tietze";
        let pats = patterns(st.cycles.len(), m, excess as usize);
        let mut plans = vec![];
        for p in &pats { if let Some(pl) = tietze(&st, p) {
            let mut w = reps.len() as f64 * (all.len() as f64).powi(pl.free.len().max(1) as i32 - 1);
            if pl.free.is_empty() { w = 1.0; }
            for (_, ct) in &pl.classes { w *= by_type[ct].len() as f64; }
            if !has("--no-solve") { if solvable(&pl).is_some() { nsolv += 1; } }
            work_est += w; maxfree = maxfree.max(pl.free.len()); nchecks_max = nchecks_max.max(pl.checks.len());
            plans.push((p.clone(), pl)); } }
        npat = plans.len();
        eprintln!("tietze: {} patterns ({} feasible, {} solved by conjugator), max free S_m vars {}, max checks {}, {:.3e} assignments before solving", pats.len(), npat, nsolv, maxfree, nchecks_max, work_est);
        if has("--dump") { for (p, pl) in plans.iter().take(6) {
            let f = |w: &Word| w.iter().map(|(v, i)| format!("{}{}", v, if *i { "'" } else { "" })).collect::<Vec<_>>().join(" ");
            eprintln!("pattern {:?}\n  free {:?} classes {:?}\n  checks {:?}", p, pl.free, pl.classes, pl.checks.iter().map(|(w, d)| (f(w), d.clone())).collect::<Vec<_>>());
            for (u, e) in &pl.exprs { eprintln!("    {} = {}", u, f(e)); } } }
        if has("--plan") { return; }
        let next = AtomicUsize::new(0);
        std::thread::scope(|s| { for _ in 0..threads { s.spawn(|| { loop {
            let j = next.fetch_add(1, Ordering::Relaxed); if j >= plans.len() { break; }
            run_pattern(&ctx, &st, &plans[j].0, &plans[j].1, &all, &reps, &by_type, &cent);
            if j % 200 == 0 { eprintln!("pattern {}/{} {:.0}s cand {} pre {} closures {} found {}", j, plans.len(), t0.elapsed().as_secs_f64(), ctx.stats.cand.load(Ordering::Relaxed), ctx.stats.pre.load(Ordering::Relaxed), ctx.stats.closures.load(Ordering::Relaxed), ctx.stats.found.load(Ordering::Relaxed)); }
        } }); } });
    }
    if let (Some(s), Some(path)) = (&ctx.seeds, get("--seeds")) {
        let s = s.lock().unwrap(); let mut v: Vec<&Key> = s.iter().collect(); v.sort();
        let mut f = std::fs::File::create(path).unwrap(); for k in v { writeln!(f, "{}", k.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")).unwrap(); }
    }
    let s = &ctx.stats; let ld = |a: &AtomicU64| a.load(Ordering::Relaxed);
    writeln!(ctx.out.lock().unwrap(), "{{\"summary\":true,\"mode\":\"{}\",\"base\":{:?},\"b\":{},\"base_genus\":{},\"m\":{},\"d\":{},\"genus\":{},\"excess\":{},\"patterns\":{},\"assign_est\":{:.0},\"candidates\":{},\"checks_ok\":{},\"genus_ok\":{},\"transitive\":{},\"three_cusp_pass\":{},\"skipped_known\":{},\"skipped_bad\":{},\"probe_rejects\":{},\"probe\":{},\"bad_states\":{},\"closures\":{},\"failed_states\":{},\"fail_hist_1_3_10_30_100_1k_10k_more\":{:?},\"found\":{},\"capped\":{},\"cap\":{},\"seconds\":{:.1},\"complete\":true}}",
        mode, bs, b, h, m, m * b, g, excess, npat, work_est, ld(&s.cand), ld(&s.checks_ok), ld(&s.genus_ok), ld(&s.trans), ld(&s.pre), ld(&s.skipped_known), ld(&s.skipped_bad), ld(&s.probe_rejects), ctx.probe, ld(&ctx.bad.n), ld(&s.closures), ld(&s.failed_states), s.fail_hist.iter().map(|a| a.load(Ordering::Relaxed)).collect::<Vec<_>>(), ld(&s.found), ld(&s.capped), cap, t0.elapsed().as_secs_f64()).unwrap();
}
