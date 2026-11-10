//! elevel: exact E-level decider for ONE moving branch point of a cover C -> E (genus-1 quotient lane).
//!
//! State: (A, B, c_1, ..., c_n) in S_k with [A,B] = A B A^-1 B^-1 = c_1 c_2 ... c_n (product x*y = x o y,
//! i.e. (x*y)[i] = x[y[i]]), one c_mv is the MOVING branch point a, the others are FIXED points P_i.
//! Moves are mapping classes of the n-punctured torus (automorphisms of pi_1 preserving peripheral classes):
//!   HA: B -> B A          HB: A -> A B          (Dehn twists about the handle curves)
//!   PC: (A,B,c_1,c_j) -> (A, c_1 B, A c_1 A^-1, c_1 c_j c_1^-1)   (push of puncture 1 through the alpha handle)
//!   BR_i: (c_i, c_{i+1}) -> (c_i c_{i+1} c_i^-1, c_i)            (braid of punctures i, i+1)
//! and their inverses. By Lemma D (frames_20261030/README.md) the collision frames of a over the point-pushing
//! orbit are exactly the frames over the whole mapping-class-group orbit, so:
//!   the family is all-smooth  <=>  in EVERY state, every pair (c_i, c_{i+1}) containing a is additive:
//!   ind(c_i c_{i+1}) = ind c_i + ind c_{i+1}.
//! A non-additive state found by literal moves is a rigorous NEGATIVE (moves are genuine mapping classes).
//! An orbit with no non-additive state under these generators is only a CANDIDATE: it must be confirmed by the
//! degree-3k engine (the generators are not proved to contain the whole push group). With r >= 1 fixed points a
//! confirmed candidate is non-isotrivial by de Franchis.

use std::collections::{HashMap, HashSet};

const K: usize = 32;
type P = [u8; K];

#[derive(Clone, Copy)]
struct Ctx { k: usize }
impl Ctx {
    fn id(&self) -> P { let mut p = [0u8; K]; for i in 0..self.k { p[i] = i as u8; } p }
    fn mul(&self, x: &P, y: &P) -> P { let mut r = [0u8; K]; for i in 0..self.k { r[i] = x[y[i] as usize]; } r }
    fn inv(&self, x: &P) -> P { let mut r = [0u8; K]; for i in 0..self.k { r[x[i] as usize] = i as u8; } r }
    fn conj(&self, g: &P, x: &P) -> P { self.mul(g, &self.mul(x, &self.inv(g))) } // g x g^-1
    fn cyc(&self, x: &P) -> Vec<usize> {
        let mut seen = [false; K]; let mut c = Vec::new();
        for s in 0..self.k { if seen[s] { continue; } let mut l = 0; let mut y = s; while !seen[y] { seen[y] = true; y = x[y] as usize; l += 1; } c.push(l); }
        c.sort_unstable_by(|a, b| b.cmp(a)); c
    }
    fn ind(&self, x: &P) -> usize { self.k - self.cyc(x).len() }
    fn comm(&self, a: &P, b: &P) -> P { self.mul(a, &self.mul(b, &self.mul(&self.inv(a), &self.inv(b)))) }
}

#[derive(Clone)]
struct St { a: P, b: P, c: Vec<P>, mv: usize }

fn transitive(cx: &Ctx, gens: &[&P]) -> bool {
    let mut seen = [false; K]; seen[0] = true; let mut q = vec![0usize]; let mut h = 0;
    while h < q.len() { let x = q[h]; h += 1; for g in gens { let y = g[x] as usize; if !seen[y] { seen[y] = true; q.push(y); } } }
    q.len() == cx.k
}
fn primitive(cx: &Ctx, gens: &[&P]) -> bool {
    let k = cx.k;
    for a in 1..k {
        let mut par: Vec<usize> = (0..k).collect();
        fn f(p: &mut Vec<usize>, x: usize) -> usize { let mut x = x; while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
        let (r0, ra) = (f(&mut par, 0), f(&mut par, a)); par[ra] = r0;
        let mut ch = true;
        while ch { ch = false;
            for g in gens { for x in 0..k { for y in (x + 1)..k { if f(&mut par, x) == f(&mut par, y) {
                let (rx, ry) = (f(&mut par, g[x] as usize), f(&mut par, g[y] as usize)); if rx != ry { par[rx] = ry; ch = true; } } } } } }
        let r = f(&mut par, 0); if (0..k).filter(|&x| f(&mut par, x) == r).count() < k { return false; }
    }
    true
}

/// canonical key modulo simultaneous conjugation (the group is transitive)
fn key(cx: &Ctx, s: &St) -> Vec<u8> {
    let k = cx.k; let n = s.c.len();
    let mut gens: Vec<P> = vec![s.a, s.b]; for i in 0..n - 1 { gens.push(s.c[i]); }
    let ng = gens.len(); for i in 0..ng { let iv = cx.inv(&gens[i]); gens.push(iv); }
    let mut best: Vec<u8> = Vec::new();
    for r in 0..k {
        let mut lab = [u8::MAX; K]; let mut order = [0usize; K]; lab[r] = 0; order[0] = r; let mut len = 1; let mut h = 0;
        while h < len { let x = order[h]; h += 1; for g in gens.iter() { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y; len += 1; } } }
        if len < k { return Vec::new(); }
        let mut key = Vec::with_capacity(ng * k + 1); key.push(s.mv as u8);
        let mut state = if best.is_empty() { 1 } else { 0 };
        'o: for g in gens[..ng].iter() { for &x in order[..k].iter() {
            let v = lab[g[x] as usize];
            if state == 0 { let b = best[key.len()]; if v < b { state = 1; } else if v > b { state = 2; break 'o; } }
            key.push(v);
        } }
        if state == 1 { best = key; }
    }
    best
}

fn check_rel(cx: &Ctx, s: &St) -> bool { let mut p = cx.id(); for c in s.c.iter() { p = cx.mul(&p, c); } p == cx.comm(&s.a, &s.b) }

/// move m: 0 HA, 1 HA^-1, 2 HB, 3 HB^-1, 4 PC, 5 PC^-1, 6+2i BR_i, 7+2i BR_i^-1
fn apply(cx: &Ctx, s: &St, m: usize) -> St {
    let mut t = s.clone(); let n = s.c.len();
    match m {
        0 => t.b = cx.mul(&s.b, &s.a),
        1 => t.b = cx.mul(&s.b, &cx.inv(&s.a)),
        2 => t.a = cx.mul(&s.a, &s.b),
        3 => t.a = cx.mul(&s.a, &cx.inv(&s.b)),
        4 => { let c1 = s.c[0]; t.b = cx.mul(&c1, &s.b); t.c[0] = cx.conj(&s.a, &c1); for j in 1..n { t.c[j] = cx.conj(&c1, &s.c[j]); } }
        5 => { let c1 = cx.conj(&cx.inv(&s.a), &s.c[0]); t.c[0] = c1; t.b = cx.mul(&cx.inv(&c1), &s.b); let ic = cx.inv(&c1); for j in 1..n { t.c[j] = cx.conj(&ic, &s.c[j]); } }
        _ => { let i = (m - 6) / 2; let (x, y) = (s.c[i], s.c[i + 1]);
            if (m - 6) % 2 == 0 { t.c[i] = cx.conj(&x, &y); t.c[i + 1] = x; if s.mv == i { t.mv = i + 1 } else if s.mv == i + 1 { t.mv = i } }
            else { t.c[i] = y; t.c[i + 1] = cx.conj(&cx.inv(&y), &x); if s.mv == i { t.mv = i + 1 } else if s.mv == i + 1 { t.mv = i } } }
    }
    t
}
fn nmoves(n: usize) -> usize { 6 + 2 * (n - 1) }
fn bad(cx: &Ctx, s: &St) -> bool {
    let n = s.c.len(); let m = s.mv;
    let nonadd = |i: usize| cx.ind(&cx.mul(&s.c[i], &s.c[i + 1])) != cx.ind(&s.c[i]) + cx.ind(&s.c[i + 1]);
    (m + 1 < n && nonadd(m)) || (m >= 1 && nonadd(m - 1))
}

enum Res { Bad(usize), Good(usize), Capped(usize) }
/// BFS until a non-additive state; `known` holds keys already proven to lie in a bad orbit.
fn resolve(cx: &Ctx, s0: &St, known: &mut HashSet<Vec<u8>>, cap: usize) -> Res {
    let n = s0.c.len(); let k0 = key(cx, s0);
    if known.contains(&k0) { return Res::Bad(0); }
    let mut seen: HashSet<Vec<u8>> = HashSet::new(); seen.insert(k0.clone());
    let mut q = vec![s0.clone()]; let mut keys = vec![k0]; let mut h = 0;
    while h < q.len() {
        let s = q[h].clone(); h += 1;
        debug_assert!(check_rel(cx, &s));
        if bad(cx, &s) { for kk in keys { known.insert(kk); } return Res::Bad(h); }
        for m in 0..nmoves(n) { let t = apply(cx, &s, m); let kt = key(cx, &t);
            if known.contains(&kt) { for kk in keys { known.insert(kk); } return Res::Bad(h); }
            if seen.insert(kt.clone()) { q.push(t); keys.push(kt); } }
        if q.len() > cap { return Res::Capped(q.len()); }
    }
    Res::Good(q.len())
}

// ---------- permutation utilities ----------
fn from_cycles(k: usize, ct: &[usize]) -> P { let mut p = [0u8; K]; for i in 0..k { p[i] = i as u8; } let mut s = 0; for &l in ct { for j in 0..l { p[s + j] = (s + (j + 1) % l) as u8; } s += l; } p }
fn all_perms(k: usize) -> Vec<P> {
    let mut out = Vec::new(); let mut a: Vec<u8> = (0..k as u8).collect(); let mut c = vec![0usize; k];
    let push = |a: &Vec<u8>, out: &mut Vec<P>| { let mut p = [0u8; K]; p[..k].copy_from_slice(a); out.push(p); };
    push(&a, &mut out); let mut i = 0;
    while i < k { if c[i] < i { if i % 2 == 0 { a.swap(0, i) } else { a.swap(c[i], i) } push(&a, &mut out); c[i] += 1; i = 0; } else { c[i] = 0; i += 1; } }
    out
}
/// one conjugator y with y x y^-1 = z (same cycle type assumed)
fn conjugator(cx: &Ctx, x: &P, z: &P) -> P {
    let cyc = |p: &P| { let mut seen = [false; K]; let mut cs: Vec<Vec<u8>> = Vec::new();
        for s in 0..cx.k { if seen[s] { continue; } let mut c = Vec::new(); let mut y = s; while !seen[y] { seen[y] = true; c.push(y as u8); y = p[y] as usize; } cs.push(c); }
        cs.sort_by(|a, b| b.len().cmp(&a.len())); cs };
    let (cxs, czs) = (cyc(x), cyc(z)); let mut y = [0u8; K];
    for (a, b) in cxs.iter().zip(czs.iter()) { for (u, v) in a.iter().zip(b.iter()) { y[*u as usize] = *v; } }
    y
}
/// all elements of the centralizer of x
fn centralizer(cx: &Ctx, x: &P) -> Vec<P> {
    // y commutes with x iff y maps cycles to cycles of equal length compatibly
    let mut seen = [false; K]; let mut cs: Vec<Vec<u8>> = Vec::new();
    for s in 0..cx.k { if seen[s] { continue; } let mut c = Vec::new(); let mut y = s; while !seen[y] { seen[y] = true; c.push(y as u8); y = x[y] as usize; } cs.push(c); }
    let mut out = Vec::new(); let mut used = vec![false; cs.len()]; let mut y = [0u8; K];
    fn rec(i: usize, cs: &Vec<Vec<u8>>, used: &mut Vec<bool>, y: &mut P, out: &mut Vec<P>) {
        if i == cs.len() { out.push(*y); return; }
        for j in 0..cs.len() { if used[j] || cs[j].len() != cs[i].len() { continue; } used[j] = true; let l = cs[i].len();
            for r in 0..l { for t in 0..l { y[cs[i][t] as usize] = cs[j][(t + r) % l]; } rec(i + 1, cs, used, y, out); } used[j] = false; }
    }
    rec(0, &cs, &mut used, &mut y, &mut out); out
}
fn parse_type(s: &str) -> Vec<usize> { s.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect() }

struct Rng(u64);
impl Rng { fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0; z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
    fn perm(&mut self, k: usize) -> P { let mut p = [0u8; K]; for i in 0..k { p[i] = i as u8; } for i in (1..k).rev() { let j = self.below(i + 1); p.swap(i, j); } p } }

/// Exhaustive r = 1: every transitive (A, B, sigma, tau) with sigma in class S, tau in class T and an
/// additive frame (ind(sigma tau) = ind sigma + ind tau), modulo simultaneous conjugation.
fn exhaust_r1(k: usize, st: &[usize], tt: &[usize], cap: usize, only_prim: bool) {
    let cx = Ctx { k }; let t0 = std::time::Instant::now();
    let sig = from_cycles(k, st); let all = all_perms(k);
    let csig = centralizer(&cx, &sig);
    let tsig = { let mut v = tt.to_vec(); v.retain(|&l| l > 1); v };
    // tau reps up to C(sigma) with additive frame
    let tau0 = from_cycles(k, &tsig); let mut taus: Vec<P> = Vec::new(); let mut seen_t: HashSet<P> = HashSet::new();
    let ct = cx.cyc(&tau0);
    for g in all.iter() { let t = cx.conj(g, &tau0); if seen_t.contains(&t) { continue; }
        for c in csig.iter() { seen_t.insert(cx.conj(c, &t)); }
        let _ = &ct; if cx.ind(&cx.mul(&sig, &t)) == cx.ind(&sig) + cx.ind(&t) { taus.push(t); } }
    let target = cx.ind(&sig) + cx.ind(&tau0);
    println!("k={k} sigma{:?} tau{:?}: {} additive (sigma,tau) classes, total index {}", st, tt, taus.len(), target);
    let mut known: HashSet<Vec<u8>> = HashSet::new();
    let (mut starts, mut nbad, mut ngood, mut ncap, mut nprim, mut maxdepth) = (0u64, 0u64, 0u64, 0u64, 0u64, 0usize);
    let mut goods: Vec<String> = Vec::new();
    let mut start_keys: HashSet<Vec<u8>> = HashSet::new();
    for tau in taus.iter() {
        let g = cx.mul(&sig, tau);
        for a in all.iter() {
            let x = cx.inv(a); let y = cx.mul(&x, &g); if cx.cyc(&x) != cx.cyc(&y) { continue; }
            let b0 = conjugator(&cx, &x, &y);
            for c in centralizer(&cx, &x).iter() {
                let b = cx.mul(&b0, c);
                if !transitive(&cx, &[a, &b, &sig]) { continue; }
                let s = St { a: *a, b, c: vec![sig, *tau], mv: 0 };
                debug_assert!(check_rel(&cx, &s));
                let ks = key(&cx, &s); if !start_keys.insert(ks) { continue; }
                let prim = primitive(&cx, &[a, &b, &sig]); if prim { nprim += 1; } else if only_prim { continue; }
                starts += 1;
                match resolve(&cx, &s, &mut known, cap) {
                    Res::Bad(d) => { nbad += 1; maxdepth = maxdepth.max(d); }
                    Res::Good(nst) => { ngood += 1; goods.push(format!("A={:?} B={:?} sigma={:?} tau={:?} orbit_states={} primitive={}", &a[..k], &b[..k], &sig[..k], &tau[..k], nst, prim)); }
                    Res::Capped(_) => { ncap += 1; }
                }
            }
        }
    }
    println!("  starts {starts} (primitive {nprim})  NEGATIVE {nbad}  ALL-ADDITIVE ORBITS {ngood}  capped {ncap}  max BFS states before hit {maxdepth}  known-bad keys {}  {:.1}s",
        known.len(), t0.elapsed().as_secs_f64());
    for gline in goods.iter().take(20) { println!("  CANDIDATE {gline}"); }
}

/// Random general-n search: sigma moving at position 0, fixed classes given; tuples sampled at random.
fn random_search(k: usize, types: &[Vec<usize>], trials: usize, cap: usize, seed: u64, genus_index: usize) {
    let cx = Ctx { k }; let mut r = Rng(seed); let n = types.len();
    let (mut tried, mut nbad, mut ngood, mut ncap) = (0, 0, 0, 0); let mut known: HashSet<Vec<u8>> = HashSet::new();
    let mut goods = Vec::new();
    let mut attempts = 0u64;
    while tried < trials && attempts < 200_000_000 { attempts += 1;
        let a = r.perm(k); let b = r.perm(k);
        let mut c: Vec<P> = Vec::new();
        for t in types[..n - 1].iter() { let g = r.perm(k); c.push(cx.conj(&g, &from_cycles(k, t))); }
        let mut pre = cx.id(); for x in c.iter() { pre = cx.mul(&pre, x); }
        let last = cx.mul(&cx.inv(&pre), &cx.comm(&a, &b));
        let mut want = types[n - 1].clone(); want.retain(|&l| l > 1);
        let mut have = cx.cyc(&last); have.retain(|&l| l > 1);
        if have != want { continue; }
        c.push(last);
        let tot: usize = c.iter().map(|x| cx.ind(x)).sum(); if tot != genus_index { continue; }
        let mut gens: Vec<&P> = vec![&a, &b]; for x in c.iter() { gens.push(x); }
        if !transitive(&cx, &gens) { continue; }
        let s = St { a, b, c: c.clone(), mv: 0 }; assert!(check_rel(&cx, &s));
        tried += 1;
        match resolve(&cx, &s, &mut known, cap) { Res::Bad(_) => nbad += 1, Res::Good(nst) => { ngood += 1; goods.push((nst, primitive(&cx, &gens))); }, Res::Capped(_) => ncap += 1 }
    }
    println!("random k={k} types {:?}: tuples {tried} NEGATIVE {nbad} ALL-ADDITIVE {ngood} capped {ncap}  good orbits (size, primitive) {:?}", types, &goods[..goods.len().min(10)]);
}

/// Positive control (genus 4): sigma = two 3-cycles permuting 3 blocks of size 2, tau1, tau2 within-block
/// transpositions, A, B block-preserving. C -> D -> E with D genus 2 constant, all frames additive.
fn control_g4(trials: usize, seed: u64, cap: usize) {
    let k = 6; let cx = Ctx { k }; let mut r = Rng(seed);
    let blocks = [[0u8, 1], [2, 3], [4, 5]];
    let sig: P = { let mut p = cx.id(); p[0] = 2; p[2] = 4; p[4] = 0; p[1] = 3; p[3] = 5; p[5] = 1; p };
    let wreath: Vec<P> = all_perms(k).into_iter().filter(|g| blocks.iter().all(|bl| { let b0 = g[bl[0] as usize] / 2; g[bl[1] as usize] / 2 == b0 })).collect();
    let (mut tried, mut good, mut badn) = (0, 0, 0); let mut known = HashSet::new();
    while tried < trials {
        let a = wreath[r.below(wreath.len())]; let b = wreath[r.below(wreath.len())];
        let bl = r.below(3); let t1: P = { let mut p = cx.id(); p[2 * bl] = (2 * bl + 1) as u8; p[2 * bl + 1] = (2 * bl) as u8; p };
        let t2 = cx.mul(&cx.inv(&cx.mul(&sig, &t1)), &cx.comm(&a, &b));
        if cx.ind(&t2) != 1 { continue; }
        if !transitive(&cx, &[&a, &b, &sig, &t1]) { continue; }
        tried += 1;
        let s = St { a, b, c: vec![sig, t1, t2], mv: 0 };
        match resolve(&cx, &s, &mut known, cap) { Res::Good(_) => good += 1, Res::Bad(_) => badn += 1, _ => {} }
    }
    println!("genus-4 positive control (C->D->E, D genus 2): tuples {tried} all-additive orbits {good} negative {badn}");
}


/// BFS depth (number of moves) to the first non-additive state; None if not found within cap.
fn depth_to_bad(cx: &Ctx, s0: &St, cap: usize, moves: &[usize]) -> Option<usize> {
    let mut seen: HashSet<Vec<u8>> = HashSet::new(); seen.insert(key(cx, s0));
    let mut layer = vec![s0.clone()]; let mut d = 0;
    while !layer.is_empty() {
        for s in layer.iter() { if bad(cx, s) { return Some(d); } }
        let mut next = Vec::new();
        for s in layer.iter() { for &m in moves { let t = apply(cx, s, m); if seen.insert(key(cx, &t)) { next.push(t); } } }
        if seen.len() > cap { return None; }
        layer = next; d += 1;
    }
    None
}
/// Lemma test: is every additive r=1 start refuted by "Nielsen-twist the handles, then push a around
/// the first handle n times", i.e. frames ((s^-1 A')^n s (s^-1 A')^-n, t), A' over first coordinates of the
/// SL2(Z)-orbit of (A,B)?  Reports starts NOT refuted this way.
fn push_lemma(k: usize, st: &[usize], tt: &[usize], samples: usize, seed: u64) {
    let cx = Ctx { k }; let mut r = Rng(seed);
    let sig = from_cycles(k, st); let mut tt2 = tt.to_vec(); tt2.retain(|&l| l > 1);
    let (mut tested, mut refuted, mut depthmax, mut hist) = (0, 0, 0usize, std::collections::BTreeMap::new());
    let mut attempts = 0u64;
    while tested < samples && attempts < 500_000_000 { attempts += 1;
        let a = r.perm(k); let g = r.perm(k); let tau = cx.conj(&g, &from_cycles(k, &tt2));
        if cx.ind(&cx.mul(&sig, &tau)) != cx.ind(&sig) + cx.ind(&tau) { continue; }
        let gg = cx.mul(&sig, &tau); let x = cx.inv(&a); let y = cx.mul(&x, &gg);
        if cx.cyc(&x) != cx.cyc(&y) { continue; }
        let cents = centralizer(&cx, &x); let b = cx.mul(&conjugator(&cx, &x, &y), &cents[r.below(cents.len())]);
        if !transitive(&cx, &[&a, &b, &sig]) { continue; }
        tested += 1;
        // SL2(Z)-orbit of (A,B) (finite), collect first coordinates
        let mut seen: HashSet<(P, P)> = HashSet::new(); let mut q = vec![(a, b)]; seen.insert((a, b)); let mut h = 0;
        while h < q.len() && q.len() < 200000 { let (x1, y1) = q[h]; h += 1;
            for (x2, y2) in [(x1, cx.mul(&y1, &x1)), (x1, cx.mul(&y1, &cx.inv(&x1))), (cx.mul(&x1, &y1), y1), (cx.mul(&x1, &cx.inv(&y1)), y1)] {
                if seen.insert((x2, y2)) { q.push((x2, y2)); } } }
        let mut ok = false;
        'outer: for (ap, _) in q.iter() { let u = cx.mul(&cx.inv(&sig), ap); let mut un = u; 
            for n in 1..=420 { let sx = cx.conj(&un, &sig); if cx.ind(&cx.mul(&sx, &tau)) != cx.ind(&sx) + cx.ind(&tau) { ok = true; *hist.entry(n).or_insert(0) += 1; depthmax = depthmax.max(n); break 'outer; }
                un = cx.mul(&un, &u); if un == cx.id() { break; } } }
        if ok { refuted += 1; }
    }
    println!("push-lemma k={k} sigma{:?} tau{:?}: tested {tested} refuted-by-twist+push {refuted}  push power histogram {:?}", st, tt, hist);
}

/// Exhaustive r = 2: (A, B, sigma, tau1, tau2), [A,B] = sigma tau1 tau2, sigma moving (position 1),
/// starts with (sigma, tau1) additive (every state of an all-smooth orbit is), tau1 up to C(sigma),
/// tau2 over its whole class, all A, B solved.
fn exhaust_r2(k: usize, st: &[usize], t1: &[usize], t2: &[usize], cap: usize) {
    let cx = Ctx { k }; let t0 = std::time::Instant::now();
    let sig = from_cycles(k, st); let all = all_perms(k); let csig = centralizer(&cx, &sig);
    let strip = |v: &[usize]| { let mut w = v.to_vec(); w.retain(|&l| l > 1); w };
    let tau10 = from_cycles(k, &strip(t1)); let tau20 = from_cycles(k, &strip(t2));
    let mut tau1s: Vec<P> = Vec::new(); let mut seen_t: HashSet<P> = HashSet::new();
    for g in all.iter() { let t = cx.conj(g, &tau10); if seen_t.contains(&t) { continue; }
        for c in csig.iter() { seen_t.insert(cx.conj(c, &t)); }
        if cx.ind(&cx.mul(&sig, &t)) == cx.ind(&sig) + cx.ind(&t) { tau1s.push(t); } }
    let mut class2: HashSet<P> = HashSet::new(); for g in all.iter() { class2.insert(cx.conj(g, &tau20)); }
    let class2: Vec<P> = class2.into_iter().collect();
    println!("k={k} sigma{:?} tau1{:?} tau2{:?}: {} additive (sigma,tau1) classes, |class tau2| {}", st, t1, t2, tau1s.len(), class2.len());
    let mut known: HashSet<Vec<u8>> = HashSet::new(); let mut start_keys: HashSet<Vec<u8>> = HashSet::new();
    let (mut starts, mut nbad, mut ngood, mut ncap, mut nprim) = (0u64, 0u64, 0u64, 0u64, 0u64); let mut goods = Vec::new();
    // index all perms by cycle type for fast A filtering
    for tau1 in tau1s.iter() { for tau2 in class2.iter() {
        let g = cx.mul(&sig, &cx.mul(tau1, tau2));
        for a in all.iter() {
            let x = cx.inv(a); let y = cx.mul(&x, &g); if cx.cyc(&x) != cx.cyc(&y) { continue; }
            let b0 = conjugator(&cx, &x, &y);
            for c in centralizer(&cx, &x).iter() {
                let b = cx.mul(&b0, c);
                if !transitive(&cx, &[a, &b, &sig, tau1]) { continue; }
                let s = St { a: *a, b, c: vec![sig, *tau1, *tau2], mv: 0 };
                let ks = key(&cx, &s); if !start_keys.insert(ks) { continue; }
                if primitive(&cx, &[a, &b, &sig, tau1]) { nprim += 1; }
                starts += 1;
                match resolve(&cx, &s, &mut known, cap) { Res::Bad(_) => nbad += 1,
                    Res::Good(nst) => { ngood += 1; goods.push(format!("A={:?} B={:?} sigma={:?} tau1={:?} tau2={:?} states={}", &a[..k], &b[..k], &sig[..k], &tau1[..k], &tau2[..k], nst)); },
                    Res::Capped(_) => ncap += 1 }
            }
        }
    } }
    println!("  starts {starts} (primitive {nprim})  NEGATIVE {nbad}  ALL-ADDITIVE ORBITS {ngood}  capped {ncap}  known-bad {}  {:.1}s", known.len(), t0.elapsed().as_secs_f64());
    for gl in goods.iter().take(20) { println!("  CANDIDATE {gl}"); }
}

/// Exhaustive general r (= types.len()-1 fixed points): sigma = types[0] moving at position 1, tau1 up to
/// C(sigma) with (sigma, tau1) additive, tau2..tau_{r-1} over full classes, tau_r and B solved.
fn exhaust_gen(k: usize, types: &[Vec<usize>], cap: usize) {
    let cx = Ctx { k }; let t0 = std::time::Instant::now(); let n = types.len(); assert!(n >= 3);
    let strip = |v: &[usize]| { let mut w = v.to_vec(); w.retain(|&l| l > 1); w };
    let sig = from_cycles(k, &strip(&types[0])); let all = all_perms(k); let csig = centralizer(&cx, &sig);
    let reps: Vec<P> = types.iter().map(|t| from_cycles(k, &strip(t))).collect();
    let mut tau1s: Vec<P> = Vec::new(); let mut seen_t: HashSet<P> = HashSet::new();
    for g in all.iter() { let t = cx.conj(g, &reps[1]); if seen_t.contains(&t) { continue; }
        for c in csig.iter() { seen_t.insert(cx.conj(c, &t)); }
        if cx.ind(&cx.mul(&sig, &t)) == cx.ind(&sig) + cx.ind(&t) { tau1s.push(t); } }
    let classes: Vec<Vec<P>> = (2..n).map(|i| { let mut h: HashSet<P> = HashSet::new(); for g in all.iter() { h.insert(cx.conj(g, &reps[i])); } h.into_iter().collect() }).collect();
    let last_ct = { let mut c = cx.cyc(&reps[n - 1]); c.sort(); c };
    let total: usize = reps.iter().map(|x| cx.ind(x)).sum();
    println!("k={k} types {:?}: {} additive (sigma,tau1) classes, middle class sizes {:?}, total index {total}", types, tau1s.len(), classes.iter().map(|c| c.len()).collect::<Vec<_>>());
    let mut known: HashSet<Vec<u8>> = HashSet::new(); let mut start_keys: HashSet<Vec<u8>> = HashSet::new();
    let (mut starts, mut nbad, mut ngood, mut ncap, mut nprim) = (0u64, 0u64, 0u64, 0u64, 0u64); let mut goods = Vec::new();
    // iterate middle choices
    let mut idx = vec![0usize; classes.len()];
    loop {
        for tau1 in tau1s.iter() {
            let mut mid: Vec<P> = vec![sig, *tau1]; for (j, c) in classes.iter().enumerate() { mid.push(c[idx[j]]); }
            let mut pre = cx.id(); for x in mid.iter() { pre = cx.mul(&pre, x); }
            { let _ = &last_ct;
                let g = pre;
                // pick A by scanning; to keep cost bounded only A with matching cycle type are used
                for a in all.iter() {
                    let x = cx.inv(a); let y = cx.mul(&x, &g); if cx.cyc(&x) != cx.cyc(&y) { continue; }
                    let b0 = conjugator(&cx, &x, &y);
                    for c in centralizer(&cx, &x).iter() {
                        let b = cx.mul(&b0, c);
                        let cs = mid.clone();
                        let mut gens: Vec<&P> = vec![a, &b]; for z in cs.iter() { gens.push(z); }
                        if !transitive(&cx, &gens) { continue; }
                        let s = St { a: *a, b, c: cs.clone(), mv: 0 };
                        let ks = key(&cx, &s); if !start_keys.insert(ks) { continue; }
                        if primitive(&cx, &gens) { nprim += 1; }
                        starts += 1;
                        match resolve(&cx, &s, &mut known, cap) { Res::Bad(_) => nbad += 1,
                            Res::Good(nst) => { ngood += 1; goods.push(format!("A={:?} B={:?} c={:?} states={}", &a[..k], &b[..k], cs.iter().map(|z| z[..k].to_vec()).collect::<Vec<_>>(), nst)); },
                            Res::Capped(_) => ncap += 1 }
                    }
                }
            }
        }
        // advance odometer
        let mut j = 0; while j < idx.len() { idx[j] += 1; if idx[j] < classes[j].len() { break; } idx[j] = 0; j += 1; }
        if j == idx.len() { break; }
    }
    println!("  starts {starts} (primitive {nprim})  NEGATIVE {nbad}  ALL-ADDITIVE ORBITS {ngood}  capped {ncap}  known-bad {}  {:.1}s", known.len(), t0.elapsed().as_secs_f64());
    for gl in goods.iter().take(10) { println!("  CANDIDATE {gl}"); }
}

/// Random additive starts (r = 1 or 2), report the exact move-distance to the first non-additive state.
fn depth_stats(k: usize, types: &[Vec<usize>], samples: usize, seed: u64, cap: usize) {
    let cx = Ctx { k }; let mut r = Rng(seed); let n = types.len();
    let strip = |v: &[usize]| { let mut w = v.to_vec(); w.retain(|&l| l > 1); w };
    let reps: Vec<P> = types.iter().map(|t| from_cycles(k, &strip(t))).collect();
    let mut hist: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new(); let mut none = 0; let mut tested = 0; let mut att = 0u64;
    let all_moves: Vec<usize> = (0..nmoves(n)).collect();
    while tested < samples && att < 2_000_000_000 { att += 1;
        let mut cs: Vec<P> = vec![reps[0]];
        for i in 1..n { let g = r.perm(k); cs.push(cx.conj(&g, &reps[i])); }
        if cx.ind(&cx.mul(&cs[0], &cs[1])) != cx.ind(&cs[0]) + cx.ind(&cs[1]) { continue; }
        let mut pre = cx.id(); for x in cs.iter() { pre = cx.mul(&pre, x); }
        let a = r.perm(k); let x = cx.inv(&a); let y = cx.mul(&x, &pre); if cx.cyc(&x) != cx.cyc(&y) { continue; }
        let cents = centralizer(&cx, &x); let b = cx.mul(&conjugator(&cx, &x, &y), &cents[r.below(cents.len())]);
        let mut gens: Vec<&P> = vec![&a, &b]; for z in cs.iter() { gens.push(z); } if !transitive(&cx, &gens) { continue; }
        tested += 1;
        let s = St { a, b, c: cs.clone(), mv: 0 };
        match depth_to_bad(&cx, &s, cap, &all_moves) { Some(d) => *hist.entry(d).or_insert(0) += 1, None => none += 1 }
    }
    println!("depth k={k} types {:?}: samples {tested} move-distance histogram {:?} unresolved {none}", types, hist);
}

// ---------------- genus-0 quotient lane: P^1-level decider ----------------
// State (c_0..c_r) in S_k, c_0 c_1 ... c_r = 1, one entry is the moving point a. Moves: braids of adjacent
// entries (sphere). A collision of a with an adjacent fixed point is smooth iff the pair is additive OR the
// complementary disk (all other entries) is additive: ind(c_i c_{i+1}) == ind c_i + ind c_{i+1} or
// == sum of the other indices. r >= 3 fixed points (Lemma D needs chi < 0).
#[derive(Clone)]
struct S0 { c: Vec<P>, mv: usize }
fn key0(cx: &Ctx, s: &S0) -> Vec<u8> {
    let k = cx.k; let n = s.c.len();
    let mut gens: Vec<P> = s.c[..n - 1].to_vec(); let ng = gens.len(); for i in 0..ng { let iv = cx.inv(&gens[i]); gens.push(iv); }
    let mut best: Vec<u8> = Vec::new();
    for r in 0..k {
        let mut lab = [u8::MAX; K]; let mut order = [0usize; K]; lab[r] = 0; order[0] = r; let mut len = 1; let mut h = 0;
        while h < len { let x = order[h]; h += 1; for g in gens.iter() { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y; len += 1; } } }
        if len < k { return Vec::new(); }
        let mut key = Vec::with_capacity(ng * k + 1); key.push(s.mv as u8); let mut state = if best.is_empty() { 1 } else { 0 };
        'o: for g in gens[..ng].iter() { for &x in order[..k].iter() { let v = lab[g[x] as usize];
            if state == 0 { let b = best[key.len()]; if v < b { state = 1; } else if v > b { state = 2; break 'o; } } key.push(v); } }
        if state == 1 { best = key; }
    }
    best
}
fn apply0(cx: &Ctx, s: &S0, m: usize) -> S0 {
    let mut t = s.clone(); let i = m / 2; let (x, y) = (s.c[i], s.c[i + 1]);
    if m % 2 == 0 { t.c[i] = cx.conj(&x, &y); t.c[i + 1] = x; } else { t.c[i] = y; t.c[i + 1] = cx.conj(&cx.inv(&y), &x); }
    if s.mv == i { t.mv = i + 1 } else if s.mv == i + 1 { t.mv = i }
    t
}
fn bad0(cx: &Ctx, s: &S0) -> bool {
    let n = s.c.len(); let inds: Vec<usize> = s.c.iter().map(|x| cx.ind(x)).collect(); let tot: usize = inds.iter().sum();
    let chk = |i: usize, j: usize| { let z = cx.ind(&cx.mul(&s.c[i], &s.c[j])); !(z == inds[i] + inds[j] || z == tot - inds[i] - inds[j]) };
    let m = s.mv;
    (m + 1 < n && chk(m, m + 1)) || (m >= 1 && chk(m - 1, m)) || (m == 0 && chk(n - 1, 0)) || (m == n - 1 && chk(n - 1, 0))
}
fn resolve0(cx: &Ctx, s0: &S0, known: &mut HashSet<Vec<u8>>, cap: usize) -> Res {
    let n = s0.c.len(); let k0 = key0(cx, s0); if known.contains(&k0) { return Res::Bad(0); }
    let mut seen: HashSet<Vec<u8>> = HashSet::new(); seen.insert(k0.clone()); let mut q = vec![s0.clone()]; let mut keys = vec![k0]; let mut h = 0;
    while h < q.len() { let s = q[h].clone(); h += 1;
        if bad0(cx, &s) { for kk in keys { known.insert(kk); } return Res::Bad(h); }
        for m in 0..2 * (n - 1) { let t = apply0(cx, &s, m); let kt = key0(cx, &t);
            if known.contains(&kt) { for kk in keys { known.insert(kk); } return Res::Bad(h); }
            if seen.insert(kt.clone()) { q.push(t); keys.push(kt); } }
        if q.len() > cap { return Res::Capped(q.len()); } }
    Res::Good(q.len())
}
/// random genus-0-lane search: r fixed points, total index 2k+4, random classes of bounded index
fn p1_random(k: usize, r: usize, trials: usize, seed: u64, cap: usize) {
    let cx = Ctx { k }; let mut rng = Rng(seed); let n = r + 1; let target = 2 * k + 4;
    let (mut tried, mut nbad, mut ngood, mut ncap, mut att) = (0, 0, 0, 0, 0u64); let mut known = HashSet::new(); let mut goods = Vec::new();
    while tried < trials && att < 3_000_000_000 { att += 1;
        // index-constrained generation: ind c_0 = s in {1,2}, every fixed entry index <= k-1-s
        let sidx = 1 + rng.below(2); let cap_i = k - 1 - sidx;
        let mut c: Vec<P> = Vec::new();
        loop { let x = rng.perm(k); if cx.ind(&x) == sidx { c.push(x); break; } }
        for _ in 1..n - 1 { loop { let x = rng.perm(k); let i = cx.ind(&x); if i >= 1 && i <= cap_i { c.push(x); break; } } }
        let mut pre = cx.id(); for x in c.iter() { pre = cx.mul(&pre, x); } c.push(cx.inv(&pre));
        if c.iter().any(|x| cx.ind(x) == 0) || cx.ind(&c[n - 1]) > cap_i { continue; }
        let tot: usize = c.iter().map(|x| cx.ind(x)).sum(); if tot != target { continue; }
        let gens: Vec<&P> = c.iter().collect(); if !transitive(&cx, &gens) { continue; }
        let s = S0 { c: c.clone(), mv: 0 }; if bad0(&cx, &s) { continue; }
        tried += 1;
        match resolve0(&cx, &s, &mut known, cap) { Res::Bad(_) => nbad += 1, Res::Capped(_) => ncap += 1,
            Res::Good(nst) => { ngood += 1; goods.push(format!("types {:?} prim={} states={} c={:?}", c.iter().map(|x| cx.cyc(x).into_iter().filter(|&l| l > 1).collect::<Vec<_>>()).collect::<Vec<_>>(), primitive(&cx, &gens), nst, c.iter().map(|x| x[..k].to_vec()).collect::<Vec<_>>())); } }
    }
    println!("p1 k={k} r={r} (genus 3, total index {target}): smooth-start tuples {tried} NEGATIVE {nbad} SURVIVING ORBITS {ngood} capped {ncap} (attempts {att})");
    for g in goods.iter().take(15) { println!("  SURVIVOR {g}"); }
}

fn parse_cycles(k: usize, s: &str) -> P {
    let mut p = [0u8; K]; for i in 0..k { p[i] = i as u8; }
    for cyc in s.split(')').filter(|c| c.contains('(')) {
        let v: Vec<usize> = cyc.trim().trim_start_matches('(').split_whitespace().map(|x| x.parse::<usize>().unwrap() - 1).collect();
        for i in 0..v.len() { p[v[i]] = v[(i + 1) % v.len()] as u8; }
    }
    p
}
/// positive/negative controls for the P^1-level decider: a known 4-tuple (any entry may be "a" for r = 3)
fn p1_check(k: usize, tuple: &str, cap: usize) {
    let cx = Ctx { k }; let c: Vec<P> = tuple.split(';').map(|t| parse_cycles(k, t)).collect();
    let mut pre = cx.id(); for x in c.iter() { pre = cx.mul(&pre, x); }
    let c = if pre != cx.id() { // the engine convention is g1 g2 g3 g4 = 1 with compose(p,q)=p o q, same as ours
        println!("note: product != 1 under this convention; using reversed order"); c.into_iter().rev().collect() } else { c };
    let mut known = HashSet::new();
    for mv in 0..c.len() { let s = S0 { c: c.clone(), mv };
        let r = match resolve0(&cx, &s, &mut known, cap) { Res::Good(n) => format!("ALL-SMOOTH orbit, {n} states"), Res::Bad(h) => format!("NEGATIVE after {h}"), Res::Capped(n) => format!("capped {n}") };
        println!("p1-check k={k} moving entry {mv}: {r}"); known.clear(); }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |f: &str, d: &str| -> String { args.iter().position(|a| a == f).map(|i| args[i + 1].clone()).unwrap_or(d.to_string()) };
    let k: usize = get("--k", "5").parse().unwrap(); let cap: usize = get("--cap", "2000000").parse().unwrap();
    if args.iter().any(|a| a == "--push-lemma") { push_lemma(k, &parse_type(&get("--sig", "3")), &parse_type(&get("--tau", "3")), get("--trials", "200").parse().unwrap(), get("--seed", "1").parse().unwrap()); return; }
    if args.iter().any(|a| a == "--exhaust-r2") { exhaust_r2(k, &parse_type(&get("--sig", "3")), &parse_type(&get("--tau1", "2")), &parse_type(&get("--tau2", "2")), cap); return; }
    if args.iter().any(|a| a == "--exhaust-gen") { let types: Vec<Vec<usize>> = get("--types", "2;2;2;2").split(';').map(parse_type).collect(); exhaust_gen(k, &types, cap); return; }
    if args.iter().any(|a| a == "--depth") { let types: Vec<Vec<usize>> = get("--types", "3;3").split(';').map(parse_type).collect(); depth_stats(k, &types, get("--trials", "200").parse().unwrap(), get("--seed", "1").parse().unwrap(), cap); return; }
    if args.iter().any(|a| a == "--p1") { p1_random(k, get("--r", "4").parse().unwrap(), get("--trials", "200").parse().unwrap(), get("--seed", "1").parse().unwrap(), cap); return; }
    if args.iter().any(|a| a == "--p1-check") { p1_check(k, &get("--tuple", ""), cap); return; }
    if args.iter().any(|a| a == "--control") { control_g4(get("--trials", "200").parse().unwrap(), get("--seed", "1").parse().unwrap(), cap); return; }
    if args.iter().any(|a| a == "--exhaust-r1") {
        exhaust_r1(k, &parse_type(&get("--sig", "3")), &parse_type(&get("--tau", "3")), cap, args.iter().any(|a| a == "--primitive-only"));
        return;
    }
    if args.iter().any(|a| a == "--random") {
        let types: Vec<Vec<usize>> = get("--types", "3;2,2").split(';').map(parse_type).collect();
        random_search(k, &types, get("--trials", "100").parse().unwrap(), cap, get("--seed", "1").parse().unwrap(), get("--total-index", "4").parse().unwrap());
        return;
    }
    let _ = HashMap::<u8, u8>::new();
    println!("modes: --exhaust-r1 --k K --sig 3 --tau 2,2 | --random --k K --types '3;2;1' --total-index 4 | --control");
}
