//! frames: test the "frames conjecture" for single-moving-point covers C -> E (genus-1 quotient lane).
//!
//! Model: degree d = 3k four-tuples on the m = 3 genus-1 dessin x = y = z = (0 1 2) (blocks of size k),
//! t1 = (id; sigma on block 0), t2 = (x; f), t3 = (y; g), t4 = (t1 t2 t3)^-1, source genus 3.
//! The full braid orbit (six moves, modulo simultaneous S_d conjugation) is the point-pushing orbit of
//! the moving branch point a of C -> E. At each state, the entry acting trivially on blocks is s, its
//! nontrivial block j carries sigma', and a collision of s with an adjacent entry t is the frame
//! (sigma', tau) with tau = t^L restricted to block j (L = length of the block cycle of t through j):
//! tau is the fixed branching of C -> E at the face point that a reaches.
//! Observed frames F_obs (up to simultaneous S_k conjugation) are compared with
//! F_full = {(x sigma' x^-1, tau) : x in G_E}, G_E = stabiliser of block j restricted to block j
//! (the monodromy group of C -> E). Conjecture: F_full is contained in F_obs.
//! A collision is smooth (additive) iff ind(sigma' tau) = ind sigma' + ind tau.
//! Conventions: compose(p,q)[x] = p[q[x]] (rightmost first); t1 t2 t3 t4 = 1.

use std::collections::{HashMap, HashSet};

const MX: usize = 32;
type A = [u8; MX];
type P = Vec<u8>;

fn acomp(p: &A, q: &A, n: usize) -> A { let mut r = [0u8; MX]; for i in 0..n { r[i] = p[q[i] as usize]; } r }
fn ainv(p: &A, n: usize) -> A { let mut r = [0u8; MX]; for i in 0..n { r[p[i] as usize] = i as u8; } r }
fn aindex(p: &[u8]) -> usize {
    let n = p.len(); let mut seen = vec![false; n]; let mut c = 0;
    for s in 0..n { if seen[s] { continue; } c += 1; let mut x = s; while !seen[x] { seen[x] = true; x = p[x] as usize; } }
    n - c
}
fn comp(p: &[u8], q: &[u8]) -> P { q.iter().map(|&x| p[x as usize]).collect() }
fn inv(p: &[u8]) -> P { let mut r = vec![0u8; p.len()]; for (i, &x) in p.iter().enumerate() { r[x as usize] = i as u8; } r }
fn ctype(p: &[u8]) -> Vec<usize> {
    let n = p.len(); let mut seen = vec![false; n]; let mut c = Vec::new();
    for s in 0..n { if seen[s] { continue; } let mut l = 0; let mut x = s; while !seen[x] { seen[x] = true; x = p[x] as usize; l += 1; } if l > 1 { c.push(l); } }
    c.sort_unstable_by(|a, b| b.cmp(a)); c
}

fn mv(g: &mut [A; 4], m: usize, n: usize) {
    let j = m % 3; let a = g[j]; let b = g[j + 1];
    if m < 3 { g[j] = acomp(&a, &acomp(&b, &ainv(&a, n), n), n); g[j + 1] = a; }
    else { g[j] = b; g[j + 1] = acomp(&ainv(&b, n), &acomp(&a, &b, n), n); }
}

fn canonical(g: &[A; 4], n: usize) -> Vec<u8> {
    let iv = [ainv(&g[0], n), ainv(&g[1], n), ainv(&g[2], n), ainv(&g[3], n)];
    let mut best: Vec<u8> = Vec::new(); let mut lab = [u8::MAX; MX]; let mut order = [0usize; MX];
    let mut key = Vec::with_capacity(3 * n);
    for r in 0..n {
        for v in lab[..n].iter_mut() { *v = u8::MAX; } lab[r] = 0; order[0] = r; let mut len = 1; let mut h = 0;
        while h < len { let x = order[h]; h += 1;
            for p in g.iter().chain(iv.iter()) { let y = p[x] as usize; if lab[y] == u8::MAX { lab[y] = len as u8; order[len] = y; len += 1; } } }
        key.clear(); let mut state = if best.is_empty() { 1 } else { 0 };
        'o: for p in g[..3].iter() { for &x in order[..n].iter() {
            let v = lab[p[x] as usize];
            if state == 0 { let b = best[key.len()]; if v < b { state = 1; } else if v > b { state = 2; break 'o; } }
            key.push(v);
        } }
        if state == 1 { best = key.clone(); }
    }
    best
}

/// canonical form of a pair (a, b) of permutations of k points up to simultaneous conjugation,
/// valid for intransitive pairs: per-orbit minimal BFS code, orbit codes sorted.
fn canon2(a: &[u8], b: &[u8]) -> Vec<u8> {
    let k = a.len(); let ai = inv(a); let bi = inv(b); let gens = [a, b, &ai[..], &bi[..]];
    let mut orb = vec![usize::MAX; k]; let mut orbits: Vec<Vec<usize>> = Vec::new();
    for s in 0..k { if orb[s] != usize::MAX { continue; } let id = orbits.len(); let mut o = vec![s]; orb[s] = id; let mut h = 0;
        while h < o.len() { let x = o[h]; h += 1; for gg in gens.iter() { let y = gg[x] as usize; if orb[y] == usize::MAX { orb[y] = id; o.push(y); } } }
        orbits.push(o); }
    let mut codes: Vec<Vec<u8>> = Vec::new();
    for o in orbits.iter() {
        let mut best: Option<Vec<u8>> = None;
        for &r in o.iter() {
            let mut lab: HashMap<usize, u8> = HashMap::new(); let mut order = vec![r]; lab.insert(r, 0); let mut h = 0;
            while h < order.len() { let x = order[h]; h += 1; for gg in gens.iter() { let y = gg[x] as usize; if !lab.contains_key(&y) { let l = order.len() as u8; lab.insert(y, l); order.push(y); } } }
            let mut code = vec![o.len() as u8];
            for gg in [a, b].iter() { for &x in order.iter() { code.push(lab[&(gg[x] as usize)]); } }
            if best.as_ref().map_or(true, |bb| code < *bb) { best = Some(code); }
        }
        codes.push(best.unwrap());
    }
    codes.sort(); codes.concat()
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0; z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
    fn perm(&mut self, k: usize) -> P { let mut p: P = (0..k as u8).collect(); for i in (1..k).rev() { let j = self.below(i + 1); p.swap(i, j); } p }
}

fn group_closure(gens: &[P], k: usize, cap: usize) -> Vec<P> {
    let id: P = (0..k as u8).collect(); let mut seen: HashSet<P> = HashSet::new(); seen.insert(id.clone());
    let mut q = vec![id]; let mut h = 0;
    while h < q.len() { let x = q[h].clone(); h += 1; for g in gens { let y = comp(g, &x); if !seen.contains(&y) { seen.insert(y.clone()); q.push(y); if q.len() > cap { return q; } } } }
    q
}
fn primitive(gens: &[P], k: usize) -> bool {
    for a in 1..k {
        // minimal block containing {0, a} via union-find closure
        let mut par: Vec<usize> = (0..k).collect();
        fn f(p: &mut Vec<usize>, x: usize) -> usize { let mut x = x; while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
        let r0 = f(&mut par, 0); let ra = f(&mut par, a); par[ra] = r0;
        let mut changed = true;
        while changed { changed = false;
            for g in gens { for x in 0..k { for y in 0..k { if f(&mut par, x) == f(&mut par, y) {
                let gx = g[x] as usize; let gy = g[y] as usize; let (rx, ry) = (f(&mut par, gx), f(&mut par, gy));
                if rx != ry { par[rx] = ry; changed = true; } } } } } }
        let r = f(&mut par, 0); let sz = (0..k).filter(|&x| f(&mut par, x) == r).count();
        if sz < k { return false; }
    }
    true
}

fn blockperm(t: &A, k: usize) -> [usize; 3] { [t[0] as usize / k, t[k] as usize / k, t[2 * k] as usize / k] }
fn restrict(t: &A, j: usize, k: usize) -> P { (0..k).map(|i| t[j * k + i] - (j * k) as u8).collect() }

/// stabiliser of block j in <g>, restricted to block j (Schreier generators)
fn ge_gens(g: &[A; 4], j: usize, k: usize) -> Vec<P> {
    let n = 3 * k; let mut id = [0u8; MX]; for i in 0..n { id[i] = i as u8; }
    let mut u: Vec<Option<A>> = vec![None; 3]; u[j] = Some(id); let mut q = vec![j]; let mut h = 0;
    while h < q.len() { let b = q[h]; h += 1; for t in g.iter() { let c = blockperm(t, k)[b]; if u[c].is_none() { u[c] = Some(acomp(t, &u[b].unwrap(), n)); q.push(c); } } }
    let mut out = Vec::new();
    for b in 0..3 { for t in g.iter() { let c = blockperm(t, k)[b];
        let s = acomp(&ainv(&u[c].unwrap(), n), &acomp(t, &u[b].unwrap(), n), n);
        let r = restrict(&s, j, k); if !out.contains(&r) { out.push(r); } } }
    out
}

struct Frame { sig: P, tau: P, j: usize, gi: usize }
fn frames_of(g: &[A; 4], k: usize) -> Vec<Frame> {
    let mut out = Vec::new();
    let si = (0..4).filter(|&i| blockperm(&g[i], k) == [0, 1, 2]).collect::<Vec<_>>();
    assert!(si.len() == 1, "expected exactly one block-trivial entry");
    let i = si[0]; let s = g[i];
    let js: Vec<usize> = (0..3).filter(|&b| (0..k).any(|a| s[b * k + a] as usize != b * k + a)).collect();
    assert!(js.len() == 1, "single moving block expected"); let j = js[0];
    let sig = restrict(&s, j, k);
    for &nb in &[(i + 1) % 4, (i + 3) % 4] {
        let t = g[nb]; let bp = blockperm(&t, k); let mut l = 1; let mut b = bp[j]; while b != j { b = bp[b]; l += 1; }
        let mut tl = t; for _ in 1..l { tl = acomp(&t, &tl, 3 * k); }
        let tau = restrict(&tl, j, k);
        if aindex(&tau) > 0 { out.push(Frame { sig: sig.clone(), tau, j, gi: nb }); }
    }
    out
}
fn additive(a: &[u8], b: &[u8]) -> bool { aindex(&comp(a, b)) == aindex(a) + aindex(b) }

fn gen_tuple(k: usize, r: &mut Rng, want_r: usize) -> Option<[A; 4]> {
    let n = 3 * k;
    for _ in 0..2_000_000 {
        let sig = r.perm(k); let is = aindex(&sig); if is == 0 || is > 3 { continue; }
        let tau2 = if r.below(2) == 0 { (0..k as u8).collect::<P>() } else { r.perm(k) };
        let i2 = aindex(&tau2); if is + i2 > 4 { continue; }
        let tau3 = if r.below(2) == 0 { (0..k as u8).collect::<P>() } else { r.perm(k) };
        let i3 = aindex(&tau3); if is + i2 + i3 > 4 { continue; }
        let g0 = r.perm(k); let g1 = r.perm(k); let g2 = comp(&tau3, &inv(&comp(&g1, &g0)));
        let id: P = (0..k as u8).collect(); let f = [id.clone(), id.clone(), tau2.clone()]; let gg = [g0, g1, g2];
        let mut t1 = [0u8; MX]; let mut t2 = [0u8; MX]; let mut t3 = [0u8; MX];
        for b in 0..3 { for a in 0..k {
            t1[b * k + a] = (b * k) as u8 + if b == 0 { sig[a] } else { a as u8 };
            t2[b * k + a] = (((b + 1) % 3) * k) as u8 + f[b][a];
            t3[b * k + a] = (((b + 1) % 3) * k) as u8 + gg[b][a];
        } }
        let t4 = ainv(&acomp(&t1, &acomp(&t2, &t3, n), n), n);
        let g = [t1, t2, t3, t4];
        let tot: usize = g.iter().map(|t| aindex(&t[..n])).sum();
        if tot != 2 * n + 4 { continue; }
        // r = number of fixed branch points of C -> E (nontrivial holonomies of t2, t3, t4)
        let mut hol = 0; for t in [&t2, &t3, &t4] { let mut tl = *t; for _ in 1..3 { tl = acomp(t, &tl, n); } if aindex(&restrict(&tl, 0, k)) > 0 { hol += 1; } }
        if hol != want_r { continue; }
        // transitivity on 3k points
        let mut seen = vec![false; n]; seen[0] = true; let mut q = vec![0usize]; let mut h = 0;
        while h < q.len() { let x = q[h]; h += 1; for t in g.iter() { let y = t[x] as usize; if !seen[y] { seen[y] = true; q.push(y); } } }
        if q.len() < n { continue; }
        return Some(g);
    }
    None
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |f: &str, d: &str| -> String { args.iter().position(|a| a == f).map(|i| args[i + 1].clone()).unwrap_or(d.to_string()) };
    let k: usize = get("--k", "5").parse().unwrap(); let trials: usize = get("--trials", "20").parse().unwrap();
    let cap: usize = get("--cap", "1500000").parse().unwrap(); let seed: u64 = get("--seed", "1").parse().unwrap();
    let want_r: usize = get("--r", "1").parse().unwrap(); let verbose = args.iter().any(|a| a == "--verbose");
    let n = 3 * k; let mut rng = Rng(seed);
    for trial in 0..trials {
        let Some(g0) = gen_tuple(k, &mut rng, want_r) else { println!("no tuple"); continue; };
        // full orbit with original labels
        let mut seen: HashSet<Vec<u8>> = HashSet::new(); seen.insert(canonical(&g0, n));
        let mut reps = vec![g0]; let mut h = 0; let mut capped = false; let mut allsmooth = true;
        while h < reps.len() { let u = reps[h]; h += 1;
            let z = aindex(&acomp(&u[0], &u[1], n)[..n]);
            if z != aindex(&u[0][..n]) + aindex(&u[1][..n]) && z != aindex(&u[2][..n]) + aindex(&u[3][..n]) { allsmooth = false; }
            for m in 0..6 { let mut w = u; mv(&mut w, m, n); let key = canonical(&w, n); if seen.insert(key) { reps.push(w); } }
            if reps.len() > cap { capped = true; break; } }
        if capped { println!("trial {trial}: CAPPED > {cap}"); continue; }
        // observed frames, keyed by tau cycle type
        let mut fobs: HashMap<Vec<usize>, HashSet<Vec<u8>>> = HashMap::new();
        let mut samples: HashMap<Vec<usize>, Vec<(P, P, usize, usize)>> = HashMap::new();
        let mut obs_nonadd = false;
        for (si, u) in reps.iter().enumerate() { for fr in frames_of(u, k) {
            let ct = ctype(&fr.tau);
            if !additive(&fr.sig, &fr.tau) { obs_nonadd = true; }
            fobs.entry(ct.clone()).or_default().insert(canon2(&fr.sig, &fr.tau));
            let e = samples.entry(ct).or_default(); if e.len() < 12 && (si * 7919 + fr.gi) % 5 == 0 || e.is_empty() { e.push((fr.sig, fr.tau, fr.j, si)); }
        } }
        // G_E from the initial state (block 0) for order/primitivity; per sample its own block
        let ge0 = ge_gens(&g0, 0, k); let ge_el = group_closure(&ge0, k, 50000); let prim = primitive(&ge0, k);
        let mut full_ok = true; let mut full_nonadd = false; let mut nfull = 0; let mut nobs = 0;
        let mut detail = String::new();
        for (ct, ss) in samples.iter() {
            let mut fu: HashSet<Vec<u8>> = HashSet::new();
            for (sig, tau, j, si) in ss.iter() {
                let el = group_closure(&ge_gens(&reps[*si], *j, k), k, 50000);
                for x in el.iter() { let sx = comp(x, &comp(sig, &inv(x))); if !additive(&sx, tau) { full_nonadd = true; } fu.insert(canon2(&sx, tau)); }
            }
            let ob = &fobs[ct]; let miss = fu.iter().filter(|c| !ob.contains(*c)).count();
            if miss > 0 && verbose {
                for (sig, tau, j, si) in ss.iter().take(1) {
                    let el = group_closure(&ge_gens(&reps[*si], *j, k), k, 50000); let mut shown: HashSet<Vec<u8>> = HashSet::new();
                    for x in el.iter() { let sx = comp(x, &comp(sig, &inv(x))); let c = canon2(&sx, tau);
                        if shown.insert(c.clone()) { println!("    {} pair sig'={:?} tau={:?} ind(sig'tau)={} (sum {}) prod_type={:?}",
                            if ob.contains(&c) {"OBS "} else {"MISS"}, sx, tau, aindex(&comp(&sx, tau)), aindex(&sx)+aindex(tau), ctype(&comp(&sx, tau))); } }
                }
            }
            if miss > 0 { full_ok = false; }
            nfull += fu.len(); nobs += ob.len();
            detail += &format!(" tau{:?}:obs{}/full{}/miss{}", ct, ob.len(), fu.len(), miss);
        }
        println!("trial {trial}: k{k} sig{:?} N={} |G_E|={} prim={} allsmooth={} frames obs={} full={} FULL_SUBSET_OBS={} obs_nonadd={} full_nonadd={} |{}",
            ctype(&restrict(&g0[0], 0, k)), reps.len(), ge_el.len(), prim, allsmooth, nobs, nfull, full_ok, obs_nonadd, full_nonadd, detail);
    }
}
