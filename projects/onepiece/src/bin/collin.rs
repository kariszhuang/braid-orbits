//! Search for all-zero braid orbits among 4-tuples whose first three entries are central
//! collineations (homologies or elations) of PG(2,q); the fourth entry g4 = (g1g2g3)^-1 may be
//! anything. g1 runs over one representative per PGL(3,q)-conjugacy type (tuples are taken up
//! to conjugation, and canonical forms are modulo S_n anyway), g2 and g3 over all central
//! collineations. Any orbit containing at least three central collineations has such a state.
//! Input: coll_q{q}.txt from make_collineations.py. Points up to 128.
//! Conventions as src/main.rs: rightmost-first composition, g1 g2 g3 g4 = 1.
use std::collections::HashSet;
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

const NM: usize = 128;
type P = [u8; NM];
type T = [P; 4];

#[derive(Clone, Copy)]
struct Cx { n: usize }
impl Cx {
    fn id(&self) -> P { let mut p = [0u8; NM]; for i in 0..NM { p[i] = i as u8; } p }
    #[inline] fn mul(&self, a: &P, b: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[i] = a[b[i] as usize]; } c }
    #[inline] fn inv(&self, a: &P) -> P { let mut c = self.id(); for i in 0..self.n { c[a[i] as usize] = i as u8; } c }
    #[inline] fn is_id(&self, a: &P) -> bool { (0..self.n).all(|i| a[i] as usize == i) }
    #[inline] fn ncyc(&self, a: &P) -> usize { let mut seen: u128 = 0; let mut k = 0;
        for s in 0..self.n { if seen >> s & 1 == 1 { continue; } k += 1; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u128 << x; x = a[x] as usize; } } k }
    #[inline] fn ind(&self, a: &P) -> usize { self.n - self.ncyc(a) }
    #[inline] fn norb(&self, gens: &[&P]) -> usize {
        let mut seen: u128 = 0; let mut k = 0; let mut st = [0u8; NM];
        for s in 0..self.n { if seen >> s & 1 == 1 { continue; } k += 1; seen |= 1u128 << s; let mut sp = 1; st[0] = s as u8;
            while sp > 0 { sp -= 1; let x = st[sp] as usize; for g in gens { let y = g[x] as usize; if seen >> y & 1 == 0 { seen |= 1u128 << y; st[sp] = y as u8; sp += 1; } } } }
        k
    }
    fn labels(&self, gens: &[&P]) -> ([u8; NM], usize) {
        let mut lab = [u8::MAX; NM]; let mut k = 0u8; let mut st = [0u8; NM];
        for s in 0..self.n { if lab[s] != u8::MAX { continue; } lab[s] = k; let mut sp = 1; st[0] = s as u8;
            while sp > 0 { sp -= 1; let x = st[sp] as usize; for g in gens { let y = g[x] as usize; if lab[y] == u8::MAX { lab[y] = k; st[sp] = y as u8; sp += 1; } } }
            k += 1; }
        (lab, k as usize)
    }
    fn fmt(&self, p: &P) -> String {
        let mut seen = [false; NM]; let mut s = String::new();
        for i in 0..self.n { if seen[i] { continue; } let mut c = vec![]; let mut x = i;
            while !seen[x] { seen[x] = true; c.push(x + 1); x = p[x] as usize; }
            if c.len() > 1 { s.push('('); s.push_str(&c.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ")); s.push(')'); } }
        if s.is_empty() { "()".into() } else { s }
    }
}
#[inline] fn defect(cx: &Cx, t: &T) -> i64 { let p = cx.mul(&t[0], &t[1]); cx.ncyc(&p) as i64 - cx.norb(&[&t[0], &t[1]]) as i64 - cx.norb(&[&t[2], &t[3]]) as i64 + 1 }
/// one-piece at {12|34}: defect 0 and at most one positive-genus piece
fn one_piece(cx: &Cx, t: &T) -> bool {
    if defect(cx, t) != 0 { return false; }
    let p = cx.mul(&t[0], &t[1]); let pi = cx.inv(&p); let mut pos = 0;
    for (a, b, q) in [(&t[0], &t[1], &p), (&t[2], &t[3], &pi)] {
        let (lab, k) = cx.labels(&[a, b]);
        let mut c = vec![0i64; k]; for i in 0..cx.n { c[lab[i] as usize] += 1; } for v in 0..k { c[v] += 2; }
        for g in [a, b, q] { let mut seen: u128 = 0; for s in 0..cx.n { if seen >> s & 1 == 1 { continue; } c[lab[s] as usize] -= 1; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u128 << x; x = g[x] as usize; } } }
        pos += c.iter().filter(|&&v| v > 0).count();
    }
    pos <= 1
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
    for j in 0..4 { let a = &t[j]; let mut seen: u128 = 0;
        for s in 0..n { if seen >> s & 1 == 1 { continue; } let mut l = 0u8; let mut x = s; while seen >> x & 1 == 0 { seen |= 1u128 << x; x = a[x] as usize; l += 1; }
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
fn unpack(cx: &Cx, k: &[u8]) -> T { let n = cx.n; let mut t = [cx.id(); 4]; for j in 0..3 { for i in 0..n { t[j][i] = k[j * n + i]; } } t[3] = cx.inv(&cx.mul(&t[0], &cx.mul(&t[1], &t[2]))); t }
fn h128(k: &[u8]) -> u128 { let (mut a, mut b) = (0xcbf2_9ce4_8422_2325u64, 0x9e37_79b9_7f4a_7c15u64);
    for &x in k { a = (a ^ x as u64).wrapping_mul(0x0000_0100_0000_01b3); b = (b.rotate_left(7) ^ x as u64).wrapping_mul(0x517c_c1b7_2722_0a95); } ((a as u128) << 64) | (b as u128) }
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
/// closure with early exit on positive defect; returns (all_zero, state hashes, states, one-piece states)
fn close(cx: &Cx, t0: &T, cap: usize) -> (bool, Vec<u128>, usize, usize) {
    let root = canon(cx, t0); let mut idx: HashSet<u128> = HashSet::new(); idx.insert(h128(&root));
    let mut st = vec![root]; let mut i = 0; let mut onep = 0;
    while i < st.len() {
        let u = unpack(cx, &st[i]); i += 1;
        if defect(cx, &u) != 0 { return (false, st.iter().map(|k| h128(k)).collect(), st.len(), 0); }
        if one_piece(cx, &u) { onep += 1; }
        for m in 0..6 { let v = mv(cx, &u, m); let k = canon(cx, &v); if idx.insert(h128(&k)) { st.push(k); if st.len() > cap { return (false, vec![], st.len(), 0); } } }
    }
    let n = st.len(); (true, st.iter().map(|k| h128(k)).collect(), n, onep)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    let file = get("--coll").expect("--coll coll_qQ.txt");
    let threads: usize = get("--threads").map(|x| x.parse().unwrap()).unwrap_or(8);
    let min_genus: usize = get("--min-genus").map(|x| x.parse().unwrap()).unwrap_or(1);
    let cap: usize = get("--cap").map(|x| x.parse().unwrap()).unwrap_or(1_000_000);
    let out = Mutex::new(std::fs::File::create(get("--out").expect("--out FILE")).unwrap());
    let txt = std::fs::read_to_string(&file).unwrap(); let mut lines = txt.lines();
    let hdr: Vec<usize> = lines.next().unwrap().split(' ').map(|x| x.parse().unwrap()).collect(); let (q, n) = (hdr[0], hdr[1]); assert!(n <= NM);
    let reps: Vec<usize> = lines.next().unwrap().split(' ').map(|x| x.parse().unwrap()).collect();
    let cx = Cx { n };
    let perms: Vec<P> = lines.map(|l| { let mut p = cx.id(); for (i, x) in l.split(' ').enumerate() { p[i] = x.parse().unwrap(); } p }).collect();
    let ind: Vec<usize> = perms.iter().map(|p| cx.ind(p)).collect();
    let t0 = std::time::Instant::now();
    eprintln!("q={} n={} central collineations {} types {}", q, n, perms.len(), reps.len());
    let np = perms.len(); let jobs: Vec<(usize, usize, usize)> = reps.iter().flat_map(|&r| (0..np).step_by(256).map(move |s| (r, s, (s + 256).min(np)))).collect();
    let next = AtomicUsize::new(0);
    let (tested, d0, pass3, probe_rej, closures, found) = (AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0));
    let done: Mutex<HashSet<u128>> = Mutex::new(HashSet::new());
    let perms_r = &perms; let ind_r = &ind;
    std::thread::scope(|s| { for _ in 0..threads { s.spawn(|| { loop {
        let j = next.fetch_add(1, Ordering::Relaxed); if j >= jobs.len() { break; }
        if j % 500 == 0 { eprintln!("job {}/{} {:.0}s tested {} d0 {} pass3 {} closures {} found {}", j, jobs.len(), t0.elapsed().as_secs_f64(), tested.load(Ordering::Relaxed), d0.load(Ordering::Relaxed), pass3.load(Ordering::Relaxed), closures.load(Ordering::Relaxed), found.load(Ordering::Relaxed)); }
        let (r, a, b) = jobs[j]; let g1 = perms_r[r]; let (mut lt, mut ld, mut lp) = (0u64, 0u64, 0u64);
        for k2 in a..b {
            let g2 = perms_r[k2]; let p12 = cx.mul(&g1, &g2); let c12 = cx.ncyc(&p12) as i64; let v12 = cx.norb(&[&g1, &g2]) as i64;
            for k3 in 0..perms_r.len() {
                let g3 = &perms_r[k3]; let qq = cx.mul(&p12, g3); let i4 = cx.ind(&qq); if i4 == 0 { continue; }
                if ind_r[r] + ind_r[k2] + ind_r[k3] + i4 < 2 * n - 2 + 2 * min_genus { continue; }
                lt += 1;
                if c12 - v12 - cx.norb(&[g3, &p12]) as i64 + 1 != 0 { continue; }
                if cx.norb(&[&g1, &g2, g3]) != 1 { continue; }
                ld += 1;
                let t: T = [g1, g2, *g3, cx.inv(&qq)];
                if defect(&cx, &rot(&t)) != 0 || defect(&cx, &mv(&cx, &t, 1)) != 0 || defect(&cx, &mv(&cx, &t, 4)) != 0 { continue; }
                lp += 1;
                { let mut rr = { let mut bb = Vec::with_capacity(3 * n); for x in &t[..3] { bb.extend_from_slice(&x[..n]); } h128(&bb) as u64 | 1 }; let mut u = t; let mut bad = false;
                  for _ in 0..64 { rr ^= rr << 13; rr ^= rr >> 7; rr ^= rr << 17; u = mv(&cx, &u, (rr % 6) as usize); if defect(&cx, &u) != 0 { bad = true; break; } }
                  if bad { probe_rej.fetch_add(1, Ordering::Relaxed); continue; } }
                let key = h128(&canon(&cx, &t));
                if done.lock().unwrap().contains(&key) { continue; }
                closures.fetch_add(1, Ordering::Relaxed);
                let (ok, hs, ns, onep) = close(&cx, &t, cap);
                { let mut dn = done.lock().unwrap(); for h in &hs { dn.insert(*h); } }
                if ns > cap { writeln!(out.lock().unwrap(), "{{\"q\":{},\"d\":{},\"capped\":true,\"tuple\":{:?}}}", q, n, t.iter().map(|p| cx.fmt(p)).collect::<Vec<_>>()).unwrap(); continue; }
                if !ok { continue; }
                found.fetch_add(1, Ordering::Relaxed);
                let tv: Vec<P> = t.to_vec(); let pr = primitive(&cx, &tv);
                let ab = (0..4).all(|a| (0..4).all(|b| cx.mul(&t[a], &t[b]) == cx.mul(&t[b], &t[a])));
                let genus = (t.iter().map(|p| cx.ind(p)).sum::<usize>() as i64 - 2 * n as i64 + 2) / 2;
                let line = format!("{{\"q\":{},\"d\":{},\"genus\":{},\"states\":{},\"one_piece_states\":{},\"primitive\":{},\"abelian\":{},\"tuple\":{:?}}}", q, n, genus, ns, onep, pr, ab, t.iter().map(|p| cx.fmt(p)).collect::<Vec<_>>());
                if pr && !ab { eprintln!("!!! PRIMITIVE NON-ABELIAN ALL-ZERO: {}", line); }
                writeln!(out.lock().unwrap(), "{}", line).unwrap();
            }
        }
        tested.fetch_add(lt, Ordering::Relaxed); d0.fetch_add(ld, Ordering::Relaxed); pass3.fetch_add(lp, Ordering::Relaxed);
    } }); } });
    let ld = |a: &AtomicU64| a.load(Ordering::Relaxed);
    let s = format!("{{\"summary\":true,\"q\":{},\"d\":{},\"central_collineations\":{},\"types\":{},\"min_genus\":{},\"tuples_tested\":{},\"defect0_first_cusp\":{},\"defect0_three_cusps\":{},\"probe_rejects\":{},\"closures\":{},\"allzero_orbits\":{},\"seconds\":{:.1},\"complete\":true}}",
        q, n, perms.len(), reps.len(), min_genus, ld(&tested), ld(&d0), ld(&pass3), ld(&probe_rej), ld(&closures), ld(&found), t0.elapsed().as_secs_f64());
    eprintln!("{}", s); writeln!(out.lock().unwrap(), "{}", s).unwrap();
}
