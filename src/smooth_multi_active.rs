// Multi-active-point unramified lifts (evidence/multi_active_20261105/MATH.md).
//
// Base: transitive (x1, x2, x3) on m points, x1 x2 x3 = 1, group H = <x1, x2>; as a
// 4-tuple (1, x1, x2, x3). Source on n = m k points (beta, s) -> beta*k + s:
//   b(beta, s) = (x1 beta, s)                (gauge: b has trivial sheet labels)
//   c(beta, s) = (x2 beta, nu(beta) s)       (unramified: trivial return on x2-cycles)
//   a(beta, s) = (beta, sigma_beta s)        (sigma_beta != 1 exactly on the active set S)
//   e = (a b c)^-1                           (unramified: trivial return on y-cycles, y = x1 x2)
// e^-1 = abc has base permutation y and sheet label sigma_{y beta} nu(beta).
// Theorem (MATH.md): if S is independent in Gamma = pairs lying in one cycle of some
// H-conjugate of some x_i, the whole braid orbit is smooth. The census certifies
// every component by literal six-move closure anyway; a non-smooth result for an
// independent S is reported as THEOREM_FAILURE.
// AtomicUsize comes from smooth_star_census.rs

fn ma_perm_from(img: &[u8], n: usize) -> Perm {
    let mut p = Perm::identity(n);
    for i in 0..n { p.img[i] = img[i]; }
    p
}

fn ma_cycles(p: Perm, n: usize) -> Vec<Vec<u8>> {
    let mut seen = [false; MAX_D];
    let mut out = Vec::new();
    for i in 0..n {
        if seen[i] { continue; }
        let mut c = Vec::new();
        let mut x = i;
        while !seen[x] { seen[x] = true; c.push(x as u8); x = p.img[x] as usize; }
        out.push(c);
    }
    out
}

fn ma_transitive(gens: &[Perm], n: usize) -> bool {
    let mut par: Vec<usize> = (0..n).collect();
    fn f(par: &mut [usize], mut x: usize) -> usize { while par[x] != x { par[x] = par[par[x]]; x = par[x]; } x }
    let mut comps = n;
    for g in gens {
        for x in 0..n {
            let (u, v) = (f(&mut par, x), f(&mut par, g.img[x] as usize));
            if u != v { par[u] = v; comps -= 1; }
        }
    }
    comps == 1
}

fn ma_pair(u: usize, v: usize, m: usize) -> u32 { let (a, b) = if u < v { (u, v) } else { (v, u) }; (a * m + b) as u32 }

// Gamma as an adjacency bitmask per point; returns None if Gamma is complete.
fn ma_gamma(x: &[Perm; 3], m: usize) -> Option<[u16; 16]> {
    let all: u32 = (m * (m - 1) / 2) as u32;
    let mut have = vec![false; m * m];
    let mut stack: Vec<(u8, u8)> = Vec::new();
    let mut count = 0u32;
    for p in x {
        for c in ma_cycles(*p, m) {
            for i in 0..c.len() { for j in i + 1..c.len() {
                let id = ma_pair(c[i] as usize, c[j] as usize, m) as usize;
                if !have[id] { have[id] = true; count += 1; stack.push((c[i], c[j])); }
            } }
        }
    }
    let gens = [x[0], x[1], x[0].inverse(m), x[1].inverse(m)];
    while let Some((u, v)) = stack.pop() {
        if count == all { return None; }
        for g in &gens {
            let (a, b) = (g.img[u as usize], g.img[v as usize]);
            let id = ma_pair(a as usize, b as usize, m) as usize;
            if !have[id] { have[id] = true; count += 1; stack.push((a, b)); }
        }
    }
    if count == all { return None; }
    let mut adj = [0u16; 16];
    for u in 0..m { for v in 0..m { if u != v && have[ma_pair(u, v, m) as usize] { adj[u] |= 1 << v; } } }
    Some(adj)
}

// canonical key of a triple under simultaneous S_m relabelling (rooted traversal)
fn ma_triple_key(x: &[Perm; 3], m: usize) -> u128 {
    let mut best = u128::MAX;
    for r in 0..m {
        let mut lab = [u8::MAX; 16];
        let mut order = [0u8; 16];
        lab[r] = 0; order[0] = r as u8;
        let (mut len, mut head) = (1usize, 0usize);
        while head < len {
            let y = order[head] as usize; head += 1;
            for p in x { let z = p.img[y] as usize; if lab[z] == u8::MAX { lab[z] = len as u8; order[len] = z as u8; len += 1; } }
        }
        let mut key = 0u128;
        for p in x { for i in 0..m { key = (key << 4) | lab[p.img[order[i] as usize] as usize] as u128; } }
        if key < best { best = key; }
    }
    best
}

fn ma_partitions(n: usize, max: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if n == 0 { out.push(cur.clone()); return; }
    for p in (1..=n.min(max)).rev() { cur.push(p); ma_partitions(n - p, p, cur, out); cur.pop(); }
}

struct MaBase { x: [Perm; 3], m: usize, h: usize, adj: [u16; 16], key: u128 }

// all transitive identity-free product-one triples of degree m with incomplete Gamma
fn ma_bases(m: usize, all_bases: bool) -> Vec<MaBase> {
    let mut parts = Vec::new();
    ma_partitions(m, m, &mut Vec::new(), &mut parts);
    let mut seen: HashSet<u128> = HashSet::new();
    let mut out = Vec::new();
    for lam in parts.iter().filter(|l| l[0] > 1) {
        let mut x1 = Perm::identity(m);
        let mut off = 0;
        for &l in lam { for i in 0..l { x1.img[off + i] = (off + (i + 1) % l) as u8; } off += l; }
        // Heap's algorithm over all x2
        let mut a: Vec<u8> = (0..m as u8).collect();
        let mut c = vec![0usize; m];
        let mut visit = |a: &[u8]| {
            let x2 = ma_perm_from(a, m);
            if x2.is_identity(m) { return; }
            let x3 = x1.compose(x2, m).inverse(m);
            if x3.is_identity(m) || !ma_transitive(&[x1, x2], m) { return; }
            let x = [x1, x2, x3];
            let gamma = ma_gamma(&x, m).or_else(|| if all_bases {
                let mut full = [0u16; 16];
                for u in 0..m { full[u] = (((1u32 << m) - 1) as u16) & !(1u16 << u); }
                Some(full)
            } else { None });
            if let Some(adj) = gamma {
                let key = ma_triple_key(&x, m);
                if seen.insert(key) {
                    let ind: usize = x.iter().map(|p| sm_index(*p, m)).sum();
                    out.push(MaBase { x, m, h: (ind + 2 - 2 * m) / 2, adj, key });
                }
            }
        };
        visit(&a);
        let mut i = 0;
        while i < m {
            if c[i] < i {
                if i % 2 == 0 { a.swap(0, i); } else { a.swap(c[i], i); }
                visit(&a);
                c[i] += 1; i = 0;
            } else { c[i] = 0; i += 1; }
        }
    }
    out
}

// active sets up to the action of H: independent in Gamma (or any, if `any`)
fn ma_active_sets(b: &MaBase, lo: usize, hi: usize, any: bool) -> Vec<u16> {
    let m = b.m;
    let gens = [b.x[0], b.x[1]];
    let img = |s: u16, g: &Perm| -> u16 { let mut t = 0u16; for i in 0..m { if s >> i & 1 == 1 { t |= 1 << g.img[i]; } } t };
    let mut out = Vec::new();
    for s in 1u16..(1u16 << m) {
        let k = s.count_ones() as usize;
        if k < lo || k > hi { continue; }
        if !any && (0..m).any(|i| s >> i & 1 == 1 && b.adj[i] & s != 0) { continue; }
        // canonical = minimum of the H-orbit of s
        let mut orbit = HashSet::from([s]);
        let mut st = vec![s];
        let mut minimal = true;
        while let Some(t) = st.pop() {
            for g in &gens { let u = img(t, g); if u < s { minimal = false; break; } if orbit.insert(u) { st.push(u); } }
            if !minimal { break; }
        }
        if minimal { out.push(s); }
    }
    out
}

fn ma_sk(k: usize) -> Vec<[u8; 8]> {
    let mut out = Vec::new();
    let mut a: Vec<u8> = (0..k as u8).collect();
    fn rec(a: &mut Vec<u8>, i: usize, out: &mut Vec<[u8; 8]>) {
        if i == a.len() { let mut p = [0u8; 8]; for j in 0..a.len() { p[j] = a[j]; } out.push(p); return; }
        for j in i..a.len() { a.swap(i, j); rec(a, i + 1, out); a.swap(i, j); }
    }
    rec(&mut a, 0, &mut out);
    out.sort();
    out
}
#[inline] fn sk_mul(p: &[u8; 8], q: &[u8; 8], k: usize) -> [u8; 8] { let mut r = [0u8; 8]; for i in 0..k { r[i] = p[q[i] as usize]; } r }
#[inline] fn sk_inv(p: &[u8; 8], k: usize) -> [u8; 8] { let mut r = [0u8; 8]; for i in 0..k { r[p[i] as usize] = i as u8; } r }
fn sk_id(k: usize) -> [u8; 8] { let mut r = [0u8; 8]; for i in 0..k { r[i] = i as u8; } r }

// The voltage system: unknowns nu(beta); one constraint per x2-cycle (return of c)
// and per y-cycle (return of abc with label sigma_{y beta} nu(beta)). Constraint graph:
// vertices = x2-cycles then y-cycles, edge beta joins its two cycles. Spanning tree from
// vertex 0; non-tree edges are free; tree edges are solved leaf-first; the root's
// equation is checked.
struct MaSystem {
    m: usize,
    cyc: Vec<Vec<u8>>,       // constraint cycles (points in cycle order)
    is_y: Vec<bool>,
    vert_of: [[usize; 16]; 2], // [0][beta] = x2-cycle vertex, [1][beta] = y-cycle vertex
    free: Vec<u8>,           // non-tree edges (points)
    order: Vec<(usize, u8)>, // (vertex, parent edge) in leaf-first order
    root: usize,
    ybase: Perm,
}

fn ma_system(x1: Perm, x2: Perm, m: usize) -> MaSystem {
    let y = x1.compose(x2, m);
    let mut cyc = Vec::new();
    let mut is_y = Vec::new();
    let mut vert_of = [[0usize; 16]; 2];
    for (t, p) in [x2, y].iter().enumerate() {
        for c in ma_cycles(*p, m) {
            for &b in &c { vert_of[t][b as usize] = cyc.len(); }
            cyc.push(c); is_y.push(t == 1);
        }
    }
    let v = cyc.len();
    // BFS spanning tree
    let mut seen = vec![false; v];
    let mut tree_edge = [false; 16];
    let mut parent_edge = vec![u8::MAX; v];
    let mut bfs = vec![0usize];
    seen[0] = true;
    let mut h = 0;
    while h < bfs.len() {
        let u = bfs[h]; h += 1;
        for &b in &cyc[u] {
            let w = if is_y[u] { vert_of[0][b as usize] } else { vert_of[1][b as usize] };
            if !seen[w] { seen[w] = true; tree_edge[b as usize] = true; parent_edge[w] = b; bfs.push(w); }
        }
    }
    assert!(seen.iter().all(|&s| s), "constraint graph connected (H transitive)");
    let free: Vec<u8> = (0..m as u8).filter(|&b| !tree_edge[b as usize]).collect();
    let order: Vec<(usize, u8)> = bfs.iter().rev().filter(|&&u| u != 0).map(|&u| (u, parent_edge[u])).collect();
    MaSystem { m, cyc, is_y, vert_of, free, order, root: 0, ybase: y }
}

impl MaSystem {
    // label of edge beta inside constraint u: c-label nu, or abc-label sigma_{y beta} nu
    fn label(&self, u: usize, b: u8, nu: &[[u8; 8]; 16], sig: &[[u8; 8]; 16], k: usize) -> [u8; 8] {
        if self.is_y[u] { sk_mul(&sig[self.ybase.img[b as usize] as usize], &nu[b as usize], k) } else { nu[b as usize] }
    }
    // return product L = w_{l-1} ... w_0 of constraint u
    fn ret(&self, u: usize, nu: &[[u8; 8]; 16], sig: &[[u8; 8]; 16], k: usize) -> [u8; 8] {
        let mut r = sk_id(k);
        for &b in &self.cyc[u] { r = sk_mul(&self.label(u, b, nu, sig, k), &r, k); }
        r
    }
    // solve the one unknown edge b of constraint u
    fn solve(&self, u: usize, b: u8, nu: &mut [[u8; 8]; 16], sig: &[[u8; 8]; 16], k: usize) {
        let c = &self.cyc[u];
        let i = c.iter().position(|&z| z == b).unwrap();
        let mut rr = sk_id(k); // R = w_{i-1} ... w_0
        for &z in &c[..i] { rr = sk_mul(&self.label(u, z, nu, sig, k), &rr, k); }
        let mut ll = sk_id(k); // L = w_{l-1} ... w_{i+1}
        for &z in &c[i + 1..] { ll = sk_mul(&self.label(u, z, nu, sig, k), &ll, k); }
        let w = sk_mul(&sk_inv(&ll, k), &sk_inv(&rr, k), k); // L w R = 1
        nu[b as usize] = if self.is_y[u] { sk_mul(&sk_inv(&sig[self.ybase.img[b as usize] as usize], k), &w, k) } else { w };
    }
}

struct MaFound {
    base: usize,
    s: u16,
    k: usize,
    sig_types: String,
    t: [Perm; 4],
}

fn ma_build(b: &MaBase, k: usize, nu: &[[u8; 8]; 16], sig: &[[u8; 8]; 16]) -> [Perm; 4] {
    let (m, n) = (b.m, b.m * k);
    let (mut pa, mut pb, mut pc) = (Perm::identity(n), Perm::identity(n), Perm::identity(n));
    for be in 0..m {
        for s in 0..k {
            let p = be * k + s;
            pa.img[p] = (be * k + sig[be][s] as usize) as u8;
            pb.img[p] = (b.x[0].img[be] as usize * k + s) as u8;
            pc.img[p] = (b.x[1].img[be] as usize * k + nu[be][s] as usize) as u8;
        }
    }
    let pe = pa.compose(pb, n).compose(pc, n).inverse(n);
    [pa, pb, pc, pe]
}

// returns of y over each cycle of the block action are trivial (blocks = fibres of `blk`)
fn ma_unramified(y: Perm, n: usize, blk: &[u8], nb: usize) -> bool {
    let mut seen = vec![false; nb];
    let mut rep = vec![usize::MAX; nb];
    for x in 0..n { if rep[blk[x] as usize] == usize::MAX { rep[blk[x] as usize] = x; } }
    for b0 in 0..nb {
        if seen[b0] { continue; }
        // walk the block cycle starting at block b0 from each of its points
        let mut len = 0;
        let mut bb = b0;
        loop { seen[bb] = true; len += 1; bb = blk[y.img[rep[bb]] as usize] as usize; if bb == b0 { break; } }
        for x in 0..n {
            if blk[x] as usize != b0 { continue; }
            let mut z = x;
            for _ in 0..len { z = y.img[z] as usize; }
            if z != x { return false; }
        }
    }
    true
}

// all block systems (as block labelings) of <t> on n points
fn ma_block_systems(t: &[Perm; 4], n: usize) -> Vec<(Vec<u8>, usize)> {
    let gens: Vec<Perm> = t[..3].to_vec();
    let closure = |seed: u64| -> u64 {
        // minimal block containing the set `seed` (Atkinson)
        let mut par: Vec<usize> = (0..n).collect();
        fn f(par: &mut [usize], mut x: usize) -> usize { while par[x] != x { par[x] = par[par[x]]; x = par[x]; } x }
        let first = seed.trailing_zeros() as usize;
        for x in 0..n { if seed >> x & 1 == 1 { let (u, v) = (f(&mut par, x), f(&mut par, first)); par[u] = v; } }
        loop {
            let mut changed = false;
            for g in &gens {
                for x in 0..n {
                    let r = f(&mut par, x);
                    if r == x { continue; }
                    let (u, v) = (f(&mut par, g.img[x] as usize), f(&mut par, g.img[r] as usize));
                    if u != v { par[u] = v; changed = true; }
                }
            }
            if !changed { break; }
        }
        let r0 = f(&mut par, first);
        let mut mask = 0u64;
        for x in 0..n { if f(&mut par, x) == r0 { mask |= 1 << x; } }
        mask
    };
    let full: u64 = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
    let mut blocks: Vec<u64> = Vec::new();
    let mut st = vec![1u64];
    let mut seen = HashSet::from([1u64]);
    while let Some(bk) = st.pop() {
        for x in 0..n {
            if bk >> x & 1 == 1 { continue; }
            let nb = closure(bk | 1 << x);
            if nb != full && seen.insert(nb) { blocks.push(nb); st.push(nb); }
        }
    }
    let mut out = Vec::new();
    for bk in blocks {
        // labeling: images of bk under the group = blocks of the system
        let size = bk.count_ones() as usize;
        let mut lab = vec![u8::MAX; n];
        let mut nbk = 0usize;
        let mut queue = vec![bk];
        let mut seenb = HashSet::from([bk]);
        while let Some(c) = queue.pop() {
            for x in 0..n { if c >> x & 1 == 1 { lab[x] = nbk as u8; } }
            nbk += 1;
            for g in &gens {
                let mut im = 0u64;
                for x in 0..n { if c >> x & 1 == 1 { im |= 1 << g.img[x]; } }
                if seenb.insert(im) { queue.push(im); }
            }
        }
        debug_assert_eq!(nbk * size, n);
        out.push((lab, nbk));
    }
    out
}

pub fn multi_active_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: &str| args.windows(2).find(|w| w[0] == name).map(|w| w[1].clone()).unwrap_or(default.into());
    let ms: Vec<usize> = arg("--m", "6").split(',').map(|x| x.parse().map_err(|_| "m")).collect::<Result<_, _>>()?;
    let ks: Vec<usize> = arg("--k", "2").split(',').map(|x| x.parse().map_err(|_| "k")).collect::<Result<_, _>>()?;
    let lo: usize = arg("--min-active", "2").parse().map_err(|_| "min-active")?;
    let hi: usize = arg("--max-active", "16").parse().map_err(|_| "max-active")?;
    let any = args.iter().any(|x| x == "--any-active");
    let cap: usize = arg("--cap", "2000000").parse().map_err(|_| "cap")?;
    let secs: u64 = arg("--closure-seconds", "600").parse().map_err(|_| "secs")?;
    let assign_cap: u64 = arg("--assign-cap", "2000000").parse().map_err(|_| "assign-cap")?;
    let workers: usize = arg("--workers", "2").parse().map_err(|_| "workers")?;
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() { return Err("--out required".into()); }
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let mut summary = BufWriter::new(OpenOptions::new().create(true).append(true).open(out.join("COMPONENTS.tsv")).map_err(|e| e.to_string())?);
    writeln!(summary, "m\tk\tbase_h\tactive\tsigma_types\tn\tgenus\tstates\tstatus\tprimitive\tblock_systems\tclass\ta\tb\tc\te").map_err(|e| e.to_string())?;
    for &m in &ms {
        let t0 = Instant::now();
        let bases = ma_bases(m, args.iter().any(|x| x == "--all-bases"));
        println!("m={m} bases_with_incomplete_gamma={} seconds={:.1}", bases.len(), t0.elapsed().as_secs_f64());
        let mut bl = BufWriter::new(File::create(out.join(format!("BASES_m{m}.tsv"))).map_err(|e| e.to_string())?);
        writeln!(bl, "base\th\tx1\tx2\tx3\tgamma_edges").map_err(|e| e.to_string())?;
        for (i, b) in bases.iter().enumerate() {
            let edges: u32 = (0..m).map(|u| b.adj[u].count_ones()).sum::<u32>() / 2;
            writeln!(bl, "{i}\t{}\t{}\t{}\t{}\t{edges}", b.h, b.x[0].cycle_string(m), b.x[1].cycle_string(m), b.x[2].cycle_string(m)).map_err(|e| e.to_string())?;
        }
        bl.flush().map_err(|e| e.to_string())?;
        for &k in &ks {
            let n = m * k;
            if n > 25 || k > 8 || m > 16 { println!("skip m={m} k={k} (n={n} > 25)"); continue; }
            let sk = ma_sk(k);
            let nonid: Vec<[u8; 8]> = sk.iter().copied().filter(|p| *p != sk_id(k)).collect();
            let known: Mutex<HashSet<SmoothKey>> = Mutex::new(HashSet::new());
            let found: Mutex<Vec<(MaFound, SmoothScan)>> = Mutex::new(Vec::new());
            let next = AtomicUsize::new(0);
            let stats = Mutex::new([0u64; 6]); // candidates, transitive, smooth_start, closures, theorem_failures, assign_capped
            std::thread::scope(|scope| {
                for _ in 0..workers {
                    scope.spawn(|| {
                        loop {
                            let bi = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if bi >= bases.len() { break; }
                            let b = &bases[bi];
                            let sys = ma_system(b.x[0], b.x[1], m);
                            for s in ma_active_sets(b, lo, hi, any) {
                                let act: Vec<usize> = (0..m).filter(|&i| s >> i & 1 == 1).collect();
                                let independent = act.iter().all(|&i| b.adj[i] & s == 0);
                                let total = (nonid.len() as u64).pow(act.len() as u32) * (sk.len() as u64).pow(sys.free.len() as u32);
                                let mut local = [0u64; 6];
                                if total > assign_cap { local[5] += 1; let mut st = stats.lock().unwrap(); for j in 0..6 { st[j] += local[j]; } continue; }
                                let mut sig = [sk_id(k); 16];
                                let mut nu = [sk_id(k); 16];
                                let mut idx = vec![0usize; act.len() + sys.free.len()];
                                'assign: loop {
                                    for (j, &be) in act.iter().enumerate() { sig[be] = nonid[idx[j]]; }
                                    for (j, &be) in sys.free.iter().enumerate() { nu[be as usize] = sk[idx[act.len() + j]]; }
                                    for &(u, be) in &sys.order { sys.solve(u, be, &mut nu, &sig, k); }
                                    local[0] += 1;
                                    if sys.ret(sys.root, &nu, &sig, k) == sk_id(k) {
                                        let t = ma_build(b, k, &nu, &sig);
                                        if t.iter().all(|p| !p.is_identity(n)) && ma_transitive(&t[..3], n) {
                                            local[1] += 1;
                                            let v = sm_values(&t, n);
                                            if v.0 {
                                                local[2] += 1;
                                                let key = sm_key(&t, n);
                                                let new = !known.lock().unwrap().contains(&key);
                                                if new {
                                                    let scan = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(secs));
                                                    local[3] += 1;
                                                    if scan.status == "nonsmooth" && independent { local[4] += 1; }
                                                    let mut kn = known.lock().unwrap();
                                                    if !kn.contains(&key) {
                                                        for st in &scan.states { kn.insert(*st); }
                                                        drop(kn);
                                                        let types = act.iter().map(|&be| {
                                                            let mut c: Vec<usize> = ma_cycles(ma_perm_from(&sig[be][..k], k), k).iter().map(|x| x.len()).filter(|&l| l > 1).collect();
                                                            c.sort_unstable_by(|x, y| y.cmp(x));
                                                            c.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(".")
                                                        }).collect::<Vec<_>>().join("|");
                                                        found.lock().unwrap().push((MaFound { base: bi, s, k, sig_types: types, t }, scan));
                                                    }
                                                }
                                            } else if independent { local[4] += 1; }
                                        }
                                    }
                                    // odometer
                                    let mut j = 0;
                                    loop {
                                        if j == idx.len() { break 'assign; }
                                        idx[j] += 1;
                                        let lim = if j < act.len() { nonid.len() } else { sk.len() };
                                        if idx[j] < lim { break; }
                                        idx[j] = 0; j += 1;
                                    }
                                }
                                let mut st = stats.lock().unwrap();
                                for j in 0..6 { st[j] += local[j]; }
                            }
                        }
                    });
                }
            });
            let st = *stats.lock().unwrap();
            let found = found.into_inner().unwrap();
            println!("m={m} k={k} n={n} candidates={} transitive={} smooth_start={} closures={} theorem_failures={} assign_capped_sets={} components={}",
                st[0], st[1], st[2], st[3], st[4], st[5], found.len());
            for (fi, (f, scan)) in found.iter().enumerate() {
                let b = &bases[f.base];
                let t = f.t;
                let g = sm_values(&t, n).2;
                let systems = ma_block_systems(&t, n);
                let mut desc = Vec::new();
                let mut localized = false;
                let mut fibres = false;
                for (lab, nb) in &systems {
                    let qa_trivial = (0..n).all(|x| lab[t[0].img[x] as usize] == lab[x]);
                    let mut active = vec![false; *nb];
                    for x in 0..n { if t[0].img[x] as usize != x { active[lab[x] as usize] = true; } }
                    let na = active.iter().filter(|&&z| z).count();
                    let unram = qa_trivial && (1..4).all(|i| ma_unramified(t[i], n, lab, *nb));
                    if qa_trivial && na == 1 && unram { localized = true; }
                    if *nb == m && (0..n).all(|x| lab[x] == lab[(x / k) * k]) { fibres = true; }
                    desc.push(format!("{}x{}:{}{}{}", nb, n / nb, if qa_trivial { "a1" } else { "a*" }, na, if unram { "u" } else { "r" }));
                }
                desc.sort();
                let class = if scan.status != "closed_smooth" { "not_all_smooth" } else if localized { "localized_old" } else { "NEW_outside_localized" };
                let prim = systems.is_empty();
                writeln!(summary, "{m}\t{}\t{}\t{}\t{}\t{n}\t{g}\t{}\t{}\t{prim}\t{}\t{class}\t{}\t{}\t{}\t{}", f.k, b.h, f.s.count_ones(), f.sig_types,
                    scan.states.len(), scan.status, desc.join(","), t[0].cycle_string(n), t[1].cycle_string(n), t[2].cycle_string(n), t[3].cycle_string(n)).map_err(|e| e.to_string())?;
                if class != "localized_old" && scan.status == "closed_smooth" {
                    let dir = out.join(format!("new_m{m}_k{k}_{fi:04}"));
                    if !dir.exists() { sm_save(&t, n, scan, &dir)?; }
                }
                let _ = fibres;
            }
            summary.flush().map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
