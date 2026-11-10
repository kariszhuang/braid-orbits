// Primitive-by-construction STAR sampler over many degrees.
//
// Each proposal is a STAR-regime start (smallest-index entry first, index-
// additive first pair, allowed sorted index list). It is kept only with a
// primitivity certificate: d prime, or some sampled group element has a cycle
// of prime length p > d/2 (any block system would force p | block count or
// p <= block size, both impossible). Then literal braid exploration decides it.
// A sampled negative is not coverage; closed and capped seeds are saved.

fn ps_is_prime(p: usize) -> bool { p >= 2 && (2..p).take_while(|q| q * q <= p).all(|q| p % q != 0) }

fn ps_cycles(p: Perm, n: usize) -> Vec<usize> {
    let mut seen = [false; MAX_D];
    let mut out = Vec::new();
    for i in 0..n {
        if seen[i] { continue; }
        let mut j = i; let mut l = 0;
        while !seen[j] { seen[j] = true; j = p.img[j] as usize; l += 1; }
        out.push(l);
    }
    out
}

fn ps_big_prime_cycle(p: Perm, n: usize) -> bool {
    ps_cycles(p, n).into_iter().any(|l| 2 * l > n && ps_is_prime(l))
}

// Uniform permutation of 0..n with exactly k cycles (Feller coupling).
fn ps_with_cycles(n: usize, k: usize, st: &[Vec<f64>], rng: &mut SplitMix64) -> Perm {
    let mut fresh = vec![false; n];
    let mut kk = k;
    for m in (1..=n).rev() {
        let p_new = if kk == 0 { 0.0 } else { st[m - 1][kk - 1] / st[m][kk] };
        let u = (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        if u < p_new { fresh[m - 1] = true; kk -= 1; }
    }
    let mut p = Perm::identity(n);
    for i in 0..n {
        if !fresh[i] {
            let j = rng.below(i);
            p.img[i] = p.img[j];
            p.img[j] = i as u8;
        }
    }
    // random relabel keeps uniformity and removes the Feller ordering bias in labels
    let mut lab: Vec<u8> = (0..n as u8).collect();
    for i in (1..n).rev() { let j = rng.below(i + 1); lab.swap(i, j); }
    let mut q = Perm::identity(n);
    for i in 0..n { q.img[lab[i] as usize] = lab[p.img[i] as usize]; }
    q
}

fn ps_stirling(nmax: usize) -> Vec<Vec<f64>> {
    let mut s = vec![vec![0.0f64; nmax + 2]; nmax + 2];
    s[0][0] = 1.0;
    for n in 1..=nmax { for k in 1..=n { s[n][k] = s[n - 1][k - 1] + (n - 1) as f64 * s[n - 1][k]; } }
    s
}

fn ps_patterns(d: usize) -> Vec<[usize; 4]> {
    let mut out = Vec::new();
    for g in 1..d {
        // genus 3 below degree 14 is already completely excluded (C209/C229/C230 and degree12)
        if g == 3 && d <= 13 { continue; }
        let s = 2 * d + 2 * g - 2;
        for a in 2..d { for b in a + 1..d { for c in b..d {
            if s < a + b + c { continue; }
            let e = s - a - b - c;
            if c <= e && e < d && (b - a) as i64 - (e - c) as i64 >= 2 * g as i64 { out.push([a, b, c, e]); }
        }}}
    }
    out
}

pub fn prim_sampler_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: &str| args.windows(2).find(|w| w[0] == name)
        .map(|w| w[1].clone()).unwrap_or(default.into());
    let degrees: Vec<usize> = arg("--degrees", "12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34")
        .split(',').map(|x| x.parse().map_err(|_| "degrees")).collect::<Result<_, _>>()?;
    let seed: u64 = arg("--seed", "1").parse().map_err(|_| "seed")?;
    let seconds: u64 = arg("--seconds", "31536000").parse().map_err(|_| "seconds")?;
    let cap: usize = arg("--cap", "300000").parse().map_err(|_| "cap")?;
    let closure_seconds: u64 = arg("--closure-seconds", "60").parse().map_err(|_| "closure")?;
    let no_cert = args.iter().any(|a| a == "--no-cert");
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() { return Err("--out required".into()); }
    if degrees.iter().any(|&d| d < 9 || d > MAX_D) { return Err("degree range".into()); }
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let st = ps_stirling(MAX_D);
    let pats: Vec<Vec<[usize; 4]>> = degrees.iter().map(|&d| ps_patterns(d)).collect();
    let mut rng = SplitMix64::new(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ 0x5EED);
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    // per degree: proposals, index_ok, transitive, primitive_cert, short_reject, scan_reject, closed, unresolved
    let mut stats = vec![[0u64; 8]; degrees.len()];
    let mut last = Instant::now();
    let mut seen_closed: HashMap<SmoothKey, ()> = HashMap::new();
    let mut log = BufWriter::new(OpenOptions::new().create(true).append(true)
        .open(out.join("PROGRESS.jsonl")).map_err(|e| e.to_string())?);
    let mut found = 0usize;
    while Instant::now() < deadline {
        let di = rng.below(degrees.len());
        let n = degrees[di];
        if pats[di].is_empty() { continue; }
        let q = pats[di][rng.below(pats[di].len())];
        stats[di][0] += 1;
        // position 1 = smallest; random order of the other three
        let mut rest = [q[1], q[2], q[3]];
        for i in (1..3).rev() { let j = rng.below(i + 1); rest.swap(i, j); }
        let (ia, ib, ic, id) = (q[0], rest[0], rest[1], rest[2]);
        let a = ps_with_cycles(n, n - ia, &st, &mut rng);
        // b = product of ib transpositions each merging two cycles of the running product
        let mut cur = a;
        let mut b = Perm::identity(n);
        let mut ok = true;
        for _ in 0..ib {
            let mut lab = [usize::MAX; MAX_D];
            let mut ncyc = 0;
            for i in 0..n { if lab[i] == usize::MAX { let mut j = i; while lab[j] == usize::MAX { lab[j] = ncyc; j = cur.img[j] as usize; } ncyc += 1; } }
            if ncyc < 2 { ok = false; break; }
            let x = rng.below(n);
            let mut y = rng.below(n);
            let mut tries = 0;
            while lab[y] == lab[x] { y = rng.below(n); tries += 1; if tries > 1000 { break; } }
            if lab[y] == lab[x] { ok = false; break; }
            let mut t = Perm::identity(n); t.img[x] = y as u8; t.img[y] = x as u8;
            cur = cur.compose(t, n);
            b = b.compose(t, n);
        }
        if !ok || sm_index(b, n) != ib { continue; }
        let p = a.compose(b, n);
        debug_assert_eq!(sm_index(p, n), ia + ib);
        // c with index ic whose fourth entry has index id: a few tries per pair
        let pinv = p.inverse(n);
        let mut hit = None;
        for _ in 0..64 {
            let c = ps_with_cycles(n, n - ic, &st, &mut rng);
            let d = c.inverse(n).compose(pinv, n);
            if sm_index(d, n) == id { hit = Some((c, d)); break; }
        }
        let Some((c, d)) = hit else { continue };
        stats[di][1] += 1;
        let t = [a, b, c, d];
        if !generators_transitive(n, &t[..3]) { continue; }
        stats[di][2] += 1;
        let cert = ps_is_prime(n) || t.iter().any(|&x| ps_big_prime_cycle(x, n))
            || [p, b.compose(c, n), a.compose(c, n), c.compose(d, n), a.compose(c, n).compose(b, n)]
                .iter().any(|&x| ps_big_prime_cycle(x, n));
        if !cert && !no_cert { continue; }
        stats[di][3] += 1;
        if sc_word(t, n).is_some() { stats[di][4] += 1; }
        else {
            let scan = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(closure_seconds));
            match scan.status {
                "nonsmooth" => stats[di][5] += 1,
                "closed_smooth" => {
                    stats[di][6] += 1;
                    let k = sm_key(&t, n);
                    if !seen_closed.contains_key(&k) {
                        sm_save(&t, n, &scan, &out.join(format!("CLOSED_d{n}_{found:04}")))?;
                        for s in &scan.states { seen_closed.insert(*s, ()); }
                        found += 1;
                        writeln!(log, "{{\"event\":\"CLOSED\",\"degree\":{n},\"pattern\":{:?},\"states\":{},\"tuple\":\"{}\"}}",
                            q, scan.states.len(), tuple_string_dynamic(n, &t)).map_err(|e| e.to_string())?;
                        log.flush().map_err(|e| e.to_string())?;
                    }
                }
                _ => {
                    stats[di][7] += 1;
                    sm_save(&t, n, &scan, &out.join(format!("UNRESOLVED_d{n}_{}", stats[di][7])))?;
                }
            }
        }
        if last.elapsed() > Duration::from_secs(60) {
            last = Instant::now();
            let rows: Vec<String> = degrees.iter().zip(&stats).map(|(d, s)| format!("\"{d}\":{:?}", s)).collect();
            writeln!(log, "{{\"event\":\"progress\",\"secs\":{:.0},\"cols\":\"proposals,index_ok,transitive,primitive_cert,short_reject,scan_reject,closed,unresolved\",\"by_degree\":{{{}}}}}",
                start.elapsed().as_secs_f64(), rows.join(",")).map_err(|e| e.to_string())?;
            log.flush().map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
