// Exact STAR-regime census for prescribed sorted index lists.
//
// In the STAR (low-index) regime every additive pair contains the unique
// smallest-index entry. Any orbit has a state with that entry in position 1,
// and there smoothness at the (1,2) cusp forces (pos1,pos2) to be additive.
// So it suffices to take every canonical index-additive first pair (a,b) with
// ind a = pattern minimum and every third entry c in S_n whose indices complete
// an allowed pattern. Each eligible transitive start gets a literal negative
// word, a full closure, or an explicit unresolved record.
use std::sync::atomic::AtomicUsize;

#[derive(Default)]
struct StarProofs {
    closed: HashMap<SmoothKey, usize>,
    bad: HashMap<SmoothKey, Vec<usize>>,
    components: usize,
}
// enumerated, index_eligible, transitive, nonsmooth, closed_occurrences, unresolved
const ST_COLUMNS: &str = "enumerated\tindex_eligible\ttransitive\tnonsmooth\tclosed_occurrences\tunresolved";

fn st_parse_patterns(raw: &str) -> Result<Vec<[usize; 4]>, String> {
    let mut out = Vec::new();
    for part in raw.split(';').filter(|s| !s.trim().is_empty()) {
        let v: Vec<usize> = part.split(',').map(|x| x.trim().parse().map_err(|_| "pattern"))
            .collect::<Result<_, _>>()?;
        if v.len() != 4 { return Err("pattern needs four indices".into()); }
        let mut s = [v[0], v[1], v[2], v[3]];
        s.sort_unstable();
        if s[0] == s[1] { return Err("STAR pattern needs a unique minimum".into()); }
        out.push(s);
    }
    if out.is_empty() { return Err("no patterns".into()); }
    Ok(out)
}

fn st_matches(patterns: &[[usize; 4]], ia: usize, ib: usize, ic: usize, id: usize) -> bool {
    let mut s = [ia, ib, ic, id];
    s.sort_unstable();
    patterns.iter().any(|p| *p == s && p[0] == ia)
}

fn st_frame(
    n: usize, fi: usize, key: [u64; 2], perms: &[(Perm, usize, Perm)], patterns: &[[usize; 4]],
    out: &Path, proofs: &Mutex<StarProofs>, cache_limit: usize, cap: usize, closure_seconds: u64,
) -> Result<[u128; 6], String> {
    let folder = out.join("frames").join(format!("frame_{fi:06}"));
    fs::create_dir(&folder).map_err(|e| e.to_string())?;
    let (a, b) = sc_unpack(key, n);
    let p = a.compose(b, n);
    let pinv = p.inverse(n);
    let (ia, ib) = (sm_index(a, n), sm_index(b, n));
    assert_eq!(ia + ib, sm_index(p, n));
    let mut words: HashMap<Vec<usize>, u64> = HashMap::new();
    let mut count = [0u128; 6];
    for &(c, ic, cinv) in perms {
        count[0] += 1;
        // cheap index-only filter before computing the fourth entry's cycles
        if !patterns.iter().any(|q| q[0] == ia) { continue; }
        let d = cinv.compose(pinv, n);
        let id = sm_index(d, n);
        if !st_matches(patterns, ia, ib, ic, id) { continue; }
        count[1] += 1;
        let t = [a, b, c, d];
        if !generators_transitive(n, &t[..3]) { continue; }
        count[2] += 1;
        let mut word = sc_word(t, n);
        let mut status = "nonsmooth";
        if word.is_none() {
            let sk = sm_key(&t, n);
            let known = {
                let cache = proofs.lock().map_err(|_| "poisoned proof cache")?;
                if cache.closed.contains_key(&sk) { Some(None) } else { cache.bad.get(&sk).cloned().map(Some) }
            };
            match known {
                Some(None) => status = "closed_smooth",
                Some(Some(w)) => word = Some(w),
                None => {
                    // closure runs without the lock; only the bookkeeping is serialized
                    let scan = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(closure_seconds));
                    status = scan.status;
                    let mut cache = proofs.lock().map_err(|_| "poisoned proof cache")?;
                    match status {
                        "nonsmooth" => {
                            word = Some(scan.word.clone());
                            if cache.bad.len() < cache_limit { cache.bad.insert(sk, scan.word); }
                        }
                        "closed_smooth" => {
                            if !cache.closed.contains_key(&sk) {
                                let id = cache.components;
                                sm_save(&t, n, &scan, &out.join(format!("component_{id:04}")))?;
                                cache.components += 1;
                                for k in scan.states { cache.closed.insert(k, id); }
                            }
                        }
                        _ => {
                            sm_save(&t, n, &scan, &out.join(format!("unresolved_{fi}_{}", count[0])))?;
                        }
                    }
                }
            }
        }
        match status {
            "nonsmooth" => {
                count[3] += 1;
                *words.entry(word.expect("negative word")).or_insert(0) += 1;
            }
            "closed_smooth" => count[4] += 1,
            _ => count[5] += 1,
        }
    }
    let mut w = BufWriter::new(File::create(folder.join("WORDS.tsv")).map_err(|e| e.to_string())?);
    writeln!(w, "word\toccurrences").map_err(|e| e.to_string())?;
    let mut ws: Vec<_> = words.into_iter().collect();
    ws.sort();
    for (word, k) in ws {
        writeln!(w, "{}\t{k}", word.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(","))
            .map_err(|e| e.to_string())?;
    }
    w.flush().map_err(|e| e.to_string())?;
    fs::write(folder.join("COUNTS.tsv"), format!("frame\t{ST_COLUMNS}\n{fi}\t{}\n",
        count.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"))).map_err(|e| e.to_string())?;
    if count[0] == perms.len() as u128 {
        fs::write(folder.join("COMPLETE"), format!("{fi}\t{}\n", count[0])).map_err(|e| e.to_string())?;
    }
    Ok(count)
}

pub fn star_census_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: &str| args.windows(2).find(|w| w[0] == name)
        .map(|w| w[1].clone()).unwrap_or(default.into());
    let n: usize = arg("--d", "9").parse().map_err(|_| "degree")?;
    let workers: usize = arg("--workers", "6").parse().map_err(|_| "workers")?;
    let cap: usize = arg("--cap", "1000000").parse().map_err(|_| "cap")?;
    let cache_limit: usize = arg("--cache-entries", "2000000").parse().map_err(|_| "cache")?;
    let closure_seconds: u64 = arg("--closure-seconds", "120").parse().map_err(|_| "closure seconds")?;
    let patterns = st_parse_patterns(&arg("--patterns", ""))?;
    if !(3..=11).contains(&n) || !(1..=10).contains(&workers) {
        return Err("degree3..11 and workers1..10 required".into());
    }
    for q in &patterns {
        if q.iter().any(|&x| x == 0 || x >= n) { return Err("indices must lie in 1..n-1".into()); }
        let s: usize = q.iter().sum();
        if s % 2 != 0 || s < 2 * n { return Err("pattern must have even sum >= 2n (positive genus)".into()); }
    }
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() { return Err("fresh --out required".into()); }
    fs::create_dir(&out).map_err(|e| e.to_string())?;
    fs::create_dir(out.join("frames")).map_err(|e| e.to_string())?;
    fs::write(out.join("PATTERNS.tsv"), patterns.iter()
        .map(|q| format!("{}\t{}\t{}\t{}\tgenus={}\n", q[0], q[1], q[2], q[3], (q.iter().sum::<usize>() + 2 - 2 * n) / 2))
        .collect::<String>()).map_err(|e| e.to_string())?;
    let start = Instant::now();
    let original_perms = sc_perms(n);
    let all_frames = sc_frames(n, &original_perms);
    let perms: Vec<_> = original_perms.iter().map(|&(p, i)| (p, i, p.inverse(n))).collect();
    drop(original_perms);
    // keep frames whose first entry is a pattern minimum and second index occurs in that pattern
    let frames: Vec<_> = all_frames.iter().copied().filter(|&(key, _)| {
        let (a, b) = sc_unpack(key, n);
        let (ia, ib) = (sm_index(a, n), sm_index(b, n));
        patterns.iter().any(|q| q[0] == ia && q[1..].contains(&ib))
    }).collect();
    let mut table = BufWriter::new(File::create(out.join("FRAMES.tsv")).map_err(|e| e.to_string())?);
    writeln!(table, "frame\tpair_orbit\ta\tb").map_err(|e| e.to_string())?;
    for (fi, &(key, z)) in frames.iter().enumerate() {
        let (a, b) = sc_unpack(key, n);
        writeln!(table, "{fi}\t{}\t{}\t{}", factorial_u128(n) / z, a.cycle_string(n), b.cycle_string(n))
            .map_err(|e| e.to_string())?;
    }
    table.flush().map_err(|e| e.to_string())?;
    println!("degree={n} workers={workers} all_frames={} kept_frames={} setup_seconds={:.3}",
        all_frames.len(), frames.len(), start.elapsed().as_secs_f64());
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let proofs = Mutex::new(StarProofs::default());
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers).map(|_| scope.spawn(|| -> Result<Vec<(usize, [u128; 6])>, String> {
            let mut result = Vec::new();
            loop {
                let fi = next.fetch_add(1, Ordering::Relaxed);
                if fi >= frames.len() { break; }
                let count = st_frame(n, fi, frames[fi].0, &perms, &patterns, &out, &proofs, cache_limit, cap, closure_seconds)?;
                result.push((fi, count));
                let completed = done.fetch_add(1, Ordering::Relaxed) + 1;
                if completed % 10 == 0 || completed == frames.len() {
                    println!("frames_processed={completed}/{} elapsed_seconds={:.3}", frames.len(), start.elapsed().as_secs_f64());
                }
            }
            Ok(result)
        })).collect();
        handles.into_iter().map(|h| h.join().map_err(|_| "worker panicked".to_string()).and_then(|r| r))
            .collect::<Result<Vec<_>, _>>()
    })?;
    let mut total = [0u128; 6];
    let mut completed_frames = 0;
    for (_, count) in results.into_iter().flatten() {
        for j in 0..6 { total[j] += count[j]; }
        if count[0] == perms.len() as u128 { completed_frames += 1; }
    }
    let complete = completed_frames == frames.len();
    let components = proofs.into_inner().map_err(|_| "poisoned proof cache")?.components;
    fs::write(out.join("COVERAGE.tsv"), format!("degree\tkept_frames\tcompleted_frames\t{ST_COLUMNS}\tclosed_components\tenumeration_complete\tfully_resolved\telapsed_seconds\n{n}\t{}\t{completed_frames}\t{}\t{components}\t{complete}\t{}\t{:.3}\n",
        frames.len(), total.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"), complete && total[5] == 0, start.elapsed().as_secs_f64()))
        .map_err(|e| e.to_string())?;
    println!("complete={complete} degree={n} counts={total:?} components={components} seconds={:.3}", start.elapsed().as_secs_f64());
    Ok(())
}
