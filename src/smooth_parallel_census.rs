// Disjoint dynamic frame ownership; immutable permutation tables are shared by
// scoped threads. Only literal negative words and fully closed components enter
// the shared proof cache. Caps/timeouts never become exclusions.
use std::sync::atomic::AtomicUsize;

#[derive(Default)]
struct ParallelProofs {
    closed: HashMap<SmoothKey, usize>,
    bad: HashMap<SmoothKey, Vec<usize>>,
    components: usize,
}
const PC_COLUMNS: &str = "enumerated\ttransitive\tboth_geodesic\tnonidentity\tidentity\tnonsmooth_orbit\tclosed_occurrences\tunresolved";

fn pc_frame(
    n: usize, fi: usize, key: [u64; 2], perms: &[(Perm, usize, Perm)],
    out: &Path, proofs: &Mutex<ParallelProofs>, cache_limit: usize,
    cap: usize, closure_seconds: u64, deadline: Instant, positive_only: bool,
) -> Result<[u128; 8], String> {
    let folder = out.join("frames").join(format!("frame_{fi:05}"));
    fs::create_dir(&folder).map_err(|e| e.to_string())?;
    let (a, b) = sc_unpack(key, n);
    let p = a.compose(b, n);
    let pinv = p.inverse(n);
    let ip = sm_index(p, n);
    assert_eq!(sm_index(a, n) + sm_index(b, n), ip);
    let left = point_orbit_count_perms(n, a, b);
    assert_eq!(left, n - ip);
    let mut records = BufWriter::with_capacity(16 << 20,
        File::create(folder.join("OUTCOMES.bin")).map_err(|e| e.to_string())?);
    let mut dictionary = BufWriter::new(
        File::create(folder.join("WITNESSES.tsv")).map_err(|e| e.to_string())?);
    writeln!(dictionary, "token\tstatus\tword_or_component").map_err(|e| e.to_string())?;
    let mut tokens: HashMap<(String, String), u32> = HashMap::new();
    let mut count = [0u128; 8];
    let mut omitted = [0u128; 2];
    for (rank, &(c, ic, cinv)) in perms.iter().enumerate() {
        if rank % 256 == 0 && Instant::now() >= deadline { break; }
        count[0] += 1;
        // (p*c)^-1 = c^-1*p^-1: table inverses eliminate repeated inversion.
        let d = cinv.compose(pinv, n);
        let t = [a, b, c, d];
        if positive_only {
            if t.iter().any(|q| q.is_identity(n)) { omitted[0] += 1; continue; }
            if ip + ic + sm_index(d, n) <= 2*n-2 { omitted[1] += 1; continue; }
        }
        if !generators_transitive(n, &t[..3]) { continue; }
        count[1] += 1;
        let right = point_orbit_count_perms(n, c, d);
        if right != 1 || (n-ip) as i64 - left as i64 - right as i64 + 1 != 0 {
            fs::write(folder.join("IMPLICATION_COUNTEREXAMPLE.txt"), tuple_string_dynamic(n, &t))
                .map_err(|e| e.to_string())?;
            return Err("smooth-implies-zero counterexample".into());
        }
        if ic + sm_index(d, n) == ip { count[2] += 1; }
        if t.iter().any(|q| q.is_identity(n)) { count[4] += 1; continue; }
        count[3] += 1;
        let mut word = sc_word(t, n);
        let mut status = "nonsmooth";
        let mut component = String::new();
        if word.is_none() {
            let sk = sm_key(&t, n);
            // Short-witness hot path never locks. Serializing rare closure/cache
            // misses prevents duplicate certificates under concurrent discovery.
            let mut cache = proofs.lock().map_err(|_| "poisoned proof cache")?;
            if let Some(&id) = cache.closed.get(&sk) {
                status = "closed_smooth";
                component = format!("component_{id:04}");
            } else if let Some(w) = cache.bad.get(&sk) {
                word = Some(w.clone());
            } else {
                let scan = sm_scan(&t, n, cap,
                    deadline.min(Instant::now() + Duration::from_secs(closure_seconds)));
                status = scan.status;
                match status {
                    "nonsmooth" => {
                        word = Some(scan.word.clone());
                        if cache.bad.len() < cache_limit { cache.bad.insert(sk, scan.word); }
                    }
                    "closed_smooth" => {
                        let id = cache.components;
                        component = format!("component_{id:04}");
                        sm_save(&t, n, &scan, &out.join(&component))?;
                        cache.components += 1;
                        for k in scan.states { cache.closed.insert(k, id); }
                    }
                    _ => {
                        component = format!("unresolved_{fi}_{rank}");
                        sm_save(&t, n, &scan, &out.join(&component))?;
                    }
                }
            }
        }
        let data = word.map(|w| w.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(","))
            .unwrap_or(component);
        match status {
            "nonsmooth" => count[5] += 1,
            "closed_smooth" => count[6] += 1,
            _ => count[7] += 1,
        }
        let signature = (status.to_string(), data);
        let token = if let Some(&id) = tokens.get(&signature) { id } else {
            let id = u32::try_from(tokens.len()).map_err(|_| "dictionary overflow")?;
            writeln!(dictionary, "{id}\t{}\t{}", signature.0, signature.1).map_err(|e| e.to_string())?;
            tokens.insert(signature, id);
            id
        };
        records.write_all(&(rank as u32).to_le_bytes()).map_err(|e| e.to_string())?;
        records.write_all(&token.to_le_bytes()).map_err(|e| e.to_string())?;
    }
    records.flush().map_err(|e| e.to_string())?;
    dictionary.flush().map_err(|e| e.to_string())?;
    if positive_only {
        fs::write(folder.join("OMITTED.tsv"), format!("identity\tnonpositive_genus\n{}\t{}\n", omitted[0], omitted[1]))
            .map_err(|e| e.to_string())?;
    }
    fs::write(folder.join("COUNTS.tsv"), format!("frame\t{PC_COLUMNS}\n{fi}\t{}\n",
        count.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"))).map_err(|e| e.to_string())?;
    // COMPLETE is last and only written after a whole frame; partial evidence
    // remains durable but cannot be mistaken for gap-free enumeration.
    if count[0] == perms.len() as u128 {
        fs::write(folder.join("COMPLETE"), format!("{fi}\t{}\n", count[0])).map_err(|e| e.to_string())?;
    }
    Ok(count)
}

pub fn parallel_smooth_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: &str| args.windows(2).find(|w| w[0] == name)
        .map(|w| w[1].clone()).unwrap_or(default.into());
    let n: usize = arg("--d", "9").parse().map_err(|_| "degree")?;
    let workers: usize = arg("--workers", "6").parse().map_err(|_| "workers")?;
    let seconds: u64 = arg("--seconds", "31536000").parse().map_err(|_| "seconds")?;
    let cap: usize = arg("--cap", "100000").parse().map_err(|_| "cap")?;
    let cache_limit: usize = arg("--cache-entries", "1000000").parse().map_err(|_| "cache")?;
    let closure_seconds: u64 = arg("--closure-seconds", "30").parse().map_err(|_| "closure seconds")?;
    let positive_only = args.iter().any(|a| a == "--positive-only");
    if !(3..=10).contains(&n) || !(1..=10).contains(&workers) || seconds == 0 || cap == 0 || closure_seconds == 0 {
        return Err("degree3..10, workers1..10 and positive limits required".into());
    }
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() { return Err("fresh --out required".into()); }
    fs::create_dir(&out).map_err(|e| e.to_string())?;
    fs::create_dir(out.join("frames")).map_err(|e| e.to_string())?;
    if positive_only { fs::write(out.join("POSITIVE_ONLY"), "Identity-free, genus>0. Counts except enumerated are restricted to this scope.\n")
        .map_err(|e|e.to_string())?; }
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    let original_perms = sc_perms(n);
    let frames = sc_frames(n, &original_perms);
    let perms: Vec<_> = original_perms.iter().map(|&(p, i)| (p, i, p.inverse(n))).collect();
    drop(original_perms);
    let mut table = BufWriter::new(File::create(out.join("FRAMES.tsv")).map_err(|e| e.to_string())?);
    writeln!(table, "frame\tpair_orbit\ta\tb").map_err(|e| e.to_string())?;
    for (fi, &(key, z)) in frames.iter().enumerate() {
        let (a,b) = sc_unpack(key,n);
        writeln!(table, "{fi}\t{}\t{}\t{}", factorial_u128(n)/z, a.cycle_string(n), b.cycle_string(n))
            .map_err(|e| e.to_string())?;
    }
    table.flush().map_err(|e| e.to_string())?;
    println!("degree={n} workers={workers} frames={} third_entries={} shared_table_bytes={} setup_seconds={:.3}",
        frames.len(), frames.len()*perms.len(), perms.len()*std::mem::size_of::<(Perm,usize,Perm)>(), start.elapsed().as_secs_f64());
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let proofs = Mutex::new(ParallelProofs::default());
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers).map(|_| scope.spawn(|| -> Result<Vec<(usize,[u128;8])>,String> {
            let mut result = Vec::new();
            while Instant::now() < deadline {
                let fi = next.fetch_add(1,Ordering::Relaxed);
                if fi >= frames.len() { break; }
                let count = pc_frame(n,fi,frames[fi].0,&perms,&out,&proofs,cache_limit,cap,closure_seconds,deadline,positive_only)?;
                result.push((fi,count));
                let completed = done.fetch_add(1,Ordering::Relaxed)+1;
                if completed % 10 == 0 || completed == frames.len() {
                    println!("frames_processed={completed}/{} elapsed_seconds={:.3}",frames.len(),start.elapsed().as_secs_f64());
                }
            }
            Ok(result)
        })).collect();
        handles.into_iter().map(|h| h.join().map_err(|_| "worker panicked".to_string()).and_then(|r|r))
            .collect::<Result<Vec<_>,_>>()
    })?;
    let mut results: Vec<_> = results.into_iter().flatten().collect();
    results.sort_unstable_by_key(|r|r.0);
    let mut total = [0u128;8];
    let mut weighted = [0u128;2];
    let mut completed_frames = 0;
    let mut progress = BufWriter::new(File::create(out.join("FRAME_COUNTS.tsv")).map_err(|e|e.to_string())?);
    writeln!(progress,"frame\t{PC_COLUMNS}").map_err(|e|e.to_string())?;
    for (fi,count) in results {
        for j in 0..8 { total[j] += count[j]; }
        let mass = factorial_u128(n)/frames[fi].1;
        weighted[0] += mass*count[1]; weighted[1] += mass*count[2];
        if count[0] == perms.len() as u128 { completed_frames += 1; }
        writeln!(progress,"{fi}\t{}",count.iter().map(|x|x.to_string()).collect::<Vec<_>>().join("\t"))
            .map_err(|e|e.to_string())?;
    }
    progress.flush().map_err(|e|e.to_string())?;
    let complete = completed_frames == frames.len();
    if complete { assert_eq!(total[0], frames.len() as u128*factorial_u128(n)); }
    let cache = proofs.into_inner().map_err(|_|"poisoned proof cache")?;
    let components = cache.components;
    fs::write(out.join("CACHE.tsv"),format!("negative_entries\tclosed_state_entries\tworkers\n{}\t{}\t{workers}\n",cache.bad.len(),cache.closed.len()))
        .map_err(|e|e.to_string())?;
    fs::write(out.join("COVERAGE.tsv"),format!("degree\tframes\tcompleted_frames\texpected\tenumerated\ttransitive_smooth_occurrences\tboth_geodesic_occurrences\tnonidentity_occurrences\tidentity_occurrences\tnonsmooth_orbit_occurrences\tclosed_occurrences\tunresolved\tclosed_components\timplication_counterexamples\tenumeration_complete\tfully_resolved\tlabeled_left_smooth\tlabeled_both_smooth\tlabeled_smooth_tuples\telapsed_seconds\n{n}\t{}\t{completed_frames}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{components}\t0\t{complete}\t{}\t{}\t{}\t{}\t{}\n",frames.len(),frames.len() as u128*factorial_u128(n),total[0],total[1],total[2],total[3],total[4],total[5],total[6],total[7],complete&&total[7]==0,weighted[0],weighted[1],2*weighted[0]-weighted[1],start.elapsed().as_secs_f64()))
        .map_err(|e|e.to_string())?;
    println!("complete={complete} degree={n} counts={total:?} components={components} seconds={:.3}",start.elapsed().as_secs_f64());
    Ok(())
}

#[cfg(test)]
mod parallel_census_tests {
    use super::*;
    #[test]
    fn inverse_table_matches_literal_fourth_entry() {
        for n in 3..=5 {
            let perms = sc_perms(n);
            for (key,_) in sc_frames(n,&perms) {
                let (a,b)=sc_unpack(key,n); let p=a.compose(b,n);
                for &(c,_) in &perms {
                    assert_eq!(c.inverse(n).compose(p.inverse(n),n),p.compose(c,n).inverse(n));
                }
            }
        }
    }
}
