// All passports, all transitive generated groups. The only reduction is exact
// simultaneous conjugation of an index-additive first pair; c ranges over S_n.
fn sc_perms(n: usize) -> Vec<(Perm, usize)> {
    let mut a: Vec<u8> = (0..n as u8).collect();
    let mut result = Vec::new();
    loop {
        let mut p = Perm::identity(n);
        p.img[..n].copy_from_slice(&a);
        result.push((p, sm_index(p, n)));
        let Some(i) = (0..n - 1).rev().find(|&i| a[i] < a[i + 1]) else {
            break;
        };
        let j = (i + 1..n).rev().find(|&j| a[i] < a[j]).unwrap();
        a.swap(i, j);
        a[i + 1..].reverse();
    }
    assert_eq!(result.len() as u128, factorial_u128(n));
    result
}
fn sc_partitions(rem: usize, min: usize, parts: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if rem == 0 {
        out.push(parts.clone());
    } else {
        for k in min..=rem {
            parts.push(k);
            sc_partitions(rem - k, k, parts, out);
            parts.pop();
        }
    }
}
// Canonicalize each connected pair action by all rooted equivariant traversals,
// then sort its components. Count automorphisms and interchangeable components.
fn sc_pair(a: Perm, b: Perm, n: usize) -> ([u64; 2], u128) {
    let mut remaining = (1u64 << n) - 1;
    let mut components = Vec::new();
    while remaining != 0 {
        let first = remaining.trailing_zeros() as usize;
        let mut mask = 0;
        let mut best = (usize::MAX, [u64::MAX; 2]);
        let mut aut = 0u128;
        let mut roots = remaining;
        // The first traversal identifies this component; subsequent roots
        // are only its points, never roots in other components.
        let mut initial = true;
        while roots != 0 {
            let root = if initial {
                first
            } else {
                roots.trailing_zeros() as usize
            };
            let mut labels = [u8::MAX; MAX_D];
            let mut order = [0usize; MAX_D];
            labels[root] = 0;
            order[0] = root;
            let mut head = 0;
            let mut len = 1;
            while head < len {
                let x = order[head];
                head += 1;
                for p in [a, b] {
                    let y = p.img[x] as usize;
                    if labels[y] == u8::MAX {
                        labels[y] = len as u8;
                        order[len] = y;
                        len += 1;
                    }
                }
            }
            if initial {
                for &x in &order[..len] {
                    mask |= 1 << x;
                }
                roots = mask;
                initial = false;
            }
            roots &= !(1 << root);
            let mut key = [0u64; 2];
            for (j, p) in [a, b].iter().enumerate() {
                for i in 0..len {
                    key[j] |= (labels[p.img[order[i]] as usize] as u64) << (4 * i);
                }
            }
            let value = (len, key);
            if value < best {
                best = value;
                aut = 1;
            } else if value == best {
                aut += 1;
            }
        }
        remaining &= !mask;
        components.push((best.0, best.1, aut));
    }
    components.sort_unstable();
    let mut key = [0u64; 2];
    let mut offset = 0;
    let mut centralizer = 1;
    let mut multiplicity = 0;
    let mut previous = None;
    for (len, local, aut) in components {
        let signature = (len, local);
        if previous == Some(signature) {
            multiplicity += 1;
        } else {
            multiplicity = 1;
            previous = Some(signature);
        }
        centralizer *= aut * multiplicity;
        for j in 0..2 {
            for i in 0..len {
                key[j] |= (((local[j] >> (4 * i)) & 15) + offset as u64) << (4 * (offset + i));
            }
        }
        offset += len;
    }
    assert_eq!(offset, n);
    (key, centralizer)
}
fn sc_unpack(k: [u64; 2], n: usize) -> (Perm, Perm) {
    let mut t = [Perm::identity(n); 2];
    for j in 0..2 {
        for i in 0..n {
            t[j].img[i] = ((k[j] >> (4 * i)) & 15) as u8;
        }
    }
    (t[0], t[1])
}
fn sc_frames(n: usize, perms: &[(Perm, usize)]) -> Vec<([u64; 2], u128)> {
    let mut parts = Vec::new();
    sc_partitions(n, 1, &mut Vec::new(), &mut parts);
    let mut frames = std::collections::BTreeMap::new();
    let mut labeled_mass = 0u128;
    for lengths in parts {
        let mut p = Perm::identity(n);
        let mut offset = 0;
        let mut multiplicities = HashMap::new();
        let mut z = 1u128;
        for &len in &lengths {
            for i in 0..len {
                p.img[offset + i] = (offset + (i + 1) % len) as u8;
            }
            offset += len;
            let m = multiplicities.entry(len).or_insert(0u128);
            *m += 1;
            z *= len as u128 * *m;
        }
        let ip = n - lengths.len();
        let mut raw = 0;
        for &(a, ia) in perms {
            let b = a.inverse(n).compose(p, n);
            if ia + sm_index(b, n) != ip {
                continue;
            }
            raw += 1;
            let (key, centralizer) = sc_pair(a, b, n);
            if let Some(old) = frames.insert(key, centralizer) {
                assert_eq!(old, centralizer);
            }
        }
        // A second mass check: geodesic factors on disjoint product cycles
        // have product Catalan counts; this is a check, not a pruning rule.
        let expected: u128 = lengths
            .iter()
            .map(|&l| {
                factorial_u128(2 * l) / (factorial_u128(l) * factorial_u128(l) * (l + 1) as u128)
            })
            .product();
        assert_eq!(raw, expected, "geodesic factorization mass");
        labeled_mass += raw * factorial_u128(n) / z;
    }
    assert_eq!(
        frames.values().map(|z| factorial_u128(n) / z).sum::<u128>(),
        labeled_mass
    );
    frames.into_iter().collect()
}
fn sc_word(t: [Perm; 4], n: usize) -> Option<Vec<usize>> {
    // Purely literal witnesses, no genus/passport exclusion heuristic.
    for first in [None, Some(0), Some(3), Some(2), Some(5)] {
        let u = first.map(|m| hurwitz_move_dynamic(n, t, m)).unwrap_or(t);
        for middle in [1, 4] {
            let v = hurwitz_move_dynamic(n, u, middle);
            let ip = sm_index(v[0].compose(v[1], n), n);
            let indices = v.map(|p| sm_index(p, n));
            if ip != indices[0] + indices[1] && ip != indices[2] + indices[3] {
                let mut word: Vec<_> = first.into_iter().collect();
                word.push(middle);
                return Some(word);
            }
        }
    }
    None
}
pub fn small_smooth_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: &str| {
        args.windows(2)
            .find(|w| w[0] == name)
            .map(|w| w[1].clone())
            .unwrap_or(default.into())
    };
    let n: usize = arg("--d", "0").parse().map_err(|_| "degree")?;
    if !(3..=9).contains(&n) {
        return Err("degree 3..9 required".into());
    }
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() {
        return Err("fresh --out required".into());
    }
    let seconds: u64 = arg("--seconds", "3600").parse().map_err(|_| "seconds")?;
    let cap: usize = arg("--cap", "50000").parse().map_err(|_| "cap")?;
    if seconds == 0 || cap == 0 {
        return Err("positive limits required".into());
    }
    fs::create_dir(&out).map_err(|e| e.to_string())?;
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    let perms = sc_perms(n);
    let frames = sc_frames(n, &perms);
    let mut table =
        BufWriter::new(File::create(out.join("FRAMES.tsv")).map_err(|e| e.to_string())?);
    writeln!(table, "frame\tpair_orbit\ta\tb").map_err(|e| e.to_string())?;
    for (i, &(key, z)) in frames.iter().enumerate() {
        let (a, b) = sc_unpack(key, n);
        writeln!(
            table,
            "{i}\t{}\t{}\t{}",
            factorial_u128(n) / z,
            a.cycle_string(n),
            b.cycle_string(n)
        )
        .map_err(|e| e.to_string())?;
    }
    table.flush().map_err(|e| e.to_string())?;
    println!(
        "degree={n} frames={} third_entries={} setup_seconds={:.3}",
        frames.len(),
        frames.len() * perms.len(),
        start.elapsed().as_secs_f64()
    );
    let mut records = BufWriter::with_capacity(
        1 << 20,
        File::create(out.join("OUTCOMES.tsv")).map_err(|e| e.to_string())?,
    );
    writeln!(records, "frame\tthird_rank\tstatus\tword_or_component").map_err(|e| e.to_string())?;
    let mut progress =
        BufWriter::new(File::create(out.join("FRAME_COUNTS.tsv")).map_err(|e| e.to_string())?);
    writeln!(progress, "frame\tenumerated\ttransitive\tboth_geodesic\tnonidentity\tidentity\tnonsmooth_orbit\tclosed_occurrences\tunresolved").map_err(|e| e.to_string())?;
    let mut total = [0u128; 8];
    let mut weighted = [0u128; 3]; // left-smooth, both-smooth, nonidentity left-smooth
    let mut closed: HashMap<SmoothKey, String> = HashMap::new();
    let mut bad: HashMap<SmoothKey, Vec<usize>> = HashMap::new();
    let mut components = 0usize;
    let mut completed_frames = 0;
    for (fi, &(key, z)) in frames.iter().enumerate() {
        if Instant::now() >= deadline {
            break;
        }
        let (a, b) = sc_unpack(key, n);
        let p = a.compose(b, n);
        let ip = sm_index(p, n);
        assert_eq!(sm_index(a, n) + sm_index(b, n), ip);
        assert_eq!(point_orbit_count_perms(n, a, b), n - ip);
        let mut count = [0u128; 8];
        for (rank, &(c, ic)) in perms.iter().enumerate() {
            if rank % 256 == 0 && Instant::now() >= deadline {
                break;
            }
            count[0] += 1;
            let d = p.compose(c, n).inverse(n);
            let t = [a, b, c, d];
            if !generators_transitive(n, &t[..3]) {
                continue;
            }
            count[1] += 1;
            // Independently count the complementary pair's point orbits for
            // EVERY transitive normal-form occurrence, including identity cases.
            let r = point_orbit_count_perms(n, c, d);
            if r != 1 {
                fs::write(
                    out.join("IMPLICATION_COUNTEREXAMPLE.txt"),
                    tuple_string_dynamic(n, &t),
                )
                .map_err(|e| e.to_string())?;
                return Err("smooth-implies-zero counterexample".into());
            }
            assert_eq!(
                (n - ip) as i64 - point_orbit_count_perms(n, a, b) as i64 - r as i64 + 1,
                0
            );
            if ic + sm_index(d, n) == ip {
                count[2] += 1;
            }
            if t.iter().any(|q| q.is_identity(n)) {
                count[4] += 1;
                continue;
            }
            count[3] += 1;
            let mut word = sc_word(t, n);
            let mut status = "nonsmooth";
            let mut component = String::new();
            if word.is_none() {
                let sk = sm_key(&t, n);
                if let Some(name) = closed.get(&sk) {
                    status = "closed_smooth";
                    component = name.clone();
                } else if let Some(w) = bad.get(&sk) {
                    word = Some(w.clone());
                } else {
                    let scan = sm_scan(
                        &t,
                        n,
                        cap,
                        deadline.min(Instant::now() + Duration::from_secs(10)),
                    );
                    status = scan.status;
                    match status {
                        "nonsmooth" => {
                            word = Some(scan.word.clone());
                            if bad.len() < 200000 {
                                bad.insert(sk, scan.word);
                            }
                        }
                        "closed_smooth" => {
                            component = format!("component_{components:04}");
                            components += 1;
                            sm_save(&t, n, &scan, &out.join(&component))?;
                            for k in scan.states {
                                closed.insert(k, component.clone());
                            }
                        }
                        _ => {
                            component = format!("unresolved_{fi}_{rank}");
                            sm_save(&t, n, &scan, &out.join(&component))?;
                        }
                    }
                }
            }
            let data = if let Some(w) = word {
                w.iter()
                    .map(|m| m.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            } else {
                component
            };
            match status {
                "nonsmooth" => count[5] += 1,
                "closed_smooth" => count[6] += 1,
                _ => count[7] += 1,
            }
            writeln!(records, "{fi}\t{rank}\t{status}\t{data}").map_err(|e| e.to_string())?;
        }
        let mass = factorial_u128(n) / z;
        weighted[0] += mass * count[1];
        weighted[1] += mass * count[2];
        weighted[2] += mass * count[3];
        for j in 0..8 {
            total[j] += count[j];
        }
        writeln!(
            progress,
            "{fi}\t{}",
            count
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join("\t")
        )
        .map_err(|e| e.to_string())?;
        progress.flush().map_err(|e| e.to_string())?;
        records.flush().map_err(|e| e.to_string())?;
        if count[0] != perms.len() as u128 {
            break;
        }
        completed_frames += 1;
        if fi % 10 == 0 {
            println!(
                "frame={}/{} checked={} transitive={} unresolved={} elapsed={:.2}",
                fi + 1,
                frames.len(),
                total[0],
                total[1],
                total[7],
                start.elapsed().as_secs_f64()
            );
        }
    }
    let complete = completed_frames == frames.len();
    if complete {
        assert_eq!(total[0], frames.len() as u128 * factorial_u128(n));
    }
    fs::write(out.join("COVERAGE.tsv"), format!("degree\tframes\tcompleted_frames\texpected\tenumerated\ttransitive_smooth_occurrences\tboth_geodesic_occurrences\tnonidentity_occurrences\tidentity_occurrences\tnonsmooth_orbit_occurrences\tclosed_occurrences\tunresolved\tclosed_components\timplication_counterexamples\tenumeration_complete\tfully_resolved\tlabeled_left_smooth\tlabeled_both_smooth\tlabeled_smooth_tuples\telapsed_seconds\n{n}\t{}\t{completed_frames}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{components}\t0\t{complete}\t{}\t{}\t{}\t{}\t{}\n", frames.len(), frames.len() as u128 * factorial_u128(n), total[0], total[1], total[2], total[3], total[4], total[5], total[6], total[7], complete && total[7] == 0, weighted[0], weighted[1], 2 * weighted[0] - weighted[1], start.elapsed().as_secs_f64())).map_err(|e| e.to_string())?;
    println!(
        "complete={complete} degree={n} counts={total:?} components={components} seconds={:.3}",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
#[cfg(test)]
mod small_census_tests {
    use super::*;
    #[test]
    fn pair_quotient_and_centralizers_match_literal_conjugation() {
        let n = 4;
        let perms = sc_perms(n);
        let mut seen = HashSet::new();
        let mut mass = 0;
        for &(a, ia) in &perms {
            for &(b, ib) in &perms {
                if ia + ib != sm_index(a.compose(b, n), n) {
                    continue;
                }
                let (key, z) = sc_pair(a, b, n);
                let literal = perms
                    .iter()
                    .filter(|&&(h, _)| a.conjugate_by(h, n) == a && b.conjugate_by(h, n) == b)
                    .count();
                assert_eq!(z, literal as u128);
                for &(h, _) in &perms {
                    assert_eq!(
                        sc_pair(a.conjugate_by(h, n), b.conjugate_by(h, n), n),
                        (key, z)
                    );
                }
                seen.insert(key);
                mass += 1;
            }
        }
        let frames = sc_frames(n, &perms);
        assert_eq!(seen.len(), frames.len());
        assert_eq!(
            mass,
            frames
                .iter()
                .map(|(_, z)| factorial_u128(n) / z)
                .sum::<u128>()
        );
    }
    #[test]
    fn weighted_reduction_matches_every_literal_s4_tuple() {
        let n = 4;
        let perms = sc_perms(n);
        let mut left = 0u128;
        let mut both = 0u128;
        let mut all = 0u128;
        for &(a, _) in &perms {
            for &(b, _) in &perms {
                for &(c, _) in &perms {
                    let t = [a, b, c, a.compose(b, n).compose(c, n).inverse(n)];
                    if !generators_transitive(n, &t[..3]) {
                        continue;
                    }
                    let ip = sm_index(a.compose(b, n), n);
                    let l = sm_index(a, n) + sm_index(b, n) == ip;
                    let r = sm_index(c, n) + sm_index(t[3], n) == ip;
                    left += l as u128;
                    both += (l && r) as u128;
                    all += (l || r) as u128;
                    if l || r {
                        assert_eq!(sm_values(&t, n).1, 0);
                    }
                }
            }
        }
        let mut reduced_left = 0;
        let mut reduced_both = 0;
        for (key, z) in sc_frames(n, &perms) {
            let (a, b) = sc_unpack(key, n);
            let p = a.compose(b, n);
            for &(c, ic) in &perms {
                let d = p.compose(c, n).inverse(n);
                if !generators_transitive(n, &[a, b, c]) {
                    continue;
                }
                reduced_left += factorial_u128(n) / z;
                if ic + sm_index(d, n) == sm_index(p, n) {
                    reduced_both += factorial_u128(n) / z;
                }
            }
        }
        assert_eq!(
            (left, both, all),
            (reduced_left, reduced_both, 2 * reduced_left - reduced_both)
        );
    }
}
