// Exact finite end-fixed family, not a census of all S18 or A18 tuples.
// Uses the production exact group-order, closure and certificate verifier.

fn lr_cycles(p: Perm, d: usize) -> Vec<Vec<usize>> {
    let mut seen = [false; MAX_D];
    let mut cycles = Vec::new();
    for x in 0..d {
        if seen[x] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut y = x;
        while !seen[y] {
            seen[y] = true;
            cycle.push(y);
            y = p.img[y] as usize;
        }
        cycles.push(cycle);
    }
    cycles
}

// Every square root either preserves an odd cycle or interleaves two equal
// cycles. These disjoint choices exhaust roots, including fixed-point pairs.
fn lr_roots(q: Perm, d: usize) -> Vec<Perm> {
    fn visit(cycles: &[Vec<usize>], used: u64, r: Perm, out: &mut Vec<Perm>) {
        let Some(i) = (0..cycles.len()).find(|&i| used & (1 << i) == 0) else {
            out.push(r);
            return;
        };
        let c = &cycles[i];
        let l = c.len();
        let used = used | (1 << i);
        if l % 2 == 1 {
            let mut next = r;
            for j in 0..l {
                next.img[c[j]] = c[(j + (l + 1) / 2) % l] as u8;
            }
            visit(cycles, used, next, out);
        }
        for j in i + 1..cycles.len() {
            if used & (1 << j) != 0 || cycles[j].len() != l {
                continue;
            }
            for offset in 0..l {
                let mut next = r;
                for z in 0..l {
                    let y = cycles[j][(z + offset) % l];
                    next.img[c[z]] = y as u8;
                    next.img[y] = c[(z + 1) % l] as u8;
                }
                visit(cycles, used | (1 << j), next, out);
            }
        }
    }
    let mut roots = Vec::new();
    visit(&lr_cycles(q, d), 0, Perm::identity(d), &mut roots);
    roots.sort_unstable_by_key(|r| r.img);
    assert!(roots.windows(2).all(|w| w[0] != w[1]));
    assert!(roots.iter().all(|r| r.compose(*r, d) == q));
    roots
}

fn lr_frame(text: &str, d: usize) -> Result<[Perm; 3], String> {
    let mut out = [Perm::identity(d); 3];
    let lines: Vec<_> = text.lines().filter(|s| !s.is_empty()).collect();
    if lines.len() != 3 {
        return Err("frame requires h,k,p rows (one-based images)".into());
    }
    for (i, name) in ["h", "k", "p"].iter().enumerate() {
        let (label, values) = lines[i].split_once('\t').ok_or("missing frame tab")?;
        if label != *name {
            return Err("wrong frame row order".into());
        }
        let images: Vec<usize> = values
            .split(',')
            .map(|s| s.parse::<usize>().map_err(|e| e.to_string()))
            .collect::<Result<_, _>>()?;
        if images.len() != d {
            return Err("frame dimension mismatch".into());
        }
        let mut sorted = images.clone();
        sorted.sort_unstable();
        if sorted != (1..=d).collect::<Vec<_>>() {
            return Err("invalid frame permutation".into());
        }
        for (x, y) in images.into_iter().enumerate() {
            out[i].img[x] = (y - 1) as u8;
        }
    }
    Ok(out)
}

fn lr_centralizer(h: Perm, k: Perm, d: usize) -> Result<Vec<Perm>, String> {
    // Two regular C3^2 orbits: coordinates (h^i k^j root) determine the
    // four independent translations and the canonical exchange of the blocks.
    let mut covered = [false; MAX_D];
    let mut blocks = Vec::new();
    for x in 0..d {
        if covered[x] {
            continue;
        }
        let mut block = Vec::new();
        let mut z = x;
        for _ in 0..3 {
            let mut y = z;
            for _ in 0..3 {
                block.push(y);
                y = k.img[y] as usize;
            }
            z = h.img[z] as usize;
        }
        let unique: HashSet<_> = block.iter().copied().collect();
        if unique.len() != 9 || block.iter().any(|&y| covered[y]) {
            return Err("frame is not two regular C3^2 blocks".into());
        }
        for &y in &block {
            covered[y] = true;
        }
        blocks.push(block);
    }
    if blocks.len() != 2 {
        return Err("frame must have two regular blocks".into());
    }
    let mut gens = Vec::new();
    for block in &blocks {
        for p in [h, k] {
            let mut g = Perm::identity(d);
            for &x in block {
                g.img[x] = p.img[x];
            }
            gens.push(g);
        }
    }
    let mut swap = Perm::identity(d);
    for i in 0..9 {
        swap.img[blocks[0][i]] = blocks[1][i] as u8;
        swap.img[blocks[1][i]] = blocks[0][i] as u8;
    }
    gens.push(swap);
    if gens
        .iter()
        .any(|g| h.conjugate_by(*g, d) != h || k.conjugate_by(*g, d) != k)
        || smart::subgroup_order(d, &gens) != 162
    {
        return Err("frame centralizer validation failed".into());
    }
    Ok(gens)
}

fn lr_pair_orbits(left: &[Perm], right: &[Perm], gens: &[Perm], d: usize) -> Vec<(usize, usize)> {
    let index = |roots: &[Perm]| {
        roots
            .iter()
            .enumerate()
            .map(|(i, &p)| (p, i))
            .collect::<HashMap<_, _>>()
    };
    let li = index(left);
    let ri = index(right);
    let actions: Vec<_> = gens
        .iter()
        .map(|&g| {
            (
                left.iter()
                    .map(|&p| li[&p.conjugate_by(g, d)])
                    .collect::<Vec<_>>(),
                right
                    .iter()
                    .map(|&p| ri[&p.conjugate_by(g, d)])
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    let m = right.len();
    let mut seen = vec![false; left.len() * m];
    let mut reps = Vec::new();
    let mut queue = Vec::with_capacity(162);
    for root in 0..seen.len() {
        if seen[root] {
            continue;
        }
        queue.clear();
        queue.push(root);
        seen[root] = true;
        let mut head = 0;
        while head < queue.len() {
            let x = queue[head];
            head += 1;
            for (a, b) in &actions {
                let y = a[x / m] * m + b[x % m];
                if !seen[y] {
                    seen[y] = true;
                    queue.push(y);
                }
            }
        }
        assert!(162 % queue.len() == 0);
        reps.push((root, queue.len()));
    }
    assert_eq!(reps.iter().map(|(_, mass)| mass).sum::<usize>(), seen.len());
    reps
}

fn lr_tuple(h: Perm, k: Perm, u: Perm, v: Perm, d: usize) -> [Perm; 4] {
    [
        u.compose(h, d),
        h.inverse(d).compose(u, d).compose(h.compose(h, d), d),
        v.compose(k, d),
        k.inverse(d).compose(v, d).compose(k.compose(k, d), d),
    ]
}

pub fn linked_roots_main(args: &[String]) -> Result<(), String> {
    let d = 18;
    let frame_path = get_arg(args, "--frame", None)?;
    let frame_text = fs::read_to_string(&frame_path).map_err(|e| e.to_string())?;
    let [h, k, p] = lr_frame(&frame_text, d)?;
    if h.compose(k, d) != k.compose(h, d)
        || !h.compose(h, d).compose(h, d).is_identity(d)
        || !k.compose(k, d).compose(k, d).is_identity(d)
        || p != h.compose(h, d).compose(k, d)
    {
        return Err("unsupported frame algebra".into());
    }
    let gens = lr_centralizer(h, k, d)?;
    let left = lr_roots(k, d);
    let right = lr_roots(h, d);
    if left.len() != 856 || right.len() != 856 {
        return Err("expected 856 roots on each side".into());
    }
    let reps = lr_pair_orbits(&left, &right, &gens, d);
    let out = PathBuf::from(get_arg(args, "--out", None)?);
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let create = |name: &str| {
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join(name))
            .map_err(|e| e.to_string())
    };
    let mut table = BufWriter::new(create("FRAME_CLASSES.tsv")?);
    let mut replay_a = BufWriter::new(create("REPLAY_A.txt")?);
    let mut replay_s = BufWriter::new(create("REPLAY_S.txt")?);
    writeln!(
        table,
        "class\tleft_root\tright_root\tlabeled_mass\tstatus\tgroup_order\tgenus\tdefect\tseed"
    )
    .map_err(|e| e.to_string())?;
    let mut counts: HashMap<String, (usize, usize)> = HashMap::new();
    let started = Instant::now();
    eprintln!(
        "exact linked roots: {} pairs -> {} frame classes",
        left.len() * right.len(),
        reps.len()
    );
    for (i, &(root, mass)) in reps.iter().enumerate() {
        let t = lr_tuple(h, k, left[root / right.len()], right[root % right.len()], d);
        if !tuple_product_is_identity(d, &t)
            || t[0].compose(t[1], d) != p
            || [0, 2].iter().any(|&mv| {
                let g = if mv == 0 { h } else { k };
                let moved = hurwitz_move_dynamic(d, t, mv);
                (0..4).any(|j| t[j].conjugate_by(g, d) != moved[j])
            })
        {
            return Err(format!("literal construction failure at class {i}"));
        }
        let genus = 1 - d as i64
            + t.iter()
                .map(|&g| (d - lr_cycles(g, d).len()) as i64)
                .sum::<i64>()
                / 2;
        let defect = lr_cycles(p, d).len() as i64
            - point_orbit_count_perms(d, t[0], t[1]) as i64
            - point_orbit_count_perms(d, t[2], t[3]) as i64
            + 1;
        let mut order = 0;
        let status = if t.iter().any(|g| g.is_identity(d)) {
            "identity_triple"
        } else if !generators_transitive(d, &t[..3]) {
            "intransitive"
        } else {
            order = smart::subgroup_order(d, &t[..3]);
            let ambient = if t.iter().all(|g| g.is_even(d)) {
                Ambient::Alternating
            } else {
                Ambient::Symmetric
            };
            if order != ambient.target_order(d) {
                "proper_transitive"
            } else if generators_transitive(d, &t[..2]) {
                "literal_first_pair_positive"
            } else {
                if ambient == Ambient::Alternating {
                    writeln!(replay_a, "{}", tuple_string_dynamic(d, &t))
                } else {
                    writeln!(replay_s, "{}", tuple_string_dynamic(d, &t))
                }
                .map_err(|e| e.to_string())?;
                if ambient == Ambient::Alternating {
                    "closure_A"
                } else {
                    "closure_S"
                }
            }
        };
        let entry = counts.entry(status.into()).or_default();
        entry.0 += 1;
        entry.1 += mass;
        writeln!(
            table,
            "{i}\t{}\t{}\t{mass}\t{status}\t{order}\t{genus}\t{defect}\t{}",
            root / right.len(),
            root % right.len(),
            tuple_string_dynamic(d, &t)
        )
        .map_err(|e| e.to_string())?;
        if i % 500 == 0 {
            eprintln!(
                "classified {}/{} in {:.1}s",
                i + 1,
                reps.len(),
                started.elapsed().as_secs_f64()
            );
        }
    }
    table.flush().map_err(|e| e.to_string())?;
    replay_a.flush().map_err(|e| e.to_string())?;
    replay_s.flush().map_err(|e| e.to_string())?;
    let mut summary = create("CLASSIFICATION.txt")?;
    writeln!(summary,"roots_left = {}\nroots_right = {}\nlabeled_pairs = {}\nframe_classes = {}\nframe_centralizer_order = 162\nclassification_complete = true\nelapsed_seconds = {:.3}",left.len(),right.len(),left.len()*right.len(),reps.len(),started.elapsed().as_secs_f64()).map_err(|e|e.to_string())?;
    let mut keys: Vec<_> = counts.keys().collect();
    keys.sort();
    for key in keys {
        let (classes, mass) = counts[key];
        writeln!(
            summary,
            "{key}_classes = {classes}\n{key}_labeled_pairs = {mass}"
        )
        .map_err(|e| e.to_string())?;
    }
    summary.sync_all().map_err(|e| e.to_string())?;
    if has_flag(args, "--classify-only") {
        return Ok(());
    }
    let threads = get_arg(args, "--threads", Some("7"))?;
    let seconds = get_arg(args, "--seconds", Some("3600"))?;
    let candidate_seconds = get_arg(args, "--candidate-seconds", Some("30"))?;
    let cap = get_arg(args, "--cap", Some("200000"))?;
    for group in ["A", "S"] {
        if counts.get(&format!("closure_{group}")).map_or(0, |c| c.0) == 0 {
            continue;
        }
        let closure = out.join(format!("closure_{group}"));
        let command: Vec<String> = [
            "linked-roots",
            "smart-hunt",
            "--group",
            group,
            "--d",
            "18",
            "--family",
            "low-index",
            "--replay-only",
            "--weak-only",
            "--cache-positive-witnesses",
            "--threads",
            &threads,
            "--seconds",
            &seconds,
            "--candidate-seconds",
            &candidate_seconds,
            "--cap",
            &cap,
            "--known-certificates",
            &get_arg(args, "--known-certificates", None)?,
            "--replay",
            &out.join(format!("REPLAY_{group}.txt")).to_string_lossy(),
            "--out",
            &closure.to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        smart::run(&command)?;
        for file in fs::read_dir(&closure).map_err(|e| e.to_string())? {
            let path = file.map_err(|e| e.to_string())?.path();
            if path
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("CERTIFIED_")
            {
                smart::run(&[
                    "linked-roots".into(),
                    "verify".into(),
                    "--certificate".into(),
                    path.to_string_lossy().into(),
                ])?;
            }
        }
    }
    // The launcher reconciles replay counts and hard files before claiming a
    // fully resolved population. Reaching this line alone is not that claim.
    Ok(())
}

#[cfg(test)]
mod linked_roots_tests {
    use super::*;
    #[test]
    fn square_roots_exhaust_all_s5() {
        fn permutations(d: usize, depth: usize, p: Perm, out: &mut Vec<Perm>) {
            if depth == d {
                out.push(p);
                return;
            }
            for i in depth..d {
                let mut q = p;
                q.img.swap(depth, i);
                permutations(d, depth + 1, q, out);
            }
        }
        let mut perms = Vec::new();
        permutations(5, 0, Perm::identity(5), &mut perms);
        for q in &perms {
            let exact: HashSet<_> = perms
                .iter()
                .copied()
                .filter(|r| r.compose(*r, 5) == *q)
                .collect();
            assert_eq!(lr_roots(*q, 5).into_iter().collect::<HashSet<_>>(), exact);
        }
    }
    #[test]
    fn regular_frame_roots_and_masses() {
        let (mut h, mut k) = (Perm::identity(18), Perm::identity(18));
        for b in 0..2 {
            for i in 0..3 {
                for j in 0..3 {
                    let x = b * 9 + i * 3 + j;
                    h.img[x] = (b * 9 + ((i + 1) % 3) * 3 + j) as u8;
                    k.img[x] = (b * 9 + i * 3 + (j + 1) % 3) as u8;
                }
            }
        }
        let (left, right) = (lr_roots(k, 18), lr_roots(h, 18));
        assert_eq!((left.len(), right.len()), (856, 856));
        assert_eq!(left.iter().filter(|r| !r.is_even(18)).count(), 450);
        let reps = lr_pair_orbits(&left, &right, &lr_centralizer(h, k, 18).unwrap(), 18);
        assert_eq!(reps.iter().map(|r| r.1).sum::<usize>(), 732736);
        for (root, _) in reps {
            let t = lr_tuple(h, k, left[root / 856], right[root % 856], 18);
            assert!(tuple_product_is_identity(18, &t));
            for (mv, g) in [(0, h), (2, k)] {
                let moved = hurwitz_move_dynamic(18, t, mv);
                assert!((0..4).all(|i| t[i].conjugate_by(g, 18) == moved[i]));
            }
        }
    }
}
