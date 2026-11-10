// Literal smoothness is a separate objective from weak pair intransitivity.
// Keys use rooted equivariant traversal; only transitive product-one tuples enter.
// Keep the original compact representation for existing binaries. The wide
// build is explicit; small-degree canonical order and exported identities stay
// unchanged, while high labels never shift beyond a u128 boundary.
const SM_KEY_LIMBS:usize=if cfg!(smooth_wide_keys){3}else{1};
const SM_KEY_WORDS:usize=3*SM_KEY_LIMBS;
const SM_MAX_DEGREE:usize=if cfg!(smooth_wide_keys){MAX_D}else{25};
type SmoothKey = [u128; SM_KEY_WORDS];
fn sm_packed_word(p:Perm,labels:&[u8;MAX_D],order:&[usize;MAX_D],n:usize)->[u128;SM_KEY_LIMBS]{
    let mut word=[0u128;SM_KEY_LIMBS];let bits=if n<=25{5}else{6};
    for i in 0..n{let value=labels[p.img[order[i]]as usize]as u128;let shift=bits*i;
        let limb=SM_KEY_LIMBS-1-shift/128;let offset=shift%128;
        word[limb]|=value<<offset;
        if offset+bits>128{word[limb-1]|=value>>(128-offset);}
    }word
}

fn sm_index(p: Perm, n: usize) -> usize {
    let mut seen = 0u64;
    let mut c = 0;
    for root in 0..n {
        if seen & (1 << root) != 0 {
            continue;
        }
        c += 1;
        let mut x = root;
        while seen & (1 << x) == 0 {
            seen |= 1 << x;
            x = p.img[x] as usize;
        }
    }
    n - c
}
fn sm_values(t: &[Perm; 4], n: usize) -> (bool, i64, i64, usize, usize, usize) {
    let indices = t.map(|p| sm_index(p, n));
    let ip = sm_index(t[0].compose(t[1], n), n);
    let l = point_orbit_count_perms(n, t[0], t[1]);
    let r = point_orbit_count_perms(n, t[2], t[3]);
    let delta = (n - ip) as i64 - l as i64 - r as i64 + 1;
    let g = 1 - n as i64 + indices.iter().sum::<usize>() as i64 / 2;
    (
        ip == indices[0] + indices[1] || ip == indices[2] + indices[3],
        delta,
        g,
        ip,
        l,
        r,
    )
}
fn sm_key(t: &[Perm; 4], n: usize) -> SmoothKey {
    assert!(n <= SM_MAX_DEGREE);
    let mut best = [u128::MAX; SM_KEY_WORDS];
    for root in 0..n {
        let mut labels = [u8::MAX; MAX_D];
        let mut order = [0usize; MAX_D];
        order[0] = root;
        labels[root] = 0;
        let mut len = 1;
        let mut head = 0;
        while head < len {
            let x = order[head];
            head += 1;
            for p in &t[..3] {
                let y = p.img[x] as usize;
                if labels[y] == u8::MAX {
                    labels[y] = len as u8;
                    order[len] = y;
                    len += 1;
                }
            }
        }
        assert_eq!(len, n);
        let mut key = [0u128; SM_KEY_WORDS];
        for j in 0..3 {
            key[j*SM_KEY_LIMBS..(j+1)*SM_KEY_LIMBS]
                .copy_from_slice(&sm_packed_word(t[j],&labels,&order,n));
        }
        best = best.min(key);
    }
    best
}
fn sm_unpack(key: SmoothKey, n: usize) -> [Perm; 4] {
    assert!(n<=SM_MAX_DEGREE);let bits=if n<=25{5}else{6};let mask=(1u128<<bits)-1;
    let mut t = [Perm::identity(n); 4];
    for j in 0..3 {
        for i in 0..n {
            let shift=bits*i;let limb=(j+1)*SM_KEY_LIMBS-1-shift/128;let offset=shift%128;
            let mut value=key[limb]>>offset;
            if offset+bits>128{value|=key[limb-1]<<(128-offset);}
            t[j].img[i]=(value&mask)as u8;
        }
    }
    t[3] = t[0].compose(t[1], n).compose(t[2], n).inverse(n);
    t
}
fn sm_parse(raw: &str, n: usize) -> Result<[Perm; 4], String> {
    let text = raw
        .trim()
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .ok_or("expected [four permutations]")?;
    let parts: Vec<_> = text.split(',').collect();
    if parts.len() != 4 {
        return Err("four entries required".into());
    }
    let mut t = [Perm::identity(n); 4];
    for (i, part) in parts.iter().enumerate() {
        let mut rest = part.trim();
        let mut used = 0u64;
        if rest.is_empty() {
            return Err("empty entry".into());
        }
        while !rest.is_empty() {
            rest = rest.strip_prefix('(').ok_or("expected (")?;
            let (cyc, after) = rest.split_once(')').ok_or("missing )")?;
            let points = cyc
                .split_whitespace()
                .map(|x| x.parse::<usize>().map_err(|_| "invalid point".to_string()))
                .collect::<Result<Vec<_>, _>>()?;
            for &x in &points {
                if x == 0 || x > n || used & (1 << (x - 1)) != 0 {
                    return Err("duplicate/out of range point".into());
                }
                used |= 1 << (x - 1);
            }
            for j in 0..points.len() {
                t[i].img[points[j] - 1] = (points[(j + 1) % points.len()] - 1) as u8;
            }
            rest = after.trim();
        }
    }
    if !tuple_product_is_identity(n, &t) {
        return Err("not product one".into());
    }
    if !generators_transitive(n, &t[..3]) {
        return Err("not transitive".into());
    }
    Ok(t)
}
struct SmoothScan {
    status: &'static str,
    states: Vec<SmoothKey>,
    word: Vec<usize>,
}
fn sm_scan(seed: &[Perm; 4], n: usize, cap: usize, deadline: Instant) -> SmoothScan {
    let root = sm_key(seed, n);
    let mut states = vec![root];
    let mut seen = HashMap::from([(root, 0usize)]);
    let mut parent = vec![(0usize, 0usize)];
    let mut head = 0;
    while head < states.len() {
        if Instant::now() >= deadline {
            return SmoothScan {
                status: "timeout",
                states,
                word: vec![],
            };
        }
        let t = sm_unpack(states[head], n);
        let values = sm_values(&t, n);
        assert!(values.1 >= 0 && values.1 <= values.2);
        if !values.0 {
            let mut word = vec![];
            let mut i = head;
            while i > 0 {
                let (p, m) = parent[i];
                word.push(m);
                i = p;
            }
            word.reverse();
            return SmoothScan {
                status: "nonsmooth",
                states,
                word,
            };
        }
        assert_eq!(values.1, 0, "smooth-implies-zero regression");
        for mv in 0..6 {
            let k = sm_key(&hurwitz_move_dynamic(n, t, mv), n);
            if seen.contains_key(&k) {
                continue;
            }
            if states.len() == cap {
                return SmoothScan {
                    status: "capped",
                    states,
                    word: vec![],
                };
            }
            seen.insert(k, states.len());
            parent.push((head, mv));
            states.push(k);
        }
        head += 1;
    }
    SmoothScan {
        status: "closed_smooth",
        states,
        word: vec![],
    }
}
fn sm_save(seed: &[Perm; 4], n: usize, result: &SmoothScan, out: &Path) -> Result<(), String> {
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let ambient = if seed.iter().all(|p| p.is_even(n)) {
        "A"
    } else {
        "S"
    };
    // u128 cannot represent arbitrary subgroup orders past 34!. Closure does
    // not need an order. Exact arbitrary-precision geometry profiling is later.
    let (order,full_group)=if n<=34{
        let order=smart::subgroup_order(n,&seed[..3]);
        let expected=factorial_u128(n)/if ambient=="A"{2}else{1};
        (order.to_string(),(order==expected).to_string())
    }else{("unknown".to_string(),"unknown".to_string())};
    let mut f = BufWriter::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("STATES.tsv"))
            .map_err(|e| e.to_string())?,
    );
    writeln!(
        f,
        "state\tsmooth\tdefect\tgenus\tproduct_index\tleft_orbits\tright_orbits\tedges\ttuple"
    )
    .map_err(|e| e.to_string())?;
    let index: HashMap<_, _> = result
        .states
        .iter()
        .enumerate()
        .map(|(i, k)| (*k, i))
        .collect();
    for (i, &key) in result.states.iter().enumerate() {
        let t = sm_unpack(key, n);
        let v = sm_values(&t, n);
        let edges = (0..6)
            .map(|m| {
                index
                    .get(&sm_key(&hurwitz_move_dynamic(n, t, m), n))
                    .map(|j| j.to_string())
                    .unwrap_or_else(|| "-".into())
            })
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            f,
            "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{edges}\t{}",
            v.0,
            v.1,
            v.2,
            v.3,
            v.4,
            v.5,
            tuple_string_dynamic(n, &t)
        )
        .map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    fs::write(out.join("SUMMARY.tsv"),format!("status\tdegree\tgenus\tstates\tgroup_order\tambient\tfull_group\tword\tseed\n{}\t{n}\t{}\t{}\t{order}\t{ambient}\t{full_group}\t{}\t{}\n",result.status,sm_values(seed,n).2,result.states.len(),result.word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(","),tuple_string_dynamic(n,seed))).map_err(|e|e.to_string())?;
    Ok(())
}
fn sm_type(p: Perm, n: usize, m: usize, f: usize) -> bool {
    let mut seen = 0u64;
    let mut triples = 0;
    let mut singles = 0;
    for root in 0..n {
        if seen & (1 << root) != 0 {
            continue;
        }
        let mut x = root;
        let mut len = 0;
        while seen & (1 << x) == 0 {
            seen |= 1 << x;
            len += 1;
            x = p.img[x] as usize;
        }
        match len {
            1 => singles += 1,
            3 => triples += 1,
            _ => return false,
        }
    }
    triples == m && singles == f
}
// Each partition/cycle orientation occurs once: least remaining point roots
// the next cycle, unordered other two points, then both directed orientations.
fn sm_triples<F: FnMut(Perm) -> bool>(
    n: usize,
    mask: u64,
    m: usize,
    f: usize,
    p: Perm,
    visit: &mut F,
) -> bool {
    if mask == 0 {
        assert_eq!((m, f), (0, 0));
        return visit(p);
    }
    let x = mask.trailing_zeros() as usize;
    let rem = mask & !(1 << x);
    if f > 0 && !sm_triples(n, rem, m, f - 1, p, visit) {
        return false;
    }
    if m > 0 {
        let mut ys = rem;
        while ys != 0 {
            let y = ys.trailing_zeros() as usize;
            ys &= ys - 1;
            let mut zs = ys;
            while zs != 0 {
                let z = zs.trailing_zeros() as usize;
                zs &= zs - 1;
                for (a, b) in [(y, z), (z, y)] {
                    let mut q = p;
                    q.img[x] = a as u8;
                    q.img[a] = b as u8;
                    q.img[b] = x as u8;
                    if !sm_triples(n, rem & !(1 << y) & !(1 << z), m - 1, f, q, visit) {
                        return false;
                    }
                }
            }
        }
    }
    true
}
fn sm_frame(m: usize, f: usize, shape: &str) -> Result<(Perm, Perm), String> {
    let n = 3 * m + f;
    let mut a = Perm::identity(n);
    a.img[0] = 1;
    a.img[1] = 2;
    a.img[2] = 0;
    let lens = match shape {
        "333" if m >= 3 => [3, 3, 3],
        "331" if m >= 2 && f == 1 => [3, 3, 1],
        _ => return Err("invalid frame shape".into()),
    };
    let mut b = Perm::identity(n);
    let mut next = 3;
    for root in 0..3 {
        if lens[root] == 3 {
            b.img[root] = next as u8;
            b.img[next] = (next + 1) as u8;
            b.img[next + 1] = root as u8;
            next += 2;
        }
    }
    let used = lens.iter().filter(|&&v| v == 3).count();
    for _ in used..m {
        b.img[next] = (next + 1) as u8;
        b.img[next + 1] = (next + 2) as u8;
        b.img[next + 2] = next as u8;
        next += 3;
    }
    assert!(sm_type(b, n, m, f));
    assert_eq!(sm_index(a.compose(b, n), n), 2 + 2 * m);
    Ok((a, b))
}
pub fn smooth_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: Option<&str>| -> Result<String, String> {
        args.windows(2)
            .find(|w| w[0] == name)
            .map(|w| w[1].clone())
            .or_else(|| default.map(str::to_string))
            .ok_or_else(|| format!("missing {name}"))
    };
    let out = PathBuf::from(arg("--out", None)?);
    let mode = arg("--mode", Some("profile"))?;
    let seconds: u64 = arg("--seconds", Some("300"))?
        .parse()
        .map_err(|_| "seconds")?;
    let cap: usize = arg("--cap", Some("50000"))?.parse().map_err(|_| "cap")?;
    if cap == 0 || seconds == 0 {
        return Err("positive limits required".into());
    }
    if mode == "profile" {
        let n: usize = arg("--d", None)?.parse().map_err(|_| "degree")?;
        if !(3..=25).contains(&n) {
            return Err("degree 3..25".into());
        }
        let t = sm_parse(&arg("--tuple", None)?, n)?;
        let result = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(seconds));
        sm_save(&t, n, &result, &out)?;
        println!("{}: {} states", result.status, result.states.len());
        return Ok(());
    }
    if mode == "block-family" {
        let k: usize = arg("--block-size", None)?
            .parse()
            .map_err(|_| "block size")?;
        if !(3..=8).contains(&k) {
            return Err("block size 3..8".into());
        }
        let n = 3 * k;
        let mut x = Perm::identity(k);
        let mut y = Perm::identity(k);
        for i in 0..k {
            x.img[i] = ((i + 1) % k) as u8;
        }
        y.img[0] = 1;
        y.img[1] = 0;
        let z = y.compose(x, k).inverse(k);
        let tau = x
            .inverse(k)
            .compose(y, k)
            .compose(x, k)
            .compose(y.inverse(k), k);
        let mut a = Perm::identity(n);
        let mut b = a;
        let mut c = a;
        for i in 0..k {
            a.img[i] = tau.img[i];
        }
        for block in 0..3 {
            for i in 0..k {
                b.img[block * k + i] = (((block + 1) % 3) * k + i) as u8;
                c.img[block * k + i] =
                    (((block + 1) % 3) * k + [x, y, z][block].img[i] as usize) as u8;
            }
        }
        let t = [a, b, c, a.compose(b, n).compose(c, n).inverse(n)];
        assert!(tuple_product_is_identity(n, &t));
        assert!(generators_transitive(n, &t));
        assert_eq!(sm_index(a, n), 2);
        assert!(t[1..].iter().all(|p| sm_type(*p, n, k, 0)));
        let result = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(seconds));
        sm_save(&t, n, &result, &out)?;
        println!(
            "block size={k}: {} {} states",
            result.status,
            result.states.len()
        );
        return Ok(());
    }
    if mode != "enumerate" {
        return Err("mode profile|enumerate|block-family".into());
    }
    let m: usize = arg("--m", None)?.parse().map_err(|_| "m")?;
    let f: usize = arg("--fixed", Some("0"))?.parse().map_err(|_| "fixed")?;
    let n = 3 * m + f;
    if !(3..=25).contains(&n) || f > 1 {
        return Err("invalid family".into());
    }
    let shape = arg("--shape", Some("333"))?;
    let (a, b) = sm_frame(m, f, &shape)?;
    fs::create_dir(&out).map_err(|e| e.to_string())?;
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    let mut count = [0usize; 7]; // enumerated,type,transitive,nonsmooth,closed,proper,capped
    let mut closed_seen = HashSet::new();
    let mut components = 0;
    let mut leads = BufWriter::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("LEADS.tsv"))
            .map_err(|e| e.to_string())?,
    );
    writeln!(leads, "candidate\tseed").map_err(|e| e.to_string())?;
    let mut resolutions = BufWriter::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("RESOLUTIONS.tsv"))
            .map_err(|e| e.to_string())?,
    );
    writeln!(resolutions, "candidate\tstatus\tword").map_err(|e| e.to_string())?;
    let mut failure = None;
    let complete = sm_triples(n, (1 << n) - 1, m, f, Perm::identity(n), &mut |c| {
        if Instant::now() >= deadline {
            return false;
        }
        count[0] += 1;
        let d = a.compose(b, n).compose(c, n).inverse(n);
        if !sm_type(d, n, m, f) {
            return true;
        }
        count[1] += 1;
        let t = [a, b, c, d];
        if !generators_transitive(n, &t[..3]) {
            return true;
        }
        count[2] += 1;
        let key = sm_key(&t, n);
        if closed_seen.contains(&key) {
            return true;
        }
        // Save before expensive closure so interruption never loses a lead.
        if writeln!(leads, "{}\t{}", count[0], tuple_string_dynamic(n, &t))
            .and_then(|_| leads.flush())
            .is_err()
        {
            failure = Some("lead write failed".into());
            return false;
        }
        let result = sm_scan(
            &t,
            n,
            cap,
            deadline.min(Instant::now() + Duration::from_secs(10)),
        );
        if writeln!(
            resolutions,
            "{}\t{}\t{}",
            count[0],
            result.status,
            result
                .word
                .iter()
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
        .and_then(|_| resolutions.flush())
        .is_err()
        {
            failure = Some("resolution write failed".into());
            return false;
        }
        match result.status {
            "nonsmooth" => count[3] += 1,
            "closed_smooth" => {
                count[4] += 1;
                closed_seen.extend(result.states.iter().copied());
                if smart::subgroup_order(n, &t[..3]) != factorial_u128(n) / 2 {
                    count[5] += 1;
                }
                let dir = out.join(format!("component_{components:04}"));
                components += 1;
                if let Err(e) = sm_save(&t, n, &result, &dir) {
                    failure = Some(e);
                    return false;
                }
            }
            _ => {
                count[6] += 1;
                let dir = out.join(format!("unresolved_{:06}", count[0]));
                if let Err(e) = sm_save(&t, n, &result, &dir) {
                    failure = Some(e);
                    return false;
                }
            }
        }
        if count[0] % 100000 == 0 {
            println!("{count:?} elapsed={:.2}", start.elapsed().as_secs_f64());
        }
        true
    });
    if let Some(e) = failure {
        return Err(e);
    }
    let expected =
        factorial_u128(n) / (3u128.pow(m as u32) * factorial_u128(m) * factorial_u128(f));
    if complete && count[0] as u128 != expected {
        return Err("enumeration mass mismatch".into());
    }
    fs::write(out.join("COVERAGE.tsv"),format!("degree\tm\tfixed\tshape\texpected\tenumerated\tmatching_fourth\ttransitive\tnonsmooth\tclosed_components\tproper_components\tunresolved\tenumeration_complete\tfully_resolved\telapsed_seconds\n{n}\t{m}\t{f}\t{shape}\t{expected}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{complete}\t{}\t{}\n",count[0],count[1],count[2],count[3],count[4],count[5],count[6],complete && count[6]==0,start.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
    println!("complete={complete} counts={count:?} expected={expected}");
    Ok(())
}
#[cfg(test)]
mod smooth_tests {
    use super::*;
    #[cfg(smooth_wide_keys)]
    #[test]fn compact_keys_preserve_old_order_and_roundtrip_all_supported_labels(){
        let mut rng=SplitMix64::new(386761);
        for n in 3..=SM_MAX_DEGREE{for _ in 0..100{
            let a=gd_shuffle(n,&mut rng);let b=gd_shuffle(n,&mut rng);let c=gd_shuffle(n,&mut rng);
            let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];if !generators_transitive(n,&t){continue;}
            let key=sm_key(&t,n);let fast=ff_key(&t,n);assert_eq!(key,fast);
            let u=sm_unpack(key,n);assert_eq!(sm_key(&u,n),key);
            assert_eq!(t.map(|p|sp_type(p,n)),u.map(|p|sp_type(p,n)));
            let relabel=gd_shuffle(n,&mut rng);assert_eq!(sm_key(&t.map(|p|p.conjugate_by(relabel,n)),n),key);
            if n<=25 {
                let mut old=[u128::MAX;3];
                for root in 0..n {
                    let mut map=[u8::MAX;MAX_D];let mut q=vec![root];map[root]=0;let mut head=0;
                    while head<q.len(){let x=q[head];head+=1;for p in &t[..3]{let y=p.img[x]as usize;
                        if map[y]==u8::MAX{map[y]=q.len()as u8;q.push(y);}
                    }}
                    let mut k=[0u128;3];for j in 0..3{for i in 0..n{k[j]|=(map[t[j].img[q[i]]as usize]as u128)<<(5*i);}}old=old.min(k);
                }
                let observed=std::array::from_fn(|j|key[(j+1)*SM_KEY_LIMBS-1]);
                assert_eq!(old,observed);
            }
        }}
    }
    #[test]
    fn triple_partition_mass() {
        for (m, f) in [(1, 0), (1, 1), (2, 0), (2, 1), (3, 0)] {
            let n = 3 * m + f;
            let mut unique = HashSet::new();
            assert!(sm_triples(
                n,
                (1 << n) - 1,
                m,
                f,
                Perm::identity(n),
                &mut |p| {
                    assert!(sm_type(p, n, m, f));
                    assert!(unique.insert(p));
                    true
                }
            ));
            assert_eq!(
                unique.len() as u128,
                factorial_u128(n) / (3u128.pow(m as u32) * factorial_u128(m) * factorial_u128(f))
            );
        }
    }
    #[test]
    fn reported_orbit() {
        let t = sm_parse(
            "[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]",
            9,
        )
        .unwrap();
        let r = sm_scan(&t, 9, 100, Instant::now() + Duration::from_secs(5));
        assert_eq!(r.status, "closed_smooth");
        assert_eq!(r.states.len(), 12);
        assert_eq!(sm_values(&t, 9).2, 2);
        assert_eq!(smart::subgroup_order(9, &t[..3]), 324);
    }
    #[test]
    fn frames_are_geodesic() {
        for (m, f, s) in [(2, 1, "331"), (3, 0, "333"), (3, 1, "331"), (4, 1, "333")] {
            let (a, b) = sm_frame(m, f, s).unwrap();
            assert_eq!(sm_index(a.compose(b, 3 * m + f), 3 * m + f), 2 + 2 * m);
        }
    }
}
