// Proof-directed transversal constructions. No weak-positive filters, random
// sampling or assumption that the generated group is the full ambient group.
fn sp_gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 { (a, b) = (b, a % b); }
    a
}
fn sp_labels(p: Perm, n: usize) -> ([usize; MAX_D], usize) {
    let mut labels = [usize::MAX; MAX_D];
    let mut count = 0;
    for root in 0..n {
        if labels[root] != usize::MAX { continue; }
        let mut x = root;
        while labels[x] == usize::MAX {
            labels[x] = count;
            x = p.img[x] as usize;
        }
        count += 1;
    }
    (labels, count)
}
// Exact alternative to multiplying permutations: a pair is index-additive
// precisely when its cycle-incidence multigraph is a forest.
fn sp_forest(a: Perm, b: Perm, n: usize) -> bool {
    let (la, ca) = sp_labels(a, n);
    let (lb, _) = sp_labels(b, n);
    let mut parent = [0usize; 2 * MAX_D];
    for (i, p) in parent.iter_mut().enumerate() { *p = i; }
    fn root(parent: &mut [usize; 2 * MAX_D], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for x in 0..n {
        let u = root(&mut parent, la[x]);
        let v = root(&mut parent, ca + lb[x]);
        if u == v { return false; }
        parent[u] = v;
    }
    true
}
fn sp_passport_possible(mut indices: [usize; 4], n: usize) -> bool {
    indices.sort_unstable();
    (indices[0] + indices[3]).min(indices[1] + indices[2]) <= n - 1
}
fn sp_type(p: Perm, n: usize) -> Vec<usize> {
    let mut used = 0u64;
    let mut lengths = Vec::new();
    for root in 0..n {
        if used & (1 << root) != 0 { continue; }
        let mut x = root;
        let mut length = 0;
        while used & (1 << x) == 0 {
            used |= 1 << x;
            length += 1;
            x = p.img[x] as usize;
        }
        if length > 1 { lengths.push(length); }
    }
    lengths.sort_unstable();
    lengths
}
// Checks a braid-stable sufficient condition in the fixed block coordinates.
// The one block-local entry need not be a single cycle, nor have order three.
fn sp_hypothesis(t: &[Perm; 4], r: usize, k: usize) -> bool {
    let n = r * k;
    let mut local = 0;
    for p in t {
        let shift = p.img[0] as usize / k;
        if (0..n).any(|x| p.img[x] as usize / k != (x / k + shift) % r) {
            return false;
        }
        if shift == 0 {
            local += 1;
            let blocks: HashSet<_> = (0..n).filter(|&x| p.img[x] as usize != x).map(|x| x / k).collect();
            if blocks.len() != 1 { return false; }
        } else if sp_gcd(shift, r) != 1 || sp_type(*p, n) != vec![r; k] {
            return false;
        }
    }
    local == 1
}
// Canonical state keys relabel points, so the original contiguous blocks are
// no longer visible. Recover a C_r-equivariant coloring, never test canonical
// states against the original labels or treat a failed guessed coloring as
// a general mathematical exclusion.
fn sp_recover_hypothesis(t: &[Perm; 4], r: usize, k: usize) -> bool {
    let n=r*k;
    let small: Vec<_>=(0..4).filter(|&j|sp_type(t[j],n)!=vec![r;k]).collect();
    if small.len()!=1 { return false; }
    let local=small[0];
    let large: Vec<_>=(0..4).filter(|&j|j!=local).collect();
    for q in 1..r {
        if sp_gcd(q,r)!=1 || sp_gcd(q+1,r)!=1 { continue; }
        let values=[1,q,(r-(q+1)%r)%r];
        for order in [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]] {
            let mut shifts=[0;4];
            for j in 0..3 { shifts[large[j]]=values[order[j]]; }
            let mut colors=[usize::MAX;MAX_D];
            let mut queue=[0usize;MAX_D];
            colors[0]=0;let mut len=1;let mut head=0;let mut valid=true;
            while head<len && valid {
                let x=queue[head];head+=1;
                for j in 0..4 {
                    let y=t[j].img[x] as usize;let expected=(colors[x]+shifts[j])%r;
                    if colors[y]==usize::MAX { colors[y]=expected;queue[len]=y;len+=1; }
                    else if colors[y]!=expected { valid=false;break; }
                }
            }
            if !valid || len!=n { continue; }
            if (0..r).any(|c|colors[..n].iter().filter(|&&x|x==c).count()!=k) { continue; }
            let support: HashSet<_>=(0..n).filter(|&x|t[local].img[x] as usize!=x).map(|x|colors[x]).collect();
            if support.len()==1 { return true; }
        }
    }
    false
}
fn sp_construct(r: usize, k: usize, q: usize, y: Perm) -> [Perm; 4] {
    assert!(r >= 3 && r % 2 == 1 && r * k <= 25);
    assert!(q > 0 && q < r && sp_gcd(q, r) == 1 && sp_gcd(q + 1, r) == 1);
    let n = r * k;
    let mut x = Perm::identity(k);
    for i in 0..k { x.img[i] = ((i + 1) % k) as u8; }
    let z = y.compose(x, k).inverse(k);
    let mut voltage = vec![Perm::identity(k); r];
    voltage[0] = x;
    voltage[q] = y;
    voltage[r - q] = z;
    // Return map for b*c on block zero. Placing its inverse there makes
    // a*b*c have order r and hence d consist of k transversal r-cycles.
    let mut return_map = Perm::identity(k);
    let mut block = 0;
    for _ in 0..r {
        return_map = voltage[block].compose(return_map, k);
        block = (block + q + 1) % r;
    }
    assert_eq!(block, 0);
    let tau = return_map.inverse(k);
    let mut a = Perm::identity(n);
    a.img[..k].copy_from_slice(&tau.img[..k]);
    let mut b = Perm::identity(n);
    let mut c = b;
    for block in 0..r {
        for i in 0..k {
            b.img[block * k + i] = (((block + 1) % r) * k + i) as u8;
            c.img[block * k + i] = (((block + q) % r) * k + voltage[block].img[i] as usize) as u8;
        }
    }
    [a, b, c, a.compose(b, n).compose(c, n).inverse(n)]
}
fn sp_next(a: &mut [u8]) -> bool {
    let Some(i) = (0..a.len() - 1).rev().find(|&i| a[i] < a[i + 1]) else { return false; };
    let j = (i + 1..a.len()).rev().find(|&j| a[i] < a[j]).unwrap();
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}
pub fn special_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: Option<&str>| -> Result<String, String> {
        args.windows(2).find(|w| w[0] == name).map(|w| w[1].clone())
            .or_else(|| default.map(str::to_string)).ok_or_else(|| format!("missing {name}"))
    };
    let r: usize = arg("--blocks", None)?.parse().map_err(|_| "blocks")?;
    let k: usize = arg("--size", None)?.parse().map_err(|_| "size")?;
    let seconds: u64 = arg("--seconds", Some("600"))?.parse().map_err(|_| "seconds")?;
    let closure_seconds: u64 = arg("--closure-seconds", Some("3"))?.parse().map_err(|_| "closure seconds")?;
    let cap: usize = arg("--cap", Some("10000"))?.parse().map_err(|_| "cap")?;
    let representatives: usize = arg("--representatives", Some("2"))?.parse().map_err(|_| "representatives")?;
    if r < 3 || r % 2 == 0 || !(3..=8).contains(&k) || r*k > 25 || seconds == 0 || closure_seconds == 0 || cap == 0 || representatives == 0 {
        return Err("odd blocks >=3, size3..8, degree<=25 and positive limits required".into());
    }
    let out = PathBuf::from(arg("--out", None)?);
    fs::create_dir(&out).map_err(|e| e.to_string())?;
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    let n = r*k;
    let mut candidates = BufWriter::new(OpenOptions::new().create_new(true).write(true).open(out.join("CANDIDATES.tsv")).map_err(|e| e.to_string())?);
    writeln!(candidates, "candidate\tq\ty_rank\tsmall_type\tgenus\tseed").map_err(|e| e.to_string())?;
    let mut selected = Vec::new();
    let mut histogram = std::collections::BTreeMap::<(usize, Vec<usize>), (usize, usize)>::new();
    let mut enumerated = 0usize;
    let mut identity = 0usize;
    let mut old_template = 0usize;
    let mut complete = true;
    let qs: Vec<_> = (1..r).filter(|&q| sp_gcd(q,r)==1 && sp_gcd(q+1,r)==1).collect();
    // Classify the entire inexpensive population FIRST. Closure cannot steal
    // the construction deadline or lose a candidate before durable export.
    'outer: for &q in &qs {
        let mut images: Vec<_> = (0..k as u8).collect();
        let mut rank = 0;
        loop {
            if enumerated % 256 == 0 && Instant::now() >= deadline { complete = false; break 'outer; }
            let mut y = Perm::identity(k);
            y.img[..k].copy_from_slice(&images);
            let t = sp_construct(r,k,q,y);
            enumerated += 1;
            let kind = sp_type(t[0],n);
            if kind.is_empty() { identity += 1; }
            else if r == 3 && kind == [3] { old_template += 1; }
            else {
                let bucket = histogram.entry((q,kind.clone())).or_default();
                bucket.0 += 1;
                if bucket.1 < representatives {
                    assert!(tuple_product_is_identity(n,&t) && generators_transitive(n,&t));
                    assert!(sp_hypothesis(&t,r,k));
                    let value = sm_values(&t,n);
                    assert!(value.0 && value.1 == 0);
                    assert_eq!(value.2,1 + ((k*(r-3)+sm_index(t[0],n))/2) as i64);
                    let id = selected.len();
                    writeln!(candidates,"{id}\t{q}\t{rank}\t{kind:?}\t{}\t{}",value.2,tuple_string_dynamic(n,&t)).map_err(|e| e.to_string())?;
                    candidates.flush().map_err(|e| e.to_string())?;
                    selected.push((id,t));
                    bucket.1 += 1;
                }
            }
            rank += 1;
            if !sp_next(&mut images) { break; }
        }
    }
    let expected = qs.len() as u128 * factorial_u128(k);
    if complete { assert_eq!(enumerated as u128,expected); }
    let mut types = BufWriter::new(File::create(out.join("TYPES.tsv")).map_err(|e| e.to_string())?);
    writeln!(types,"q\tsmall_type\toccurrences\tselected").map_err(|e| e.to_string())?;
    for ((q,kind),(count,promoted)) in &histogram {
        writeln!(types,"{q}\t{kind:?}\t{count}\t{promoted}").map_err(|e| e.to_string())?;
    }
    types.flush().map_err(|e| e.to_string())?;
    let mut results = BufWriter::new(File::create(out.join("RESOLUTIONS.tsv")).map_err(|e| e.to_string())?);
    writeln!(results,"candidate\tstatus\tfolder\tstates").map_err(|e| e.to_string())?;
    let mut known = HashMap::<SmoothKey, String>::new();
    let mut closed = 0;
    let mut partial = 0;
    let mut aliases = 0;
    for (id,t) in &selected {
        let key = sm_key(t,n);
        if let Some(folder) = known.get(&key) {
            aliases += 1;
            writeln!(results,"{id}\talias\t{folder}\t-").map_err(|e| e.to_string())?;
        } else if Instant::now() >= deadline {
            partial += 1;
            writeln!(results,"{id}\tdeferred\t-\t0").map_err(|e| e.to_string())?;
        } else {
            let scan = sm_scan(t,n,cap,deadline.min(Instant::now()+Duration::from_secs(closure_seconds)));
            if scan.status == "nonsmooth" { return Err("transversal theorem contradicted by literal closure".into()); }
            let name = if scan.status == "closed_smooth" { format!("component_{closed:04}") } else { format!("partial_{id:04}") };
            for &s in &scan.states {
                assert!(sp_recover_hypothesis(&sm_unpack(s,n),r,k),"braid-stable hypothesis regression");
            }
            sm_save(t,n,&scan,&out.join(&name))?;
            if scan.status == "closed_smooth" {
                closed += 1;
                for &s in &scan.states { known.insert(s,name.clone()); }
            } else { partial += 1; }
            writeln!(results,"{id}\t{}\t{name}\t{}",scan.status,scan.states.len()).map_err(|e| e.to_string())?;
        }
        results.flush().map_err(|e| e.to_string())?;
    }
    let qualifying: usize = histogram.values().map(|x| x.0).sum();
    assert_eq!(enumerated,identity+old_template+qualifying);
    fs::write(out.join("COVERAGE.tsv"),format!("degree\tblocks\tsize\texpected\tenumerated\tidentity\told_template\tqualifying_occurrences\tselected\tclosed_components\taliases\tpartial_state_counts\tenumeration_complete\telapsed_seconds\n{n}\t{r}\t{k}\t{expected}\t{enumerated}\t{identity}\t{old_template}\t{qualifying}\t{}\t{closed}\t{aliases}\t{partial}\t{complete}\t{}\n",selected.len(),start.elapsed().as_secs_f64())).map_err(|e| e.to_string())?;
    println!("r={r} k={k} complete={complete} occurrences={enumerated} qualifying={qualifying} closed={closed} partial={partial}");
    Ok(())
}
#[cfg(test)]
mod special_tests {
    use super::*;
    #[test]
    fn forest_and_commuting_reductions_all_s5_pairs() {
        let mut permutations = Vec::new();
        let mut images: Vec<_> = (0..5).collect();
        loop {
            let mut p=Perm::identity(5);p.img[..5].copy_from_slice(&images);
            permutations.push(p); if !sp_next(&mut images) { break; }
        }
        for &a in &permutations { for &b in &permutations {
            let additive=sm_index(a,5)+sm_index(b,5)==sm_index(a.compose(b,5),5);
            assert_eq!(sp_forest(a,b,5),additive);
            if a.compose(b,5)==b.compose(a,5) {
                let disjoint=(0..5).all(|x|a.img[x] as usize==x || b.img[x] as usize==x);
                assert_eq!(additive,disjoint);
            }
        }}
    }
    #[test]
    fn passport_reduction_all_index_vectors() {
        for a in 0..8 { for b in 0..8 { for c in 0..8 { for d in 0..8 {
            let explicit=(a+b).min(c+d)<=7 && (a+c).min(b+d)<=7 && (a+d).min(b+c)<=7;
            assert_eq!(sp_passport_possible([a,b,c,d],8),explicit);
        }}}}
    }
    #[test]
    fn every_small_voltage_and_six_moves() {
        for r in [3,5,7] { for k in 3..=4 {
            if r*k>25 { continue; }
            let mut images: Vec<_>=(0..k as u8).collect();
            loop {
                let mut y=Perm::identity(k);y.img[..k].copy_from_slice(&images);
                for q in 1..r {
                    if sp_gcd(q,r)!=1 || sp_gcd(q+1,r)!=1 { continue; }
                    let t=sp_construct(r,k,q,y);
                    if t[0].is_identity(r*k) { continue; }
                    assert!(tuple_product_is_identity(r*k,&t));
                    assert!(generators_transitive(r*k,&t));
                    assert!(sp_hypothesis(&t,r,k));
                    for mv in 0..6 {
                        let h=hurwitz_move_dynamic(r*k,t,mv);
                        assert!(sp_hypothesis(&h,r,k));
                        assert!(sm_values(&h,r*k).0 && sm_values(&h,r*k).1==0);
                    }
                }
                if !sp_next(&mut images) { break; }
            }
        }}
    }
}
