// Stronger biased constructions. Neither failed construction nor finite-word
// survival is an exclusion/certificate. Share all permutation and BFS kernels.
fn dd_cycle_labels(p: Perm, n: usize) -> [u8; MAX_D] {
    let mut labels = [u8::MAX; MAX_D];
    for root in 0..n {
        if labels[root] != u8::MAX { continue; }
        let mut x = root;
        while labels[x] == u8::MAX {
            labels[x] = root as u8; x = p.img[x] as usize;
        }
    }
    labels
}

// At most C(MAX_D,2) unordered pairs. Byte endpoints suffice (MAX_D=34),
// reducing per-step clearing from 18 KiB to 1.1 KiB without unsafe memory or
// changing enumeration order, uniform selection, or the random stream.
const DD_PAIR_CAPACITY: usize = MAX_D * (MAX_D - 1) / 2;

fn dd_joint_third(frame: &ForestFrame, n: usize, target: usize,
                  rng: &mut SplitMix64) -> Option<Perm> {
    let mut c = Perm::identity(n);
    let mut products = frame.factors;
    for _ in 0..target {
        let mut labels = [[0u8; MAX_D]; 6];
        labels[0] = dd_cycle_labels(c, n);
        for j in 0..5 { labels[j+1] = dd_cycle_labels(products[j], n); }
        let mut choices = [(0u8, 0u8); DD_PAIR_CAPACITY]; let mut size = 0;
        for x in 0..n { for y in x+1..n {
            if labels.iter().all(|l| l[x] != l[y]) {
                choices[size] = (x as u8,y as u8); size += 1;
            }
        }}
        if size == 0 { return None; }
        let (x,y) = choices[rng.below(size)];
        let (x,y) = (x as usize,y as usize);
        // Right multiplication by (x y) swaps images. Connecting distinct
        // cycles increases each index by exactly one, including ind(c).
        c.img.swap(x,y);
        for p in &mut products { p.img.swap(x,y); }
    }
    Some(c)
}

fn dd_full_cycle(n: usize, rng: &mut SplitMix64) -> Perm {
    let mut order: Vec<_> = (0..n).collect(); rng.shuffle(&mut order);
    let mut p = Perm::identity(n);
    for j in 0..n { p.img[order[j]] = order[(j+1)%n] as u8; }
    p
}

fn dd_solved_shapes(n: usize) -> Vec<(usize,usize)> {
    (3..=n/2).filter(|r| r%2==1 && n%r==0 && (n/r)>=3)
        // Complete arbitrary S4-pair population C141; do not rediscover it.
        .map(|r|(r,n/r)).filter(|&(r,k)| (r,k)!=(3,4)).collect()
}

fn dd_solved_tuple(n: usize, rng: &mut SplitMix64) -> Option<([Perm;4],usize,usize)> {
    let shapes = dd_solved_shapes(n); if shapes.is_empty() { return None; }
    let (r,k) = shapes[rng.below(shapes.len())];
    let mut beta = Perm::identity(r);
    for j in 0..r { beta.img[j] = ((j+1)%r) as u8; }
    // Random full gamma; conditioning beta*gamma to be a full cycle. A bounded
    // failed draw is not a mathematical exclusion of that quotient action.
    let mut gamma = None;
    for _ in 0..128 {
        let g = dd_full_cycle(r,rng);
        if sm_index(beta.compose(g,r),r) != r-1 { continue; }
        if (r,k)==(5,3) {
            // Retire only the exact marked quotient of C140, modulo the
            // centralizer of the fixed beta. Same group order is NOT enough.
            let (_,old_gamma)=nc_quotients();let mut h=Perm::identity(r);
            let mut old=false;
            for _ in 0..r {
                if h.compose(old_gamma,r).compose(h.inverse(r),r)==g {old=true;break;}
                h=beta.compose(h,r);
            }
            if old {continue;}
        }
        gamma = Some(g); break;
    }
    let gamma = gamma?;
    let free: Vec<_> = (0..r-1).map(|_|{
        let mut p = Perm::identity(k); rng.shuffle(&mut p.img[..k]); p
    }).collect();
    let tuple = nc_construct(beta,gamma,r,k,&free);
    Some((tuple,r,k))
}

fn dd_read_negative(path: &Path, n: usize) -> Result<HashMap<SmoothKey,Vec<usize>>,String> {
    let raw = fs::read_to_string(path).map_err(|e|e.to_string())?;
    let mut cache = HashMap::new();
    for line in raw.lines() {
        let f: Vec<_> = line.splitn(2,'\t').collect();
        if f.len()!=2 { return Err("negative bank tuple,word".into()); }
        let t = sm_parse(f[0],n)?;
        let word:Vec<usize> = f[1].split(',').map(|s|s.parse().map_err(|_|"bank move"))
            .collect::<Result<_,_>>()?;
        if word.is_empty() || word.len()>64 || word.iter().any(|&m|m>=6)
            || !generators_transitive(n,&t) || !tuple_product_is_identity(n,&t)
            || !ff_replay_negative(t,n,&word) { return Err("invalid literal negative bank".into()); }
        if cache.len()>=65536 { return Err("negative bank exceeds bound".into()); }
        cache.insert(ff_key(&t,n),word);
    }
    Ok(cache)
}

#[cfg(test)]
mod discovery_design_tests {
    use super::*;
    #[test]
    fn joint_products_are_geodesic_not_a_heuristic() {
        let mut rng=SplitMix64::new(72921); let mut successes=0;
        for n in 9..=18 { for _ in 0..100 {
            let f=ForestFrame::new(gd_forest(n,&mut rng),n);
            let target=1+rng.below(n-1);
            if let Some(c)=dd_joint_third(&f,n,target,&mut rng) {
                successes+=1;assert_eq!(sm_index(c,n),target);
                for p in f.factors { assert_eq!(ff_product_index(p,c,n),sm_index(p,n)+target); }
            }
        }}
        assert!(successes>0);
    }
    #[test]
    fn generalized_voltage_constraints_and_all_neighbors() {
        let mut rng=SplitMix64::new(47321);let mut qualifying=0;
        for n in [18,20,21,24,25] { for _ in 0..100 {
            let (t,r,k)=dd_solved_tuple(n,&mut rng).unwrap();
            assert!(tuple_product_is_identity(n,&t));
            assert!(t[1..].iter().all(|p|sp_type(*p,n)==vec![r;k]));
            if t[0].is_identity(n) {continue;}
            assert!(nc_hypothesis(&t,r,k));
            if !generators_transitive(n,&t) {continue;}
            qualifying+=1;
            assert!(sm_values(&t,n).0);assert_eq!(sm_values(&t,n).1,0);
            for m in 0..6 {let u=hurwitz_move_dynamic(n,t,m);
                assert!(nc_hypothesis(&u,r,k));assert!(sm_values(&u,n).0);
            }
        }}
        assert!(qualifying>0);
        assert!(dd_solved_shapes(12).is_empty());
        assert!(dd_solved_tuple(15,&mut rng).is_some());
    }
}
