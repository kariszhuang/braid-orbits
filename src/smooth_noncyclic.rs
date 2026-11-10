// Gauge-fixed finite population: beta=(12345), gamma=(12453), k=3.
// Composition acts rightmost first. No cyclic-translation assumption.
fn nc_quotients() -> (Perm, Perm) {
    let mut beta = Perm::identity(5);
    beta.img[..5].copy_from_slice(&[1, 2, 3, 4, 0]);
    let mut gamma = Perm::identity(5);
    gamma.img[..5].copy_from_slice(&[1, 3, 0, 4, 2]);
    (beta, gamma)
}
fn nc_orbit(p: Perm, r: usize) -> Vec<usize> {
    let mut block = 0;
    let mut order = Vec::with_capacity(r);
    for _ in 0..r {
        assert!(!order.contains(&block));
        order.push(block);
        block = p.img[block] as usize;
    }
    assert_eq!(block, 0);
    order
}
fn nc_lift(beta: Perm, gamma: Perm, k: usize, voltage: &[Perm], tau: Perm) -> [Perm; 4] {
    let n = voltage.len() * k;
    let mut a = Perm::identity(n);
    a.img[..k].copy_from_slice(&tau.img[..k]);
    let mut b = Perm::identity(n);
    let mut c = b;
    for j in 0..voltage.len() { for u in 0..k {
        b.img[j*k+u] = (beta.img[j] as usize*k+u) as u8;
        c.img[j*k+u] = (gamma.img[j] as usize*k+voltage[j].img[u] as usize) as u8;
    }}
    [a, b, c, a.compose(b,n).compose(c,n).inverse(n)]
}
fn nc_construct(beta: Perm, gamma: Perm, r: usize, k: usize, free: &[Perm]) -> [Perm; 4] {
    assert_eq!(free.len(), r-1);
    let order = nc_orbit(gamma,r);
    let mut voltage = vec![Perm::identity(k);r];
    let mut return_map = Perm::identity(k);
    for (j,&v) in free.iter().enumerate() {
        voltage[order[j]] = v;
        return_map = v.compose(return_map,k);
    }
    voltage[order[r-1]] = return_map.inverse(k);
    let mut return_map = Perm::identity(k);
    for j in nc_orbit(beta.compose(gamma,r),r) {
        return_map = voltage[j].compose(return_map,k);
    }
    nc_lift(beta,gamma,k,&voltage,return_map.inverse(k))
}
// Exact sufficient condition in the given coordinates, preserved by every
// literal Hurwitz move. Block permutations may generate any group.
fn nc_hypothesis(t: &[Perm;4], r: usize, k: usize) -> bool {
    let n=r*k;
    let mut small=0;
    for p in t {
        let mut q=Perm::identity(r);
        for j in 0..r {
            q.img[j]=(p.img[j*k] as usize/k) as u8;
            if (0..k).any(|u|p.img[j*k+u] as usize/k!=q.img[j] as usize) { return false; }
        }
        if q.is_identity(r) {
            small+=1;
            let blocks: HashSet<_>=(0..n).filter(|&x|p.img[x] as usize!=x).map(|x|x/k).collect();
            if blocks.len()!=1 { return false; }
        } else if sp_type(q,r)!=[r] || sp_type(*p,n)!=vec![r;k] { return false; }
    }
    small==1
}
fn nc_local_permutations(k: usize) -> Vec<Perm> {
    let mut out=Vec::new();let mut images:Vec<_>=(0..k as u8).collect();
    loop {
        let mut p=Perm::identity(k);p.img[..k].copy_from_slice(&images);out.push(p);
        if !sp_next(&mut images) { break; }
    }
    out
}
pub fn noncyclic_main(args: &[String]) -> Result<(), String> {
    let arg=|name:&str,default:Option<&str>|->Result<String,String>{
        args.windows(2).find(|w|w[0]==name).map(|w|w[1].clone())
            .or_else(||default.map(str::to_string)).ok_or_else(||format!("missing {name}"))
    };
    let cap:usize=arg("--cap",Some("10000"))?.parse().map_err(|_|"cap")?;
    let seconds:u64=arg("--closure-seconds",Some("5"))?.parse().map_err(|_|"seconds")?;
    let limit:usize=arg("--closures",Some("4"))?.parse().map_err(|_|"closures")?;
    if cap==0 || seconds==0 || limit==0 { return Err("positive limits required".into()); }
    let out=PathBuf::from(arg("--out",None)?);
    fs::create_dir(&out).map_err(|e|e.to_string())?;
    let start=Instant::now();let (beta,gamma)=nc_quotients();
    let r=5;let k=3;let n=r*k;
    let local=nc_local_permutations(k);let expected=local.len().pow((r-1) as u32);
    let mut frames=BufWriter::new(File::create(out.join("FRAMES.tsv")).map_err(|e|e.to_string())?);
    let mut seeds=BufWriter::new(File::create(out.join("SEEDS.tsv")).map_err(|e|e.to_string())?);
    writeln!(frames,"rank\tstatus\tcandidate").map_err(|e|e.to_string())?;
    writeln!(seeds,"candidate\trank\tgenus\tseed").map_err(|e|e.to_string())?;
    let mut index=HashMap::new();let mut selected=Vec::new();let mut identity=0;let mut intransitive=0;let mut qualifying=0;
    for rank in 0..expected {
        let mut value=rank;let mut free=Vec::with_capacity(r-1);
        for _ in 0..r-1 { free.push(local[value%local.len()]);value/=local.len(); }
        let t=nc_construct(beta,gamma,r,k,&free);
        assert!(tuple_product_is_identity(n,&t));
        assert!(t[1..].iter().all(|p|sp_type(*p,n)==vec![r;k]));
        if t[0].is_identity(n) {
            identity+=1;writeln!(frames,"{rank}\tidentity\t-").map_err(|e|e.to_string())?;continue;
        }
        if !generators_transitive(n,&t) {
            intransitive+=1;writeln!(frames,"{rank}\tintransitive\t-").map_err(|e|e.to_string())?;continue;
        }
        assert!(nc_hypothesis(&t,r,k));
        let v=sm_values(&t,n);assert!(v.0 && v.1==0 && v.2==5);
        // Check the definitions and theorem preservation on all six neighbors
        // of EVERY qualifying occurrence, not only closure representatives.
        for m in 0..6 {
            let h=hurwitz_move_dynamic(n,t,m);assert!(nc_hypothesis(&h,r,k));
            assert_eq!(sm_values(&h,n).0,true);assert_eq!(sm_values(&h,n).1,0);
        }
        qualifying+=1;let key=sm_key(&t,n);
        let id=if let Some(&id)=index.get(&key) { id } else {
            let id=selected.len();index.insert(key,id);selected.push(t);
            writeln!(seeds,"{id}\t{rank}\t5\t{}",tuple_string_dynamic(n,&t)).map_err(|e|e.to_string())?;id
        };
        writeln!(frames,"{rank}\tqualifying\t{id}").map_err(|e|e.to_string())?;
    }
    frames.flush().map_err(|e|e.to_string())?;seeds.flush().map_err(|e|e.to_string())?;
    fs::write(out.join("CONSTRUCTION.tsv"),format!("expected\tenumerated\tidentity\tintransitive\tqualifying\tunique_seed_states\tcomplete\telapsed_seconds\n{expected}\t{expected}\t{identity}\t{intransitive}\t{qualifying}\t{}\ttrue\t{}\n",selected.len(),start.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
    println!("construction complete: {expected} choices, {qualifying} qualifying, {} seed states",selected.len());
    let mut results=BufWriter::new(File::create(out.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?);
    writeln!(results,"candidate\tstatus\tfolder\tstates").map_err(|e|e.to_string())?;
    let mut known=HashMap::<SmoothKey,String>::new();let mut scans=0;
    for (id,t) in selected.iter().enumerate() {
        if let Some(folder)=known.get(&sm_key(t,n)) {
            writeln!(results,"{id}\talias\t{folder}\t-").map_err(|e|e.to_string())?;
        } else if scans==limit {
            writeln!(results,"{id}\tdeferred\t-\t0").map_err(|e|e.to_string())?;
        } else {
            scans+=1;
            let scan=sm_scan(t,n,cap,Instant::now()+Duration::from_secs(seconds));
            if scan.status=="nonsmooth" { return Err("literal contradiction to block theorem".into()); }
            let name=if scan.status=="closed_smooth" {format!("component_{id:04}")} else {format!("partial_{id:04}")};
            sm_save(t,n,&scan,&out.join(&name))?;
            if scan.status=="closed_smooth" { for &s in &scan.states {known.insert(s,name.clone());} }
            writeln!(results,"{id}\t{}\t{name}\t{}",scan.status,scan.states.len()).map_err(|e|e.to_string())?;
            println!("candidate {id}: {} {} states",scan.status,scan.states.len());
        }
        results.flush().map_err(|e|e.to_string())?;
    }
    Ok(())
}
#[cfg(test)]
mod noncyclic_tests {
    use super::*;
    #[test]
    fn quotient_is_a5_not_cyclic() {
        let (b,c)=nc_quotients();
        assert_eq!(sp_type(b,5),[5]);assert_eq!(sp_type(c,5),[5]);
        assert_eq!(sp_type(b.compose(c,5),5),[5]);
        assert_eq!(smart::subgroup_order(5,&[b,c]),60);
    }
    #[test]
    fn exhaustive_constraint_solver_equals_direct_definitions() {
        let (b,c)=nc_quotients();let local=nc_local_permutations(3);
        // Literal image tuples, NOT canonical states: verifies occurrence
        // coverage too. The direct enumerator never solves a return equation.
        let literal=|t:[Perm;4]| t.map(|p|p.img[..15].to_vec());
        let mut solved=HashSet::new();
        for rank in 0usize..1296 {
            let mut value=rank;let mut free=Vec::new();
            for _ in 0..4 {free.push(local[value%6]);value/=6;}
            solved.insert(literal(nc_construct(b,c,5,3,&free)));
        }
        let mut direct=HashSet::new();let mut frames=0;
        for rank in 0usize..7776 {
            let mut value=rank;let mut v=Vec::new();
            for _ in 0..5 {v.push(local[value%6]);value/=6;}
            for &a in &local {
                frames+=1;let t=nc_lift(b,c,3,&v,a);
                if sp_type(t[2],15)==[5,5,5] && sp_type(t[3],15)==[5,5,5] {direct.insert(literal(t));}
            }
        }
        assert_eq!(frames,46656);assert_eq!(solved.len(),1296);assert_eq!(solved,direct);
    }
    #[test]
    fn all_six_moves_and_relabeling() {
        let (b,c)=nc_quotients();let local=nc_local_permutations(3);let mut accepted=0;
        for rank in 0usize..1296 {
            let mut value=rank;let mut free=Vec::new();
            for _ in 0..4 {free.push(local[value%6]);value/=6;}
            let t=nc_construct(b,c,5,3,&free);
            if t[0].is_identity(15) || !generators_transitive(15,&t) {continue;}
            accepted+=1;assert!(nc_hypothesis(&t,5,3));assert_eq!(sm_values(&t,15).2,5);
            for m in 0..6 {
                let h=hurwitz_move_dynamic(15,t,m);assert!(nc_hypothesis(&h,5,3));
                assert!(sm_values(&h,15).0);assert_eq!(sm_values(&h,15).1,0);
            }
            let mut relabel=Perm::identity(15);
            for x in 0..15 {relabel.img[x]=(14-x) as u8;}
            let h=t.map(|p|relabel.compose(p,15).compose(relabel.inverse(15),15));
            assert_eq!(sm_key(&h,15),sm_key(&t,15));
        }
        assert!(accepted>0);
    }
}
