// Complete gauge-fixed degree12 family with arbitrary S4 local voltages.
// Existing exact permutation, smoothness, canonicalization and closure engines
// are reused; the fixed-full-cycle restriction of C138 is deliberately absent.
fn smaller_quotient() -> Perm {
    let mut p=Perm::identity(3);
    p.img[..3].copy_from_slice(&[1,2,0]);
    p
}
pub fn smaller_main(args: &[String]) -> Result<(), String> {
    let arg=|name:&str,default:Option<&str>|->Result<String,String>{
        args.windows(2).find(|w|w[0]==name).map(|w|w[1].clone())
            .or_else(||default.map(str::to_string)).ok_or_else(||format!("missing {name}"))
    };
    let cap:usize=arg("--cap",Some("10000"))?.parse().map_err(|_|"cap")?;
    let seconds:u64=arg("--closure-seconds",Some("5"))?.parse().map_err(|_|"seconds")?;
    let limit:usize=arg("--closures",Some("100"))?.parse().map_err(|_|"closures")?;
    if cap==0 || seconds==0 || limit==0 { return Err("positive limits required".into()); }
    let out=PathBuf::from(arg("--out",None)?);
    fs::create_dir(&out).map_err(|e|e.to_string())?;
    let start=Instant::now();let beta=smaller_quotient();let gamma=beta;
    let r=3;let k=4;let n=r*k;
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
        let v=sm_values(&t,n);assert!(v.0 && v.1==0 && v.2==2);
        // Check the definitions and theorem preservation on all six neighbors
        // of EVERY qualifying occurrence, not only closure representatives.
        for m in 0..6 {
            let h=hurwitz_move_dynamic(n,t,m);assert!(nc_hypothesis(&h,r,k));
            assert_eq!(sm_values(&h,n).0,true);assert_eq!(sm_values(&h,n).1,0);
        }
        qualifying+=1;let key=sm_key(&t,n);
        let id=if let Some(&id)=index.get(&key) { id } else {
            let id=selected.len();index.insert(key,id);selected.push(t);
            writeln!(seeds,"{id}\t{rank}\t2\t{}",tuple_string_dynamic(n,&t)).map_err(|e|e.to_string())?;id
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
mod smaller_tests {
    use super::*;
    #[test]
    fn exact_solver_matches_every_direct_frame() {
        let b=smaller_quotient();let local=nc_local_permutations(4);
        let literal=|t:[Perm;4]| t.map(|p|p.img[..12].to_vec());
        let mut solved=HashSet::new();
        for &x in &local { for &y in &local {
            solved.insert(literal(nc_construct(b,b,3,4,&[x,y])));
        }}
        let mut direct=HashSet::new();let mut frames=0;
        for &x in &local {for &y in &local {for &z in &local {for &a in &local {
            frames+=1;let t=nc_lift(b,b,4,&[x,y,z],a);
            if sp_type(t[2],12)==[3;4] && sp_type(t[3],12)==[3;4] {
                direct.insert(literal(t));
            }
        }}}}
        assert_eq!(frames,331776);assert_eq!(solved.len(),576);assert_eq!(solved,direct);
    }
    #[test]
    fn all_accepted_frames_and_literal_moves() {
        let b=smaller_quotient();let local=nc_local_permutations(4);
        let mut accepted=0;
        for &x in &local {for &y in &local {
            let t=nc_construct(b,b,3,4,&[x,y]);
            if t[0].is_identity(12) || !generators_transitive(12,&t) {continue;}
            accepted+=1;
            for m in 0..6 {
                let h=hurwitz_move_dynamic(12,t,m);
                assert!(nc_hypothesis(&h,3,4));
                let v=sm_values(&h,12);assert!(v.0 && v.1==0 && v.2==2);
            }
            let mut relabel=Perm::identity(12);
            for x in 0..12 {relabel.img[x]=(11-x) as u8;}
            let h=t.map(|p|relabel.compose(p,12).compose(relabel.inverse(12),12));
            assert_eq!(sm_key(&h,12),sm_key(&t,12));
        }}
        assert!(accepted>0);
    }
}
