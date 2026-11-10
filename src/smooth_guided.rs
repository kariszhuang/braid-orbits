// Two biased construction lanes, not exhaustive coverage. All final decisions
// use literal definitions and the same six-move closure as the exact census.
fn gd_find(parent: &mut [usize; MAX_D], mut x: usize) -> usize {
    while parent[x] != x { parent[x] = parent[parent[x]]; x = parent[x]; }
    x
}
fn gd_join(parent: &mut [usize; MAX_D], x: usize, y: usize) -> bool {
    let x=gd_find(parent,x); let y=gd_find(parent,y);
    if x==y { false } else { parent[x]=y; true }
}
fn gd_transpose(p: &mut Perm, x: usize, y: usize) { p.img.swap(x,y); }
fn gd_shuffle(n:usize,rng:&mut SplitMix64)->Perm {
    let mut p=Perm::identity(n); rng.shuffle(&mut p.img[..n]); p
}
fn gd_forest(n:usize,rng:&mut SplitMix64)->[Perm;2] {
    let mut parent=std::array::from_fn(|i|i);
    let mut pair=[Perm::identity(n);2];
    let edges=2+rng.below(n-2);
    let unequal=rng.below(2)==0;
    for j in 0..edges {
        loop {
            let x=rng.below(n); let y=rng.below(n);
            if gd_join(&mut parent,x,y) {
                let entry=if j<2 {j} else if unequal {1} else {rng.below(2)};
                gd_transpose(&mut pair[entry],x,y); break;
            }
        }
    }
    debug_assert_eq!(sm_index(pair[0],n)+sm_index(pair[1],n),sm_index(pair[0].compose(pair[1],n),n));
    pair
}
// Random partition with prescribed index: disjoint cycles, including fixed
// points, chosen before the expensive braid calculation.
fn gd_index_perm(n:usize,index:usize,rng:&mut SplitMix64)->Perm {
    let mut points:Vec<_>=(0..n).collect();rng.shuffle(&mut points);
    let mut p=Perm::identity(n); let mut used=0; let mut remaining=index;
    while remaining>0 {
        // Leave enough points for the still-unallocated index.
        let slack=n-used-remaining;
        let contribution=if slack<=1 {remaining} else {1+rng.below(remaining)};
        let len=contribution+1;
        for j in 0..len {p.img[points[used+j]]=points[used+(j+1)%len] as u8;}
        used+=len;remaining-=contribution;
    }
    debug_assert_eq!(sm_index(p,n),index); p
}
fn gd_block_pair(r:usize,k:usize,rng:&mut SplitMix64)->[Perm;2] {
    let n=r*k; let mut pair=[Perm::identity(n);2];
    let mut parent=std::array::from_fn(|i|i);
    let mut blocks:Vec<_>=(0..r).collect();rng.shuffle(&mut blocks);
    // A tree of block matchings: every point edge joins different components.
    // Both first entries move blocks, unlike the exhausted one-local-entry form.
    let unequal=rng.below(2)==0;
    for j in 1..r {
        let x=blocks[j];let y=blocks[rng.below(j)];let voltage=gd_shuffle(k,rng);
        let entry=if j<=2 {j-1} else if unequal {1} else {rng.below(2)};
        for u in 0..k {
            let a=x*k+u;let b=y*k+voltage.img[u] as usize;
            assert!(gd_join(&mut parent,a,b));gd_transpose(&mut pair[entry],a,b);
        }
    }
    // Local forest extensions create independent internal motions without
    // sacrificing the exact index-additive first pair or the block system.
    for _ in 0..rng.below(k) {
        for _ in 0..16 {
            let block=rng.below(r);let x=block*k+rng.below(k);let y=block*k+rng.below(k);
            if gd_join(&mut parent,x,y) {gd_transpose(&mut pair[rng.below(2)],x,y);break;}
        }
    }
    debug_assert_eq!(sm_index(pair[0],n)+sm_index(pair[1],n),sm_index(pair[0].compose(pair[1],n),n));
    pair
}
fn gd_block_perm(r:usize,k:usize,rng:&mut SplitMix64)->Perm {
    let n=r*k;let q=gd_shuffle(r,rng);let mut p=Perm::identity(n);
    for j in 0..r {let v=gd_shuffle(k,rng);for u in 0..k {
        p.img[j*k+u]=(q.img[j] as usize*k+v.img[u] as usize) as u8;
    }} p
}
fn gd_passport(t:&[Perm;4],n:usize)->bool {
    let mut i=t.map(|p|sm_index(p,n));i.sort_unstable();
    i[0]+i[3]<=n-1 || i[1]+i[2]<=n-1
}
pub fn guided_main(args:&[String])->Result<(),String> {
    let arg=|s:&str,d:&str|args.windows(2).find(|w|w[0]==s).map(|w|w[1].clone()).unwrap_or(d.into());
    let n:usize=arg("--d","11").parse().map_err(|_|"degree")?;
    let seconds:u64=arg("--seconds","60").parse().map_err(|_|"seconds")?;
    let seed:u64=arg("--seed","1").parse().map_err(|_|"seed")?;
    // Finite matched benchmark only; default preserves historical behavior.
    let attempts:u64=arg("--attempts","0").parse().map_err(|_|"attempts")?;
    let mode=arg("--lane","forest");
    if !(9..=14).contains(&n) || seconds==0 || !["forest","blocks"].contains(&mode.as_str()) {return Err("degree9..14, positive seconds, forest|blocks".into());}
    let shapes:Vec<_>=(3..=n/2).filter(|r|n%r==0).map(|r|(r,n/r)).collect();
    if mode=="blocks" && shapes.is_empty() {return Err("no admissible at-least-three-block shape".into());}
    let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh output required".into());}
    fs::create_dir(&out).map_err(|e|e.to_string())?;
    let mut rng=SplitMix64::new(seed);let start=Instant::now();let deadline=start+Duration::from_secs(seconds);
    let mut count=[0u64;8]; // constructed, scope, transitive, passport, short-reject, duplicate, closed, unresolved
    let mut seen=HashSet::new(); let mut closed=HashMap::new();
    let mut resolutions=BufWriter::new(File::create(out.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?);
    writeln!(resolutions,"candidate\tstatus\tword\tfolder\tseed").map_err(|e|e.to_string())?;
    while Instant::now()<deadline && (attempts==0 || count[0]<attempts) {
        let (pair,c)=if mode=="forest" {
            let pair=gd_forest(n,&mut rng);let index=1+rng.below(n-1);
            (pair,gd_index_perm(n,index,&mut rng))
        } else {
            let(r,k)=shapes[rng.below(shapes.len())];
            (gd_block_pair(r,k,&mut rng),gd_block_perm(r,k,&mut rng))
        };
        let t=[pair[0],pair[1],c,pair[0].compose(pair[1],n).compose(c,n).inverse(n)];count[0]+=1;
        if t.iter().any(|p|p.is_identity(n)) || t.iter().map(|p|sm_index(*p,n)).sum::<usize>()<=2*n-2 {continue;}count[1]+=1;
        if !generators_transitive(n,&t[..3]) {continue;}count[2]+=1;
        if !gd_passport(&t,n) {continue;}count[3]+=1;
        if sc_word(t,n).is_some() {count[4]+=1;continue;}
        let key=sm_key(&t,n);if closed.contains_key(&key) || seen.contains(&key) {count[5]+=1;continue;}
        if seen.len()<1000000 {seen.insert(key);}
        let id=count[0];
        // Durable seed precedes closure. A crashed last row remains a lead.
        let mut lead=OpenOptions::new().create(true).append(true).open(out.join("LEADS.tsv")).map_err(|e|e.to_string())?;
        writeln!(lead,"{id}\t{}",tuple_string_dynamic(n,&t)).map_err(|e|e.to_string())?;lead.flush().map_err(|e|e.to_string())?;
        let scan=sm_scan(&t,n,100000,deadline.min(Instant::now()+Duration::from_secs(10)));
        let mut folder=String::new();
        if scan.status!="nonsmooth" {
            folder=format!("{}_{id:09}",if scan.status=="closed_smooth" {"component"} else {"unresolved"});
            sm_save(&t,n,&scan,&out.join(&folder))?;
            if scan.status=="closed_smooth" {count[6]+=1;for &s in &scan.states {closed.insert(s,folder.clone());}} else {count[7]+=1;}
        }
        writeln!(resolutions,"{id}\t{}\t{}\t{folder}\t{}",scan.status,scan.word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(","),tuple_string_dynamic(n,&t)).map_err(|e|e.to_string())?;
        resolutions.flush().map_err(|e|e.to_string())?;
        println!("degree={n} lane={mode} candidate={id} status={} states={}",scan.status,scan.states.len());
    }
    fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\t{mode}\t{seed}\t{}\t{}\n",count.iter().map(|v|v.to_string()).collect::<Vec<_>>().join("\t"),start.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
    Ok(())
}
#[cfg(test)]
mod guided_tests {
    use super::*;
    #[test] fn constructors_preserve_geodesy_and_requested_index() {
        let mut rng=SplitMix64::new(41);
        for n in 9..=14 {for _ in 0..200 {
            let p=gd_forest(n,&mut rng);assert_eq!(sm_index(p[0],n)+sm_index(p[1],n),sm_index(p[0].compose(p[1],n),n));
            for index in 0..n {assert_eq!(sm_index(gd_index_perm(n,index,&mut rng),n),index);}
        }}
        for (r,k) in [(3,3),(6,2),(4,3),(7,2)] {for _ in 0..200 {
            let p=gd_block_pair(r,k,&mut rng);let n=r*k;
            assert_eq!(sm_index(p[0],n)+sm_index(p[1],n),sm_index(p[0].compose(p[1],n),n));
            for q in p {for j in 0..r {assert!((0..k).all(|u|q.img[j*k+u] as usize/k==q.img[j*k] as usize/k));}}
        }}
    }
    #[test] fn passport_filter_keeps_known_control() {
        let t=sm_parse("[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]",9).unwrap();
        assert!(gd_passport(&t,9));
    }
}
