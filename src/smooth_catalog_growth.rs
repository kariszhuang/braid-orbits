// Single-threaded workers. Independent audits run synchronously in place of
// the producer, so a two-worker controller has at most two heavy CPU jobs.
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn cg_epoch() -> f64 { SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64() }
fn cg_write(path: &Path, text: &str) -> Result<(), String> {
    let tmp=path.with_extension("tmp");
    { let mut f=File::create(&tmp).map_err(|e|e.to_string())?;
      f.write_all(text.as_bytes()).map_err(|e|e.to_string())?;
      f.sync_all().map_err(|e|e.to_string())?; }
    fs::rename(tmp,path).map_err(|e|e.to_string())
}
fn cg_profile(p: Perm,n:usize) -> Vec<usize> {
    let mut s:Vec<_>=qf_cycles(p,n).iter().map(Vec::len).collect(); s.sort_unstable(); s
}
fn cg_basic(t:&[Perm;4],n:usize)->bool {
    tuple_product_is_identity(n,t) && t.iter().all(|p|!p.is_identity(n))
        && generators_transitive(n,&t[..3]) && sm_values(t,n).2>=2
}
struct CgScan { states:Vec<SmoothKey>, edges:Vec<[u32;6]>, head:usize,
    status:&'static str, word:Vec<usize> }
fn cg_scan(t:&[Perm;4],n:usize,cap:usize,seconds:f64,stop:&Path)->CgScan {
    let start=Instant::now(); let root=ff_key(t,n);
    let mut states=vec![root]; let mut index=HashMap::from([(root,0u32)]);
    let mut edges=vec![[u32::MAX;6]];let mut parents=vec![(0usize,0usize)];let mut head=0;
    let mut status="closed_smooth";let mut word=vec![];
    'scan: while head<states.len() {
        if head%16==0 {
            if stop.exists(){status="stopped";break;}
            if start.elapsed().as_secs_f64()>=seconds {status="timeout";break;}
        }
        let u=sm_unpack(states[head],n);
        if !sm_values(&u,n).0 {
            status="nonsmooth";let mut at=head;
            while at!=0 {let (p,m)=parents[at];word.push(m);at=p;} word.reverse();break;
        }
        for m in 0..6 {
            let k=ff_key(&hurwitz_move_dynamic(n,u,m),n);
            let next=if let Some(&j)=index.get(&k){j}else{
                if states.len()>=cap{status="capped";break 'scan;}
                let j=states.len() as u32;index.insert(k,j);states.push(k);
                edges.push([u32::MAX;6]);parents.push((head,m));j
            };
            edges[head][m]=next;
        }
        head+=1;
    }
    CgScan{states,edges,head,status,word}
}
fn cg_save(folder:&Path,t:&[Perm;4],n:usize,s:&CgScan,source:&str)->Result<(),String>{
    fs::create_dir(folder).map_err(|e|e.to_string())?;
    let seed=tuple_string_dynamic(n,t);
    let genus=sm_values(t,n).2;
    // Order calculation is deliberately deferred: it is irrelevant to closure
    // and the legacy u128 group-order engine must not be called at large d.
    let ambient=if t.iter().all(|p|p.is_even(n)){"A"}else{"S"};
    fs::write(folder.join("SUMMARY.tsv"),format!("status\tdegree\tgenus\tstates\tgroup_order\tambient\tfull_group\tword\tseed\n{}\t{n}\t{genus}\t{}\tunknown\t{ambient}\tunknown\t\t{seed}\n",s.status,s.states.len())).map_err(|e|e.to_string())?;
    cg_write(&folder.join("INPUT.json"),&format!("{{\"degree\":{n},\"genus\":{genus},\"source\":\"{source}\",\"status\":\"{}\",\"states\":{},\"expanded\":{},\"epoch\":{}}}\n",s.status,s.states.len(),s.head,cg_epoch()))?;
    if s.status=="closed_smooth" {
        let mut out=BufWriter::new(File::create(folder.join("STATES.tsv")).map_err(|e|e.to_string())?);
        writeln!(out,"state\tsmooth\tdefect\tgenus\tproduct_index\tleft_orbits\tright_orbits\tedges\ttuple").map_err(|e|e.to_string())?;
        for (i,&key) in s.states.iter().enumerate(){
            let u=sm_unpack(key,n);let v=sm_values(&u,n);
            let e=s.edges[i].iter().map(|j|j.to_string()).collect::<Vec<_>>().join(",");
            writeln!(out,"{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{e}\t{}",v.0,v.1,v.2,v.3,v.4,v.5,tuple_string_dynamic(n,&u)).map_err(|e|e.to_string())?;
        }
        out.flush().map_err(|e|e.to_string())?;
    } else {
        // Retain the native frontier plus the raw root. Unfinished work is never
        // inserted into the positive cache or counted as a completed component.
        let mut out=BufWriter::new(File::create(folder.join("FRONTIER.bin")).map_err(|e|e.to_string())?);
        out.write_all(b"CGF1").map_err(|e|e.to_string())?;
        out.write_all(&(n as u64).to_le_bytes()).map_err(|e|e.to_string())?;
        out.write_all(&(s.head as u64).to_le_bytes()).map_err(|e|e.to_string())?;
        for key in &s.states {for limb in key {out.write_all(&limb.to_le_bytes()).map_err(|e|e.to_string())?;}}
        out.flush().map_err(|e|e.to_string())?;
    }
    Ok(())
}

fn cg_localized(r:usize,k:usize,rng:&mut SplitMix64,mixed:bool)->Option<[Perm;4]>{
    let n=r*k;
    let (beta,gamma)=if mixed {
        let b=if rng.below(2)==0{dd_full_cycle(r,rng)}else{let mut p=Perm::identity(r);rng.shuffle(&mut p.img[..r]);p};
        let mut c=Perm::identity(r);rng.shuffle(&mut c.img[..r]);(b,c)
    }else{
        let qs:Vec<_>=(1..r).filter(|&x|sp_gcd(x,r)==1&&sp_gcd(x+1,r)==1).collect();
        if qs.is_empty(){return None;}let q=qs[rng.below(qs.len())];
        let mut b=Perm::identity(r);let mut c=b;
        for j in 0..r {b.img[j]=((j+1)%r)as u8;c.img[j]=((j+q)%r)as u8;}(b,c)
    };
    let frame=QuotientForest::new(beta,gamma,r,k)?;
    let free:Vec<_>=frame.free_edges.iter().map(|_|{let mut p=Perm::identity(k);rng.shuffle(&mut p.img[..k]);p}).collect();
    let t=frame.solve(&free);
    if !cg_basic(&t,n){return None;}
    assert!(qf_hypothesis(&t,r,k),"unramified construction guard");
    Some(t)
}
struct CgBase { n:usize,id:String,t:[Perm;4] }
fn cg_bases(path:&Path)->Result<Vec<CgBase>,String>{
    let mut out=vec![];
    for line in fs::read_to_string(path).map_err(|e|e.to_string())?.lines(){
        let parts:Vec<_>=line.splitn(3,'\t').collect();if parts.len()!=3{return Err("base fields".into());}
        let n=parts[0].parse::<usize>().map_err(|_|"base degree")?;
        let mut t=sm_parse(parts[2],n)?;
        let a=(0..4).min_by_key(|&j|sm_index(t[j],n)).unwrap();
        for j in (0..a).rev(){t=hurwitz_move_dynamic(n,t,j+3);}
        if !sm_values(&t,n).0{return Err("base first-pair geodesy".into());}
        out.push(CgBase{n,id:parts[1].to_string(),t});
    }
    if out.is_empty(){return Err("empty double-lift bank".into());}Ok(out)
}
fn cg_double(base:&CgBase,rng:&mut SplitMix64)->Option<[Perm;4]>{
    let r=base.n;let n=2*r;let b=base.t;
    let p=b[0].compose(b[1],r);
    let roots=qf_cycles(p.compose(b[2],r),r).len();
    let frame=QuotientForest::tree(p,b[2],r,2,rng.below(roots))?;
    let branches=qf_cycles(b[2],r).len();let ram=rng.below(branches);
    let free:Vec<_>=frame.free_edges.iter().map(|_|{
        let mut p=Perm::identity(2);if rng.below(2)==0{p.img.swap(0,1);}p}).collect();
    let (voltage,_)=frame.voltages(&free,Some(ram));
    let mut t=[Perm::identity(n);4];
    for x in 0..r {for u in 0..2 {
        t[0].img[2*x+u]=(2*b[0].img[x]as usize+u)as u8;
        t[1].img[2*x+u]=(2*b[1].img[x]as usize+u)as u8;
        t[2].img[2*x+u]=(2*b[2].img[x]as usize+voltage[x].img[u]as usize)as u8;
    }}
    t[3]=t[0].compose(t[1],n).compose(t[2],n).inverse(n);
    if cg_basic(&t,n){Some(t)}else{None}
}
fn cg_mutate(mut t:[Perm;4],n:usize,rng:&mut SplitMix64)->[Perm;4]{
    for _ in 0..rng.below(4){t=hurwitz_move_dynamic(n,t,rng.below(6));}
    let j=rng.below(3);let x=rng.below(n);let mut y=rng.below(n-1);if y>=x{y+=1;}
    let mut h=Perm::identity(n);h.img.swap(x,y);
    if rng.below(2)==0 {t[j]=h.compose(t[j],n).compose(h,n);}
    else {t[j]=t[j].compose(h,n);}
    t[3]=t[0].compose(t[1],n).compose(t[2],n).inverse(n);t
}

pub fn catalog_growth_main(args:&[String])->Result<(),String>{
    let arg=|s:&str,default:&str|args.windows(2).find(|w|w[0]==s).map(|w|w[1].clone()).unwrap_or_else(||default.into());
    let run=PathBuf::from(arg("--run","")); let lane=arg("--lane","0").parse::<usize>().map_err(|_|"lane")?;
    let mut rng=SplitMix64::new(arg("--seed","8675309").parse().map_err(|_|"seed")?);
    if lane>3||SM_MAX_DEGREE<50{return Err("requires wide keys and lanes0..3".into());}
    let cap=arg("--cap","20000").parse::<usize>().map_err(|_|"cap")?;
    let seconds=arg("--closure-seconds","4").parse::<f64>().map_err(|_|"seconds")?;
    if let Some(i)=args.iter().position(|a|a=="--profile-seed") {
        let n=arg("--degree","9").parse::<usize>().map_err(|_|"profile degree")?;
        let t=sm_parse(args.get(i+1).ok_or("profile seed")?,n)?;
        let s=cg_scan(&t,n,cap,seconds,&run.join("STOP"));
        cg_save(&run,&t,n,&s,"profile_control")?;
        println!("{} {}",s.status,s.states.len());return Ok(());
    }
    let max_attempts=arg("--attempts","0").parse::<u64>().map_err(|_|"attempts")?;
    let python=arg("--python","python3");let audit=arg("--audit","");let auditor=arg("--auditor","");
    let draw_offset=arg("--draw-offset","0").parse::<u64>().map_err(|_|"draw offset")?;
    let bases=cg_bases(&PathBuf::from(arg("--bases","")))?;
    fs::create_dir_all(run.join(format!("lane_{lane}"))).map_err(|e|e.to_string())?;
    let state=run.join(format!("worker_{lane}.json"));let stop=run.join("STOP");
    // Extend the original degree25 catalog construction. A warmed membership
    // bank and the controller's exact IDs handle already published components.
    let shapes0=[(9,3),(7,4),(5,6),(3,9),(7,5),(9,4),(11,3),(5,7),(3,10),(5,8),(7,6),(9,5),(5,9),(3,11),(5,10),(7,7),(11,4),(13,3),(9,6),(3,12)];
    let shapes1=[(7,4),(8,4),(9,3),(9,4),(10,3),(6,5),(6,6),(12,3),(8,5),(10,4),(12,4),(9,5)];
    let mut known=HashSet::<(usize,SmoothKey)>::new();let mut attempted=HashSet::<(usize,SmoothKey)>::new();
    let mut pending_seen=HashSet::<(usize,SmoothKey)>::new();
    let mut successes:Vec<(usize,[Perm;4])>=vec![];let mut words:Vec<Vec<usize>>=vec![];
    if let Some(i)=args.iter().position(|a|a=="--warm-bank") {
        use std::io::{BufRead,BufReader};
        let path=args.get(i+1).ok_or("warm bank path")?;
        cg_write(&state,&format!("{{\"lane\":{lane},\"pid\":{},\"epoch\":{},\"stage\":\"warming_verified_cache\"}}\n",std::process::id(),cg_epoch()))?;
        let mut component_ids=HashSet::new();let file=BufReader::new(File::open(path).map_err(|e|e.to_string())?);
        for line in file.lines(){
            if stop.exists(){return Ok(());}
            let line=line.map_err(|e|e.to_string())?;let v:Vec<_>=line.splitn(3,'\t').collect();
            if v.len()!=3{return Err("warm-bank fields".into());}
            let n=v[0].parse::<usize>().map_err(|_|"warm degree")?;
            let t=sm_parse(v[2],n)?;
            if known.len()<500000{known.insert((n,ff_key(&t,n)));}
            if successes.len()<4096&&component_ids.insert(v[1].to_string()){successes.push((n,t));}
        }
        println!("warm cache lane={lane} verified_states={} mutation_seeds={}",known.len(),successes.len());
    }
    // Keep the trusted predecessor bank immutable when the bounded hot cache
    // rotates. Rotation may cause redundant work, never erase old membership.
    let warm_known=std::mem::take(&mut known);
    let (mut draws,mut qualified,mut duplicates,mut scans,mut verified,mut pending,mut negative,mut aliases)=(0u64,0u64,0u64,0u64,0u64,0u64,0u64,0u64);
    let mut max_states=0usize;let mut last_tick=Instant::now();let started=Instant::now();let mut stage="starting";
    let save=|stage:&str,rng:u64,draws:u64,qualified:u64,duplicates:u64,scans:u64,verified:u64,pending:u64,negative:u64,aliases:u64,max_states:usize,cache:usize|cg_write(&state,&format!("{{\"lane\":{lane},\"pid\":{},\"epoch\":{},\"stage\":\"{stage}\",\"rng\":{rng},\"draws\":{draws},\"qualified\":{qualified},\"duplicate_attempts\":{duplicates},\"closure_attempts\":{scans},\"verified_closures\":{verified},\"unresolved\":{pending},\"negative\":{negative},\"aliases\":{aliases},\"max_states\":{max_states},\"cached_states\":{cache},\"seconds\":{}}}\n",std::process::id(),cg_epoch(),started.elapsed().as_secs_f64()));
    let mut ready=OpenOptions::new().create(true).append(true).open(run.join(format!("READY_{lane}.jsonl"))).map_err(|e|e.to_string())?;
    while !stop.exists() && (max_attempts==0||draws<max_attempts){
        if last_tick.elapsed().as_secs_f64()>=5.0 {
            save(stage,rng.state,draws,qualified,duplicates,scans,verified,pending,negative,aliases,max_states,known.len()+warm_known.len())?;last_tick=Instant::now();
        }
        draws+=1;stage="constructing";
        let draw_number=draw_offset+draws;
        let (n,t,source)=if lane==0 {
            // Measured predecessor: 5408/5415 cyclic positives came from r=3.
            // Prioritize that productive family and extend its block sizes,
            // retaining a minority of other quotient shapes for exploration.
            let (r,k)=if draw_number%5!=0{(3,9+((draw_number/5)%12)as usize)}
                else{shapes0[(draw_number as usize/5)%shapes0.len()]};
            let Some(t)=cg_localized(r,k,&mut rng,false)else{continue;};
            (r*k,t,format!("cyclic_unramified_r{r}_k{k}"))
        }else if lane==2 {
            // A second unramified lane explores the full measured shape bank
            // without lane 0's heavy r=3 / k=9..20 weighting.
            let (r,k)=shapes0[rng.below(shapes0.len())];let Some(t)=cg_localized(r,k,&mut rng,false)else{continue;};
            (r*k,t,format!("cyclic_unramified_balanced_r{r}_k{k}"))
        }else if lane==1 && draw_number%10<5 || lane==3 && draw_number%3==1 {
            let base=&bases[rng.below(bases.len())];let Some(t)=cg_double(base,&mut rng)else{continue;};
            (2*base.n,t,format!("ramified_double_{}",base.id))
        }else if lane==1 && draw_number%10==9 && !successes.is_empty() || lane==3 && draw_number%3==0 && !successes.is_empty(){
            let (n,seed)=successes[rng.below(successes.len())];
            let t=cg_mutate(seed,n,&mut rng);if !cg_basic(&t,n){continue;}
            (n,t,"mutation_of_verified".into())
        }else{
            let (r,k)=if lane==3 {shapes1[(draw_number as usize/3)%shapes1.len()]}else{shapes1[rng.below(shapes1.len())]};
            let Some(t)=cg_localized(r,k,&mut rng,true)else{continue;};
            (r*k,t,format!("mixed_unramified_r{r}_k{k}"))
        };
        if !sm_values(&t,n).0{negative+=1;continue;}
        qualified+=1;let root=ff_key(&t,n);
        if warm_known.contains(&(n,root))||known.contains(&(n,root)){aliases+=1;continue;}
        // This is an alias of connected UNRESOLVED work, never a positive or
        // negative assertion. The saved root/frontier remains on the work list.
        if pending_seen.contains(&(n,root)){duplicates+=1;continue;}
        if !attempted.insert((n,root)){duplicates+=1;continue;}
        if attempted.len()>200000{attempted.clear();}
        let mut bad=false;
        for word in &words {let mut u=t;for &m in word{u=hurwitz_move_dynamic(n,u,m);}
            if !sm_values(&u,n).0{bad=true;break;}}
        if bad{negative+=1;continue;}
        scans+=1;stage="full_six_closure";
        save(stage,rng.state,draws,qualified,duplicates,scans,verified,pending,negative,aliases,max_states,known.len()+warm_known.len())?;
        let s=cg_scan(&t,n,cap,seconds,&stop);max_states=max_states.max(s.states.len());
        if s.status=="nonsmooth" {
            if words.len()<256{words.push(s.word.clone());}
            negative+=1;
            let path=run.join(format!("lane_{lane}/negative_{scans:09}.json"));
            let arrays:Vec<Vec<u8>>=t.iter().map(|p|p.img[..n].to_vec()).collect();
            cg_write(&path,&format!("{{\"degree\":{n},\"tuple\":{arrays:?},\"word\":{:?},\"source\":\"{source}\"}}\n",s.word))?;
            continue;
        }
        let folder=run.join(format!("lane_{lane}/candidate_{scans:09}"));
        cg_save(&folder,&t,n,&s,&source)?;
        if s.status!="closed_smooth" {
            if pending_seen.len()+s.states.len()>100000{pending_seen.clear();}
            for &key in &s.states{pending_seen.insert((n,key));}
            pending+=1;continue;
        }
        stage="independent_audit";
        save(stage,rng.state,draws,qualified,duplicates,scans,verified,pending,negative,aliases,max_states,known.len()+warm_known.len())?;
        let output=Command::new(&python).arg(&audit).arg("--folder").arg(&folder).arg("--auditor").arg(&auditor).output().map_err(|e|e.to_string())?;
        fs::write(folder.join("AUDIT.log"),[output.stdout,output.stderr].concat()).map_err(|e|e.to_string())?;
        if !output.status.success(){return Err(format!("independent audit failed: {}",folder.display()));}
        let metadata=folder.join("VERIFIED.json");
        if !metadata.is_file(){return Err("missing independent verification marker".into());}
        verified+=1;
        if known.len()+s.states.len()>500000{known.clear();}
        for &key in &s.states{known.insert((n,key));}
        if successes.len()<4096 {successes.push((n,t));}
        else {successes[verified as usize%4096]=(n,t);}
        writeln!(ready,"{{\"metadata\":\"{}\"}}",metadata.display()).map_err(|e|e.to_string())?;
        ready.sync_all().map_err(|e|e.to_string())?;
        println!("verified lane={lane} degree={n} genus={} states={} source={source}",sm_values(&t,n).2,s.states.len());
    }
    save("stopped",rng.state,draws,qualified,duplicates,scans,verified,pending,negative,aliases,max_states,known.len()+warm_known.len())?;
    Ok(())
}

#[cfg(test)]
mod catalog_growth_tests {
    use super::*;
    #[test] fn new_degree_constructors_and_wide_keys(){
        let mut rng=SplitMix64::new(431002);let mut good=0;
        for (r,k) in [(3,9),(3,13),(3,17),(3,20),(5,6),(7,4),(9,3),(11,3),(13,3)] {for _ in 0..6 {
            let Some(t)=cg_localized(r,k,&mut rng,false)else{continue;};let n=r*k;
            assert!(cg_basic(&t,n));assert!(qf_hypothesis(&t,r,k));let key=sm_key(&t,n);
            assert_eq!(sm_key(&sm_unpack(key,n),n),key);
            let mut u=t;for _ in 0..20{assert!(sm_values(&u,n).0);u=hurwitz_move_dynamic(n,u,rng.below(6));}
            good+=1;
        }}assert!(good>=20);
    }
    #[test] fn general_double_control_literally_smooth(){
        let t=sm_parse("[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]",9).unwrap();
        let b=CgBase{n:9,id:"control".into(),t};let mut rng=SplitMix64::new(1492);
        for _ in 0..20{let mut u=cg_double(&b,&mut rng).unwrap();assert_eq!(sm_values(&u,18).2,4);
            for _ in 0..30{assert!(sm_values(&u,18).0);u=hurwitz_move_dynamic(18,u,rng.below(6));}}
    }
}
