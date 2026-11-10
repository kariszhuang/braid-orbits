//! Shared discovery and benchmark validation behind Metal screening.
use super::*;

pub(crate) struct Session {
    parents: Vec<(usize, Vec<[Perm;4]>)>,
    types: HashMap<usize, Vec<Vec<usize>>>,
    cache: Mutex<Cache>,
    ledger: Mutex<File>,
    pending: Mutex<File>,
    results: Mutex<File>,
    out: PathBuf,
}
struct Cache {
    seen: HashSet<(usize,Key)>,
    bad: HashSet<(usize,Key)>,
    good: HashSet<(usize,Key)>,
}

pub(crate) fn open(out:&Path,known:&Path)->Result<Session,String> {
    fs::create_dir_all(out).map_err(|e|e.to_string())?;
    let mut parents=Vec::new();let mut bad=HashSet::new();
    for name in ["CERTIFIED_S18_genus1_weak.txt","CERTIFIED_S18_XXYY_genus1_weak.txt","CERTIFIED_A25_weak.txt"] {
        let path=known.join(name);verify(&path)?;let case=parse_certificate_case(&path)?;
        let states=case.states.unwrap();bad.extend(states.iter().map(|&k|(case.d,k)));
        parents.push((case.d,states.iter().map(|&k|unpack(k,case.d)).collect()));
    }
    // Load all other curated weak components too, not just the mutation parents.
    let mut paths=Vec::new();collect_certificate_paths(known,&mut paths)?;
    for path in paths {
        let text=fs::read_to_string(&path).map_err(|e|e.to_string())?;
        if !text.contains("status=certified_weak_and_strong_failure") {continue;}
        verify(&path)?;let case=parse_certificate_case(&path)?;
        bad.extend(case.states.unwrap().into_iter().map(|k|(case.d,k)));
    }
    fn partitions(n:usize,min:usize,p:&mut Vec<usize>,out:&mut Vec<Vec<usize>>) {
        if n==0 {if p.len()>=3 {out.push(p.clone());}return;}
        for x in min..=n {p.push(x);partitions(n-x,x,p,out);p.pop();}
    }
    let mut types=HashMap::new();
    for d in 9..=17 {let mut all=Vec::new();partitions(d,1,&mut Vec::new(),&mut all);types.insert(d,all);}
    let ledger=OpenOptions::new().create_new(true).write(true).open(out.join("VALIDATION.tsv")).map_err(|e|e.to_string())?;
    let pending=OpenOptions::new().create_new(true).write(true).open(out.join("survivors.txt")).map_err(|e|e.to_string())?;
    let mut results=OpenOptions::new().create_new(true).write(true).open(out.join("RESULTS.tsv")).map_err(|e|e.to_string())?;
    writeln!(results,"status\tseed").map_err(|e|e.to_string())?;
    Ok(Session{parents,types,cache:Mutex::new(Cache{seen:HashSet::new(),bad,good:HashSet::new()}),ledger:Mutex::new(ledger),pending:Mutex::new(pending),results:Mutex::new(results),out:out.to_path_buf()})
}

pub(crate) fn generate(s:&Session,d:usize,family:u32,seed:u64,bytes:&mut[u8]) {
    let mut rng=SplitMix64::new(seed);
    for row in bytes.chunks_exact_mut(3*d) {
        let mut t=[Perm::identity(d);4];
        if family>=2 {
            let (n,states)=&s.parents[(family-2) as usize];assert_eq!(*n,d);
            t=states[rng.below(states.len())];
            // Both passport-preserving entry conjugations and passport-changing
            // entry multiplication; neither is a braid move. No injected controls.
            let changes=if rng.below(4)==0 {2} else {1};
            for _ in 0..changes {
                let entry=rng.below(3);let mut h=Perm::identity(d);
                let radius=[1,2,4,8,16][rng.below(5)];
                for _ in 0..radius {let x=rng.below(d);let y=(x+1+rng.below(d-1))%d;h.img.swap(x,y);}
                t[entry]=if rng.below(4)==0 {
                    if d==25 && !h.is_even(d) {h.img.swap(0,1);}
                    h.compose(t[entry],d)
                } else {t[entry].conjugate_by(h,d)};
            }
        } else if family==1 {
            let all=&s.types[&d];
            // Condition first-three cycle types on admitting a fourth index <=d-3.
            // No exact fourth passport or full-group claim until validation.
            loop {
                let ix:Vec<_>=(0..3).map(|_|rng.below(all.len())).collect();
                let sum:usize=ix.iter().map(|&i|d-all[i].len()).sum();
                if sum+d-3<2*d-2 {continue;}
                for j in 0..3 {t[j]=perm_of_type(&all[ix[j]],d,Some(&mut rng));}break;
            }
        } else {
            for p in &mut t[..3] {
                let mut labels=std::array::from_fn::<_,MAX_D,_>(|i|i);rng.shuffle(&mut labels[..d]);
                let k=1+rng.below(d/2);
                for pair in labels[..2*k].chunks_exact(2) {p.img.swap(pair[0],pair[1]);}
            }
        }
        for j in 0..3 {row[j*d..(j+1)*d].copy_from_slice(&t[j].img[..d]);}
    }
}

// Statistics: screened survivors, point-transitive, duplicate starts, known bad,
// proper group, new generating start, probe witness, exact witness, new verified
// failure, unresolved, explored states, cached good. Samples != components.
pub(crate) fn validate(s:&Session,d:usize,bytes:&[u8],flags:&[u32],workers:usize,cap:usize,millis:u64)->Result<[u64;12],String> {
    let next=AtomicU64::new(0);let totals:Vec<_>=(0..12).map(|_|AtomicU64::new(0)).collect();
    let errors=Mutex::new(Vec::new());
    thread::scope(|scope| {
        for _ in 0..workers {let next=&next;let totals=&totals;let errors=&errors;
            scope.spawn(move || {
                let add=|i:usize,n:u64|{totals[i].fetch_add(n,Ordering::Relaxed);};
                loop {
                    let i=next.fetch_add(1,Ordering::Relaxed) as usize;if i>=flags.len() {break;}if flags[i]!=1 {continue;}
                    add(0,1);let mut t=[Perm::identity(d);4];let row=&bytes[i*3*d..(i+1)*3*d];
                    for j in 0..3 {t[j].img[..d].copy_from_slice(&row[j*d..(j+1)*d]);}
                    t[3]=product_one_fourth(d,t[0],t[1],t[2]);
                    if !generators_transitive(d,&t[..3]) {continue;}add(1,1);
                    let key=canonical(&t,d);
                    {let mut c=s.cache.lock().unwrap();
                     if c.bad.contains(&(d,key)) {add(3,1);continue;}
                     if c.good.contains(&(d,key)) {add(11,1);continue;}
                     if !c.seen.insert((d,key)) {add(2,1);continue;}}
                    let ambient=if t[..3].iter().all(|g|g.is_even(d)) {Ambient::Alternating} else {Ambient::Symmetric};
                    if subgroup_order(d,&t[..3])!=ambient.target_order(d) {add(4,1);continue;}add(5,1);
                    let text=tuple_string_dynamic(d,&t);
                    let result=(|| -> Result<(),String> {
                        // Record BEFORE closure. The host has saved the entire
                        // batch already; amortize durable ledger flush per batch.
                        {let mut f=s.ledger.lock().unwrap();writeln!(f,"pending\t{}_{d}\tseed={text}",ambient.short()).map_err(|e|e.to_string())?;}
                        writeln!(s.pending.lock().unwrap(),"group={}_{} seed={text}",ambient.short(),d).map_err(|e|e.to_string())?;
                        let mut rng=SplitMix64::new(key[0] as u64^key[1] as u64^key[2] as u64);
                        if probe(d,t,&mut rng,64,true).0 {
                            add(6,1);s.cache.lock().unwrap().good.insert((d,key));
                            let mut f=s.ledger.lock().unwrap();writeln!(f,"probe_witness\t{}_{d}\tseed={text}",ambient.short()).map_err(|e|e.to_string())?;
                            writeln!(s.results.lock().unwrap(),"weak_found\t{text}").map_err(|e|e.to_string())?;return Ok(());
                        }
                        let r=explore_goal(d,t,cap,Instant::now()+Duration::from_millis(millis),true);add(10,r.states.len() as u64);
                        if r.status.ends_with("_found") {add(7,1);s.cache.lock().unwrap().good.extend(r.states.iter().map(|&k|(d,k)));}
                        else if r.status=="certified_weak_and_strong_failure" {
                            let representative=*r.states.iter().min().unwrap();let mut c=s.cache.lock().unwrap();
                            if c.bad.contains(&(d,representative)) {add(3,1);}
                            else {
                                let path=s.out.join(format!("CERTIFIED_{}{d}_{:032x}_{:032x}_{:032x}.txt",ambient.short(),representative[0],representative[1],representative[2]));
                                write_certificate(&path,d,ambient,&t,&r).map_err(|e|e.to_string())?;verify(&path)?;
                                c.bad.extend(r.states.iter().map(|&k|(d,k)));add(8,1);
                            }
                        } else {add(9,1);}
                        let mut f=s.ledger.lock().unwrap();writeln!(f,"{}\t{}_{d}\tstates={}\tseed={text}",r.status,ambient.short(),r.states.len()).map_err(|e|e.to_string())?;
                        writeln!(s.results.lock().unwrap(),"{}\t{text}",r.status).map_err(|e|e.to_string())?;
                        Ok(())
                    })();
                    if let Err(e)=result {errors.lock().unwrap().push(e);break;}
                }
            });
        }
    });
    s.ledger.lock().unwrap().sync_data().map_err(|e|e.to_string())?;
    s.pending.lock().unwrap().sync_data().map_err(|e|e.to_string())?;
    s.results.lock().unwrap().sync_data().map_err(|e|e.to_string())?;
    let errors=errors.into_inner().unwrap();if !errors.is_empty() {return Err(errors.join("; "));}
    Ok(std::array::from_fn(|i|totals[i].load(Ordering::Relaxed)))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hybrid_controls_and_deterministic_proposals() {
        let out=env::temp_dir().join(format!("hurwitz_hybrid_test_{}_{}",std::process::id(),std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        let s=open(&out,Path::new("evidence/certificates")).unwrap();
        for (d,states) in &s.parents {
            let t=states[0];let bytes:Vec<_>=t[..3].iter().flat_map(|p|p.img[..*d].iter().copied()).collect();
            let stats=validate(&s,*d,&bytes,&[1],2,1000,100).unwrap();assert_eq!(stats[3],1);assert_eq!(stats[8],0);
        }
        for (d,family) in [(9,0),(10,1),(11,1),(12,0),(14,1),(17,1),(18,2),(18,3),(25,4)] {
            let mut a=vec![0;3*d*64];let mut b=a.clone();generate(&s,d,family,42,&mut a);generate(&s,d,family,42,&mut b);assert_eq!(a,b);
            for row in a.chunks_exact(d) {let mut sorted=row.to_vec();sorted.sort_unstable();assert_eq!(sorted,(0..d as u8).collect::<Vec<_>>());}
        }
        drop(s);for name in ["VALIDATION.tsv","survivors.txt","RESULTS.tsv"] {fs::remove_file(out.join(name)).unwrap();}fs::remove_dir(out).unwrap();
    }
}
