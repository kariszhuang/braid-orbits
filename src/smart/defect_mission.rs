//! Bounded defect missions. A finite positive witness rejects an all-zero claim;
//! only checked six-move closure certifies it. Reuse the sparse group engine.
use super::*;

struct Scan {
    status: &'static str,
    states: Vec<Key>,
    histogram: BTreeMap<(i64, i64, i64), usize>,
    witness: Option<[Perm; 4]>,
}

fn checked_defect(t: &[Perm; 4], d: usize) -> Result<(i64, i64, i64), String> {
    let values = defect_of(t, d);
    if values.0 < 0 || values.1 < 0 || values.2 < 0
        || values.0 + values.1 + values.2 != genus_of(t, d) {
        return Err("defect/pair-genus decomposition failed".into());
    }
    // sign(p)=(-1)^index(p); this parity identity is valid in both A and S.
    let parity = (d - index(t[0].compose(t[1], d), d)) as i64
        - point_orbit_count_perms(d, t[0], t[1]) as i64
        - point_orbit_count_perms(d, t[2], t[3]) as i64 + 1;
    if values.0.rem_euclid(2) != parity.rem_euclid(2) {
        return Err("defect parity identity failed".into());
    }
    Ok(values)
}

fn scan(case: &SpectrumCase, cap: usize, deadline: Instant, screen: bool) -> Result<Scan, String> {
    let root = canonical(&case.seed, case.d);
    let mut result = Scan { status: "capped", states: vec![root],
        histogram: BTreeMap::new(), witness: None };
    let mut seen = HashSet::from([root]);
    let mut head = 0;
    while head < result.states.len() {
        if Instant::now() >= deadline { result.status = "timeout"; return Ok(result); }
        let t = unpack(result.states[head], case.d);
        head += 1;
        let values = checked_defect(&t, case.d)?;
        *result.histogram.entry(values).or_default() += 1;
        if values.0 > 0 {
            result.witness.get_or_insert(t);
            if screen { result.status = "positive_defect_witness"; return Ok(result); }
        }
        // Finite closure under positive generators implies inverse closure.
        for mv in 0..3 {
            let key = canonical(&hurwitz_move_dynamic(case.d, t, mv), case.d);
            if seen.contains(&key) { continue; }
            if seen.len() == cap { return Ok(result); }
            seen.insert(key); result.states.push(key);
        }
    }
    result.status = "closed";
    Ok(result)
}

// Every full A_d generating tuple has trivial S_d centralizer (d>=4), so its
// S_d class splits into exactly two A_d classes. Edge-conjugator sign describes
// the two-sheet cover. An odd closed voltage joins them into one braid orbit.
fn inner_graph(d: usize, ambient: Ambient, graph: &SpectrumGraph)
    -> Result<(SpectrumGraph, usize), String> {
    if ambient == Ambient::Symmetric { return Ok((graph.clone(), 1)); }
    let n = graph.states.len();
    if n > u32::MAX as usize / 2 { return Err("inner lift exceeds u32 capacity".into()); }
    let mut rows = vec![[0u32; 6]; 2*n];
    for (i, &key) in graph.states.iter().enumerate() {
        let t = unpack(key, d);
        for mv in 0..3 {
            let j = graph.transitions[i][mv] as usize;
            let moved = hurwitz_move_dynamic(d, t, mv);
            let h = tuple_conjugator(d, unpack(graph.states[j], d), moved)?;
            let odd = usize::from(!h.is_even(d));
            for sheet in 0..2 {
                let from = 2*i+sheet; let to = 2*j+(sheet^odd);
                rows[from][mv] = to as u32;
                rows[to][mv+3] = from as u32;
            }
        }
    }
    let mut labels = vec![u32::MAX; 2*n];
    let mut queue = vec![0usize]; labels[0] = 0;
    for head in 0..2*n {
        if head == queue.len() { break; }
        for &next in &rows[queue[head]] {
            let j = next as usize;
            if labels[j] == u32::MAX { labels[j] = queue.len() as u32; queue.push(j); }
        }
    }
    let size = queue.len();
    if size != n && size != 2*n { return Err("invalid inner-class covering degree".into()); }
    let transitions = queue.iter().map(|&i| rows[i].map(|j| labels[j as usize])).collect();
    Ok((SpectrumGraph { states: Vec::new(), transitions }, 2*n/size))
}

fn zero_certificate(case: &SpectrumCase, graph: &SpectrumGraph, path: &Path) -> Result<(), String> {
    let mut out = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path).map_err(|e|e.to_string())?);
    writeln!(out, "HURWITZ_ZERO_DEFECT_V1\ngroup={}\nd={}\nseed={}\nstates={}",
        case.ambient.short(), case.d, tuple_string_dynamic(case.d, &case.seed), graph.states.len()).map_err(|e|e.to_string())?;
    for key in &graph.states { writeln!(out,"STATE\t{:032x}\t{:032x}\t{:032x}",key[0],key[1],key[2]).map_err(|e|e.to_string())?; }
    out.flush().map_err(|e|e.to_string())?;
    out.get_ref().sync_all().map_err(|e|e.to_string())
}

fn verify_zero(path: &Path) -> Result<(), String> {
    let text = fs::read_to_string(path).map_err(|e|e.to_string())?;
    if text.lines().next() != Some("HURWITZ_ZERO_DEFECT_V1") {return Err("unsupported defect certificate".into());}
    let value = |name: &str| text.lines().find_map(|l|l.strip_prefix(name)).ok_or_else(||format!("missing {name}"));
    let d: usize = parse_num(value("d=")?.to_string(),"d")?;
    if !(4..=MAX_D).contains(&d) {return Err("invalid degree".into());}
    let ambient = Ambient::parse(value("group=")?)?;
    let seed = parse_tuple(value("seed=")?, d)?; validate_tuple(d, ambient, &seed)?;
    if genus_of(&seed,d)<=0 || seed.iter().any(|p|p.is_identity(d)) {return Err("certificate must be positive-genus and identity-free".into());}
    let states = text.lines().filter_map(|l|l.strip_prefix("STATE\t"))
        .map(|l|parse_hex_key(&l.split('\t').collect::<Vec<_>>())).collect::<Result<Vec<_>,_>>()?;
    let expected: usize = parse_num(value("states=")?.to_string(),"states")?;
    if states.len()!=expected || expected==0 {return Err("wrong state count".into());}
    let mut unique=HashSet::new();
    for &key in &states {
        if !unique.insert(key) {return Err("duplicate state".into());}
        let bits=if d<=16 {4} else {5};
        for &word in &key {
            if d<=25 {
                if word>>(bits*d)!=0 {return Err("out-of-range state".into());}
                let mut seen=0u64;
                for i in 0..d {
                    let x=((word>>(bits*i))&((1u128<<bits)-1)) as usize;
                    if x>=d || seen&(1u64<<x)!=0 {return Err("invalid encoded permutation".into());}
                    seen|=1u64<<x;
                }
            } else if word>=factorial_u128(d) {return Err("out-of-range permutation rank".into());}
        }
        let t=unpack(key,d);validate_tuple(d,ambient,&t)?;
        if canonical(&t,d)!=key || genus_of(&t,d)!=genus_of(&seed,d)
            || t.iter().any(|p|p.is_identity(d)) || checked_defect(&t,d)?.0!=0 {
            return Err("noncanonical, wrong-genus, identity or positive-defect state".into());
        }
    }
    let graph = graph_from_states(d, states)?;
    if expected!=graph.states.len() || graph_components(&graph)!=1
        || !graph.states.contains(&canonical(&seed,d)) {return Err("duplicate/disconnected/wrong-root states".into());}
    println!("VERIFIED zero defect {}{}: {} S_d-quotient states",ambient.short(),d,expected);
    Ok(())
}

fn write_case(case: &SpectrumCase, result: Scan, folder: &Path, state_rows: bool, spectra: bool)
    -> Result<(),String> {
    fs::create_dir(folder).map_err(|e|e.to_string())?;
    let create = |name:&str| File::create(folder.join(name)).map_err(|e|e.to_string());
    let mut hist=BufWriter::new(create("DEFECT.tsv")?);
    writeln!(hist,"defect\tleft_genus\tright_genus\tcount").map_err(|e|e.to_string())?;
    for ((delta,left,right),count) in &result.histogram {writeln!(hist,"{delta}\t{left}\t{right}\t{count}").map_err(|e|e.to_string())?;}
    hist.flush().map_err(|e|e.to_string())?;
    let mut inner_size=String::new(); let mut split=String::new();
    let mut bipartite=String::new();
    let mut lambda2=String::new(); let mut lambda_min=String::new();
    let mut graph=None;
    if result.status=="closed" {
        let quotient=graph_from_states(case.d,result.states.clone())?;
        if graph_components(&quotient)!=1 {return Err("disconnected quotient".into());}
        let (inner,splits)=inner_graph(case.d,case.ambient,&quotient)?;
        let n=inner.transitions.len(); inner_size=n.to_string(); split=splits.to_string();
        bipartite=graph_is_bipartite(&inner).to_string();
        if case.kind=="certificate" && case.input_states!=Some(quotient.states.len()) {
            return Err("certificate quotient size changed".into());
        }
        if case.kind=="exact_component" && case.input_states.is_some_and(|m|m!=n && m!=quotient.states.len()) {
            return Err(format!("{}: neither inner nor S_d size matches catalog",case.id));
        }
        if spectra && n>1 {
            let e=eigen_for_graph(&inner,512,384,case_seed(case,0x1aa3));
            lambda2=format_float(e.lambda2);lambda_min=format_float(e.lambda_min);
        }
        if quotient.states.len()<=512 {
            fs::write(folder.join("GRAPH_CODE.txt"),canonical_labeled_code(&quotient,512)).map_err(|e|e.to_string())?;
        }
        if result.witness.is_none() && genus_of(&case.seed,case.d)>0
            && !case.seed.iter().any(|p|p.is_identity(case.d)) {
            let path=folder.join("ZERO_DEFECT.txt");zero_certificate(case,&quotient,&path)?;verify_zero(&path)?;
        }
        graph=Some(quotient);
    }
    if state_rows {
        let mut rows=BufWriter::new(create("STATES.tsv")?);
        writeln!(rows,"state\tdefect\tleft_genus\tright_genus\tleft_point_orbits\tright_point_orbits\ttuple").map_err(|e|e.to_string())?;
        for (i,&key) in result.states.iter().enumerate() {
            let t=unpack(key,case.d);let (delta,left,right)=checked_defect(&t,case.d)?;
            let points=|gens:&[Perm]| point_orbits(case.d,gens).into_iter().map(|v|
                v.into_iter().map(|x|(x+1).to_string()).collect::<Vec<_>>().join(",")).collect::<Vec<_>>().join(";");
            writeln!(rows,"{i}\t{delta}\t{left}\t{right}\t{}\t{}\t{}",points(&t[..2]),points(&t[2..]),tuple_string_dynamic(case.d,&t)).map_err(|e|e.to_string())?;
        }
        rows.flush().map_err(|e|e.to_string())?;
    }
    let representative=graph.as_ref().and_then(|g|g.states.first()).unwrap_or(&result.states[0]);
    let measured=result.histogram.values().sum::<usize>();
    let histogram_unit=if result.status.starts_with("walk_") {"correlated_early_stopped_walk"}
        else if result.status=="closed" {"all_S_quotient_states"} else {"distinct_BFS_prefix"};
    let mut summary=BufWriter::new(create("SUMMARY.tsv")?);
    writeln!(summary,"case_id\tgroup\tdegree\tgenus\tidentity_entry\tstatus\thistogram_unit\tdiscovered_states\tmeasured_states\tS_quotient_states\tG_orbit_states\tinner_orbits_in_S_saturation\tinner_exact_bipartite\tinner_lambda2\tinner_lambda_min\trepresentative\twitness\tseed\tsource").map_err(|e|e.to_string())?;
    writeln!(summary,"{}\t{}\t{}\t{}\t{}\t{}\t{histogram_unit}\t{}\t{measured}\t{}\t{inner_size}\t{split}\t{bipartite}\t{lambda2}\t{lambda_min}\t{:032x}/{:032x}/{:032x}\t{}\t{}\t{}",
        case.id,case.ambient.short(),case.d,genus_of(&case.seed,case.d),case.seed.iter().any(|p|p.is_identity(case.d)),result.status,result.states.len(),
        if result.status=="closed" {result.states.len().to_string()} else {String::new()},representative[0],representative[1],representative[2],
        result.witness.as_ref().map(|t|tuple_string_dynamic(case.d,t)).unwrap_or_default(),tuple_string_dynamic(case.d,&case.seed),case.source).map_err(|e|e.to_string())?;
    summary.flush().map_err(|e|e.to_string())?;summary.get_ref().sync_all().map_err(|e|e.to_string())?;
    fs::write(folder.join("DONE.tmp"),result.status).map_err(|e|e.to_string())?;
    fs::rename(folder.join("DONE.tmp"),folder.join("DONE")).map_err(|e|e.to_string())
}

pub(super) fn run(args:&[String])->Result<(),String> {
    let mode=get_arg(args,"--mode",Some("screen"))?;
    if mode=="verify-zero" {return verify_zero(Path::new(&get_arg(args,"--certificate",None)?));}
    if !["screen","profile","walk"].contains(&mode.as_str()) {return Err("mode must be screen|profile|walk|verify-zero".into());}
    let cap:usize=parse_num(get_arg(args,"--cap",Some("250000"))?,"cap")?;
    let seconds:u64=parse_num(get_arg(args,"--seconds",Some("30"))?,"seconds")?;
    let threads:usize=parse_num(get_arg(args,"--threads",Some("5"))?,"threads")?;
    if cap==0 || seconds==0 || !(1..=128).contains(&threads) {return Err("positive cap/seconds and 1..128 threads required".into());}
    let mut cases=Vec::new();let mut roots=HashMap::new();
    if let Some(raw)=optional_arg(args,"--components") {
        for path in raw.split(',') {for case in parse_components_file(Path::new(path))? {add_spectrum_case(&mut cases,&mut roots,case);}}
    }
    if let Some(raw)=optional_arg(args,"--certificates") {
        let mut paths=Vec::new();for path in raw.split(',') {collect_certificate_paths(Path::new(path),&mut paths)?;}
        for path in paths {verify(&path)?;add_spectrum_case(&mut cases,&mut roots,parse_certificate_case(&path)?);}
    }
    if let Some(raw)=optional_arg(args,"--tuple") {
        let d:usize=parse_num(get_arg(args,"--d",None)?,"d")?;
        if !(4..=MAX_D).contains(&d) {return Err("invalid degree".into());}
        let ambient=Ambient::parse(&get_arg(args,"--group",None)?)?;
        let seed=parse_tuple(&raw,d)?;validate_tuple(d,ambient,&seed)?;
        cases.push(SpectrumCase {id:"explicit-tuple".into(),kind:"tuple".into(),source:"CLI".into(),ambient,d,seed,states:None,input_states:None,aliases:Vec::new()});
    }
    if let Some(path)=optional_arg(args,"--select-case-ids") {
        let text=fs::read_to_string(path).map_err(|e|e.to_string())?;
        let selected=text.lines().collect::<HashSet<_>>();cases.retain(|c|selected.contains(c.id.as_str()));
    }
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    if mode=="walk" {
        fs::create_dir_all(&out).map_err(|e|e.to_string())?;
        check_options(args,&out)?;
        return walks(args,&out,cap,seconds,threads);
    }
    if mode=="screen" || has_flag(args,"--positive-nonidentity") {
        cases.retain(|c|genus_of(&c.seed,c.d)>0 && !c.seed.iter().any(|p|p.is_identity(c.d)));
    }
    if cases.is_empty() {return Err("no eligible cases".into());}
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    check_options(args,&out)?;
    let next=std::sync::atomic::AtomicUsize::new(0);let errors=Mutex::new(Vec::<String>::new());
    thread::scope(|scope| {
        for _ in 0..threads.min(cases.len()) {
            let cases=&cases;let next=&next;let errors=&errors;let out=&out;let screen=mode=="screen";
            scope.spawn(move||loop {
                let i=next.fetch_add(1,Ordering::Relaxed);if i>=cases.len() {break;}
                let case=&cases[i];
                let folder=out.join(format!("case_{i:06}"));
                if folder.join("DONE").is_file() {continue;}
                // Do not overwrite an interrupted case; the runner preserves it
                // and creates a fresh attempt with unchanged inputs/binary.
                if folder.exists() {errors.lock().unwrap().push(format!("incomplete case requires fresh attempt: {}",folder.display()));continue;}
                let result=(|| {
                    let result=if let Some(states)=&case.states {
                        let graph=graph_from_states(case.d,states.clone())?;
                        let mut histogram=BTreeMap::new();let mut witness=None;
                        for &key in &graph.states {
                            let t=unpack(key,case.d);let v=checked_defect(&t,case.d)?;
                            *histogram.entry(v).or_default()+=1;
                            if v.0>0 {witness.get_or_insert(t);}
                        }
                        Scan {status:"closed",states:graph.states,histogram,witness}
                    } else {scan(case,cap,Instant::now()+Duration::from_secs(seconds),screen)?};
                    write_case(case,result,&folder,has_flag(args,"--state-rows"),has_flag(args,"--inner-spectra"))
                })();
                if let Err(e)=result {errors.lock().unwrap().push(format!("{}: {e}",case.id));}
                eprintln!("defect mission case={i}/{} id={}",cases.len(),case.id);
            });
        }
    });
    let errors=errors.into_inner().unwrap();
    if !errors.is_empty() {return Err(errors.join("; "));}
    fs::write(out.join("DONE"),format!("cases={}\n",cases.len())).map_err(|e|e.to_string())
}

fn check_options(args:&[String],out:&Path)->Result<(),String> {
    let path=out.join("OPTIONS.txt");let current=args[1..].join("\n");
    if path.exists() {
        if fs::read_to_string(path).map_err(|e|e.to_string())?!=current {return Err("saved mission options changed; use a new output".into());}
    } else {fs::write(path,current).map_err(|e|e.to_string())?;}
    Ok(())
}

fn walks(args:&[String],out:&Path,cap:usize,seconds:u64,threads:usize)->Result<(),String> {
    let d:usize=parse_num(get_arg(args,"--d",None)?,"d")?;
    if !(4..=MAX_D).contains(&d) {return Err("invalid degree".into());}
    let ambient=Ambient::parse(&get_arg(args,"--group",None)?)?;
    let draws:usize=parse_num(get_arg(args,"--draws",Some("128"))?,"draws")?;
    let steps:usize=parse_num(get_arg(args,"--walk-steps",Some("4096"))?,"steps")?;
    let master:u64=parse_num(get_arg(args,"--rng",Some("20261002"))?,"rng")?;
    let family=get_arg(args,"--family",Some("low-index"))?;
    if draws==0 || steps==0 || !["uniform","low-index","passport"].contains(&family.as_str()) {return Err("positive draws/steps and uniform|low-index|passport required".into());}
    let passport=if family=="passport" {Some(parse_passport(&get_arg(args,"--passport",None)?,d,ambient)?)} else {None};
    let next=std::sync::atomic::AtomicUsize::new(0);let errors=Mutex::new(Vec::<String>::new());
    thread::scope(|scope| {
        for _ in 0..threads.min(draws) {
            let next=&next;let errors=&errors;let family=&family;let passport=&passport;
            scope.spawn(move||loop {
                let i=next.fetch_add(1,Ordering::Relaxed);if i>=draws {break;}
                if out.join(format!("draw_{i:06}/DONE")).is_file() {continue;}
                let result=(|| {
                    let mut rng=SplitMix64::new(master^(i as u64).wrapping_mul(0x9e3779b97f4a7c15));
                    let deadline=Instant::now()+Duration::from_secs(seconds);
                    let mut seed=None;let mut proposals=0usize;
                    while Instant::now()<deadline && proposals<10000 {
                        proposals+=1;
                        let t=if let Some(p)=passport {passport_seed_filtered(d,ambient,p,false,true,&mut rng)}
                            else if family=="low-index" {low_index_seed_filtered(d,ambient,false,true,&mut rng)}
                            else {
                                let a=random_perm(d,ambient,&mut rng);let b=random_perm(d,ambient,&mut rng);let c=random_perm(d,ambient,&mut rng);
                                let t=[a,b,c,product_one_fourth(d,a,b,c)];
                                if !t.iter().any(|p|p.is_identity(d)) && genus_of(&t,d)>0
                                    && defect_of(&t,d).0==0 && generators_transitive(d,&t[..3])
                                    && tuple_generates_dynamic(d,ambient,&t) {Some(t)} else {None}
                            };
                        if let Some(t)=t {
                            if !t.iter().any(|p|p.is_identity(d)) && genus_of(&t,d)>0 && checked_defect(&t,d)?.0==0 {seed=Some(t);break;}
                        }
                    }
                    let folder=out.join(format!("draw_{i:06}"));
                    let Some(seed)=seed else {fs::create_dir(&folder).map_err(|e|e.to_string())?;fs::write(folder.join("DONE"),format!("no_eligible_seed proposal_batches={proposals}")).map_err(|e|e.to_string())?;return Ok(());};
                    let mut t=seed;let mut histogram=BTreeMap::new();let mut witness=None;
                    for _ in 0..steps {
                        if Instant::now()>=deadline {break;}
                        let v=checked_defect(&t,d)?;*histogram.entry(v).or_default()+=1;
                        if v.0>0 {witness=Some(t);break;}
                        // Half-lazy uniform six-move chain. Early-stopped frequencies
                        // are explicitly discovery observations, NOT stationary samples.
                        if rng.below(2)==1 {t=hurwitz_move_dynamic(d,t,rng.below(6));}
                    }
                    let case=SpectrumCase {id:format!("walk-{i}"),kind:"walk".into(),source:format!("rng={master}; draw={i}; family={family}; proposal_batches={proposals}"),ambient,d,seed,states:None,input_states:None,aliases:Vec::new()};
                    let result=if let Some(w)=witness {Scan {status:"walk_positive_defect",states:vec![canonical(&seed,d)],histogram,witness:Some(w)}}
                        else {scan(&case,cap,deadline,true)?};
                    write_case(&case,result,&folder,false,false)
                })();
                if let Err(e)=result {errors.lock().unwrap().push(format!("draw {i}: {e}"));}
            });
        }
    });
    let errors=errors.into_inner().unwrap();if !errors.is_empty() {return Err(errors.join("; "));}
    fs::write(out.join("DONE"),format!("draws={draws}\nSampling is biased; walk steps are correlated; no unseen positive defect is not a proof.\n")).map_err(|e|e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odd_defect_in_alternating_group_and_allocation_free_counters() {
        let t=parse_tuple("[(1 2 3 4 5), (1 2 6), (1 6 2), (1 5 4 3 2)]",6).unwrap();
        validate_tuple(6,Ambient::Alternating,&t).unwrap();
        assert_eq!(genus_of(&t,6),1);assert_eq!(checked_defect(&t,6).unwrap(),(1,0,0));
        let t=parse_tuple("[(1 2 3), (1 2 4), (1 2 3 4 5), (1 2 3 5 4)]",5).unwrap();
        validate_tuple(5,Ambient::Alternating,&t).unwrap();
        assert_eq!(genus_of(&t,5),2);assert_eq!(checked_defect(&t,5).unwrap(),(1,0,1));
        // Two elementary families cover every alternating degree d>=4;
        // symmetric d-cycle/transposition family has defect one in every d.
        for d in 4..=MAX_D {
            let t=if d%2==0 {
                let a=cycle_perm(d,&(0..d-1).collect::<Vec<_>>());let b=cycle_perm(d,&[0,1,d-1]);
                [a,b,b.inverse(d),a.inverse(d)]
            } else {
                let a=cycle_perm(d,&[0,1,2]);let b=cycle_perm(d,&[0,1,3]);let c=cycle_perm(d,&(0..d).collect::<Vec<_>>());
                [a,b,c,product_one_fourth(d,a,b,c)]
            };
            validate_tuple(d,Ambient::Alternating,&t).unwrap();
            assert_eq!(checked_defect(&t,d).unwrap().0,1);
            assert_eq!(genus_of(&t,d),if d%2==0 {1} else {2});
            let a=cycle_perm(d,&(0..d).collect::<Vec<_>>());let b=cycle_perm(d,&[0,1]);
            let t=[a,b,b.inverse(d),a.inverse(d)];validate_tuple(d,Ambient::Symmetric,&t).unwrap();
            assert_eq!(checked_defect(&t,d).unwrap(),(1,0,0));
        }
        let mut rng=SplitMix64::new(6161);
        for d in 4..=14 {
            for _ in 0..32 {
                let a=random_perm(d,Ambient::Alternating,&mut rng);
                let b=random_perm(d,Ambient::Alternating,&mut rng);
                let c=random_perm(d,Ambient::Alternating,&mut rng);
                let t=[a,b,c,product_one_fourth(d,a,b,c)];
                assert_eq!(point_orbit_count_perms(d,a,b),point_orbits(d,&[a,b]).len());
                let reference=cycle_type(a.compose(b,d),d).len() as i64
                    -point_orbits(d,&t[..2]).len() as i64-point_orbits(d,&t[2..]).len() as i64+1;
                assert_eq!(defect_of(&t,d).0,reference);
                if generators_transitive(d,&t[..3]) {checked_defect(&t,d).unwrap();}
            }
        }
    }

    #[test]
    fn inner_cover_matches_independent_even_conjugation_and_known_a25() {
        let d=4;let a=cycle_perm(d,&[0,1,2]);let b=cycle_perm(d,&[0,1,3]);
        let seed=[a,b,b.inverse(d),a.inverse(d)];validate_tuple(d,Ambient::Alternating,&seed).unwrap();
        let graph=enumerate_spectrum_graph(d,seed).unwrap();
        let (inner,components)=inner_graph(d,Ambient::Alternating,&graph).unwrap();
        let mut labels=vec![seed];let mut seen=HashSet::from([seed]);
        let mut head=0;
        while head<labels.len() {
            let t=labels[head];head+=1;
            for mv in 0..3 {let t=hurwitz_move_dynamic(d,t,mv);if seen.insert(t) {labels.push(t);}}
        }
        let mut even=Vec::new();let mut h=Perm::identity(d);
        loop {if h.is_even(d) {even.push(h);}if !next_permutation(&mut h.img[..d]) {break;}}
        let classes=labels.iter().map(|t|even.iter().map(|&h|std::array::from_fn::<_,3,_>(|j|t[j].conjugate_by(h,d).img)).min().unwrap()).collect::<HashSet<_>>();
        assert_eq!(inner.transitions.len(),classes.len());assert_eq!(graph_components(&inner),1);
        assert_eq!(components*inner.transitions.len(),2*graph.states.len());
        let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_A25_order12_genus12_weak.txt")).unwrap();
        let graph=graph_from_states(25,case.states.unwrap()).unwrap();
        let (inner,components)=inner_graph(25,Ambient::Alternating,&graph).unwrap();
        assert_eq!(components,1);assert_eq!(inner.transitions,vec![[1;6],[0;6]]);
        for name in ["CERTIFIED_A25_weak.txt","CERTIFIED_A25_order3_weak.txt","CERTIFIED_A25_order6_genus12_weak.txt"] {
            let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();
            let graph=graph_from_states(25,case.states.unwrap()).unwrap();
            let (inner,components)=inner_graph(25,Ambient::Alternating,&graph).unwrap();
            assert_eq!(components,1);assert_eq!(inner.transitions.len(),8);
            assert!(graph_is_bipartite(&inner));
        }
        let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_A25_mixed_genus0_weak.txt")).unwrap();
        let graph=graph_from_states(25,case.states.unwrap()).unwrap();
        let (inner,components)=inner_graph(25,Ambient::Alternating,&graph).unwrap();
        assert_eq!(components,1);assert_eq!(inner.transitions.len(),24);
        assert!(graph_is_bipartite(&inner));
        // Exact move-labeled isomorphism to S4 with adjacent transpositions.
        // This checks every edge, rather than inferring a graph from eigenvalues.
        let mut labels=vec![[0u8,1,2,3]];let mut positions=HashMap::from([(labels[0],0usize)]);let mut vertices=vec![0usize];
        let mut head=0;
        while head<labels.len() {
            for mv in 0..3 {
                let mut next=labels[head];next.swap(mv,mv+1);
                let target=inner.transitions[vertices[head]][mv] as usize;
                assert_eq!(target,inner.transitions[vertices[head]][mv+3] as usize);
                if let Some(&i)=positions.get(&next) {assert_eq!(vertices[i],target);} else {
                    positions.insert(next,labels.len());labels.push(next);vertices.push(target);
                }
            }
            head+=1;
        }
        assert_eq!(labels.len(),24);assert_eq!(vertices.iter().copied().collect::<HashSet<_>>().len(),24);
    }

    #[test]
    fn bounded_scans_and_certificate_input_guards() {
        let seed=parse_tuple("[(1 2 3),(1 2 4),(1 2 3 4 5),(1 2 3 5 4)]",5).unwrap();
        let case=SpectrumCase {id:"check".into(),kind:"tuple".into(),source:"test".into(),
            ambient:Ambient::Alternating,d:5,seed,states:None,input_states:None,aliases:Vec::new()};
        assert_eq!(scan(&case,100,Instant::now()+Duration::from_secs(5),true).unwrap().status,"positive_defect_witness");
        assert_eq!(scan(&case,1,Instant::now()+Duration::from_secs(5),false).unwrap().status,"capped");
        assert_eq!(scan(&case,100,Instant::now(),false).unwrap().status,"timeout");
        let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
        let folder=std::env::temp_dir().join(format!("hurwitz-zero-check-{}-{stamp}",std::process::id()));
        fs::create_dir(&folder).unwrap();
        let path=folder.join("ZERO.txt");
        fs::write(&path,format!("HURWITZ_ZERO_DEFECT_V1\ngroup=A\nd=5\nseed={}\nstates=1\nSTATE\t0\t0\t0\n",tuple_string_dynamic(5,&seed))).unwrap();
        assert!(verify_zero(&path).unwrap_err().contains("invalid encoded permutation"));
        let key=canonical(&seed,5);
        fs::write(&path,format!("HURWITZ_ZERO_DEFECT_V1\ngroup=A\nd=5\nseed={}\nstates=1\nSTATE\t{:032x}\t{:032x}\t{:032x}\n",tuple_string_dynamic(5,&seed),key[0],key[1],key[2])).unwrap();
        assert!(verify_zero(&path).unwrap_err().contains("positive-defect state"));
        let options=vec!["binary".into(),"defect-mission".into(),"--cap".into(),"10".into()];
        check_options(&options,&folder).unwrap();check_options(&options,&folder).unwrap();
        let mut changed=options;changed[3]="11".into();
        assert!(check_options(&changed,&folder).is_err());
        fs::remove_file(path).unwrap();fs::remove_file(folder.join("OPTIONS.txt")).unwrap();
        fs::remove_dir(folder).unwrap();
    }
}
