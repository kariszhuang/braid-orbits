// Independent optimized adapter. Historical census/search implementations stay
// available as literal oracles. Keys and six-move certificate formats unchanged.
#[inline]
fn ff_product_index(a: Perm, b: Perm, n: usize) -> usize {
    let mut seen = 0u64;
    let mut cycles = 0;
    for root in 0..n {
        if seen & (1u64 << root) != 0 { continue; }
        cycles += 1;
        let mut x = root;
        while seen & (1u64 << x) == 0 {
            seen |= 1u64 << x;
            x = a.img[b.img[x] as usize] as usize;
        }
    }
    n - cycles
}

// A positive-genus passport has at most one eligible additive side. Stop as
// soon as the cycle count cannot equal its target; no heuristic rejection.
#[inline]
fn ff_product_matches(a: Perm, b: Perm, n: usize, left: usize, right: usize) -> bool {
    let target = match (left < n, right < n) {
        (false, false) => return false,
        (true, true) if left != right => {
            let i = ff_product_index(a, b, n);
            return i == left || i == right;
        }
        (true, _) => left,
        _ => right,
    };
    let expected = n - target;
    let mut remaining = (1u64 << n) - 1;
    let mut count = 0;
    while remaining != 0 {
        count += 1;
        if count > expected { return false; }
        let mut x = remaining.trailing_zeros() as usize;
        while remaining & (1u64 << x) != 0 {
            remaining &= !(1u64 << x);
            x = a.img[b.img[x] as usize] as usize;
        }
        if (remaining.count_ones() as usize) < expected - count { return false; }
    }
    count == expected
}

// Feasibility depends on three small integers, not on the labelled pair.
// Cache it once per (degree, first index, second index), process-wide. Entries
// are immutable, lazily built and shared by every frame; no candidate is saved
// or excluded merely because it appears in this cache.
struct ForestIndexPolicy {
    allowed_v: [u64; MAX_D],
    third_indices: Vec<usize>,
}
static FF_INDEX_POLICIES: [std::sync::OnceLock<Vec<std::sync::OnceLock<ForestIndexPolicy>>>; MAX_D + 1] =
    [const { std::sync::OnceLock::new() }; MAX_D + 1];

fn ff_index_policy(n: usize, ia: usize, ib: usize) -> &'static ForestIndexPolicy {
    assert!(n <= MAX_D && ia < n && ib < n);
    let policies = FF_INDEX_POLICIES[n].get_or_init(||
        (0..n*n).map(|_| std::sync::OnceLock::new()).collect());
    policies[ia*n+ib].get_or_init(|| {
        let mut allowed_v = [0u64; MAX_D];
        let mut third_indices = Vec::new();
        let e = ia+ib;
        for t in 1..n {
            for v in 1..n {
                if e+t+v >= 2*n && (e+t+v)%2 == 0
                    && (ia+t < n || ib+v < n)
                    && (ia+v < n || ib+t < n) {
                    allowed_v[t] |= 1u64<<v;
                }
            }
            if allowed_v[t] != 0 { third_indices.push(t); }
        }
        ForestIndexPolicy { allowed_v, third_indices }
    })
}

struct ForestFrame {
    pair: [Perm; 2],
    p: Perm,
    factors: [Perm; 5],
    indices: [usize; 2],
    policy: &'static ForestIndexPolicy,
    b_inverse: Perm,
    twists: std::cell::OnceCell<Vec<TwistFactor>>,
    point_orbits: std::cell::OnceCell<[u64; MAX_D]>,
}
struct TwistFactor { factor: Perm, first_index: usize, word: Vec<usize> }
#[derive(Clone,Copy,Debug)]
enum DeepFailure { Middle(usize,usize), Twist(usize), Word(usize) }
const FF_DEEP_WORDS: [&[usize];6] = [&[1,3,1],&[0,4,4],&[4,0,4],&[3,1,1],
                                    &[5,1,1,1],&[5,3,1,1,1]];
impl ForestFrame {
    fn new(pair: [Perm; 2], n: usize) -> Self {
        let [a,b] = pair;
        let p = a.compose(b,n);
        let indices = [sm_index(a,n),sm_index(b,n)];
        let factors = [b.inverse(n).compose(a,n).compose(b,n), a, b,
                       a.compose(b,n).compose(a.inverse(n),n),
                       p.inverse(n).compose(b,n).compose(p,n)];
        let [ia,ib] = indices;
        Self {pair,p,factors,indices,policy:ff_index_policy(n,ia,ib),b_inverse:b.inverse(n),
              twists:std::cell::OnceCell::new(),point_orbits:std::cell::OnceCell::new()}
    }
    // Equivalent to the old ten words, including the first failing witness.
    // d=(abc)^-1; cycle count is invariant under conjugation and inversion.
    fn short_word(&self, c: Perm, ic: usize, id: usize, n: usize) -> Option<usize> {
        let [ia,ib] = self.indices;
        for (j,&f) in self.factors.iter().enumerate() {
            let (left,right) = if j<2 {(ia+ic,ib+id)} else {(ib+ic,ia+id)};
            if left >= n && right >= n { return Some(j); }
            if !ff_product_matches(f,c,n,left,right) { return Some(j); }
        }
        None
    }
    fn tuple(&self,c:Perm,n:usize)->[Perm;4] {
        [self.pair[0],self.pair[1],c,self.p.compose(c,n).inverse(n)]
    }
    // Geodesy makes p-cycles exactly the <a,b> orbits. Cache their masks only
    // when a frame actually reaches this stage; no cost on discarded frames.
    fn transitive(&self,c:Perm,n:usize)->bool {
        let masks=self.point_orbits.get_or_init(||{
            let mut masks=[0u64;MAX_D];let mut remaining=(1u64<<n)-1;
            while remaining!=0 {
                let root=remaining.trailing_zeros() as usize;let mut x=root;let mut cycle=0u64;
                loop {cycle|=1u64<<x;x=self.p.img[x] as usize;if x==root{break;}}
                remaining&=!cycle;let mut points=cycle;
                while points!=0 {let x=points.trailing_zeros() as usize;points&=points-1;masks[x]=cycle;}
            }
            masks
        });
        let full=(1u64<<n)-1;let mut seen=masks[0];let mut todo=seen;
        while seen!=full && todo!=0 {
            let x=todo.trailing_zeros() as usize;todo&=todo-1;
            let next=masks[c.img[x] as usize]&!seen;seen|=next;todo|=next;
        }
        seen==full
    }
    fn twist_factors(&self,n:usize)->&[TwistFactor] {
        self.twists.get_or_init(||{
            let mut result:Vec<TwistFactor>=vec![];
            for end in [0,3] {
                let mut t=[self.pair[0],self.pair[1],Perm::identity(n),Perm::identity(n)];
                for power in 1..=4 {
                    t=hurwitz_move_dynamic(n,t,end);
                    if power<2{continue;}
                    let first_index=self.indices[power%2];
                    for middle in [4,1] {
                        let factor=if middle==4 {t[0]} else {t[1].inverse(n).compose(t[0],n).compose(t[1],n)};
                        if self.factors.iter().enumerate().any(|(j,&f)|f==factor&&first_index==self.indices[usize::from(j>=2)])
                            || result.iter().any(|f|f.factor==factor&&f.first_index==first_index){continue;}
                        let mut word=vec![end;power];word.push(middle);
                        result.push(TwistFactor{factor,first_index,word});
                    }
                }
            }
            result
        })
    }
    // Negative-only sieve: the seed has passed the five original tests. No
    // canonicalization, group-order computation, hashing or allocation here.
    fn deep_failure(&self,c:Perm,ic:usize,id:usize,n:usize,level:usize)->Option<DeepFailure> {
        if level==0{return None;}
        let [a,b]=self.pair;let ci=c.inverse(n);let bi=self.b_inverse;
        let mut forward=Perm::identity(n);let mut backward=forward;
        for x in 0..n {
            forward.img[x]=b.img[c.img[b.img[ci.img[bi.img[x] as usize] as usize] as usize] as usize];
            backward.img[x]=ci.img[b.img[c.img[x] as usize] as usize];
        }
        let [ia,ib]=self.indices;
        for (j,f) in [forward,backward].into_iter().enumerate() {
            if !ff_product_matches(a,f,n,ia+ib,ic+id){return Some(DeepFailure::Middle(if j==0{1}else{4},2));}
        }
        if level==1{return None;}
        // Odd middle powers interchange only b,c; indices are already known.
        let t=self.tuple(c,n);let mut ix=[ia,ib,ic,id];ix.swap(1,2);
        for middle in [1,4] {
            let mut u=t;for _ in 0..3{u=hurwitz_move_dynamic(n,u,middle);}
            if !ff_product_matches(u[0],u[1],n,ix[0]+ix[1],ix[2]+ix[3]){return Some(DeepFailure::Middle(middle,3));}
        }
        for (j,f) in self.twist_factors(n).iter().enumerate() {
            let other=ia+ib-f.first_index;
            if !ff_product_matches(f.factor,c,n,f.first_index+ic,other+id){return Some(DeepFailure::Twist(j));}
        }
        if level==2{return None;}
        for (j,word) in FF_DEEP_WORDS.iter().enumerate() {
            let mut u=t;let mut ix=[ia,ib,ic,id];
            for &m in *word{u=hurwitz_move_dynamic(n,u,m);let k=m%3;ix.swap(k,k+1);}
            if !ff_product_matches(u[0],u[1],n,ix[0]+ix[1],ix[2]+ix[3]){return Some(DeepFailure::Word(j));}
        }
        None
    }
    fn deep_word(&self,failure:DeepFailure,n:usize)->Vec<usize> {
        match failure {
            DeepFailure::Middle(m,k)=>vec![m;k],
            DeepFailure::Twist(j)=>self.twist_factors(n)[j].word.clone(),
            DeepFailure::Word(j)=>FF_DEEP_WORDS[j].to_vec(),
        }
    }
}
const FF_WORDS: [&[usize];5] = [&[1],&[4],&[0,1],&[0,4],&[3,1]];

// Same RNG sequence and proposal law as gd_index_perm; only the allocation
// changes. The rest of MAX_D is initialized, so no uninitialized memory/unsafe.
fn ff_index_perm(n:usize,index:usize,rng:&mut SplitMix64)->Perm {
    let mut points:[usize;MAX_D] = std::array::from_fn(|i|i);
    rng.shuffle(&mut points[..n]);
    let mut p=Perm::identity(n);let mut used=0;let mut remaining=index;
    while remaining>0 {
        let slack=n-used-remaining;
        let contribution=if slack<=1 {remaining} else {1+rng.below(remaining)};
        let len=contribution+1;
        for j in 0..len {p.img[points[used+j]]=points[used+(j+1)%len] as u8;}
        used+=len;remaining-=contribution;
    }
    p
}

// Compatible canonicalization: compare completed leading words before packing
// later ones. All roots still participate; no heuristic root exclusion.
fn ff_key(t:&[Perm;4],n:usize)->SmoothKey {
    assert!(n<=SM_MAX_DEGREE);
    let mut best=[u128::MAX;SM_KEY_WORDS];
    for root in 0..n {
        let mut labels=[u8::MAX;MAX_D];let mut order=[0usize;MAX_D];
        order[0]=root;labels[root]=0;let mut len=1;let mut head=0;
        while head<len {
            let x=order[head];head+=1;
            for p in &t[..3] {
                let y=p.img[x] as usize;
                if labels[y]==u8::MAX {labels[y]=len as u8;order[len]=y;len+=1;}
            }
        }
        assert_eq!(len,n);
        let mut key=[0u128;SM_KEY_WORDS];let mut ordering=std::cmp::Ordering::Equal;
        for j in 0..3 {
            let start=j*SM_KEY_LIMBS;let end=start+SM_KEY_LIMBS;
            key[start..end].copy_from_slice(&sm_packed_word(t[j],&labels,&order,n));
            if ordering==std::cmp::Ordering::Equal{ordering=key[start..end].cmp(&best[start..end]);}
            if ordering==std::cmp::Ordering::Greater {break;}
        }
        if ordering!=std::cmp::Ordering::Greater {best=key;}
    }
    best
}
fn ff_path(parent:&[(usize,usize)],mut i:usize)->Vec<usize> {
    let mut word=vec![];
    while i>0 {let(p,m)=parent[i];word.push(m);i=p;}
    word.reverse();word
}
fn ff_smooth(t:&[Perm;4],ix:[usize;4],n:usize)->bool {
    let ip=ff_product_index(t[0],t[1],n);
    ip==ix[0]+ix[1] || ip==ix[2]+ix[3]
}
fn ff_scan(seed:&[Perm;4],n:usize,cap:usize,deadline:Instant)->SmoothScan {
    ff_scan_cached(seed,n,cap,deadline,&mut HashMap::new(),0)
}
fn ff_replay_negative(seed:[Perm;4],n:usize,word:&[usize])->bool {
    let mut t=seed;for &m in word{t=hurwitz_move_dynamic(n,t,m);}
    let ix=t.map(|p|sm_index(p,n));!ff_smooth(&t,ix,n)
}
fn ff_cache_negative(states:&[SmoothKey],parent:&[(usize,usize)],word:&[usize],n:usize,
                     cache:&mut HashMap<SmoothKey,Vec<usize>>,limit:usize) {
    // A state's inverse root path followed by the root witness is a proof for
    // that state. Store only bounded short words, and literal-replay each one.
    for (i,&key) in states.iter().take(256).enumerate() {
        if cache.len()>=limit{break;}
        if cache.contains_key(&key){continue;}
        let path=ff_path(parent,i);
        let mut proof:Vec<_>=path.into_iter().rev().map(|m|(m+3)%6).collect();
        proof.extend_from_slice(word);
        if proof.len()>64{continue;}
        assert!(ff_replay_negative(sm_unpack(key,n),n,&proof));
        cache.insert(key,proof);
    }
}
fn ff_scan_cached(seed:&[Perm;4],n:usize,cap:usize,deadline:Instant,
                  cache:&mut HashMap<SmoothKey,Vec<usize>>,cache_limit:usize)->SmoothScan {
    assert!(cap>0);
    let ix=seed.map(|p|sm_index(p,n));
    if !ff_smooth(seed,ix,n) {
        return SmoothScan{status:"nonsmooth",states:vec![],word:vec![]};
    }
    let root=ff_key(seed,n);let mut states=vec![root];
    let mut seen=HashMap::from([(root,0usize)]);
    let mut indices=vec![ix];let mut parent=vec![(0usize,0usize)];let mut head=0;
    while head<states.len() {
        if Instant::now()>=deadline {return SmoothScan{status:"timeout",states,word:vec![]};}
        let t=sm_unpack(states[head],n);
        if let Some(proof)=cache.get(&states[head]) {
            assert!(ff_replay_negative(t,n,proof));
            let mut word=ff_path(&parent,head);word.extend_from_slice(proof);
            ff_cache_negative(&states,&parent,&word,n,cache,cache_limit);
            return SmoothScan{status:"nonsmooth",states,word};
        }
        // End moves preserve smoothness. Check both middle successors before
        // doing ANY of the six costly canonicalizations, but keep all generators.
        let middle=[hurwitz_move_dynamic(n,t,1),hurwitz_move_dynamic(n,t,4)];
        let mut mid_ix=indices[head];mid_ix.swap(1,2);
        for (j,u) in middle.iter().enumerate() {
            if !ff_smooth(u,mid_ix,n) {
                let mut word=ff_path(&parent,head);word.push(if j==0 {1} else {4});
                ff_cache_negative(&states,&parent,&word,n,cache,cache_limit);
                return SmoothScan{status:"nonsmooth",states,word};
            }
        }
        for mv in 0..6 {
            let u=if mv==1 {middle[0]} else if mv==4 {middle[1]} else {hurwitz_move_dynamic(n,t,mv)};
            let k=ff_key(&u,n);
            if let std::collections::hash_map::Entry::Vacant(entry)=seen.entry(k) {
                if states.len()==cap {return SmoothScan{status:"capped",states,word:vec![]};}
                let mut next=indices[head];let j=mv%3;next.swap(j,j+1);
                entry.insert(states.len());parent.push((head,mv));indices.push(next);states.push(k);
            }
        }
        head+=1;
    }
    SmoothScan{status:"closed_smooth",states,word:vec![]}
}

pub fn forest_fast_main(args:&[String])->Result<(),String> {
    let arg=|s:&str,d:&str|args.windows(2).find(|w|w[0]==s).map(|w|w[1].clone()).unwrap_or(d.into());
    if args.iter().any(|s|s=="--io-bench") {return ff_io_bench(Path::new(&arg("--out","")));}
    if args.iter().any(|s|s=="--check-controls") {
        let input=fs::read_to_string(arg("--check-controls","")).map_err(|e|e.to_string())?;
        let mut export=if args.iter().any(|s|s=="--export-known") {
            Some(BufWriter::new(OpenOptions::new().create_new(true).write(true)
                .open(arg("--export-known","")).map_err(|e|e.to_string())?))
        }else{None};
        let mut count=0;let mut states=0;
        for line in input.lines(){
            let parts:Vec<_>=line.splitn(3,'\t').collect();if parts.len()!=3{return Err("control degree,states,tuple".into());}
            let n:usize=parts[0].parse().map_err(|_|"control degree")?;let expected:usize=parts[1].parse().map_err(|_|"control states")?;
            let t=sm_parse(parts[2],n)?;let f=ForestFrame::new([t[0],t[1]],n);
            let ic=sm_index(t[2],n);let id=sm_index(t[3],n);
            assert_eq!(f.short_word(t[2],ic,id,n),None);assert!(f.policy.allowed_v[ic]&(1<<id)!=0);
            assert!(f.deep_failure(t[2],ic,id,n,3).is_none());
            let scan=ff_scan(&t,n,1_000_000,Instant::now()+Duration::from_secs(120));
            assert_eq!(scan.status,"closed_smooth");assert_eq!(scan.states.len(),expected);
            for &key in &scan.states{
                let u=sm_unpack(key,n);assert_eq!(ff_key(&u,n),sm_key(&u,n));assert_eq!(sm_values(&u,n).1,0);
                if let Some(file)=export.as_mut(){writeln!(file,"{n}\t{}",tuple_string_dynamic(n,&u)).map_err(|e|e.to_string())?;}
                let frame=ForestFrame::new([u[0],u[1]],n);
                assert!(frame.deep_failure(u[2],sm_index(u[2],n),sm_index(u[3],n),n,3).is_none());
            }
            count+=1;states+=expected;
        }
        println!("controls={count} states={states} all_closed=true all_compatible_keys=true");return Ok(());
    }
    if args.iter().any(|s|s=="--gpu-middle-bench") {return ff_middle_gpu_bench(args);}
    if args.iter().any(|s|s=="--bench") {return ff_bench(args);}
    if args.iter().any(|s|s=="--screen-seeds") {
        let n:usize=arg("--d","11").parse().map_err(|_|"degree")?;
        if !(2..=25).contains(&n){return Err("screen degree2..25".into());}
        let input=fs::read_to_string(arg("--screen-seeds","")).map_err(|e|e.to_string())?;
        for line in input.lines() {
            let t=sm_parse(line,n)?;let f=ForestFrame::new([t[0],t[1]],n);
            if let Some(fail)=f.deep_failure(t[2],sm_index(t[2],n),sm_index(t[3],n),n,3) {
                let word=f.deep_word(fail,n);assert!(ff_replay_negative(t,n,&word));
                println!("{}",word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(","));
            }else{println!("survives");}
        }
        return Ok(());
    }
    let n:usize=arg("--d","11").parse().map_err(|_|"degree")?;
    let seconds:u64=arg("--seconds","60").parse().map_err(|_|"seconds")?;
    let seed:u64=arg("--seed","1").parse().map_err(|_|"seed")?;
    let reuse:usize=arg("--frame-reuse","256").parse().map_err(|_|"reuse")?;
    let attempts:u64=arg("--attempts","0").parse().map_err(|_|"attempts")?;
    let closure_seconds:u64=arg("--closure-seconds","10").parse().map_err(|_|"closure seconds")?;
    let closure_cap:usize=arg("--cap","100000").parse().map_err(|_|"cap")?;
    if closure_seconds==0||closure_cap==0{return Err("positive closure limits".into());}
    let mode=arg("--lane","forest");
    let legacy=args.iter().any(|s|s=="--legacy-proposals");
    let reference=args.iter().any(|s|s=="--reference");
    let deep:usize=arg("--deep-level","3").parse().map_err(|_|"deep-level")?;
    let orbit_cache=!args.iter().any(|s|s=="--no-orbit-cache");
    if deep>3{return Err("deep-level0..3".into());}
    let max_degree=if mode=="triangle"{SM_MAX_DEGREE}else{25};
    if n<9||n>max_degree||seconds==0||reuse==0||!["forest","blocks","joint","solved","quotient","primitive","mutate","double_lift","triangle","returns"].contains(&mode.as_str()) {return Err(format!("degree9..{max_degree}, positive seconds/reuse, valid lane"));}
    if ["primitive","mutate","double_lift","triangle","returns"].contains(&mode.as_str())&&reuse!=1{return Err("mechanism lanes require frame-reuse1".into());}
    if mode=="triangle" && n<11{return Err("completed9/10 triangle populations are retired".into());}
    if mode=="returns" && !(2..=n/3).any(|r|n%r==0){return Err("return forests need a composite block shape".into());}
    if mode=="solved" && (reuse!=1 || dd_solved_shapes(n).is_empty()) {return Err("solved lane requires frame-reuse1 and an unretired shape".into());}
    if mode=="quotient" && (reuse!=1 || qf_shapes(n).is_empty()) {return Err("quotient lane requires frame-reuse1 and a composite shape".into());}
    let shapes:Vec<_>=(3..=n/2).filter(|r|n%r==0).map(|r|(r,n/r)).collect();
    if mode=="blocks"&&shapes.is_empty(){return Err("no block shape".into());}
    let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh --out required".into());}
    fs::create_dir(&out).map_err(|e|e.to_string())?;
    if let Some(path)=args.windows(2).find(|w|w[0]=="--replay").map(|w|&w[1]){
        let raw=fs::read_to_string(path).map_err(|e|e.to_string())?;let t=sm_parse(raw.trim(),n)?;
        if t.iter().any(|p|p.is_identity(n))||!sm_values(&t,n).0||sm_values(&t,n).2<=0{return Err("replay scope".into());}
        let start=Instant::now();let text=tuple_string_dynamic(n,&t);
        let mut lead=File::create(out.join("LEADS.tsv")).map_err(|e|e.to_string())?;
        writeln!(lead,"1\t{text}").map_err(|e|e.to_string())?;lead.flush().map_err(|e|e.to_string())?;
        let scan=ff_scan(&t,n,closure_cap,start+Duration::from_secs(seconds));
        let folder=if scan.status=="nonsmooth"{""}else{"replay"};
        if !folder.is_empty(){sm_save(&t,n,&scan,&out.join(folder))?;}
        fs::write(out.join("RESOLUTIONS.tsv"),format!("candidate\tstatus\tword\tfolder\tseed\n1\t{}\t{}\t{folder}\t{text}\n",scan.status,scan.word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(","))).map_err(|e|e.to_string())?;
        fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\treplay\t{seed}\t1\t1\t1\t1\t0\t0\t{}\t{}\t{}\n",usize::from(scan.status=="closed_smooth"),usize::from(["timeout","capped"].contains(&scan.status)),start.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
        return Ok(());
    }
    let primitive_name=arg("--primitive-group","psl11");
    let primitive_cursor=arg("--primitive-cursor","0").parse::<usize>().map_err(|_|"primitive cursor")?;
    let mut primitive=if mode=="primitive"{Some(PrimitiveFrames::new(&primitive_name,n,primitive_cursor)?)}else{None};
    let mut mutation=if mode=="mutate"{Some(MutationSeeds::read(Path::new(&arg("--mutation-seeds","")),n)?)}else{None};
    let mut double=if mode=="double_lift"{Some(RamifiedDoubleFrames::new(Path::new(&arg("--double-bases","")),n,primitive_cursor)?)}else{None};
    let target_genus=arg("--target-genus","0").parse::<usize>().map_err(|_|"target genus")?;
    let triangle_min_genus=arg("--triangle-min-genus","0").parse::<usize>().map_err(|_|"triangle minimum genus")?;
    let triangle_low_min=arg("--triangle-low-min","2").parse::<usize>().map_err(|_|"triangle pair-index minimum")?;
    let triangle_low_max=arg("--triangle-low-max","63").parse::<usize>().map_err(|_|"triangle pair-index maximum")?;
    let triangle_passport_min=arg("--triangle-passport-min","6").parse::<usize>().map_err(|_|"triangle passport minimum")?;
    let triangle_passport_max=arg("--triangle-passport-max","31").parse::<usize>().map_err(|_|"triangle passport maximum")?;
    let triangle_strength=arg("--triangle-strength-divisor","8").parse::<usize>().map_err(|_|"triangle strength divisor")?;
    let triangle_method=arg("--triangle-third-method","mixed");
    let triangle_exact_raw=arg("--triangle-exact-passport","");
    let triangle_exact_profile:Option<[usize;4]>=if triangle_exact_raw.is_empty(){None}else{
        let values=triangle_exact_raw.split(',').map(|x|x.parse::<usize>().map_err(|_|"exact triangle passport"))
            .collect::<Result<Vec<_>,_>>()?;
        let profile:[usize;4]=values.try_into().map_err(|_|"exact triangle passport needs four indices")?;
        if profile.windows(2).any(|w|w[0]>w[1]){return Err("exact triangle passport must be sorted".into());}
        Some(profile)
    };
    if mode=="triangle"&&(triangle_strength==0||triangle_low_min<2||triangle_low_max<triangle_low_min
        ||triangle_passport_min<6||triangle_passport_max<triangle_passport_min
        ||!["random","joint","strong","mixed","centralizer","pairwise_low","middle_repair"].contains(&triangle_method.as_str())){return Err("invalid triangle passport shard/method".into());}
    if let Some(profile)=triangle_exact_profile {
        let sum:usize=profile.iter().sum();
        if mode!="triangle"||profile[0]<triangle_passport_min||profile[0]>triangle_passport_max
            ||profile[3]>=n||profile[1]+profile[2]>n-1||profile[0]+profile[3]<n
            ||sum<2*n||sum%2!=0||(target_genus>0&&1+sum/2!=n+target_genus) {
            return Err("exact triangle passport conflicts with job scope".into());
        }
    }
    let triangle_profiles=if mode=="triangle" {
        if let Some(profile)=triangle_exact_profile {vec![profile]} else if triangle_method=="centralizer" {
        sn_triangle_passports(n,triangle_passport_min,triangle_passport_max,triangle_min_genus.max(3),target_genus)
        }else{Vec::new()}
    }else{Vec::new()};
    let mut proposal:Option<MechanismProposal>=None;
    let mut mechanism_log=BufWriter::new(File::create(out.join("MECHANISM.tsv")).map_err(|e|e.to_string())?);
    writeln!(mechanism_log,"candidate\tsource\tarm\tpoints\tprimitive").map_err(|e|e.to_string())?;
    let mut leads=BufWriter::new(File::create(out.join("LEADS.tsv")).map_err(|e|e.to_string())?);
    let mut resolutions=BufWriter::new(File::create(out.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?);
    writeln!(resolutions,"candidate\tstatus\tword\tfolder\tseed").map_err(|e|e.to_string())?;
    let mut rng=SplitMix64::new(seed);let start=Instant::now();let deadline=start+Duration::from_secs(seconds);
    let mut count=[0u64;8];let mut seen=HashSet::new();let mut closed=HashMap::new();
    let(mut triangle_budget_reject,mut triangle_shard_reject,mut triangle_profile_reject,mut triangle_genus_reject,
        mut triangle_high_index_bound_reject,mut triangle_fifteen_bulk_reject,mut triangle_fixed_free_reject,
        mut triangle_low_involution_reject,mut triangle_two_special_reject)=(0u64,0u64,0u64,0u64,0u64,0u64,0u64,0u64,0u64);
    let mut triangle_short_rejects=[0u64;5];
    let mut triangle_short_examples:[Option<String>;5]=std::array::from_fn(|_|None);
    let mut triangle_deep_rejects=0u64;
    let mut triangle_deep_examples:Vec<(String,String,String)>=Vec::with_capacity(32);
    let mut triangle16_repair_stats=Triangle16RepairStats::default();
    let mut negative_cache=if args.iter().any(|s|s=="--negative-bank") {
        dd_read_negative(Path::new(&arg("--negative-bank","")),n)?
    }else{HashMap::new()};
    let mut known=HashSet::new();
    if args.iter().any(|s|s=="--known-bank") {
        for line in fs::read_to_string(arg("--known-bank","")).map_err(|e|e.to_string())?.lines() {
            let (degree,text)=line.split_once('\t').ok_or("known bank degree,tuple")?;
            if degree.parse::<usize>().map_err(|_|"known degree")?!=n{continue;}
            let t=sm_parse(text,n)?;
            if !tuple_product_is_identity(n,&t)||!generators_transitive(n,&t)||!sm_values(&t,n).0{return Err("invalid known state".into());}
            known.insert(ff_key(&t,n));
        }
    }
    let mut inputs=BufWriter::new(File::create(out.join("INPUT_BANKS.tsv")).map_err(|e|e.to_string())?);
    for flag in ["--known-bank","--negative-bank"] {
        if args.iter().any(|s|s==flag){writeln!(inputs,"{flag}\t{}",arg(flag,"")).map_err(|e|e.to_string())?;}
    }
    inputs.flush().map_err(|e|e.to_string())?;
    let mut aliases=BufWriter::new(File::create(out.join("KNOWN.tsv")).map_err(|e|e.to_string())?);
    let mut known_hits=0u64;let mut joint_failed=0u64;let mut solved=None;
    let mut frame:Option<ForestFrame>=None;let mut shape=(0,0);let mut uses=0;
    let mut frames=0u64;let mut feasible_draws=0u64;let mut deep_reject=0u64;let mut next_flush=start+Duration::from_secs(1);
    let mut frame_log=BufWriter::new(File::create(out.join("FRAMES.tsv")).map_err(|e|e.to_string())?);
    let compact=args.iter().any(|s|s=="--compact-frames");let mut last_frame_rng=seed;
    writeln!(frame_log,"frame\tfirst_candidate\ta\tb\tthird_indices{}",if compact{"\trng_before"}else{""}).map_err(|e|e.to_string())?;
    let cache_primitive=mode=="primitive"&&!args.iter().any(|s|s=="--uncached-primitive-frames");
    while Instant::now()<deadline && (attempts==0||count[0]<attempts) {
        if cache_primitive {
            proposal=primitive.as_mut().unwrap().next(n);
            let Some(p)=proposal.as_ref()else{break;};let pair=[p.tuple[0],p.tuple[1]];
            if frame.as_ref().is_none_or(|f|f.pair!=pair) {
                frame=Some(ForestFrame::new(pair,n));frames+=1;
                let f=frame.as_ref().unwrap();
                writeln!(frame_log,"{frames}\t{}\t{}\t{}\t{}{}",count[0]+1,pair[0].cycle_string(n),pair[1].cycle_string(n),
                    f.policy.third_indices.iter().map(|i|i.to_string()).collect::<Vec<_>>().join(","),if compact{format!("\t{}",rng.state)}else{String::new()}).map_err(|e|e.to_string())?;
            }
            uses=0;
        }
        // Legacy mode preserves the old random stream exactly. New mode rotates
        // independently generated frames; no invariant blocks in the forest lane.
        if !cache_primitive && (frame.is_none()||legacy||uses>=reuse) {
            last_frame_rng=rng.state;
            let pair=if ["primitive","mutate","double_lift","triangle","returns"].contains(&mode.as_str()) {
                    proposal=if mode=="triangle"{sn_triangle_targeted_tuple_with_method(n,&mut rng,triangle_low_min,triangle_low_max,triangle_strength,&triangle_method,&triangle_profiles,triangle_exact_profile,&mut triangle16_repair_stats)}else if mode=="returns"{sn_return_tuple(n,&mut rng)}
                        else if let Some(p)=primitive.as_mut(){p.next(n)}else if let Some(p)=double.as_mut(){p.next()}else{Some(mutation.as_mut().unwrap().next(n,&mut rng))};
                    if let Some(p)=proposal.as_ref(){[p.tuple[0],p.tuple[1]]}
                    else if ["triangle","returns"].contains(&mode.as_str()){count[0]+=1;continue;}else{break;}
                }else if mode=="forest"||mode=="joint" {gd_forest(n,&mut rng)}
                else if mode=="solved" || mode=="quotient" {
                    solved=if mode=="quotient" {qf_tuple(n,&mut rng)}else{dd_solved_tuple(n,&mut rng)};
                    if let Some((tuple,r,k))=solved {shape=(r,k);[tuple[0],tuple[1]]}
                    else{count[0]+=1;continue;}
                }else {shape=shapes[rng.below(shapes.len())];gd_block_pair(shape.0,shape.1,&mut rng)};
            frame=Some(ForestFrame::new(pair,n));uses=0;frames+=1;
            if !legacy && (!compact || frames%4096==1) {let f=frame.as_ref().unwrap();
                writeln!(frame_log,"{frames}\t{}\t{}\t{}\t{}{}",count[0]+1,f.pair[0].cycle_string(n),f.pair[1].cycle_string(n),f.policy.third_indices.iter().map(|i|i.to_string()).collect::<Vec<_>>().join(","),if compact{format!("\t{last_frame_rng}")}else{String::new()}).map_err(|e|e.to_string())?;}
        }
        uses+=1;count[0]+=1;let f=frame.as_ref().unwrap();
        let (c,ic)=if mode=="forest"||mode=="joint" {
            let index=if legacy {1+rng.below(n-1)} else {
                if f.policy.third_indices.is_empty(){uses=reuse;continue;}
                feasible_draws+=1;f.policy.third_indices[rng.below(f.policy.third_indices.len())]
            };
            // The constructor realizes the chosen index exactly; differential
            // tests cover every index. Do not traverse these cycles twice.
            if mode=="joint" {
                if let Some(c)=dd_joint_third(f,n,index,&mut rng){(c,index)}
                else{joint_failed+=1;continue;}
            }else{(ff_index_perm(n,index,&mut rng),index)}
        } else if ["primitive","mutate","double_lift","triangle","returns"].contains(&mode.as_str()) {let c=proposal.as_ref().unwrap().tuple[2];(c,sm_index(c,n))
        } else if mode=="solved" || mode=="quotient" {let tuple=solved.unwrap().0;
            if tuple[0].is_identity(n){continue;}
            assert!(qf_hypothesis(&tuple,shape.0,shape.1));
            (tuple[2],sm_index(tuple[2],n))
        } else {let c=gd_block_perm(shape.0,shape.1,&mut rng);(c,sm_index(c,n))};
        let id=ff_product_index(f.p,c,n);
        let total=f.indices[0]+f.indices[1]+ic+id;
        if f.indices.contains(&0)||id==0||ic==0||total<2*n||total%2!=0 || (target_genus>0&&1+total/2!=n+target_genus) {continue;}
        let mut t=if mode=="triangle"{Some(f.tuple(c,n))}else{None};
        if mode=="triangle" {
            let tuple=t.unwrap();let mut passport=tuple.map(|p|sm_index(p,n));passport.sort_unstable();
            if passport[1]+passport[2]>n-1||passport[0]+passport[3]<n {
                triangle_budget_reject+=1;continue;
            }
            if passport[0]<triangle_passport_min||passport[0]>triangle_passport_max {
                triangle_shard_reject+=1;continue;
            }
            if triangle_exact_profile.is_some_and(|expected|passport!=expected) {
                triangle_profile_reject+=1;continue;
            }
            // Riemann--Hurwitz: g = 1 + (sum ind(g_i) - 2d)/2.
            let genus=1+(total-2*n)/2;
            if genus<triangle_min_genus.max(3){triangle_genus_reject+=1;continue;}
            // Reviewed C12/C13 necessary inequality for a positive-genus
            // TRIANGLE component: h >= floor(d/2) + 2g + 1. It is especially
            // useful here because the discovery target is exact source genus 3.
            if passport[3]<n/2+2*genus+1 {triangle_high_index_bound_reject+=1;continue;}
            // C12 independently reviewed all-degree exclusion of the genus-3
            // distinguished 15-cycle/transposition-bulk family. Require the
            // entry to have maximum index before applying the theorem.
            if genus==3&&(0..4).any(|j|sm_index(tuple[j],n)==passport[3]&&sn_fifteen_transposition_bulk(tuple[j],n)) {
                triangle_fifteen_bulk_reject+=1;continue;
            }
            let equal_low_budget=n%2==0&&passport[0]==passport[1]&&passport[1]==passport[2]
                &&passport[0]==n/2-1;
            if equal_low_budget {
                if (0..3).any(|j|sn_is_involution(tuple[j],n)) {
                    triangle_low_involution_reject+=1;continue;
                }
                if (0..3).filter(|&j|sn_three_two_bulk_with_three_fixed(tuple[j],n)).count()>=2 {
                    triangle_two_special_reject+=1;continue;
                }
            }
            // In ties, any maximum-index entry may serve as the distinguished
            // fourth branch; the theorem requires such an entry to be FPF.
            let highest=passport[3];
            let has_fixed_point_free_high=(0..4).any(|j|sm_index(tuple[j],n)==highest
                && (0..n).all(|x|tuple[j].img[x]as usize!=x));
            if !has_fixed_point_free_high {
                triangle_fixed_free_reject+=1;continue;
            }
        }
        count[1]+=1;
        // Preserve legacy counter ordering for direct comparison. New mode
        // uses cheaper passport/short-word rejection before transitivity.
        if legacy {let tuple=f.tuple(c,n);if !generators_transitive(n,&tuple[..3]){continue;}count[2]+=1;t=Some(tuple);}
        if f.policy.allowed_v[ic]&(1u64<<id)==0 {continue;}count[3]+=1;
        let short_code=if reference {None}else{f.short_word(c,ic,id,n)};
        let word=if reference {sc_word(t.unwrap_or_else(||f.tuple(c,n)),n)} else {short_code.map(|j|FF_WORDS[j].to_vec())};
        if mode=="triangle" {
            if let Some(j)=short_code {
                triangle_short_rejects[j]+=1;
                if triangle_short_examples[j].is_none(){
                    triangle_short_examples[j]=Some(tuple_string_dynamic(n,&f.tuple(c,n)));
                }
            }
        }
        if word.is_some(){count[4]+=1;continue;}
        if let Some(failure)=f.deep_failure(c,ic,id,n,deep){
            count[4]+=1;deep_reject+=1;
            if mode=="triangle" {
                triangle_deep_rejects+=1;
                if triangle_deep_examples.len()<32 {
                    triangle_deep_examples.push((format!("{failure:?}"),
                        f.deep_word(failure,n).iter().map(|m|m.to_string()).collect::<Vec<_>>().join(","),
                        tuple_string_dynamic(n,&f.tuple(c,n))));
                }
            }
            continue;
        }
        let tuple=t.unwrap_or_else(||f.tuple(c,n));
        if !legacy {
            // For a geodesic pair, p-cycles ARE the pair's point orbits.
            if !(if orbit_cache{f.transitive(c,n)}else{generators_transitive(n,&[f.p,c])}){continue;}count[2]+=1;
        }
        let key=ff_key(&tuple,n);if closed.contains_key(&key)||seen.contains(&key){count[5]+=1;continue;}
        if known.contains(&key) {
            count[5]+=1;known_hits+=1;
            seen.insert(key); // Exact checked membership: log once per job.
            writeln!(aliases,"{}\t{}",count[0],tuple_string_dynamic(n,&tuple)).map_err(|e|e.to_string())?;
            continue;
        }
        let candidate=count[0];let raw=tuple_string_dynamic(n,&tuple);
        if let Some(p)=proposal.as_ref() {
            writeln!(mechanism_log,"{candidate}\t{}\t{}\t{}\t{}",p.source,p.arm,p.points.iter().map(|x|(x+1).to_string()).collect::<Vec<_>>().join(","),mb_primitive(&tuple,n)).map_err(|e|e.to_string())?;
        }
        writeln!(leads,"{candidate}\t{raw}").map_err(|e|e.to_string())?;
        leads.flush().map_err(|e|e.to_string())?; // Must precede closure for recovery.
        // Proposal time and admitted closure time are distinct budgets.
        let closure_deadline=Instant::now()+Duration::from_secs(closure_seconds);
        let limit=if args.iter().any(|s|s=="--full-closure-budget"){closure_deadline}else{deadline.min(closure_deadline)};
        let scan=if reference{sm_scan(&tuple,n,closure_cap,limit)}else{ff_scan_cached(&tuple,n,closure_cap,limit,&mut negative_cache,65536)};
        // Unresolved is never cached as a negative result or silently retired.
        // Avoid repeating an unresolved start within this bounded job, but
        // retain its durable lead and status. This is NOT a negative cache.
        if seen.len()<1_000_000{seen.insert(key);}
        let mut folder=String::new();
        if scan.status!="nonsmooth" {
            folder=format!("{}_{candidate:09}",if scan.status=="closed_smooth"{"component"}else{"unresolved"});
            sm_save(&tuple,n,&scan,&out.join(&folder))?;
            if scan.status=="closed_smooth" {count[6]+=1;for &s in &scan.states{closed.insert(s,folder.clone());}}else{count[7]+=1;}
        }
        writeln!(resolutions,"{candidate}\t{}\t{}\t{folder}\t{raw}",scan.status,scan.word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
        if Instant::now()>=next_flush {
            resolutions.flush().map_err(|e|e.to_string())?;frame_log.flush().map_err(|e|e.to_string())?;
            println!("degree={n} lane={mode} constructed={} closed={} unresolved={}",count[0],count[6],count[7]);
            next_flush=Instant::now()+Duration::from_secs(1);
        }
    }
    leads.flush().map_err(|e|e.to_string())?;resolutions.flush().map_err(|e|e.to_string())?;frame_log.flush().map_err(|e|e.to_string())?;
    aliases.flush().map_err(|e|e.to_string())?;
    mechanism_log.flush().map_err(|e|e.to_string())?;
    if let Some(p)=primitive.as_ref(){p.save(&out,&primitive_name,n)?;}
    if let Some(p)=double.as_ref(){p.save(&out)?;}
    if let Some(p)=mutation.as_ref(){fs::write(out.join("MUTATION_ARMS.tsv"),format!("cross_transposition\tinternal_transposition\tcross_three_cycle\ttwo_cross_transpositions\n{}\n",p.arm_counts.iter().map(|v|v.to_string()).collect::<Vec<_>>().join("\t"))).map_err(|e|e.to_string())?;}
    if args.iter().any(|s|s=="--export-negative-cache") {
        let mut file=BufWriter::new(File::create(out.join("NEGATIVE_CACHE.tsv")).map_err(|e|e.to_string())?);
        let mut keys:Vec<_>=negative_cache.keys().copied().collect();keys.sort();
        for key in keys.into_iter().take(4096) {
            writeln!(file,"{}\t{}",tuple_string_dynamic(n,&sm_unpack(key,n)),negative_cache[&key].iter().map(|m|m.to_string()).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
        }
        file.flush().map_err(|e|e.to_string())?;
    }
    fs::write(out.join("DESIGN.tsv"),format!("compact_frames\tstride\trng_final\tlast_frame_rng\tknown_hits\tjoint_failed\n{compact}\t4096\t{}\t{last_frame_rng}\t{known_hits}\t{joint_failed}\n",rng.state)).map_err(|e|e.to_string())?;
    fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\t{mode}\t{seed}\t{}\t{}\n",count.iter().map(|v|v.to_string()).collect::<Vec<_>>().join("\t"),start.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
    fs::write(out.join("POLICY.tsv"),format!("legacy_proposals\treference\tframe_reuse\tframes\tfeasible_draws\tnegative_cache\tcounter_order\tdeep_level\tdeep_reject\torbit_cache\n{legacy}\t{reference}\t{reuse}\t{frames}\t{feasible_draws}\t{}\t{}\t{deep}\t{deep_reject}\t{orbit_cache}\n",negative_cache.len(),if legacy{"scope,transitive,passport,short"}else{"scope,passport,short,transitive"})).map_err(|e|e.to_string())?;
    if mode=="triangle" {
        let exact_label=triangle_exact_profile.map_or_else(String::new,|p|p.iter().map(|x|x.to_string()).collect::<Vec<_>>().join(","));
        fs::write(out.join("TRIANGLE_FILTERS.tsv"),format!("degree\tminimum_genus\tpair_index_min\tpair_index_max\tpassport_min\tpassport_max\texact_passport\tthird_method\tpassport_profiles\tstrong_divisor\tbudget_reject\tshard_reject\tprofile_reject\tgenus_reject\thigh_index_bound_reject\tfifteen_cycle_bulk_reject\thigh_fixed_point_reject\tequal_low_involution_reject\tequal_low_two_special_reject\tshort_0\tshort_1\tshort_2\tshort_3\tshort_4\tdeep_reject\n{n}\t{}\t{triangle_low_min}\t{triangle_low_max}\t{triangle_passport_min}\t{triangle_passport_max}\t{exact_label}\t{triangle_method}\t{}\t{triangle_strength}\t{triangle_budget_reject}\t{triangle_shard_reject}\t{triangle_profile_reject}\t{triangle_genus_reject}\t{triangle_high_index_bound_reject}\t{triangle_fifteen_bulk_reject}\t{triangle_fixed_free_reject}\t{triangle_low_involution_reject}\t{triangle_two_special_reject}\t{}\t{triangle_deep_rejects}\n",triangle_min_genus.max(3),triangle_profiles.len(),triangle_short_rejects.iter().map(|x|x.to_string()).collect::<Vec<_>>().join("\t"))).map_err(|e|e.to_string())?;
        let mut examples=BufWriter::new(File::create(out.join("TRIANGLE_SHORT_EXAMPLES.tsv")).map_err(|e|e.to_string())?);
        writeln!(examples,"factor\tword\tseed").map_err(|e|e.to_string())?;
        for (j,seed) in triangle_short_examples.iter().enumerate(){if let Some(seed)=seed{
            writeln!(examples,"{j}\t{}\t{seed}",FF_WORDS[j].iter().map(|m|m.to_string()).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
        }}
        examples.flush().map_err(|e|e.to_string())?;
        let mut examples=BufWriter::new(File::create(out.join("TRIANGLE_DEEP_EXAMPLES.tsv")).map_err(|e|e.to_string())?);
        writeln!(examples,"failure\tword\tseed").map_err(|e|e.to_string())?;
        for (failure,word,seed) in &triangle_deep_examples{writeln!(examples,"{failure}\t{word}\t{seed}").map_err(|e|e.to_string())?;}
        if triangle_method=="middle_repair" {
            fs::write(out.join("TRIANGLE16_REPAIR.tsv"),format!(
                "pair_success\tpair_failure\tc_attempts\tc_generation_failure\tc_index_reject\tproduct_index_reject\tshort_0\tshort_1\tshort_2\tshort_3\tshort_4\tdeep_reject\tinvolution_reject\tspecial_pair_reject\taccepted\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                triangle16_repair_stats.pair_success,triangle16_repair_stats.pair_failure,
                triangle16_repair_stats.c_attempts,triangle16_repair_stats.c_generation_failure,
                triangle16_repair_stats.c_index_reject,triangle16_repair_stats.product_index_reject,
                triangle16_repair_stats.short_reject[0],triangle16_repair_stats.short_reject[1],
                triangle16_repair_stats.short_reject[2],triangle16_repair_stats.short_reject[3],
                triangle16_repair_stats.short_reject[4],triangle16_repair_stats.deep_reject,
                triangle16_repair_stats.involution_reject,triangle16_repair_stats.special_pair_reject,
                triangle16_repair_stats.accepted)).map_err(|e|e.to_string())?;
        }
        examples.flush().map_err(|e|e.to_string())?;
    }
    println!("FINISHED degree={n} lane={mode} constructed={} seconds={:.6}",count[0],start.elapsed().as_secs_f64());
    Ok(())
}

fn ff_middle_gpu_bench(args:&[String])->Result<(),String> {
    let arg=|s:&str,d:&str|args.windows(2).find(|w|w[0]==s).map(|w|w[1].clone()).unwrap_or(d.into());
    let n:usize=arg("--d","14").parse().map_err(|_|"degree")?;
    let size:usize=arg("--size","262144").parse().map_err(|_|"size")?;
    if !(9..=14).contains(&n)||size==0{return Err("GPU trace bounds".into());}
    let path=PathBuf::from(arg("--out",""));
    let mut file=BufWriter::new(OpenOptions::new().create_new(true).write(true).open(path).map_err(|e|e.to_string())?);
    file.write_all(&(n as u32).to_le_bytes()).map_err(|e|e.to_string())?;
    file.write_all(&(size as u32).to_le_bytes()).map_err(|e|e.to_string())?;
    let mut rng=SplitMix64::new(202610123);let frames:Vec<_>=(0..128).map(|_|ForestFrame::new(gd_forest(n,&mut rng),n)).collect();
    let mut trace=Vec::with_capacity(size);let mut proposals=0usize;
    while trace.len()<size {
        let j=rng.below(frames.len());let f=&frames[j];if f.policy.third_indices.is_empty(){continue;}
        let ic=f.policy.third_indices[rng.below(f.policy.third_indices.len())];let c=ff_index_perm(n,ic,&mut rng);let id=ff_product_index(f.p,c,n);proposals+=1;
        if f.policy.allowed_v[ic]&(1<<id)==0||f.short_word(c,ic,id,n).is_some(){continue;}
        let failure=f.deep_failure(c,ic,id,n,1);let expected=match failure{None=>0,Some(DeepFailure::Middle(1,2))=>1,Some(DeepFailure::Middle(4,2))=>2,_=>unreachable!()};
        if let Some(failure)=failure {assert!(ff_replay_negative(f.tuple(c,n),n,&f.deep_word(failure,n)));}
        for p in [f.pair[0],f.pair[1],f.b_inverse,c,c.inverse(n)] {file.write_all(&p.img[..n]).map_err(|e|e.to_string())?;}
        file.write_all(&[f.indices[0] as u8,f.indices[1] as u8,ic as u8,id as u8,expected]).map_err(|e|e.to_string())?;
        trace.push((j,c,ic,id,expected));
    }
    file.flush().map_err(|e|e.to_string())?;
    println!("stage\tdegree\trepeat\tvariant\tseconds\tchecksum\tunits");
    for repeat in 0..=7 {
        let start=Instant::now();let mut sum=0usize;
        for &(j,c,ic,id,expected) in &trace {
            let failure=frames[j].deep_failure(std::hint::black_box(c),ic,id,n,1);
            let result=match failure{None=>0,Some(DeepFailure::Middle(1,2))=>1,_=>2};
            assert_eq!(result,expected);sum+=result as usize;
        }
        if repeat>0 {println!("middle_square\t{n}\t{repeat}\trust_cpu\t{}\t{sum}\t{size}",start.elapsed().as_secs_f64());}
    }
    println!("# proposals={proposals} every_negative_literal_verified=true");Ok(())
}

// Benchmark implementation below keeps the trace in RAM and excludes generation
// unless generation is the named stage. Paired repeat order alternates.
fn ff_bench(args:&[String])->Result<(),String> {
    let arg=|s:&str,d:&str|args.windows(2).find(|w|w[0]==s).map(|w|w[1].clone()).unwrap_or(d.into());
    let n:usize=arg("--d","12").parse().map_err(|_|"degree")?;
    let size:usize=arg("--size","100000").parse().map_err(|_|"size")?;
    let rounds:usize=arg("--rounds","7").parse().map_err(|_|"rounds")?;
    if !(9..=14).contains(&n)||size==0||rounds<3{return Err("benchmark bounds".into());}
    let mut rng=SplitMix64::new(0x20261012+n as u64);
    let frames:Vec<_>=(0..64).map(|_|ForestFrame::new(gd_forest(n,&mut rng),n)).collect();
    let trace:Vec<_>=(0..size).map(|j|{
        let index=1+rng.below(n-1);let c=ff_index_perm(n,index,&mut rng);
        let f=j%frames.len();let id=ff_product_index(frames[f].p,c,n);(f,c,index,id)
    }).collect();
    for &(j,c,ic,id) in &trace {
        assert_eq!(sc_word(frames[j].tuple(c,n),n),frames[j].short_word(c,ic,id,n).map(|k|FF_WORDS[k].to_vec()));
    }
    if let Some(path)=args.windows(2).find(|w|w[0]=="--gpu-trace").map(|w|&w[1]) {
        let mut out=BufWriter::new(OpenOptions::new().create_new(true).write(true).open(path).map_err(|e|e.to_string())?);
        out.write_all(&(n as u32).to_le_bytes()).map_err(|e|e.to_string())?;
        out.write_all(&(size as u32).to_le_bytes()).map_err(|e|e.to_string())?;
        for &(j,c,ic,id) in &trace {
            let f=&frames[j];for &p in f.factors.iter().chain(std::iter::once(&c)){out.write_all(&p.img[..n]).map_err(|e|e.to_string())?;}
            out.write_all(&[f.indices[0] as u8,f.indices[1] as u8,ic as u8,id as u8,f.short_word(c,ic,id,n).map_or(0,|j|j+1) as u8]).map_err(|e|e.to_string())?;
        }
        out.flush().map_err(|e|e.to_string())?;println!("GPU_TRACE degree={n} candidates={size} literal_verified=true");return Ok(());
    }
    println!("stage\tdegree\trepeat\tvariant\tseconds\tchecksum\tunits");
    for stage in ["short_screen","preclosure","key","closure"] {
        for repeat in 0..=rounds {
            for fast in if repeat%2==0{[false,true]}else{[true,false]} {
                let start=Instant::now();let mut sum=0usize;
                let mut units=0;
                for &(j,c,ic,id) in &trace {
                    let f=&frames[j];
                    if stage=="closure" {
                        if id==0||f.indices[0]+f.indices[1]+ic+id<2*n||f.policy.allowed_v[ic]&(1<<id)==0||f.short_word(c,ic,id,n).is_some()||!generators_transitive(n,&[f.p,c]){continue;}
                        if units==256{break;}
                        units+=1;let t=f.tuple(c,n);let deadline=Instant::now()+Duration::from_secs(10);
                        let scan=if fast{ff_scan(&t,n,100000,deadline)}else{sm_scan(&t,n,100000,deadline)};
                        assert_eq!(scan.status,"nonsmooth");assert!(ff_replay_negative(t,n,&scan.word));sum+=1;
                        continue;
                    }
                    units+=1;
                    let v=match stage {
                        "short_screen" => if fast{f.short_word(c,ic,id,n).is_some()}else{sc_word(f.tuple(c,n),n).is_some()},
                        "key" => {let t=f.tuple(c,n);if !generators_transitive(n,&t[..3]){false}else{let key=if fast{ff_key(&t,n)}else{sm_key(&t,n)};std::hint::black_box(key);true}},
                        _ => {
                            let total=f.indices[0]+f.indices[1]+ic+id;
                            if id==0||total<2*n||total%2!=0{false}
                            else if fast {
                                f.policy.allowed_v[ic]&(1u64<<id)!=0&&f.short_word(c,ic,id,n).is_none()&&generators_transitive(n,&[f.p,c])
                            }else{
                                let t=f.tuple(c,n);generators_transitive(n,&t[..3])&&gd_passport(&t,n)&&sc_word(t,n).is_none()
                            }
                        }
                    };sum+=std::hint::black_box(v) as usize;
                }
                if repeat>0{println!("{stage}\t{n}\t{repeat}\t{}\t{}\t{sum}\t{units}",if fast{"fast"}else{"old"},start.elapsed().as_secs_f64());}
            }
        }
    }
    Ok(())
}

fn ff_io_bench(out:&Path)->Result<(),String>{
    if out.as_os_str().is_empty(){return Err("fresh I/O output required".into());}
    fs::create_dir(out).map_err(|e|e.to_string())?;
    println!("stage\tdegree\trepeat\tvariant\tseconds\tchecksum\tunits");
    for repeat in 0..=7{for fast in if repeat%2==0{[false,true]}else{[true,false]}{
        let prefix=out.join(format!("{repeat}_{}",if fast{"fast"}else{"old"}));
        let lead_path=prefix.with_extension("leads");let res_path=prefix.with_extension("resolutions");
        let start=Instant::now();
        let mut res=BufWriter::new(File::create(&res_path).map_err(|e|e.to_string())?);
        let mut lead=if fast{Some(BufWriter::new(File::create(&lead_path).map_err(|e|e.to_string())?))}else{None};
        for id in 0..200{
            let record=format!("{id}\t[(1 2 3), (2 4 5), (3 5 6), (1 6 4)]\n");
            if let Some(f)=lead.as_mut(){f.write_all(record.as_bytes()).map_err(|e|e.to_string())?;f.flush().map_err(|e|e.to_string())?;}
            else{let mut f=OpenOptions::new().create(true).append(true).open(&lead_path).map_err(|e|e.to_string())?;f.write_all(record.as_bytes()).map_err(|e|e.to_string())?;f.flush().map_err(|e|e.to_string())?;}
            res.write_all(record.as_bytes()).map_err(|e|e.to_string())?;
            if !fast{res.flush().map_err(|e|e.to_string())?;}
        }
        res.flush().map_err(|e|e.to_string())?;drop(lead);drop(res);
        let elapsed=start.elapsed().as_secs_f64();
        assert_eq!(fs::read(&lead_path).unwrap(),fs::read(&res_path).unwrap());
        if repeat>0{println!("t7_logging\t0\t{repeat}\t{}\t{elapsed}\t200\t200",if fast{"fast"}else{"old"});}
    }}Ok(())
}

#[cfg(test)]
mod forest_fast_tests {
    use super::*;
    #[test]fn index_policy_is_exact_and_shared_for_every_index_triple() {
        for n in 2..=25 {for ia in 0..n {for ib in 0..n {
            let policy=ff_index_policy(n,ia,ib);
            assert!(std::ptr::eq(policy,ff_index_policy(n,ia,ib)));
            let mut expected_indices=Vec::new();
            for t in 0..n {
                let mut mask=0u64;
                for v in 0..n {
                    let expected=t>0 && v>0 && ia+ib+t+v>=2*n
                        && (ia+ib+t+v)%2==0
                        && (ia+t<n || ib+v<n) && (ia+v<n || ib+t<n);
                    if expected {mask|=1u64<<v;}
                }
                assert_eq!(policy.allowed_v[t],mask);
                if mask!=0 {expected_indices.push(t);}
            }
            assert_eq!(policy.third_indices,expected_indices);
        }}}
    }
    fn permutations(n:usize)->Vec<Perm> {
        fn visit(a:&mut Perm,n:usize,i:usize,out:&mut Vec<Perm>){if i==n{out.push(*a);return;}for j in i..n{a.img.swap(i,j);visit(a,n,i+1,out);a.img.swap(i,j);}}
        let mut out=vec![];visit(&mut Perm::identity(n),n,0,&mut out);out
    }
    fn check_deep(f:&ForestFrame,c:Perm,n:usize) {
        let t=f.tuple(c,n);let ic=sm_index(c,n);let id=sm_index(t[3],n);
        let ix=t.map(|p|sm_index(p,n));
        for factor in f.twist_factors(n) {
            let mut u=t;for &m in &factor.word{u=hurwitz_move_dynamic(n,u,m);}
            assert_eq!(ff_product_index(factor.factor,c,n),sm_index(u[0].compose(u[1],n),n));
            assert_eq!(factor.first_index+ic,sm_index(u[0],n)+sm_index(u[1],n));
        }
        for level in 1..=3 {
            if let Some(failure)=f.deep_failure(c,ic,id,n,level){
                let word=f.deep_word(failure,n);assert!(ff_replay_negative(t,n,&word));
            }else{
                for middle in [1,4] {
                    let mut u=t;
                    for power in 1..=if level==1{2}else{3} {
                        u=hurwitz_move_dynamic(n,u,middle);
                        if power>=2{assert!(sm_values(&u,n).0,"missed literal middle failure {middle},{power}; indices {ix:?}");}
                    }
                }
            }
        }
        assert_eq!(f.transitive(c,n),generators_transitive(n,&[f.p,c]));
    }
    #[test]fn all_s4_deep_filters_match_literal_words_and_orbit_masks() {
        let n=4;let all=permutations(n);
        for &a in &all{for &b in &all{
            let f=ForestFrame::new([a,b],n);
            for &c in &all{check_deep(&f,c,n);}
        }}
    }
    #[test]fn larger_deep_filters_and_relabeling_match_literals() {
        let mut rng=SplitMix64::new(202610123);
        for n in 5..=25 {for _ in 0..500 {
            let f=ForestFrame::new([gd_shuffle(n,&mut rng),gd_shuffle(n,&mut rng)],n);
            let c=gd_shuffle(n,&mut rng);check_deep(&f,c,n);
            let relabel=gd_shuffle(n,&mut rng);
            let g=ForestFrame::new(f.pair.map(|p|p.conjugate_by(relabel,n)),n);
            let other=c.conjugate_by(relabel,n);let ic=sm_index(c,n);let id=ff_product_index(f.p,c,n);
            for level in 1..=3{
                assert_eq!(f.deep_failure(c,ic,id,n,level).map(|w|f.deep_word(w,n)),
                           g.deep_failure(other,ic,id,n,level).map(|w|g.deep_word(w,n)));
            }
        }}
    }
    #[test]fn all_s4_short_words_and_keys_match() {
        let n=4;let all=permutations(n);
        for &a in &all{for &b in &all{let f=ForestFrame::new([a,b],n);for &c in &all {
            let t=f.tuple(c,n);let ic=sm_index(c,n);let id=sm_index(t[3],n);
            assert_eq!(sc_word(t,n),f.short_word(c,ic,id,n).map(|j|FF_WORDS[j].to_vec()));
            assert_eq!(ff_product_index(a,b,n),sm_index(a.compose(b,n),n));
            if generators_transitive(n,&t[..3]){assert_eq!(ff_key(&t,n),sm_key(&t,n));}
        }}}
    }
    #[test]fn target_cycle_count_matches_literal_all_s4_and_larger() {
        let all=permutations(4);
        for &a in &all {for &b in &all {for l in 0..=5 {for r in 0..=5 {
            let i=sm_index(a.compose(b,4),4);
            assert_eq!(ff_product_matches(a,b,4,l,r),i==l||i==r);
        }}}}
        let mut rng=SplitMix64::new(20261012);
        for n in 9..=25 {for _ in 0..2000 {
            let a=gd_shuffle(n,&mut rng);let b=gd_shuffle(n,&mut rng);
            let l=rng.below(2*n);let r=rng.below(2*n);
            let i=sm_index(a.compose(b,n),n);
            assert_eq!(ff_product_matches(a,b,n,l,r),i==l||i==r);
        }}
    }
    #[test]fn all_geodesic_s4_passports_and_transitivity_match() {
        let n=4;let all=permutations(n);
        for &a in &all{for &b in &all{
            let f=ForestFrame::new([a,b],n);if f.indices.iter().sum::<usize>()!=sm_index(f.p,n){continue;}
            for &c in &all {
                let t=f.tuple(c,n);let ic=sm_index(c,n);let id=sm_index(t[3],n);
                let ix=t.map(|p|sm_index(p,n));
                let scope=ix.iter().all(|&i|i>0)&&ix.iter().sum::<usize>()>=2*n;
                assert_eq!(scope&&gd_passport(&t,n),scope&&f.policy.allowed_v[ic]&(1<<id)!=0);
                assert_eq!(generators_transitive(n,&t[..3]),generators_transitive(n,&[f.p,c]));
            }
        }}
    }
    #[test]fn stack_constructor_identical_rng_and_randomized_keys() {
        let mut old=SplitMix64::new(81);let mut new=SplitMix64::new(81);
        for n in 9..=14{for i in 0..n{for _ in 0..100{assert_eq!(gd_index_perm(n,i,&mut old).img,ff_index_perm(n,i,&mut new).img);}}}
        let mut targets=SplitMix64::new(20261012);
        for n in 2..=25{for i in 0..n{for _ in 0..100{
            assert_eq!(sm_index(ff_index_perm(n,i,&mut targets),n),i);
        }}}
        for n in 9..=14{for _ in 0..1000{
            let a=gd_shuffle(n,&mut new);let b=gd_shuffle(n,&mut new);let c=gd_shuffle(n,&mut new);
            let f=ForestFrame::new([a,b],n);let t=f.tuple(c,n);
            assert_eq!(sc_word(t,n),f.short_word(c,sm_index(c,n),sm_index(t[3],n),n).map(|j|FF_WORDS[j].to_vec()));
            if generators_transitive(n,&t[..3]){assert_eq!(ff_key(&t,n),sm_key(&t,n));}
        }}
    }
    #[test]fn complete_s4_scans_and_literal_negative_words() {
        let n=4;let all=permutations(n);
        for &a in &all{for &b in &all{for &c in &all{
            let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];
            if !generators_transitive(n,&t[..3])||!sm_values(&t,n).0{continue;}
            let deadline=Instant::now()+Duration::from_secs(60);
            let old=sm_scan(&t,n,10000,deadline);let new=ff_scan(&t,n,10000,deadline);
            assert_eq!(old.status,new.status);
            if new.status=="closed_smooth"{assert_eq!(old.states.iter().collect::<HashSet<_>>(),new.states.iter().collect::<HashSet<_>>());}
            else{let mut u=t;for m in new.word{u=hurwitz_move_dynamic(n,u,m);}assert!(!sm_values(&u,n).0);}
        }}}
    }
    #[test]fn known_positive_closed_control_and_caps() {
        let t=sm_parse("[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]",9).unwrap();
        let deadline=Instant::now()+Duration::from_secs(60);
        let scan=ff_scan(&t,9,1000,deadline);assert_eq!(scan.status,"closed_smooth");assert_eq!(scan.states.len(),12);
        assert_eq!(ff_scan(&t,9,1,deadline).status,"capped");
        assert_eq!(ff_scan(&t,9,1000,Instant::now()).status,"timeout");
    }
    #[test]fn negative_cache_only_contains_replayed_proofs() {
        let mut rng=SplitMix64::new(911);let mut cache=HashMap::new();
        for _ in 0..1000 {
            let f=ForestFrame::new(gd_forest(9,&mut rng),9);let c=gd_shuffle(9,&mut rng);let t=f.tuple(c,9);
            if !generators_transitive(9,&t[..3]){continue;}
            let scan=ff_scan_cached(&t,9,10000,Instant::now()+Duration::from_secs(1),&mut cache,1000);
            if scan.status=="nonsmooth"{assert!(ff_replay_negative(t,9,&scan.word));}
        }
        assert!(!cache.is_empty());
        for (key,word) in cache{assert!(ff_replay_negative(sm_unpack(key,9),9,&word));}
    }
}
