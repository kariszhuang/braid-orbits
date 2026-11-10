//! Export/import genuinely explicit two-braid representations for the rigid solver.
use super::*;

// This catalog acts on abstract quotient states, NOT on the d points of a
// four-tuple. It deliberately exports no counterexample certificate.
pub(super) fn fiber_actions(args:&[String])->Result<(),String> {
    let n:usize=parse_num(get_arg(args,"--d",Some("10"))?,"--d")?;
    if !(4..=10).contains(&n) {return Err("fiber-action census supports 4..10 abstract states".into());}
    let raw=get_arg(args,"--fiber-type",None)?;
    let mut lengths=raw.split('.').map(|s|parse_num::<usize>(s.into(),"fiber cycle length")).collect::<Result<Vec<_>,_>>()?;
    if lengths.contains(&0) || lengths.iter().sum::<usize>()!=n {return Err("fiber-type must include fixed points and sum to d".into());}
    lengths.sort_unstable();
    let u=perm_of_type(&lengths,n,None);
    let class=conjugacy_class(n,&lengths);
    let begin:usize=parse_num(get_arg(args,"--scan-start",Some("0"))?,"--scan-start")?;
    let count:usize=parse_num(get_arg(args,"--scan-count",Some("1000000"))?,"--scan-count")?;
    let seconds:u64=parse_num(get_arg(args,"--seconds",Some("60"))?,"--seconds")?;
    if begin>class.len() || count==0 || seconds==0 {return Err("invalid fiber action cursor/budget".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let mut rows=BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join("ACTIONS.tsv")).map_err(|e|e.to_string())?);
    writeln!(rows,"states\tauxiliary_genus\tgenerated_group_order\tu_zero_based\tv_zero_based\trepresentation_key").map_err(|e|e.to_string())?;
    let start=Instant::now();let mut next=begin;let mut braid=0;let mut modular=0;let mut connected=0;let mut seen=HashSet::new();
    for &v in class.iter().skip(begin).take(count) {
        if start.elapsed().as_secs()>=seconds {break;}
        next+=1;
        let uv=u.compose(v,n);let x=uv.compose(u,n);
        if x!=v.compose(u,n).compose(v,n) {continue;}braid+=1;
        if !uv.compose(uv,n).compose(uv,n).is_identity(n) {continue;}modular+=1;
        if !generators_transitive(n,&[u,v]) {continue;}connected+=1;
        let key=colored_representation_key(n,[u,v,Perm::identity(n)]);
        if !seen.insert(key.clone()) {continue;}
        let genus=1-n as i64+(index(x,n)+index(uv,n)+index(u,n)) as i64/2;
        if genus<0 {return Err("negative transitive auxiliary genus".into());}
        let images=|p:Perm|p.img[..n].iter().map(u8::to_string).collect::<Vec<_>>().join(",");
        writeln!(rows,"{n}\t{genus}\t{}\t{}\t{}\t{key}",subgroup_order(n,&[u,v]),images(u),images(v)).map_err(|e|e.to_string())?;
    }
    rows.flush().map_err(|e|e.to_string())?;
    let mut summary=OpenOptions::new().write(true).create_new(true).open(out.join("SUMMARY.txt")).map_err(|e|e.to_string())?;
    writeln!(summary,"search mode = abstract fiber action\nfiber type = {raw}\nscan start = {begin}\nscan next = {next}\nscan total = {}\nbraid compatible = {braid}\ncentral identity = {modular}\ntransitive solutions = {connected}\ninequivalent abstract actions = {}\nsearch complete = {}\nelapsed seconds = {:.3}\nThese are abstract B3 actions with (uv)^3=1. Realization by product-one generating four-tuples is NOT established. No weak counterexample is certified by this catalog.",class.len(),seen.len(),next==class.len() && begin==0,start.elapsed().as_secs_f64()).map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())
}

// Independent central lifts on mixed point blocks. On each cyclic fiber,
// x^2=y^3=z is enforced exactly, with an independently specified shift z.
// Unlike the old library, shifts need not be units and distinct blocks need
// not be actions of a single marked finite group.
fn voltage_block(n:usize,m:usize,shift:usize,rng:&mut SplitMix64)->Result<(Perm,Perm),String> {
    if m==0 || n%m!=0 || shift>=m || n>MAX_D {return Err("invalid coupled block degree/modulus/shift".into());}
    let base=n/m;let mut choices=Vec::new();
    for pairs in 0..=base/2 {for triples in 0..=base/3 {
        if pairs+2*triples>=base-1 && (base==2*pairs || (0..m).any(|s|2*s%m==shift))
            && (base==3*triples || (0..m).any(|s|3*s%m==shift)) {choices.push((pairs,triples));}
    }}
    if choices.is_empty() {return Err("incompatible fixed-fiber voltage equations".into());}
    for _ in 0..4096 {
        let (pairs,triples)=choices[rng.below(choices.len())];
        let mut xt=vec![1;base-2*pairs];xt.extend(vec![2;pairs]);
        let mut yt=vec![1;base-3*triples];yt.extend(vec![3;triples]);
        let x0=perm_of_type(&xt,base,None);let y0=perm_of_type(&yt,base,Some(rng));
        if !generators_transitive(base,&[x0,y0]) {continue;}
        let mut lift=|p:Perm,power:usize| {
            let mut q=Perm::identity(n);
            for cycle in permutation_cycles(base,p) {
                let roots=(0..m).filter(|&s|(power/cycle.len())*s%m==shift).collect::<Vec<_>>();
                let total=roots[rng.below(roots.len())];
                let mut phases=vec![0;cycle.len()];let mut sum=0;
                for phase in phases.iter_mut().take(cycle.len()-1) {*phase=rng.below(m);sum+=*phase;}
                phases[cycle.len()-1]=(total+m-sum%m)%m;
                for (i,&x) in cycle.iter().enumerate() {for j in 0..m {
                    q.img[x*m+j]=(cycle[(i+1)%cycle.len()]*m+(j+phases[i])%m) as u8;
                }}
            }q
        };
        let x=lift(x0,2);let y=lift(y0,3);
        let mut z=Perm::identity(n);for i in 0..base {for j in 0..m {z.img[i*m+j]=(i*m+(j+shift)%m) as u8;}}
        if x.compose(x,n)!=z || y.compose(y,n).compose(y,n)!=z {return Err("voltage construction invariant failed".into());}
        let h=y.inverse(n).compose(x,n);let k=x.inverse(n).compose(y.compose(y,n),n);
        if h.compose(k,n).compose(h,n)!=k.compose(h,n).compose(k,n) {return Err("coupled braid relation failed".into());}
        return Ok((h,k));
    }
    Err("bounded modular-base construction found no connected base".into())
}

pub(super) fn coupled(args:&[String])->Result<(usize,Perm,Perm,Vec<String>),String> {
    let spec=get_arg(args,"--representation-blocks",None)?;
    let seed:u64=parse_num(get_arg(args,"--representation-seed",Some("20261008"))?,"--representation-seed")?;
    let mut rng=SplitMix64::new(seed);let labels=spec.split(',').map(str::to_string).collect::<Vec<_>>();
    let mut blocks=Vec::new();let mut d=0usize;
    for label in &labels {
        let parts=label.split(':').map(|s|parse_num::<usize>(s.into(),"coupled block")).collect::<Result<Vec<_>,_>>()?;
        if parts.len()!=3 || parts[0]==0 || parts[1]==0 || parts[0]%parts[1]!=0 || parts[2]>=parts[1] {return Err("blocks must be degree:modulus:shift".into());}
        d=d.checked_add(parts[0]).ok_or("coupled degree overflow")?;blocks.push(parts);
    }
    if !(4..=MAX_D).contains(&d) {return Err("coupled total degree must be 4..34".into());}
    let (mut h,mut k,mut offset)=(Perm::identity(d),Perm::identity(d),0);
    for parts in blocks {
        let n=parts[0];let (a,b)=voltage_block(n,parts[1],parts[2],&mut rng)?;
        for i in 0..n {h.img[offset+i]=(offset+a.img[i] as usize) as u8;k.img[offset+i]=(offset+b.img[i] as usize) as u8;}
        offset+=n;
    }
    Ok((d,h,k,labels))
}

pub(super) fn read(path:&Path)->Result<(usize,Perm,Perm),String> {
    let text=fs::read_to_string(path).map_err(|e|e.to_string())?;
    let mut lines=text.lines();
    if lines.next()!=Some("HURWITZ_BRAID_REPRESENTATION_V1") {return Err("invalid braid-action header".into());}
    let mut fields=HashMap::new();
    for line in lines {if line.starts_with('#') || line.is_empty() {continue;}
        let (key,value)=line.split_once('=').ok_or("braid action fields must be key=value")?;
        if fields.insert(key,value).is_some() {return Err("duplicate braid-action field".into());}
    }
    if fields.len()!=3 {return Err("braid action needs exactly degree, h, k".into());}
    let d:usize=fields.get("degree").ok_or("missing degree")?.parse().map_err(|_|"invalid degree")?;
    if !(4..=MAX_D).contains(&d) {return Err("braid-action degree must be 4..34".into());}
    let parse=|name:&str|->Result<Perm,String> {
        let images=fields.get(name).ok_or("missing braid relabeling")?.split(',').map(|s|s.parse::<usize>().map_err(|_|"invalid image".to_string())).collect::<Result<Vec<_>,_>>()?;
        if images.len()!=d || images.iter().copied().collect::<HashSet<_>>()!=HashSet::from_iter(0..d) {return Err("braid-action images must be a permutation of 0..degree-1".into());}
        let mut p=Perm::identity(d);for (i,&x) in images.iter().enumerate() {p.img[i]=x as u8;}Ok(p)
    };
    let h=parse("h")?;let k=parse("k")?;
    if h.compose(k,d).compose(h,d)!=k.compose(h,d).compose(k,d) {return Err("external relabelings fail hkh=khk".into());}
    Ok((d,h,k))
}

fn write(path:&Path,d:usize,h:Perm,k:Perm)->Result<(),String> {
    let mut out=OpenOptions::new().write(true).create_new(true).open(path).map_err(|e|e.to_string())?;
    writeln!(out,"HURWITZ_BRAID_REPRESENTATION_V1\ndegree={d}\nh={}\nk={}",
        h.img[..d].iter().map(u8::to_string).collect::<Vec<_>>().join(","),
        k.img[..d].iter().map(u8::to_string).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
    out.sync_all().map_err(|e|e.to_string())
}

// A bounded new representation library, not a census of modular bases/lifts.
// Explicit inputs feed the SAME rigid propagation and six-move verifier.
pub(super) fn generate(args:&[String])->Result<(),String> {
    let d:usize=parse_num(get_arg(args,"--d",None)?,"--d")?;
    let m:usize=parse_num(get_arg(args,"--modulus",Some("3"))?,"--modulus")?;
    let draws:usize=parse_num(get_arg(args,"--draws",Some("32"))?,"--draws")?;
    let seed:u64=parse_num(get_arg(args,"--seed",Some("20261004"))?,"--seed")?;
    if !(4..=MAX_D).contains(&d) || m<2 || d%m!=0 || d/m<3 || !(1..=100000).contains(&draws) {
        return Err("representation library needs d 4..34, modulus>=2 dividing d, base>=3, draws 1..100000".into());
    }
    let out=PathBuf::from(get_arg(args,"--out",None)?);fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let mut rows=OpenOptions::new().write(true).create_new(true).open(out.join("REPRESENTATIONS.tsv")).map_err(|e|e.to_string())?;
    writeln!(rows,"draw\tdegree\tmodulus\tseed\tstatus\trepresentation_key\tfile").map_err(|e|e.to_string())?;
    let mut excluded=HashSet::new();
    if let Some(path)=optional_arg(args,"--exclude-representations") {
        let text=fs::read_to_string(path).map_err(|e|e.to_string())?;
        let header=text.lines().next().ok_or("empty representation exclusion table")?.split('\t').collect::<Vec<_>>();
        let column=header.iter().position(|&s|s=="representation_key").ok_or("exclusion table needs representation_key column")?;
        for line in text.lines().skip(1) {let fields=line.split('\t').collect::<Vec<_>>();
            let key=fields.get(column).ok_or("truncated representation exclusion row")?;
            if !key.is_empty() {excluded.insert(key.to_string());}
        }
    }
    let mut rng=SplitMix64::new(seed);let mut seen=HashSet::new();
    for draw in 0..draws {let local_seed=rng.next_u64();
        match modular_braid_block(d,m,local_seed) {
            Ok((h,k))=>{
                let key=braid_representation_key(d,h,k);
                let (status,file)=if excluded.contains(&key) {("excluded",String::new())}
                    else if !seen.insert(key.clone()) {("duplicate",String::new())}
                    else {let name=format!("action_{draw:06}.txt");write(&out.join(&name),d,h,k)?;("constructed",name)};
                writeln!(rows,"{draw}\t{d}\t{m}\t{local_seed}\t{status}\t{key}\t{file}").map_err(|e|e.to_string())?;
            }
            Err(_)=>{writeln!(rows,"{draw}\t{d}\t{m}\t{local_seed}\tbounded_construction_failed\t\t").map_err(|e|e.to_string())?;}
        }
    }
    rows.sync_all().map_err(|e|e.to_string())
}

// Search small STATE actions, a different object from point relabelings h,k.
// Every generator is conjugate to q1. Fix q1 to a cycle-type representative,
// then enumerate q2 in its class and q3 in its commuting class slice. These
// constraints are necessary for a product-one conjugacy-quotient action, not
// sufficient for realization by a tuple in a requested ambient group.
pub(super) fn state_actions(args:&[String])->Result<(),String> {
    let n:usize=parse_num(get_arg(args,"--d",None)?,"--d")?;
    let seconds:u64=parse_num(get_arg(args,"--seconds",Some("10"))?,"--seconds")?;
    let cap:u64=parse_num(get_arg(args,"--max-tests",Some("100000"))?,"--max-tests")?;
    let first:usize=parse_num(get_arg(args,"--start-type",Some("0"))?,"--start-type")?;
    let second:usize=parse_num(get_arg(args,"--start-second",Some("0"))?,"--start-second")?;
    let third:usize=parse_num(get_arg(args,"--start-third",Some("0"))?,"--start-third")?;
    if !(2..=8).contains(&n) || seconds==0 || cap==0 {return Err("small state actions need degree 2..8 and positive time/test limits".into());}
    fn partitions(left:usize,minimum:usize,current:&mut Vec<usize>,out:&mut Vec<Vec<usize>>) {
        if left==0 {out.push(current.clone());return;}
        for i in minimum..=left {current.push(i);partitions(left-i,i,current,out);current.pop();}
    }
    let mut types=Vec::new();partitions(n,1,&mut Vec::new(),&mut types);
    if first>=types.len() {return Err("start-type outside cycle-type list".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let mut rows=OpenOptions::new().write(true).create_new(true).open(out.join("STATE_ACTIONS.tsv")).map_err(|e|e.to_string())?;
    writeln!(rows,"state_count\tcycle_type\trepresentation_key\tfile\tno_generator_fixed_anywhere").map_err(|e|e.to_string())?;
    let mut summary=OpenOptions::new().write(true).create_new(true).open(out.join("SUMMARY.txt")).map_err(|e|e.to_string())?;
    let deadline=Instant::now()+Duration::from_secs(seconds);let mut next=(first,second,third);
    let (mut tests,mut constructed)=(0u64,0usize);let mut seen=HashSet::new();let mut complete=true;
    'types: for (ti,ty) in types.iter().enumerate().skip(first) {
        let a=perm_of_type(ty,n,None);let class=conjugacy_class(n,ty);
        let commuting=class.iter().copied().filter(|&c|a.compose(c,n)==c.compose(a,n)).collect::<Vec<_>>();
        let begin=if ti==first {second} else {0};
        if begin>=class.len() {return Err("start-second outside class".into());}
        for (bi,&b) in class.iter().enumerate().skip(begin) {
            if Instant::now()>=deadline || tests>=cap {complete=false;next=(ti,bi,if ti==first && bi==second {third} else {0});break 'types;}
            let ab=a.compose(b,n);
            if a.compose(b,n).compose(a,n)!=b.compose(a,n).compose(b,n) || !small_power(n,ab,3).is_identity(n) {
                continue;
            }
            let start=if ti==first && bi==second {third} else {0};
            if start>=commuting.len() {return Err("start-third outside commuting class".into());}
            for (ci,&c) in commuting.iter().enumerate().skip(start) {
                if Instant::now()>=deadline || tests>=cap {complete=false;next=(ti,bi,ci);break 'types;}
                tests+=1;
                if b.compose(c,n).compose(b,n)!=c.compose(b,n).compose(c,n)
                    || !small_power(n,ab.compose(c,n),4).is_identity(n)
                    || !a.compose(b,n).compose(c,n).compose(c,n).compose(b,n).compose(a,n).is_identity(n)
                    || !generators_transitive(n,&[a,b,c]) {continue;}
                let no_fixed=[a,b,c].iter().all(|p|(0..n).all(|i|p.img[i] as usize!=i));
                if has_flag(args,"--no-fixed-anywhere") && !no_fixed {continue;}
                let key=colored_representation_key(n,[a,b,c]);if !seen.insert(key.clone()) {continue;}
                let name=format!("state_action_{constructed:06}.txt");constructed+=1;
                let mut file=OpenOptions::new().write(true).create_new(true).open(out.join(&name)).map_err(|e|e.to_string())?;
                writeln!(file,"HURWITZ_BRAID_STATE_ACTION_V1\ndegree={n}").map_err(|e|e.to_string())?;
                for (i,p) in [a,b,c].iter().enumerate() {writeln!(file,"sigma{}={}",i+1,p.img[..n].iter().map(u8::to_string).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;}
                file.sync_all().map_err(|e|e.to_string())?;
                writeln!(rows,"{n}\t{ty:?}\t{key}\t{name}\t{no_fixed}").map_err(|e|e.to_string())?;
            }
        }
        next=(ti+1,0,0);
    }
    // Normalize exhausted inner cursors so a resume never points one past a list.
    if !complete && next.0<types.len() {
        let class=conjugacy_class(n,&types[next.0]);
        if next.1==class.len() {next=(next.0+1,0,0);}
        else {let a=perm_of_type(&types[next.0],n,None);
            let len=class.iter().filter(|&&c|a.compose(c,n)==c.compose(a,n)).count();
            if next.2==len {next=(next.0,next.1+1,0);if next.1==class.len() {next=(next.0+1,0,0);}}
        }
        if next.0==types.len() {complete=true;}
    }
    writeln!(summary,"search mode = small abstract braid-state actions\nstate count = {n}\nthird-generator tests = {tests}\ndistinct retained move actions = {constructed}\ntraversal exhausted = {complete}\nnext type = {}\nnext second = {}\nnext third = {}\nCoverage is finite permutation actions satisfying the checked relations, not realizable tuples or counterexamples; deduplication is per invocation.",next.0,next.1,next.2).map_err(|e|e.to_string())?;
    rows.sync_all().map_err(|e|e.to_string())?;summary.sync_all().map_err(|e|e.to_string())
}

fn small_power(d:usize,p:Perm,n:usize)->Perm {(0..n).fold(Perm::identity(d),|q,_|q.compose(p,d))}

#[cfg(test)] mod tests {
    use super::*;
    #[test] fn mixed_nonunit_voltages_obey_braid_and_are_reproducible() {
        let args=["binary","rigid-glue","--representation-blocks","9:3:2,8:2:1","--representation-seed","71"].map(str::to_string);
        let (d,h,k,labels)=coupled(&args).unwrap();let again=coupled(&args).unwrap();
        assert_eq!(d,17);assert_eq!((h,k),(again.1,again.2));assert_eq!(labels.len(),2);
        assert_eq!(h.compose(k,d).compose(h,d),k.compose(h,d).compose(k,d));
        let mut rng=SplitMix64::new(15);let (a,b)=voltage_block(12,4,2,&mut rng).unwrap();
        assert_eq!(a.compose(b,12).compose(a,12),b.compose(a,12).compose(b,12));
        assert!(voltage_block(9,3,1,&mut rng).is_ok());
        assert!(voltage_block(9,2,1,&mut rng).is_err());
    }
    #[test] fn explicit_representations_round_trip_and_match_rigid_constraints() {
        let temporary=std::env::temp_dir().join(format!("hurwitz-braid-action-{}-{}",std::process::id(),SplitMix64::new(20261004).next_u64()));
        fs::create_dir(&temporary).unwrap();let path=temporary.join("action.txt");
        let (h,k)=modular_braid_block(9,3,19).unwrap();write(&path,9,h,k).unwrap();assert_eq!(read(&path).unwrap(),(9,h,k));
        let center=power_perm(9,h.compose(k,9),3);assert!(!center.is_identity(9));
        fs::remove_file(path).unwrap();fs::remove_dir(temporary).unwrap();
    }
    fn power_perm(d:usize,p:Perm,n:usize)->Perm {(0..n).fold(Perm::identity(d),|q,_|q.compose(p,d))}
    #[test] fn state_action_shards_recover_the_exact_small_population() {
        let root=std::env::temp_dir().join(format!("hurwitz-state-actions-{}",std::process::id()));fs::create_dir(&root).unwrap();
        let mut cursor=(0usize,0usize,0usize);let mut combined=HashSet::new();let mut complete=false;
        let base=vec!["binary".into(),"braid-library".into(),"--state-actions".into(),"--d".into(),"4".into(),"--seconds".into(),"10".into()];
        for shard in 0..128 {let out=root.join(format!("shard_{shard}"));let mut args=base.clone();
            args.extend(["--max-tests".into(),"1".into(),"--start-type".into(),cursor.0.to_string(),"--start-second".into(),cursor.1.to_string(),
                "--start-third".into(),cursor.2.to_string(),"--out".into(),out.to_string_lossy().into_owned()]);super::super::run(&args).unwrap();
            for row in fs::read_to_string(out.join("STATE_ACTIONS.tsv")).unwrap().lines().skip(1) {combined.insert(row.split('\t').nth(2).unwrap().to_string());}
            let text=fs::read_to_string(out.join("SUMMARY.txt")).unwrap();let value=|label:&str|text.lines().find_map(|s|s.strip_prefix(label)).unwrap();
            if value("traversal exhausted = ")=="true" {complete=true;break;}
            cursor=(value("next type = ").parse().unwrap(),value("next second = ").parse().unwrap(),value("next third = ").parse().unwrap());
        }
        assert!(complete);let out=root.join("full");let mut args=base;args.extend(["--max-tests".into(),"100000".into(),"--out".into(),out.to_string_lossy().into_owned()]);super::super::run(&args).unwrap();
        let full=fs::read_to_string(out.join("STATE_ACTIONS.tsv")).unwrap().lines().skip(1).map(|s|s.split('\t').nth(2).unwrap().to_string()).collect::<HashSet<_>>();
        assert!(!full.is_empty());assert_eq!(combined,full);fs::remove_dir_all(root).unwrap();
    }
}
