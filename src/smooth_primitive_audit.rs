// Independent finite affine-action and omission audit. Uses only this binary's
// literal kernels, never the production group/frame/short-word implementation.
fn pa_form(x:usize,y:usize)->usize {
    ((x&1)*((y>>1)&1)+((x>>1)&1)*(y&1)+((x>>2)&1)*((y>>3)&1)+((x>>3)&1)*((y>>2)&1))&1
}
fn pa_conjugate(a:&P,h:&P,hi:&P,n:usize)->P {mul(&mul(h,a,n),hi,n)}
fn pa_cycles(s:&str,n:usize)->Result<P,String>{Ok(parse_tuple(&format!("[{s}, (), (), ()]"),n)?[0])}
fn pa_affine(p:&P,symplectic:bool)->bool {
    let t=p[0];let basis:[u8;4]=std::array::from_fn(|i|p[1<<i]^t);
    for x in 0..16 {let mut value=t;for i in 0..4{if x&(1<<i)!=0{value^=basis[i];}}if value!=p[x]{return false;}}
    if symplectic {for i in 0..4 {for j in 0..4{if pa_form(basis[i] as usize,basis[j] as usize)!=pa_form(1<<i,1<<j){return false;}}}}
    true
}
fn pa_words()->Vec<Vec<u8>> {
    let mut words=vec![vec![1],vec![4],vec![0,1],vec![0,4],vec![3,1],vec![1,1],vec![4,4],vec![1,1,1],vec![4,4,4]];
    for end in [0,3] {for length in 2..=4 {for middle in [4,1]{let mut w=vec![end;length];w.push(middle);words.push(w);}}}
    words.extend([vec![1,3,1],vec![0,4,4],vec![4,0,4],vec![3,1,1],vec![5,1,1,1],vec![5,3,1,1,1]]);words
}
fn primitive_population(folder:&Path)->Result<String,String> {
    let began=Instant::now();let n=16;
    let raw=std::fs::read_to_string(folder.join("PRIMITIVE.tsv")).map_err(|e|e.to_string())?;
    let lines:Vec<_>=raw.lines().collect();if lines.len()!=2{return Err("finite header shape".into());}
    let f:Vec<_>=lines[1].split('\t').collect();
    if f.len()!=8||!["asp16","agl16"].contains(&f[0])||f[1]!="16"{return Err("affine action declaration".into());}
    let symplectic=f[0]=="asp16";let expected=if symplectic{11520}else{322560};
    if f[2].parse::<usize>().ok()!=Some(expected){return Err("affine action order declaration".into());}
    let raw=std::fs::read_to_string(folder.join("GROUP.tsv")).map_err(|e|e.to_string())?;let mut group=vec![];
    for (i,line) in raw.lines().enumerate(){let (rank,s)=line.split_once('\t').ok_or("group shape")?;
        if rank.parse::<usize>().map_err(|_|"group rank")?!=i{return Err("group rank order".into());}group.push(pa_cycles(s,n)?);}
    if group.len()!=expected||group.windows(2).any(|w|w[0]>=w[1])||group.iter().any(|p|!pa_affine(p,symplectic)) {
        return Err("not the complete declared affine action".into());
    }
    // A symplectic basis has15*8*3*2=720 choices; all16 translations give11520.
    // Thus uniqueness, size and the defining equations give exact coverage.
    let mut orbit=[false;16];for p in &group{if p[0]==0{orbit[p[1] as usize]=true;}}
    if orbit.iter().filter(|&&v|v).count()!=15{return Err("point stabilizer not transitive; no primitive action proof".into());}
    let ids:HashMap<_,_>=group.iter().enumerate().map(|(i,&p)|(p,i)).collect();
    let inverse:Vec<_>=group.iter().map(|p|inv(p,n)).collect();let indices:Vec<_>=group.iter().map(|p|index(p,n)).collect();
    let mut seen=vec![false;group.len()];let mut frames=vec![];
    for (i,a) in group.iter().enumerate(){
        if seen[i]{continue;}
        for (j,h) in group.iter().enumerate(){let q=pa_conjugate(a,h,&inverse[j],n);seen[*ids.get(&q).ok_or("action not conjugation closed")?]=true;}
        if indices[i]==0||indices[i]>=n-1{continue;}
        let center:Vec<_>=group.iter().enumerate().filter(|(_,h)|mul(a,h,n)==mul(h,a,n)).map(|(j,_)|j).collect();
        let mut second=vec![false;group.len()];
        for (j,b) in group.iter().enumerate(){
            if second[j]||indices[j]==0||indices[i]+indices[j]>=n||index(&mul(a,b,n),n)!=indices[i]+indices[j]{continue;}
            for &l in &center{second[*ids.get(&pa_conjugate(b,&group[l],&inverse[l],n)).ok_or("centralizer action")?]=true;}
            frames.push([*a,*b]);
        }
    }
    let raw=std::fs::read_to_string(folder.join("PAIR_FRAMES.tsv")).map_err(|e|e.to_string())?;let mut actual=vec![];
    for (i,line) in raw.lines().enumerate(){let fields:Vec<_>=line.split('\t').collect();if fields.len()!=3||fields[0].parse::<usize>().ok()!=Some(i){return Err("pair frame shape".into());}
        actual.push([pa_cycles(fields[1],n)?,pa_cycles(fields[2],n)?]);}
    if actual!=frames||f[3].parse::<usize>().ok()!=Some(frames.len()){return Err("independent conjugacy/centralizer frame omission".into());}
    let start=f[4].parse::<usize>().map_err(|_|"start")?;let end=f[5].parse::<usize>().map_err(|_|"end")?;
    let total=frames.len()*group.len();if f[6].parse::<usize>().ok()!=Some(total)||start>end||end>total{return Err("finite interval".into());}
    let raw=std::fs::read_to_string(folder.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?;
    let mut journal=HashMap::new();let mut unresolved=0;
    for line in raw.lines().skip(1){let fields:Vec<_>=line.split('\t').collect();if fields.len()!=5{return Err("journal shape".into());}
        let t=parse_tuple(fields[4],n)?;let k=component_key(&t,n)?;
        let word=if fields[2].is_empty(){vec![]}else{fields[2].split(',').map(|s|s.parse::<u8>().map_err(|_|"witness move")).collect::<Result<Vec<_>,_>>()?};
        if ["capped","timeout"].contains(&fields[1]){unresolved+=1;}
        journal.insert(k,(fields[1].to_string(),word));
    }
    let mut known=std::collections::HashSet::new();
    if let Ok(raw)=std::fs::read_to_string(folder.join("KNOWN.tsv")){for line in raw.lines(){let (_,s)=line.split_once('\t').ok_or("known alias")?;known.insert(component_key(&parse_tuple(s,n)?,n)?);}}
    let words=pa_words();let mut qualifying=0;let mut negative=0;let mut accounted=0;
    for rank in start..end {
        let [a,b]=frames[rank/group.len()];let c=group[rank%group.len()];let d=inv(&mul(&mul(&a,&b,n),&c,n),n);let t=[a,b,c,d];
        let sum=t.iter().map(|p|index(p,n)).sum::<usize>();
        if t.contains(&id(n))||sum<2*n||sum%2!=0||orbits(&t,n)!=1{continue;}
        qualifying+=1;if !smooth(&t,n){return Err("constructed pair not additive".into());}
        let mut rejected=false;
        for word in &words{let mut u=t;for &m in word{braid(&mut u,m,n);}if !smooth(&u,n){rejected=true;break;}}
        if rejected{negative+=1;continue;}
        let k=component_key(&t,n)?;
        if let Some((status,word))=journal.get(&k){
            if status=="nonsmooth"{let mut u=t;for &m in word{if m>=6{return Err("bad move".into());}braid(&mut u,m,n);}if smooth(&u,n){return Err("false journal alias negative".into());}negative+=1;}
            else{accounted+=1;}
        }else if known.contains(&k){accounted+=1;}else{return Err("unaccounted eligible omitted seed".into());}
    }
    Ok(format!("{{\"passed\":true,\"group\":\"{}\",\"degree\":16,\"group_order\":{expected},\"ambient_action_primitive\":true,\"frames\":{},\"start\":{start},\"next\":{end},\"total\":{total},\"independent_complete_population\":{},\"counts\":{{\"qualifying\":{qualifying},\"literal_negative\":{negative},\"accounted_survivor\":{accounted}}},\"no_unresolved\":{},\"definition_checked\":\"affine equations and exact cardinality; point stabilizer transitive\",\"seconds\":{}}}",f[0],frames.len(),start==0&&end==total,unresolved==0,began.elapsed().as_secs_f64()))
}
