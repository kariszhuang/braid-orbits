//! Independent negative-only audit kernel. No production modules or dependencies.
//! Enumerates every third permutation in lexicographic rank order, checks scope,
//! transitivity/defect, and literally replays each persisted braid witness.
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read, Write},
    path::Path,
    time::Instant,
};
const MAX: usize = 63;
type P = [u8; MAX];
type T = [P; 4];
fn id(n: usize) -> P {
    let mut p = [0; MAX];
    for (i, x) in p[..n].iter_mut().enumerate() {
        *x = i as u8;
    }
    p
}
fn mul(a: &P, b: &P, n: usize) -> P {
    let mut q = [0; MAX];
    for i in 0..n {
        q[i] = a[b[i] as usize];
    }
    q
}
fn inv(a: &P, n: usize) -> P {
    let mut q = [0; MAX];
    for i in 0..n {
        q[a[i] as usize] = i as u8;
    }
    q
}
fn index(a: &P, n: usize) -> usize {
    let mut seen = 0u64;
    let mut cycles = 0;
    for root in 0..n {
        if seen & (1 << root) != 0 {
            continue;
        }
        cycles += 1;
        let mut x = root;
        while seen & (1 << x) == 0 {
            seen |= 1 << x;
            x = a[x] as usize;
        }
    }
    n - cycles
}
fn orbits(g: &[P], n: usize) -> usize {
    // Union of permutation edges; group orbits are exactly its components.
    let mut parent = id(n);
    fn root(p: &mut P, mut x: usize) -> usize {
        while p[x] as usize != x {
            p[x] = p[p[x] as usize];
            x = p[x] as usize;
        }
        x
    }
    let mut count = n;
    for a in g {
        for (x, &y) in a[..n].iter().enumerate() {
            let r = root(&mut parent, x);
            let s = root(&mut parent, y as usize);
            if r != s {
                parent[r] = s as u8;
                count -= 1;
            }
        }
    }
    count
}
fn braid(t: &mut T, m: u8, n: usize) {
    let i = m as usize % 3;
    let (a, b) = (t[i], t[i + 1]);
    if m < 3 {
        t[i] = mul(&mul(&a, &b, n), &inv(&a, n), n);
        t[i + 1] = a;
    } else {
        t[i] = b;
        t[i + 1] = mul(&mul(&inv(&b, n), &a, n), &b, n);
    }
}
fn smooth(t: &T, n: usize) -> bool {
    let k = index(&mul(&t[0], &t[1], n), n);
    k == index(&t[0], n) + index(&t[1], n) || k == index(&t[2], n) + index(&t[3], n)
}
fn parse_p(s: &str, n: usize) -> Result<P, String> {
    let v = s
        .split(',')
        .map(|x| x.parse::<u8>().map_err(|e| e.to_string()))
        .collect::<Result<Vec<_>, _>>()?;
    if v.len() != n {
        return Err("permutation length".into());
    }
    let mut p = [0; MAX];
    let mut seen = 0u64;
    for (i, x) in v.into_iter().enumerate() {
        if x as usize >= n || seen & (1 << x) != 0 {
            return Err("invalid permutation".into());
        }
        seen |= 1 << x;
        p[i] = x;
    }
    Ok(p)
}
fn next(p: &mut P, n: usize) -> bool {
    let Some(i) = (0..n - 1).rev().find(|&i| p[i] < p[i + 1]) else {
        return false;
    };
    let j = (i + 1..n).rev().find(|&j| p[j] > p[i]).unwrap();
    p.swap(i, j);
    p[i + 1..n].reverse();
    true
}
fn table(n: usize) -> Vec<(P, u8)> {
    let mut p = id(n);
    let mut out = Vec::with_capacity((1..=n).product());
    loop {
        out.push((p, index(&p, n) as u8));
        if !next(&mut p, n) {
            break;
        }
    }
    out
}
fn words(path: &Path) -> Result<Vec<Vec<u8>>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = text.lines();
    if lines.next() != Some("token\tstatus\tword_or_component") {
        return Err("dictionary header".into());
    }
    let mut words = Vec::new();
    for line in lines {
        let parts: Vec<_> = line.split('\t').collect();
        if parts.len() != 3 || parts[0].parse::<usize>().ok() != Some(words.len()) {
            return Err("dictionary token numbering".into());
        }
        if parts[1] != "nonsmooth" {
            return Err("unsupported nonnegative status: use independent Python audit".into());
        }
        let word = parts[2]
            .split(',')
            .map(|s| s.parse::<u8>().map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        if word.is_empty() || word.iter().any(|&m| m >= 6) {
            return Err("invalid witness move".into());
        }
        words.push(word);
    }
    Ok(words)
}
fn audit(
    folder: &Path,
    a: P,
    b: P,
    positive: bool,
    limit: usize,
    perms: &[(P, u8)],
    n: usize,
) -> Result<String, String> {
    let began = Instant::now();
    let dictionary = words(&folder.join("WITNESSES.tsv"))?;
    let mut input = BufReader::with_capacity(
        1024 * 1024,
        File::open(folder.join("OUTCOMES.bin")).map_err(|e| e.to_string())?,
    );
    let identity = id(n);
    let p = mul(&a, &b, n);
    let ip = index(&p, n);
    let left = orbits(&[a, b], n);
    if index(&a, n) + index(&b, n) != ip {
        return Err("first pair not geodesic".into());
    }
    let mut count = [0u64; 8];
    let mut omitted = [0u64; 2];
    let length = if limit == 0 {
        perms.len()
    } else {
        limit.min(perms.len())
    };
    for (rank, &(c, ic)) in perms[..length].iter().enumerate() {
        count[0] += 1;
        let d = inv(&mul(&p, &c, n), n);
        let t = [a, b, c, d];
        let ic = ic as usize;
        let jd = index(&d, n);
        if positive {
            if t.contains(&identity) {
                omitted[0] += 1;
                continue;
            }
            if ip + ic + jd <= 2 * n - 2 {
                omitted[1] += 1;
                continue;
            }
        }
        if orbits(&t[..3], n) != 1 {
            continue;
        }
        count[1] += 1;
        let right = orbits(&t[2..], n);
        if (n - ip + 1) != left + right || right != 1 {
            return Err(format!("implication at rank {rank}"));
        }
        if ic + jd == ip {
            count[2] += 1
        }
        if t.contains(&identity) {
            count[4] += 1;
            continue;
        }
        count[3] += 1;
        let mut record = [0u8; 8];
        input
            .read_exact(&mut record)
            .map_err(|e| format!("missing rank {rank}: {e}"))?;
        let rr = u32::from_le_bytes(record[..4].try_into().unwrap()) as usize;
        let token = u32::from_le_bytes(record[4..].try_into().unwrap()) as usize;
        if rr != rank || token >= dictionary.len() {
            return Err(format!("rank/token mismatch {rank}"));
        }
        let mut moved = t;
        for &m in &dictionary[token] {
            braid(&mut moved, m, n)
        }
        if smooth(&moved, n) {
            return Err(format!("false literal word at rank {rank}"));
        }
        count[5] += 1;
    }
    if length == perms.len() {
        let mut byte = [0u8; 1];
        if input.read(&mut byte).map_err(|e| e.to_string())? != 0 {
            return Err("extra outcome bytes".into());
        }
    }
    Ok(format!(
        "{{\"count\":{:?},\"omitted\":{:?},\"complete\":{},\"seconds\":{}}}",
        count,
        omitted,
        length == perms.len(),
        began.elapsed().as_secs_f64()
    ))
}

fn parse_cycles(s: &str, n: usize) -> Result<P, String> {
    let mut p = id(n);
    let mut used = 0u64;
    let mut text = s.trim();
    while !text.is_empty() {
        if !text.starts_with('(') { return Err("invalid cycle syntax".into()); }
        let end = text.find(')').ok_or("unterminated cycle")?;
        let mut first = None;
        let mut last = None;
        for token in text[1..end].split_whitespace() {
            let v: usize = token.parse().map_err(|_| "invalid cycle point")?;
            if v == 0 || v > n || used & (1u64 << (v-1)) != 0 {
                return Err("out-of-range/repeated cycle point".into());
            }
            let x = v-1;
            used |= 1u64 << x;
            if let Some(y) = last { p[y] = x as u8; } else { first = Some(x); }
            last = Some(x);
        }
        if let (Some(x), Some(y)) = (first, last) { p[y] = x as u8; }
        text = text[end+1..].trim();
    }
    Ok(p)
}
fn parse_tuple(s: &str, n: usize) -> Result<T, String> {
    let s = s.trim();
    if !s.starts_with('[') || !s.ends_with(']') { return Err("tuple brackets".into()); }
    let parts: Vec<_> = s[1..s.len()-1].split(',').collect();
    if parts.len() != 4 { return Err("four entries required".into()); }
    Ok([parse_cycles(parts[0],n)?,parse_cycles(parts[1],n)?,
        parse_cycles(parts[2],n)?,parse_cycles(parts[3],n)?])
}
fn json_string(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"), '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"), '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}",c as u32)),
            c => out.push(c),
        }
    }
    out.push('"'); out
}
fn read_row<R: BufRead>(reader: &mut R, line: &mut String) -> Result<bool,String> {
    line.clear();
    if reader.read_line(line).map_err(|e|e.to_string())? == 0 { return Ok(false); }
    while line.ends_with(['\n','\r']) { line.pop(); }
    Ok(true)
}
fn guided(folder: &Path, n: usize) -> Result<String,String> {
    let began = Instant::now();
    let text = std::fs::read_to_string(folder.join("COUNTS.tsv")).map_err(|e|e.to_string())?;
    let mut rows = text.lines();
    if rows.next() != Some("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds") {
        return Err("guided count header".into());
    }
    let fields: Vec<_> = rows.next().ok_or("missing guided counts")?.split('\t').collect();
    if fields.len()!=12 || rows.next().is_some() || fields[0].parse::<usize>().ok()!=Some(n) {
        return Err("guided degree/count shape".into());
    }
    if !["forest","blocks","joint","solved","quotient","primitive","mutate","double_lift","replay","triangle","returns","collision","legacy-types","elliptic","rational","mutation","refactor","full-fiber","partial-fiber","general-fiber","four-quotient","joint-genus3","involution-boundary","intrinsic"].contains(&fields[1]) { return Err("guided lane".into()); }
    let counts: Vec<u64> = fields[2..11].iter().map(|s|s.parse().map_err(|_|"invalid count"))
        .collect::<Result<_,_>>()?;
    let elapsed: f64 = fields[11].parse().map_err(|_|"invalid elapsed")?;
    if !elapsed.is_finite() || elapsed < 0.0 { return Err("invalid elapsed".into()); }
    let constructed=counts[1];
    let mut leads=BufReader::with_capacity(1024*1024,File::open(folder.join("LEADS.tsv")).map_err(|e|e.to_string())?);
    let mut resolutions=BufReader::with_capacity(1024*1024,File::open(folder.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?);
    let mut lead=String::new();let mut row=String::new();
    if !read_row(&mut resolutions,&mut row)? || row!="candidate\tstatus\tword\tfolder\tseed" {
        return Err("guided resolution header".into());
    }
    let mut records=0u64;let mut negative=0u64;let mut closed=0u64;let mut partial=0u64;
    let mut previous=0u64;let mut exceptions=Vec::new();
    while read_row(&mut resolutions,&mut row)? {
        if !read_row(&mut leads,&mut lead)? { return Err("missing durable lead".into()); }
        let f:Vec<_>=row.split('\t').collect();
        let l:Vec<_>=lead.splitn(2,'\t').collect();
        if f.len()!=5 || l.len()!=2 || f[0]!=l[0] || f[4]!=l[1] { return Err("reordered/mismatched durable lead".into()); }
        let candidate=f[0].parse::<u64>().map_err(|_|"invalid candidate")?;
        if candidate<=previous || candidate>constructed { return Err("candidate order/range".into()); }
        previous=candidate;records+=1;
        let mut t=parse_tuple(f[4],n)?;
        let identity=id(n);
        let mut prod=identity;for p in &t {prod=mul(&prod,p,n);}
        let ix=t.map(|p|index(&p,n));
        let sum=ix.iter().sum::<usize>();
        if prod!=identity || t.contains(&identity) || sum<2*n || sum%2!=0 || orbits(&t,n)!=1 || !smooth(&t,n) {
            return Err(format!("out-of-scope seed {candidate}"));
        }
        match f[1] {
            "nonsmooth" => {
                if f[2].is_empty() || !f[3].is_empty() { return Err("negative word/folder".into()); }
                for m in f[2].split(',') {
                    let m=m.parse::<u8>().map_err(|_|"invalid witness move")?;
                    if m>=6 { return Err("invalid witness move".into()); }
                    braid(&mut t,m,n);
                }
                if smooth(&t,n) { return Err(format!("false literal word {candidate}")); }
                negative+=1;
            }
            "closed_smooth" | "timeout" | "capped" => {
                if !f[2].is_empty() || f[3].is_empty() || !f[3].bytes().all(|c|c.is_ascii_alphanumeric()||c==b'_') {
                    return Err("unsafe/invalid exceptional folder or word".into());
                }
                if f[1]=="closed_smooth" {closed+=1;}else{partial+=1;}
                exceptions.push(format!("{{\"candidate\":{},\"status\":{},\"word\":\"\",\"folder\":{},\"seed\":{}}}",
                    json_string(f[0]),json_string(f[1]),json_string(f[3]),json_string(f[4])));
            }
            _ => return Err("unsupported resolution status".into()),
        }
    }
    if read_row(&mut leads,&mut lead)? { return Err("extra unclassified durable lead".into()); }
    if closed!=counts[7] || partial!=counts[8] { return Err("guided status counts mismatch".into()); }
    Ok(format!("{{\"records\":{records},\"nonsmooth\":{negative},\"closed\":{closed},\"unresolved\":{partial},\"exceptions\":[{}],\"seconds\":{}}}",
        exceptions.join(","),began.elapsed().as_secs_f64()))
}

// Independent lexicographic rooted canonicalization, deliberately NOT the
// production numeric packed-key ordering. Every possible point root participates.
type Key = [P;3];
fn component_key(t: &T,n: usize) -> Result<Key,String> {
    let mut best=[[u8::MAX;MAX];3];
    for root in 0..n {
        let mut labels=[u8::MAX;MAX];let mut queue=[0usize;MAX];
        queue[0]=root;labels[root]=0;let mut len=1;let mut head=0;
        while head<len {
            let x=queue[head];head+=1;
            for p in &t[..3] {let y=p[x] as usize;
                if labels[y]==u8::MAX {labels[y]=len as u8;queue[len]=y;len+=1;}
            }
        }
        if len!=n {return Err("component canonicalization requires transitivity".into());}
        let mut key=[[0;MAX];3];
        for j in 0..3 {for x in 0..n {key[j][x]=labels[t[j][queue[x]] as usize];}}
        if key<best {best=key;}
    }
    Ok(best)
}

fn component(folder: &Path,n: usize,export: &Path) -> Result<String,String> {
    let began=Instant::now();
    let text=std::fs::read_to_string(folder.join("SUMMARY.tsv")).map_err(|e|e.to_string())?;
    let mut summary=text.lines();
    if summary.next()!=Some("status\tdegree\tgenus\tstates\tgroup_order\tambient\tfull_group\tword\tseed") {return Err("component summary header".into());}
    let f:Vec<_>=summary.next().ok_or("component summary missing")?.split('\t').collect();
    if f.len()!=9 || summary.next().is_some() || f[1].parse::<usize>().ok()!=Some(n) {return Err("component summary shape".into());}
    if !["closed_smooth","capped","timeout"].contains(&f[0]) {return Err("component status".into());}
    let closed=f[0]=="closed_smooth";let expected=f[3].parse::<usize>().map_err(|_|"component count")?;
    if expected==0 || expected>1_000_000 {return Err("component state bound".into());}
    let seed=parse_tuple(f[8],n)?;
    let file=BufReader::with_capacity(1024*1024,File::open(folder.join("STATES.tsv")).map_err(|e|e.to_string())?);
    let mut lines=file.lines();
    if lines.next().transpose().map_err(|e|e.to_string())?.as_deref()!=Some("state\tsmooth\tdefect\tgenus\tproduct_index\tleft_orbits\tright_orbits\tedges\ttuple") {
        return Err("component state header".into());
    }
    let mut keys=HashMap::with_capacity(expected);let mut tuples=Vec::with_capacity(expected);
    let mut texts=Vec::with_capacity(expected);let mut edges=Vec::with_capacity(expected);let mut order=Vec::with_capacity(expected);
    for line in lines {
        let line=line.map_err(|e|e.to_string())?;let row:Vec<_>=line.split('\t').collect();
        let i=tuples.len();
        if row.len()!=9 || row[0].parse::<usize>().ok()!=Some(i) || i>=expected {return Err("component row order/shape".into());}
        let t=parse_tuple(row[8],n)?;let mut prod=id(n);
        for p in &t {prod=mul(&prod,p,n);}
        if prod!=id(n) || orbits(&t,n)!=1 {return Err("invalid component product/transitivity".into());}
        let key=component_key(&t,n)?;
        if keys.insert(key,i).is_some() {return Err("duplicate conjugacy state".into());}
        let total=t.iter().map(|p|index(p,n)).sum::<usize>() as i64;
        if total%2!=0 {return Err("component genus parity".into());}
        let ip=index(&mul(&t[0],&t[1],n),n);let left=orbits(&t[..2],n);let right=orbits(&t[2..],n);
        let defect=(n-ip) as i64-left as i64-right as i64+1;let genus=1-n as i64+total/2;
        let is_smooth=smooth(&t,n);
        if row[1]!=if is_smooth {"true"}else{"false"} {return Err("smoothness mismatch".into());}
        let vals=[defect,genus,ip as i64,left as i64,right as i64];
        for j in 0..5 {if row[j+2].parse::<i64>().ok()!=Some(vals[j]) {return Err(format!("component statistic mismatch {j}"));}}
        if is_smooth && defect!=0 {return Err("smooth-implies-zero violation".into());}
        if closed && !is_smooth {return Err("false smooth certificate".into());}
        let edge:Vec<_>=row[7].split(',').map(|x|if x=="-" {Ok(None)}else{x.parse::<usize>().map(Some).map_err(|_|"component edge")}).collect::<Result<_,_>>()?;
        if edge.len()!=6 {return Err("six component edges required".into());}
        tuples.push(t);texts.push(row[8].to_string());edges.push(edge);order.push(key);
    }
    if tuples.len()!=expected {return Err("component count mismatch".into());}
    if closed {
        for i in 0..expected {for m in 0..6 {
            let mut t=tuples[i];braid(&mut t,m as u8,n);let key=component_key(&t,n)?;
            let j=*keys.get(&key).ok_or("missing six-move neighbor")?;
            if edges[i][m]!=Some(j) {return Err("edge mismatch".into());}
        }}
        let mut reached=vec![false;expected];reached[0]=true;let mut queue=vec![0usize];let mut head=0;
        while head<queue.len() {let i=queue[head];head+=1;
            for j in &edges[i] {let j=j.ok_or("closed edge missing")?;
                if !reached[j] {reached[j]=true;queue.push(j);}
            }
        }
        if queue.len()!=expected {return Err("disconnected certificate".into());}
        if !keys.contains_key(&component_key(&seed,n)?) {return Err("wrong component root".into());}
    }
    // Export after ALL checks succeed. Caller chooses a fresh temporary path;
    // a failed audit never yields publishable verified-key records.
    let mut out=std::io::BufWriter::new(std::fs::OpenOptions::new().write(true).create_new(true).open(export).map_err(|e|e.to_string())?);
    for (key,text) in order.iter().zip(texts) {
        let key:Vec<Vec<u8>>=key.iter().map(|p|p[..n].to_vec()).collect();
        writeln!(out,"{:?}\t{}",key,text).map_err(|e|e.to_string())?;
    }
    out.flush().map_err(|e|e.to_string())?;
    Ok(format!("{{\"states\":{expected},\"closed\":{closed},\"all_rows_checked\":true,\"all_six_edges_checked\":{closed},\"seconds\":{}}}",began.elapsed().as_secs_f64()))
}
include!("../smooth_primitive_audit.rs");
fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    if args.len()==3 && args[1]=="--primitive-population" {
        println!("{}",primitive_population(Path::new(&args[2]))?);return Ok(());
    }
    if args.len()==4 && args[1]=="--known-keys" {
        let began=Instant::now();let file=BufReader::new(File::open(&args[2]).map_err(|e|e.to_string())?);
        let mut output=std::io::BufWriter::new(std::fs::OpenOptions::new().write(true).create_new(true).open(&args[3]).map_err(|e|e.to_string())?);
        let mut records=0usize;
        for line in file.lines() {
            let line=line.map_err(|e|e.to_string())?;let (degree,text)=line.split_once('\t').ok_or("known degree,tuple")?;
            let n=degree.parse::<usize>().map_err(|_|"known degree")?;
            if !(2..=MAX).contains(&n){return Err("known degree2..63".into());}
            let t=parse_tuple(text,n)?;let mut product=id(n);
            for p in &t{product=mul(&product,p,n);}
            if product!=id(n)||orbits(&t,n)!=1||!smooth(&t,n){return Err("invalid known literal state".into());}
            let key=component_key(&t,n)?;let key:Vec<Vec<u8>>=key.iter().map(|p|p[..n].to_vec()).collect();
            writeln!(output,"{n}\t{key:?}").map_err(|e|e.to_string())?;records+=1;
        }
        output.flush().map_err(|e|e.to_string())?;
        println!("{{\"records\":{records},\"passed\":true,\"seconds\":{}}}",began.elapsed().as_secs_f64());return Ok(());
    }
    if args.len()==5 && args[1]=="--component" {
        let n=args[3].parse::<usize>().map_err(|_|"component degree")?;
        if !(2..=MAX).contains(&n) {return Err("component degree2..63".into());}
        println!("{}",component(Path::new(&args[2]),n,Path::new(&args[4]))?);return Ok(());
    }
    if args.len()==4 && args[1]=="--negative-bank" {
        let n=args[3].parse::<usize>().map_err(|_|"invalid degree")?;
        if !(2..=MAX).contains(&n){return Err("bank degree2..63".into());}
        let file=BufReader::new(File::open(&args[2]).map_err(|e|e.to_string())?);
        let mut records=0;
        for line in file.lines() {
            let line=line.map_err(|e|e.to_string())?;
            let (seed,word)=line.split_once('\t').ok_or("bank tuple,word")?;
            let mut t=parse_tuple(seed,n)?;let mut prod=id(n);
            for p in &t {prod=mul(&prod,p,n);}
            if prod!=id(n)||orbits(&t,n)!=1 {return Err("invalid bank tuple".into());}
            if word.is_empty()||word.split(',').count()>64 {return Err("bank word length".into());}
            for move_id in word.split(',') {
                let m=move_id.parse::<u8>().map_err(|_|"bank move")?;
                if m>=6{return Err("bank move range".into());}braid(&mut t,m,n);
            }
            if smooth(&t,n){return Err("false bank negative".into());}
            records+=1;
            if records>65536{return Err("bank bound".into());}
        }
        println!("{{\"records\":{records},\"passed\":true}}");return Ok(());
    }
    if args.len()==4 && args[1]=="--guided" {
        let n=args[3].parse::<usize>().map_err(|_|"invalid degree")?;
        if !(2..=MAX).contains(&n) {return Err("guided degree2..63".into());}
        println!("{}",guided(Path::new(&args[2]),n)?);return Ok(());
    }
    if args.len() != 3 {
        return Err("usage: smooth_negative_audit --server|--literal-server degree".into());
    }
    let n = args[2].parse::<usize>().map_err(|e| e.to_string())?;
    if !(2..=MAX).contains(&n) || (args[1]=="--server" && n>10) {
        return Err("literal degree2..63; census degree2..10".into());
    }
    let key_server=args[1]=="--key-server";
    let literal = args[1] == "--literal-server" || key_server;
    if !literal && args[1] != "--server" {
        return Err("invalid mode".into());
    }
    let perms = if literal { Vec::new() } else { table(n) };
    println!("READY");
    io::stdout().flush().map_err(|e| e.to_string())?;
    for line in io::stdin().lock().lines() {
        let line = line.map_err(|e| e.to_string())?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5 {
            return Err("expected five tab-separated fields".into());
        }
        let answer = if literal {
            let mut t = [
                parse_p(fields[0], n)?,
                parse_p(fields[1], n)?,
                parse_p(fields[2], n)?,
                parse_p(fields[3], n)?,
            ];
            if !fields[4].is_empty() {
                for m in fields[4].split(',') {
                    let m = m.parse::<u8>().map_err(|e| e.to_string())?;
                    if m >= 6 {
                        return Err("invalid move".into());
                    }
                    braid(&mut t, m, n)
                }
            }
            if key_server {
                let key:Vec<Vec<u8>>=component_key(&t,n)?.iter().map(|p|p[..n].to_vec()).collect();
                println!("{:?}",key);io::stdout().flush().map_err(|e|e.to_string())?;continue;
            }
            let image: Vec<Vec<u8>> = t.iter().map(|p| p[..n].to_vec()).collect();
            let inds: Vec<_> = t.iter().map(|p| index(p, n)).collect();
            let pairs: Vec<_> = (0..3).map(|i| orbits(&t[i..i + 2], n)).collect();
            format!("{{\"tuple\":{:?},\"indices\":{:?},\"pair_orbits\":{:?},\"transitive\":{},\"smooth\":{}}}",image,inds,pairs,orbits(&t,n)==1,smooth(&t,n))
        } else {
            let positive = match fields[3] {
                "0" => false,
                "1" => true,
                _ => return Err("invalid scope flag".into()),
            };
            audit(
                Path::new(fields[0]),
                parse_p(fields[1], n)?,
                parse_p(fields[2], n)?,
                positive,
                fields[4].parse::<usize>().map_err(|e| e.to_string())?,
                &perms,
                n,
            )?
        };
        println!("{answer}");
        io::stdout().flush().map_err(|e| e.to_string())?;
    }
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("INDEPENDENT_AUDIT_FAILED: {error}");
        std::process::exit(1)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lex_population() {
        for n in 2..=8 {
            let ps = table(n);
            assert_eq!(ps.len(), (1..=n).product());
            assert!(ps.windows(2).all(|w| w[0].0 < w[1].0));
        }
    }
    #[test]
    fn inverse_moves() {
        for (a, _) in table(4) {
            for (b, _) in table(4) {
                let t = [a, b, id(4), id(4)];
                for m in 0..3 {
                    let mut q = t;
                    braid(&mut q, m, 4);
                    braid(&mut q, m + 3, 4);
                    assert_eq!(q, t);
                }
            }
        }
    }
    #[test]
    fn independent_component_key_conjugacy_invariance() {
        let ps=table(4);let mut good=0;
        for &(a,_) in &ps {for &(b,_) in &ps {
            let t=[a,b,id(4),inv(&mul(&a,&b,4),4)];
            if orbits(&t,4)!=1 {continue;}good+=1;
            let key=component_key(&t,4).unwrap();
            for &(h,_) in &ps {let hi=inv(&h,4);
                let u=t.map(|p|mul(&mul(&h,&p,4),&hi,4));
                assert_eq!(component_key(&u,4).unwrap(),key);
            }
        }}
        assert!(good>100);
    }
}
