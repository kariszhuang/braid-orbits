//! Exact odd-degree involution/full-cycle passport. Independent literal audit
//! and Frobenius coverage counts live in verify_genus_three_triangle.py.
use std::{fs,io::{BufWriter,Write},path::Path,time::Instant};
type P=Vec<usize>;
fn mul(a:&P,b:&P)->P{b.iter().map(|&x|a[x]).collect()}
fn inv(a:&P)->P{let mut b=vec![0;a.len()];for(i,&x)in a.iter().enumerate(){b[x]=i;}b}
fn full(a:&P)->bool{let mut x=0;for _ in 1..a.len(){x=a[x];if x==0{return false;}}a[x]==0}
fn braid(t:&mut [P;4],m:usize){let i=m%3;let(a,b)=(t[i].clone(),t[i+1].clone());if m<3{t[i]=mul(&mul(&a,&b),&inv(&a));t[i+1]=a;}else{t[i]=b.clone();t[i+1]=mul(&mul(&inv(&b),&a),&b);}}
fn index(a:&P)->usize{let mut seen=vec![false;a.len()];let mut c=0;for r in 0..a.len(){if !seen[r]{c+=1;let mut x=r;while !seen[x]{seen[x]=true;x=a[x];}}}a.len()-c}
fn smooth(t:&[P;4])->bool{let k=index(&mul(&t[0],&t[1]));k==index(&t[0])+index(&t[1])||k==index(&t[2])+index(&t[3])}
fn text(p:&P)->String{let mut seen=vec![false;p.len()];let mut out=String::new();for root in 0..p.len(){if seen[root]||p[root]==root{continue;}let mut c=vec![];let mut x=root;while !seen[x]{seen[x]=true;c.push((x+1).to_string());x=p[x];}out.push('(');out.push_str(&c.join(" "));out.push(')');}if out.is_empty(){out.push_str("()");}out}
struct Run{n:usize,a:P,b:P,h:P,rank:u64,qualified:u64,maximum:usize,leads:BufWriter<fs::File>,proofs:BufWriter<fs::File>}
impl Run{
    fn leaf(&mut self,c:&P,fixed:usize){
        self.rank+=1;let d=inv(&mul(&self.h,c));if !full(&d){return;}self.qualified+=1;
        let t=[self.a.clone(),self.b.clone(),c.clone(),d];assert!(smooth(&t));
        // a fixes0 and h is +1. Positive/negative outer twists move that fixed
        // point to c's fixed point. The shorter direction gives length<=n.
        let mut word=if fixed<=self.n-fixed{vec![0;2*fixed]}else{vec![3;2*(self.n-fixed)]};word.push(4);
        let mut u=t.clone();for &m in &word{braid(&mut u,m);}assert!(!smooth(&u));self.maximum=self.maximum.max(word.len());
        let seed=format!("[{}]",t.iter().map(text).collect::<Vec<_>>().join(", "));
        writeln!(self.leads,"{}\t{}",self.rank,seed).unwrap();
        writeln!(self.proofs,"{}\tnonsmooth\t{}\t\t{}",self.rank,word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(","),seed).unwrap();
    }
    fn enumerate(&mut self,c:&mut P,unused:u64,fixed:usize){
        if unused==0{self.leaf(c,fixed);return;}
        let x=unused.trailing_zeros()as usize;let rest=unused&!(1<<x);let mut options=rest;
        while options!=0{let y=options.trailing_zeros()as usize;options&=options-1;c[x]=y;c[y]=x;self.enumerate(c,rest&!(1<<y),fixed);}
    }
}
fn main(){let args=std::env::args().collect::<Vec<_>>();assert_eq!(args.len(),3);let n=args[1].parse::<usize>().unwrap();assert!([5,9,13].contains(&n));let out=Path::new(&args[2]);fs::create_dir(out).unwrap();let began=Instant::now();
    let a=(0..n).map(|i|(n-i)%n).collect::<P>();let b=(0..n).map(|i|(n-i-1)%n).collect::<P>();let h=(0..n).map(|i|(i+1)%n).collect::<P>();assert_eq!(mul(&a,&b),h);
    let mut run=Run{n,a,b,h,rank:0,qualified:0,maximum:0,leads:BufWriter::new(fs::File::create(out.join("LEADS.tsv")).unwrap()),proofs:BufWriter::new(fs::File::create(out.join("RESOLUTIONS.tsv")).unwrap())};
    writeln!(run.proofs,"candidate\tstatus\tword\tfolder\tseed").unwrap();
    for fixed in 0..n{let mut c=(0..n).collect::<P>();run.enumerate(&mut c,((1u64<<n)-1)&!(1<<fixed),fixed);}
    run.leads.flush().unwrap();run.proofs.flush().unwrap();
    fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\ttriangle\t0\t{}\t{}\t{}\t{}\t0\t0\t0\t0\t{}\n",run.rank,run.qualified,run.qualified,run.qualified,began.elapsed().as_secs_f64())).unwrap();
    fs::write(out.join("COVERAGE.json"),format!("{{\"degree\":{n},\"genus\":{},\"involutions\":{},\"necessary_Sn_states\":{},\"negative\":{},\"unresolved\":0,\"maximum_word\":{},\"quotient\":\"simultaneous Sn conjugation\"}}\n",(n-1)/4,run.rank,run.qualified,run.qualified,run.maximum)).unwrap();
}
