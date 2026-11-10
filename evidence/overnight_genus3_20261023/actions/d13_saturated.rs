//! Necessary saturated S13 genus-three slices; no unclosed start is a negative.
use std::{collections::BTreeSet,fs,io::Write,path::Path};
const N:usize=13;type P=[u8;N];
fn id()->P{std::array::from_fn(|i|i as u8)}
fn mul(a:P,b:P)->P{std::array::from_fn(|i|a[b[i]as usize])}
fn inv(a:P)->P{let mut r=[0;N];for i in 0..N{r[a[i]as usize]=i as u8;}r}
fn cycles(a:P)->Vec<Vec<usize>>{let(mut seen,mut out)=(0u16,vec![]);for mut i in 0..N{if seen&(1<<i)!=0{continue;}let mut c=vec![];while seen&(1<<i)==0{seen|=1<<i;c.push(i);i=a[i]as usize;}out.push(c);}out}
fn ix(a:P)->usize{N-cycles(a).len()}
fn canonical(a:P)->P{(0..N).map(|k|std::array::from_fn(|i|((a[(i+N-k)%N]as usize+k)%N)as u8)).min().unwrap()}
fn setcycle(a:&mut P,c:&[usize]){for i in 0..c.len(){a[c[i]]=c[(i+1)%c.len()]as u8;}}
fn matchings(points:&[usize],left:usize,edges:&mut Vec<(usize,usize)>,out:&mut Vec<P>){
 if left==0{let mut a=id();for &(x,y)in edges.iter(){a[x]=y as u8;a[y]=x as u8;}out.push(a);return;}
 for i in 0..points.len(){if points.len()-i<2*left{break;}for j in i+1..points.len(){let rest:Vec<_>=points[i+1..].iter().copied().filter(|&x|x!=points[j]).collect();edges.push((points[i],points[j]));matchings(&rest,left-1,edges,out);edges.pop();}}
}
fn permute(a:&mut[usize],at:usize,out:&mut Vec<Vec<usize>>){if at==a.len(){out.push(a.to_vec());return;}for j in at..a.len(){a.swap(at,j);permute(a,at+1,out);a.swap(at,j);}}
fn orders(c:&[usize])->Vec<Vec<usize>>{let mut z=c[1..].to_vec();let mut out=vec![];permute(&mut z,0,&mut out);out.into_iter().map(|v|{let mut r=vec![c[0]];r.extend(v);r}).collect()}
fn partitions(rem:u16,lens:&mut Vec<usize>,cs:&mut Vec<Vec<usize>>,out:&mut Vec<Vec<Vec<usize>>>){
 if rem==0{if lens.is_empty(){out.push(cs.clone());}return;}if lens.is_empty(){return;}let anchor=rem.trailing_zeros()as usize;let rest=rem&!(1<<anchor);let mut sub=rest;let lengths:BTreeSet<_>=lens.iter().copied().collect();
 loop{let len=(sub.count_ones()+1)as usize;if lengths.contains(&len){let pos=lens.iter().position(|&x|x==len).unwrap();lens.remove(pos);let mask=sub|(1<<anchor);cs.push((0..N).filter(|i|mask&(1<<i)!=0).collect());partitions(rem^mask,lens,cs,out);cs.pop();lens.insert(pos,len);}if sub==0{break;}sub=(sub-1)&rest;}
}
fn forest(labels:&[usize;N],a_cycles:&[Vec<usize>],k:usize)->bool{
 let mut par=[0,1,2,3];
 fn root(par:&[usize;4],mut x:usize)->usize{while par[x]!=x{x=par[x];}x}
 for c in a_cycles{if c.len()==1{continue;}let x=labels[(c[0]+k)%N];for &point in &c[1..]{let y=labels[(point+k)%N];let(rx,ry)=(root(&par,x),root(&par,y));if rx==ry{return false;}par[rx]=ry;}}true
}
fn emit_orders(j:usize,os:&[Vec<Vec<usize>>],e:P,visit:&mut impl FnMut(P)){if j==os.len(){visit(e);return;}for c in &os[j]{let mut next=e;setcycle(&mut next,c);emit_orders(j+1,os,next,visit);}}
fn step(t:[P;4],code:usize)->[P;4]{let mut u=t;let j=code/2;let(a,b)=(t[j],t[j+1]);if code%2==0{u[j]=mul(mul(a,b),inv(a));u[j+1]=a;}else{u[j]=b;u[j+1]=mul(mul(inv(b),a),b);}u}
fn smooth(t:[P;4])->bool{let r=ix(mul(t[0],t[1]));r==ix(t[0])+ix(t[1])||r==ix(t[2])+ix(t[3])}
fn main(){
 let args:Vec<_>=std::env::args().collect();let kind=args.get(1).expect("double|three-two|triple");let out=Path::new(args.get(2).expect("new output directory"));assert!(!out.exists());fs::create_dir(out).unwrap();
 let s=if kind=="double"{2}else{3};let profile:Vec<usize>=if s==2{vec![4,4,5]}else{vec![3,3,3,4]};let p:P=std::array::from_fn(|i|((i+1)%N)as u8);let mut all=vec![];
 if kind=="double"||kind=="triple"{matchings(&(0..N).collect::<Vec<_>>(),s,&mut vec![],&mut all);}else{assert_eq!(kind,"three-two");for x in 0..N{for y in x+1..N{for z in y+1..N{for order in [[x,y,z],[x,z,y]]{for u in 0..N{if order.contains(&u){continue;}for v in u+1..N{if order.contains(&v){continue;}let mut a=id();setcycle(&mut a,&order);setcycle(&mut a,&[u,v]);all.push(a);}}}}}}}
 let all_count=all.len();let mut frames=BTreeSet::new();let(mut full,mut bounded)=(0,0);for a in all{let c=mul(inv(a),p);let cc=cycles(c);if cc.len()!=s+1{continue;}full+=1;let mut sizes:Vec<_>=cc.iter().map(Vec::len).collect();sizes.sort();if sizes!=profile{continue;}bounded+=1;frames.insert(canonical(a));}
 let mut supports=vec![];partitions((1<<N)-1,&mut profile.clone(),&mut vec![],&mut supports);let(mut accepted,mut generated,mut wrong,mut negative,mut pending)=(0u64,0u64,0u64,0u64,0u64);let mut records=fs::File::create(out.join("UNRESOLVED.jsonl")).unwrap();let mut binary=fs::File::create(out.join("ENUMERATION.bin")).unwrap();let mut frame_rows=vec![];
 for &a in &frames{let c=mul(inv(a),p);let ac=cycles(a);let twists:Vec<P>=(0..N).map(|k|std::array::from_fn(|i|((a[(i+N-k)%N]as usize+k)%N)as u8)).collect();let(mut sup,mut gen,mut pass,mut pend)=(0u64,0u64,0u64,0u64);
  for blocks in &supports{let mut labels=[0;N];for (j,c)in blocks.iter().enumerate(){for &x in c{labels[x]=j;}}if !(0..N).all(|k|forest(&labels,&ac,k)){continue;}accepted+=1;sup+=1;let os:Vec<_>=blocks.iter().map(|c|orders(c)).collect();
   emit_orders(0,&os,id(),&mut |e|{generated+=1;gen+=1;binary.write_all(&a).unwrap();binary.write_all(&e).unwrap();let b=mul(inv(p),inv(e));if ix(b)!=6+s{wrong+=1;return;}pass+=1;let t=[a,c,b,e];assert_eq!(mul(mul(mul(a,c),b),e),id());assert!(smooth(t));
    if let Some(k)=(0..N).find(|&k|ix(mul(twists[k],b))!=6+2*s){let mut u=t;for _ in 0..2*k{u=step(u,0);}u=step(u,3);assert!(!smooth(u));negative+=1;}else{pending+=1;pend+=1;writeln!(records,"{{\"tuple\":{:?},\"status\":\"unresolved\"}}",t).unwrap();}
   });
  }frame_rows.push(format!("{{\"a\":{:?},\"supports\":{},\"e\":{},\"index_pass\":{},\"unresolved\":{}}}",a,sup,gen,pass,pend));
 }
 assert_eq!(generated,wrong+negative+pending);let report=format!("{{\"kind\":\"{}\",\"all_a\":{},\"full_frames_labelled\":{},\"bounded_frames_labelled\":{},\"rotation_frames\":{},\"support_partitions\":{},\"accepted_supports\":{},\"generated_e\":{},\"wrong_b_index\":{},\"literal_negative\":{},\"unresolved\":{},\"frames\":[{}]}}\n",kind,all_count,full,bounded,frames.len(),supports.len(),accepted,generated,wrong,negative,pending,frame_rows.join(","));fs::write(out.join("RESULT.json"),&report).unwrap();println!("{}",report);
}
