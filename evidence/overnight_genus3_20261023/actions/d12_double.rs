//! Exact necessary S12 double-transposition slice. Undecided starts are retained.
use std::{collections::{BTreeMap,BTreeSet},fs,io::Write,path::Path};
const N:usize=12;
type P=[u8;N];
fn id()->P{std::array::from_fn(|i|i as u8)}
fn mul(a:P,b:P)->P{std::array::from_fn(|i|a[b[i]as usize])}
fn inv(a:P)->P{let mut z=[0;N];for i in 0..N{z[a[i]as usize]=i as u8;}z}
fn cycles(a:P)->Vec<Vec<usize>>{let(mut seen,mut out)=(0u16,vec![]);for mut i in 0..N{if seen&(1<<i)!=0{continue;}let mut c=vec![];while seen&(1<<i)==0{seen|=1<<i;c.push(i);i=a[i]as usize;}out.push(c);}out}
fn ix(a:P)->usize{N-cycles(a).len()}
fn canon(edges:[(usize,usize);2])->[(usize,usize);2]{
 (0..N).map(|k|{let mut es=edges.map(|(x,y)|{let(a,b)=((x+k)%N,(y+k)%N);(a.min(b),a.max(b))});es.sort();es}).min().unwrap()
}
fn permute(a:&mut[usize],at:usize,out:&mut Vec<Vec<usize>>){if at==a.len(){out.push(a.to_vec());return;}for j in at..a.len(){a.swap(at,j);permute(a,at+1,out);a.swap(at,j);}}
fn orientations(block:&[usize])->Vec<Vec<usize>>{let mut tail=block[1..].to_vec();let mut out=vec![];permute(&mut tail,0,&mut out);out.into_iter().map(|v|{let mut c=vec![block[0]];c.extend(v);c}).collect()}
fn partitions(rem:u16,parts:&mut Vec<Vec<usize>>,out:&mut Vec<Vec<Vec<usize>>>){
 if rem==0{if parts.len()==3{out.push(parts.clone());}return;}
 let left=3-parts.len();let size=rem.count_ones()as usize;if size<3*left||size>6*left{return;}
 let anchor=rem.trailing_zeros()as usize;let rest=rem&!(1<<anchor);let mut sub=rest;
 loop{let m=sub|(1<<anchor);let len=m.count_ones()as usize;if (3..=6).contains(&len){parts.push((0..N).filter(|i|m&(1<<i)!=0).collect());partitions(rem^m,parts,out);parts.pop();}if sub==0{break;}sub=(sub-1)&rest;}
}
fn tree_partition(blocks:&[Vec<usize>],edges:[(usize,usize);2])->bool{
 let mut label=[0;N];for (j,c)in blocks.iter().enumerate(){for &x in c{label[x]=j;}}
 for k in 0..N{let mut pairs=[(0,0);2];for (j,&(x,y))in edges.iter().enumerate(){let(a,b)=(label[(x+k)%N],label[(y+k)%N]);if a==b{return false;}pairs[j]=(a.min(b),a.max(b));}if pairs[0]==pairs[1]{return false;}}
 true
}
fn step(t:[P;4],code:usize)->[P;4]{let mut u=t;let j=code/2;let(a,b)=(t[j],t[j+1]);if code%2==0{u[j]=mul(mul(a,b),inv(a));u[j+1]=a;}else{u[j]=b;u[j+1]=mul(mul(inv(b),a),b);}u}
fn smooth(t:[P;4])->bool{let r=ix(mul(t[0],t[1]));r==ix(t[0])+ix(t[1])||r==ix(t[2])+ix(t[3])}
fn main(){
 let args:Vec<_>=std::env::args().collect();let out=Path::new(args.get(1).expect("new output directory"));assert!(!out.exists());fs::create_dir(out).unwrap();
 let p:P=std::array::from_fn(|i|((i+1)%N)as u8);let pi=inv(p);let mut frames=BTreeSet::new();let(mut all,mut full,mut bounded)=(0,0,0);
 for w in 0..N{for x in w+1..N{for y in x+1..N{for z in y+1..N{for edges in [[(w,x),(y,z)],[(w,y),(x,z)],[(w,z),(x,y)]]{
  all+=1;let mut a=id();for(u,v)in edges{a[u]=v as u8;a[v]=u as u8;}let cs=cycles(mul(a,p));if cs.len()!=3{continue;}full+=1;if cs.iter().any(|c|!(3..=6).contains(&c.len())){continue;}bounded+=1;frames.insert(canon(edges));
 }}}}}
 let mut parts=vec![];partitions((1<<N)-1,&mut vec![],&mut parts);let mut qualified=fs::File::create(out.join("UNRESOLVED.jsonl")).unwrap();let mut enumeration=fs::File::create(out.join("ENUMERATION.bin")).unwrap();
 let(mut supports,mut generated,mut wrong,mut negative,mut pending)=(0u64,0u64,0u64,0u64,0u64);let mut frame_rows=vec![];let mut word_hist=BTreeMap::new();
 for edges in &frames{
  let mut a=id();for &(u,v)in edges{a[u]=v as u8;a[v]=u as u8;}let c=mul(a,p);let(mut fparts,mut fgen,mut findex,mut fsurv)=(0,0,0,0);
  let twists:Vec<P>=(0..N).map(|k|{let mut r=id();for &(u,v)in edges{r[(u+k)%N]=((v+k)%N)as u8;r[(v+k)%N]=((u+k)%N)as u8;}r}).collect();
  for blocks in &parts{if !tree_partition(blocks,*edges){continue;}supports+=1;fparts+=1;let orders:Vec<_>=blocks.iter().map(|c|orientations(c)).collect();
   for x in &orders[0]{for y in &orders[1]{for z in &orders[2]{
    let mut e=id();for cy in [x,y,z]{for j in 0..cy.len(){e[cy[j]]=cy[(j+1)%cy.len()]as u8;}}
    enumeration.write_all(&a).unwrap();enumeration.write_all(&e).unwrap();generated+=1;fgen+=1;let b=mul(pi,inv(e));if ix(b)!=8{wrong+=1;continue;}findex+=1;
    let t=[a,c,b,e];assert_eq!(mul(mul(mul(a,c),b),e),id());assert!(smooth(t));
    if let Some(k)=(0..N).find(|&k|ix(mul(twists[k],b))!=10){
     let mut u=t;for _ in 0..2*k{u=step(u,0);}u=step(u,3);assert!(!smooth(u));negative+=1;*word_hist.entry(2*k+1).or_insert(0u64)+=1;
    }else{pending+=1;fsurv+=1;writeln!(qualified,"{{\"tuple\":{:?},\"frame\":{:?},\"status\":\"unresolved\"}}",t,edges.map(|(x,y)|[x,y])).unwrap();}
   }}}
  }
  frame_rows.push(format!("{{\"edges\":{:?},\"supports\":{},\"generated\":{},\"index8\":{},\"unresolved\":{}}}",edges.map(|(x,y)|[x,y]),fparts,fgen,findex,fsurv));
 }
 assert_eq!(all,1485);assert_eq!(full,990);assert_eq!(frames.len(),16);assert_eq!(generated,wrong+negative+pending);
 let report=format!("{{\"a_class_labelled\":{},\"full_product_frames_labelled\":{},\"bounded_frames_labelled\":{},\"rotation_classes\":{},\"support_partitions\":{},\"cyclic_tree_supports\":{},\"e_generated\":{},\"wrong_b_index\":{},\"literal_negative\":{},\"unresolved\":{},\"frames\":[{}]}}\n",all,full,bounded,frames.len(),parts.len(),supports,generated,wrong,negative,pending,frame_rows.join(","));
 fs::write(out.join("RESULT.json"),&report).unwrap();println!("{}",report);eprintln!("negative word length histogram: {:?}",word_hist);
}
