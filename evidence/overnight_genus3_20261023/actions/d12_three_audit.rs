//! Independent full 4^3 class support enumeration and literal twist audit.
use std::{fs,io::Write,path::Path};
const N:usize=12;
type P=[u8;N];
type T=[P;4];
const ORD:[[usize;3];6]=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];
const MAP:[[usize;3];6]=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
fn mul(a:P,b:P)->P{std::array::from_fn(|i|a[b[i]as usize])}
fn inv(a:P)->P{let mut z=[0;N];for i in 0..N{z[a[i]as usize]=i as u8;}z}
fn idx(a:P)->usize{let(mut seen,mut c)=(0u16,0);for mut i in 0..N{if seen&(1<<i)!=0{continue;}c+=1;while seen&(1<<i)==0{seen|=1<<i;i=a[i]as usize;}}N-c}
fn rainbow_block(b:&[usize])->bool{b.iter().fold(0u8,|m,&i|m|(1<<(i%4))).count_ones()==4}
fn bad_color(b:P)->Option<usize>{let mut seen=0u16;for mut i in 0..N{if seen&(1<<i)!=0{continue;}let mut colors=0u8;while seen&(1<<i)==0{let color=i%4;if colors&(1<<color)!=0{return Some(color);}colors|=1<<color;seen|=1<<i;i=b[i]as usize;}}None}
fn step(t:T,code:usize)->T{let mut u=t;let j=code/2;let(a,b)=(t[j],t[j+1]);if code%2==0{u[j]=mul(mul(a,b),inv(a));u[j+1]=a;}else{u[j]=b;u[j+1]=mul(mul(inv(b),a),b);}u}
fn ordinal(e:P)->usize{
 let mut maps=[[0usize;3];3];let mut ors=[0usize;3];
 for i in 0..3{let anchor=4*i;let mut x=e[anchor]as usize;let mut colors=[0;3];
  for j in 0..3{let color=x%4;assert!(color>0);maps[color-1][i]=x/4;colors[j]=color;x=e[x]as usize;}
  assert_eq!(x,anchor);ors[i]=ORD.iter().position(|&v|v==colors).unwrap();
 }
 let outer=maps.iter().fold(0,|v,m|v*6+MAP.iter().position(|x|x==m).unwrap());
 outer*216+ors.iter().fold(0,|v,o|v*6+o)
}
fn main(){
 let args:Vec<_>=std::env::args().collect();let out=Path::new(args.get(1).expect("output directory"));assert!(!out.exists());fs::create_dir(out).unwrap();
 let id:P=std::array::from_fn(|i|i as u8);let p:P=std::array::from_fn(|i|((i+1)%N)as u8);
 let mut a=id;a[0]=4;a[4]=8;a[8]=0;let c=mul(inv(a),p);
 let mut enumeration=vec![[0u8;N];46656];let mut seen_ord=vec![false;46656];
 let(mut layouts,mut rainbow_layouts,mut e_count,mut wrong,mut witnesses)=(0usize,0usize,0usize,0usize,0usize);
 let mut word_counts=[0usize;4];
 for i in 1..N{for j in i+1..N{for k in j+1..N{
  let b0=[0,i,j,k];let rem:Vec<_>=(0..N).filter(|x|!b0.contains(x)).collect();let h=rem[0];
  for v in 1..rem.len(){for w in v+1..rem.len(){for z in w+1..rem.len(){
   let b1=[h,rem[v],rem[w],rem[z]];let b2:Vec<_>=rem.iter().copied().filter(|x|!b1.contains(x)).collect();assert_eq!(b2.len(),4);
   layouts+=1;if !rainbow_block(&b0)||!rainbow_block(&b1)||!rainbow_block(&b2){continue;}rainbow_layouts+=1;
   for o0 in ORD{for o1 in ORD{for o2 in ORD{
    let mut e=id;
    for (block,order) in [(&b0[..],o0),(&b1[..],o1),(&b2[..],o2)]{
     let cycle=[block[0],block[order[0]],block[order[1]],block[order[2]]];
     for r in 0..4{e[cycle[r]]=cycle[(r+1)%4]as u8;}
    }
    let order=ordinal(e);assert!(!seen_ord[order]);seen_ord[order]=true;enumeration[order]=e;e_count+=1;
    let b=mul(inv(p),inv(e));if idx(b)!=8{wrong+=1;continue;}
    let r=bad_color(b).expect("surviving degree12 candidate: preserve and investigate");
    let t=[a,c,b,e];assert_eq!(mul(mul(mul(a,c),b),e),id);
    let mut endpoint=t;for _ in 0..2*r{endpoint=step(endpoint,0);}endpoint=step(endpoint,3);
    let product_index=idx(mul(endpoint[0],endpoint[1]));
    assert_ne!(product_index,idx(endpoint[0])+idx(endpoint[1]));
    assert_ne!(product_index,idx(endpoint[2])+idx(endpoint[3]));witnesses+=1;word_counts[r]+=1;
   }}}
  }}}
 }}}
 assert_eq!(layouts,5775);assert_eq!(layouts*216,1247400);assert_eq!(rainbow_layouts,216);assert_eq!(e_count,46656);assert!(seen_ord.iter().all(|&x|x));assert_eq!(wrong+witnesses,e_count);
 let mut raw=fs::File::create(out.join("ENUMERATION.bin")).unwrap();for e in enumeration{raw.write_all(&e).unwrap();}
 let report=format!("{{\"full_cycle_support_partitions\":{},\"full_labelled_4_cubed_class\":{},\"rainbow_support_partitions\":{},\"necessary_e\":{},\"wrong_b_index\":{},\"literal_negative_candidates\":{},\"unresolved\":0,\"word_counts_by_color\":{:?},\"maximum_negative_word_length\":7}}\n",layouts,layouts*216,rainbow_layouts,e_count,wrong,witnesses,word_counts);
 fs::write(out.join("RESULT.json"),&report).unwrap();print!("{}",report);
}
