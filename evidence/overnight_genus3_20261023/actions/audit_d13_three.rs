//! Independent full4^3*1 support class, no assumed cyclic colouring construction.
use std::{fs,io::Write,path::Path};
const N:usize=13;type P=[u8;N];type T=[P;4];
const ORD:[[usize;3];6]=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];
fn mul(a:P,b:P)->P{std::array::from_fn(|i|a[b[i]as usize])}
fn inv(a:P)->P{let mut r=[0;N];for i in 0..N{r[a[i]as usize]=i as u8;}r}
fn ix(a:P)->usize{let(mut seen,mut c)=(0u16,0);for mut i in 0..N{if seen&(1<<i)!=0{continue;}c+=1;while seen&(1<<i)==0{seen|=1<<i;i=a[i]as usize;}}N-c}
fn step(t:T,code:usize)->T{let mut u=t;let j=code/2;let(a,b)=(t[j],t[j+1]);if code%2==0{u[j]=mul(mul(a,b),inv(a));u[j+1]=a;}else{u[j]=b;u[j+1]=mul(mul(inv(b),a),b);}u}
fn bad(t:T,w:&[usize])->bool{let mut u=t;for &j in w{u=step(u,j);}let r=ix(mul(u[0],u[1]));r!=ix(u[0])+ix(u[1])&&r!=ix(u[2])+ix(u[3])}
fn main(){
 let args:Vec<_>=std::env::args().collect();let out=Path::new(args.get(1).expect("newoutputdirectory"));assert!(!out.exists());fs::create_dir(out).unwrap();
 let id:P=std::array::from_fn(|i|i as u8);let p:P=std::array::from_fn(|i|((i+1)%N)as u8);let mut a=id;a[0]=4;a[4]=8;a[8]=0;let c=mul(inv(a),p);let words:Vec<Vec<usize>>=vec![vec![3],vec![2],vec![3,3],vec![0,0,3],vec![2,2],vec![0,3,3],vec![1,1,2]];
 let(mut layouts,mut good,mut wrong,mut negative,mut unresolved)=(0usize,0usize,0usize,0usize,0usize);let mut records=vec![];
 for fixed in 0..N{let rem:Vec<_>=(0..N).filter(|&x|x!=fixed).collect();let anchor=rem[0];
  for i in 1..12{for j in i+1..12{for k in j+1..12{let b0=[anchor,rem[i],rem[j],rem[k]];let tail:Vec<_>=rem.iter().copied().filter(|x|!b0.contains(x)).collect();
   for x in 1..8{for y in x+1..8{for z in y+1..8{let b1=[tail[0],tail[x],tail[y],tail[z]];let b2:Vec<_>=tail.iter().copied().filter(|x|!b1.contains(x)).collect();layouts+=1;
    let mut color=[3;N];for (j,b)in[&b0[..],&b1[..],&b2[..]].iter().enumerate(){for &v in b.iter(){color[v]=j;}}
    if !(0..N).all(|s|{let(x,y,z)=(color[s],color[(s+4)%N],color[(s+8)%N]);x!=y&&x!=z&&y!=z}){continue;}good+=1;
    for o0 in ORD{for o1 in ORD{for o2 in ORD{let mut e=id;for(block,order)in[(&b0[..],o0),(&b1[..],o1),(&b2[..],o2)]{let cy=[block[0],block[order[0]],block[order[1]],block[order[2]]];for r in 0..4{e[cy[r]]=cy[(r+1)%4]as u8;}}
     records.push(e);let b=mul(inv(p),inv(e));if ix(b)!=9{wrong+=1;continue;}let t=[a,c,b,e];assert_eq!(mul(mul(mul(a,c),b),e),id);if words.iter().any(|w|bad(t,w)){negative+=1;}else{unresolved+=1;}
    }}}
   }}}
  }}}
 }
 assert_eq!(layouts,75075);assert_eq!(good,13);assert_eq!(records.len(),2808);records.sort();assert!(records.windows(2).all(|v|v[0]!=v[1]));let mut bin=fs::File::create(out.join("SORTED_ENUMERATION.bin")).unwrap();for e in records{bin.write_all(&e).unwrap();}
 let report=format!("{{\"full_support_partitions\":{},\"full_labelled_4_cubed_1_class\":{},\"necessary_supports\":{},\"necessary_e\":2808,\"wrong_b_index\":{},\"literal_negatives\":{},\"maximum_word_length\":3,\"unresolved\":{}}}\n",layouts,layouts*216,good,wrong,negative,unresolved);fs::write(out.join("RESULT.json"),&report).unwrap();print!("{}",report);
}
