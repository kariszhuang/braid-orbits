//! Bounded simultaneous cycle-count completion. A budget stop is UNDECIDED.
//! Images c(u)=v simultaneously set (F_j c)(u)=F_j(v) for every constraint.
//! No block partition, group, smooth-orbit, or geometric assumption is made.
pub const LIMIT:usize=63;
pub const PRODUCTS:usize=7;
pub type Image=[u8;LIMIT];
#[derive(Clone,Copy)]
struct Paths {head:[[u8;LIMIT];PRODUCTS],tail:[[u8;LIMIT];PRODUCTS],closed:[u8;PRODUCTS]}
impl Paths {
    fn new()->Self {let id=std::array::from_fn(|i|i as u8);Self{head:[id;PRODUCTS],tail:[id;PRODUCTS],closed:[0;PRODUCTS]}}
    fn link(&mut self,u:usize,v:usize,left:&[Image;PRODUCTS],goals:&[u8;PRODUCTS],remaining:usize)->bool {
        for j in 0..PRODUCTS {
            let w=left[j][v]as usize;let h=self.head[j][u]as usize;
            if w==h {self.closed[j]+=1;}else{
                let t=self.tail[j][w]as usize;self.head[j][t]=h as u8;self.tail[j][h]=t as u8;
            }
            if self.closed[j]>goals[j]||self.closed[j]as usize+remaining<goals[j]as usize{return false;}
        }true
    }
}
pub struct Answer {pub image:Option<Image>,pub nodes:u64,pub exhausted:bool}
fn next(seed:&mut u64)->u64 {*seed=seed.wrapping_add(0x9e3779b97f4a7c15);let mut z=*seed;z=(z^(z>>30)).wrapping_mul(0xbf58476d1ce4e5b9);z=(z^(z>>27)).wrapping_mul(0x94d049bb133111eb);z^(z>>31)}
struct Search<'a> {n:usize,left:&'a[Image;PRODUCTS],inverse:[Image;PRODUCTS],goals:[u8;PRODUCTS],rows:Image,nodes:u64,budget:u64,seed:u64,limited:bool}
impl Search<'_> {
    fn visit(&mut self,s:Paths,depth:usize,unused:u64,c:&mut Image)->bool {
        if depth==self.n{return s.closed==self.goals;}
        if self.nodes>=self.budget {self.limited=true;return false;}
        let u=self.rows[depth]as usize;let remaining=self.n-depth;let mut allowed=unused;
        for j in 0..PRODUCTS {
            let v=self.inverse[j][s.head[j][u]as usize]as usize;
            // Once the cycle quota is met, another closure is forbidden.
            if s.closed[j]==self.goals[j]{allowed&=!(1u64<<v);}
            // If all remaining assignments must close, this edge is forced.
            if s.closed[j]as usize+remaining==self.goals[j]as usize {allowed&=1u64<<v;}
        }
        let pivot=(next(&mut self.seed)%(self.n as u64))as usize;
        let mut order=[0u8;LIMIT];let mut len=0;
        for v in pivot..self.n {if allowed&(1u64<<v)!=0{order[len]=v as u8;len+=1;}}
        for v in 0..pivot {if allowed&(1u64<<v)!=0{order[len]=v as u8;len+=1;}}
        for &v in &order[..len] {
            if self.nodes>=self.budget {self.limited=true;return false;}
            self.nodes+=1;let mut t=s;
            if t.link(u,v as usize,self.left,&self.goals,remaining-1){
                c[u]=v;if self.visit(t,depth+1,unused&!(1u64<<v),c){return true;}
            }
        }false
    }
}
pub fn solve(n:usize,left:&[Image;PRODUCTS],indices:[usize;PRODUCTS],seed:u64,budget:u64)->Answer {
    assert!((1..=LIMIT).contains(&n)&&indices.iter().all(|&k|k<n)&&budget>0);
    let mut inverse=[[0;LIMIT];PRODUCTS];
    for j in 0..PRODUCTS {let mut seen=0u64;for u in 0..n{let v=left[j][u]as usize;assert!(v<n&&seen&(1u64<<v)==0);seen|=1u64<<v;inverse[j][v]=u as u8;}}
    let mut search=Search{n,left,inverse,goals:indices.map(|k|(n-k)as u8),rows:std::array::from_fn(|i|i as u8),nodes:0,budget,seed,limited:false};
    for i in (1..n).rev(){let j=(next(&mut search.seed)%((i+1)as u64))as usize;search.rows.swap(i,j);}
    let mut c=std::array::from_fn(|i|i as u8);
    let found=search.visit(Paths::new(),0,(1u64<<n)-1,&mut c);
    Answer{image:found.then_some(c),nodes:search.nodes,exhausted:!found&&!search.limited}
}
#[cfg(test)]
mod tests {
    use super::*;
    fn index(p:Image,n:usize)->usize{let mut used=0u64;let mut k=0;for r in 0..n{if used&(1<<r)==0{k+=1;let mut x=r;while used&(1<<x)==0{used|=1<<x;x=p[x]as usize;}}}n-k}
    fn perms(n:usize)->Vec<Image>{fn rec(k:usize,n:usize,p:&mut Image,v:&mut Vec<Image>){if k==n{v.push(*p);return;}for j in k..n{p.swap(k,j);rec(k+1,n,p,v);p.swap(k,j);}}let mut out=vec![];rec(0,n,&mut std::array::from_fn(|i|i as u8),&mut out);out}
    #[test]fn all_small_index_vectors_match_literal_existence(){
        let n=4;let p=perms(n);
        for frame in 0..24 {let left=std::array::from_fn(|j|p[(frame+7*j)%24]);
            let attainable:std::collections::HashSet<_>=p.iter().map(|c|std::array::from_fn(|j|index(std::array::from_fn(|u|if u<n{left[j][c[u]as usize]}else{u as u8}),n))).collect();
            for code in 0usize..256 {
                let goals=std::array::from_fn(|j|if j<4{(code>>(2*j))&3}else{(code>>(2*(j-4)))&3});
                let answer=solve(n,&left,goals,code as u64,100000);
                assert_eq!(answer.image.is_some(),attainable.contains(&goals));assert!(answer.image.is_some()||answer.exhausted);
                if let Some(c)=answer.image{for j in 0..PRODUCTS{assert_eq!(index(std::array::from_fn(|u|if u<n{left[j][c[u]as usize]}else{u as u8}),n),goals[j]);}}
            }
        }
    }
    #[test]fn budget_is_never_an_exclusion(){let id=std::array::from_fn(|i|i as u8);let a=solve(15,&[id;PRODUCTS],[14;PRODUCTS],1,1);assert!(a.image.is_none());assert!(!a.exhausted);assert_eq!(a.nodes,1);}
}
