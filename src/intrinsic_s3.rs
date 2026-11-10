//! Elliptic S3 covers in two-handle tree--cotree coordinates.
//! This module constructs tuples, and never certifies braid smoothness.
use std::collections::VecDeque;

pub const MAX_Q:usize=21;
pub const MAX_N:usize=63;
pub const PERMS:[[u8;3];6]=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
pub const TRANS:[u8;3]=[1,2,5];
pub type Lift=[[u8;MAX_N];4];

#[derive(Clone)]
pub struct S3 {pub mul:[[u8;6];6],pub inv:[u8;6]}
impl S3 {
    pub fn new()->Self{
        let mut mul=[[0;6];6];let mut inv=[0;6];
        for a in 0..6{for b in 0..6{
            let p=std::array::from_fn(|x|PERMS[a][PERMS[b][x]as usize]);
            mul[a][b]=PERMS.iter().position(|&v|v==p).unwrap()as u8;
        }}
        for a in 0..6{inv[a]=(0..6).find(|&b|mul[a][b]==0).unwrap()as u8;}
        Self{mul,inv}
    }
}
pub fn compose(a:&[u8],b:&[u8])->Vec<u8>{b.iter().map(|&x|a[x as usize]).collect()}
pub fn inverse(p:&[u8])->Vec<u8>{let mut out=vec![0;p.len()];for(x,&y)in p.iter().enumerate(){out[y as usize]=x as u8;}out}
pub fn cycles(p:&[u8])->Vec<Vec<usize>>{
    let mut seen=0u32;let mut out=Vec::new();
    for start in 0..p.len(){if seen&(1<<start)!=0{continue;}
        let mut row=Vec::new();let mut at=start;
        while seen&(1<<at)==0{seen|=1<<at;row.push(at);at=p[at]as usize;}out.push(row);
    }out
}
struct Dsu {parent:Vec<usize>}
impl Dsu {
    fn new(n:usize)->Self{Self{parent:(0..n).collect()}}
    fn root(&mut self,mut x:usize)->usize{while self.parent[x]!=x{self.parent[x]=self.parent[self.parent[x]];x=self.parent[x];}x}
    fn join(&mut self,a:usize,b:usize)->bool{let a=self.root(a);let b=self.root(b);if a==b{return false;}self.parent[b]=a;true}
}
#[derive(Clone,Copy)]
enum Atom {Edge(usize,bool),Local(usize)}
pub struct IntrinsicFrame {
    pub beta:Vec<u8>,pub gamma:Vec<u8>,pub delta:Vec<u8>,pub primal:Vec<usize>,
    pub dual:Vec<usize>,pub handles:[usize;2],pub faces:Vec<Vec<usize>>,
    leaves:Vec<(usize,usize)>,words:Vec<Vec<Atom>>,pub s3:S3,
}
impl IntrinsicFrame {
    pub fn new(beta:Vec<u8>,gamma:Vec<u8>)->Result<Self,String>{
        let q=beta.len();if gamma.len()!=q||q<3||q>MAX_Q{return Err("intrinsic degree3..21".into());}
        for p in [&beta,&gamma]{let mut labels=p.to_vec();labels.sort_unstable();if labels!=(0..q as u8).collect::<Vec<_>>(){return Err("quotient permutation".into());}}
        let delta=compose(&beta,&gamma);let black=cycles(&beta);let white=cycles(&gamma);let faces=cycles(&delta);
        if black.len()+white.len()+faces.len()!=q{return Err("intrinsic quotient genus must be one".into());}
        let labels=|rows:&Vec<Vec<usize>>|{let mut v=vec![0;q];for(j,row)in rows.iter().enumerate(){for &x in row{v[x]=j;}}v};
        let bl=labels(&black);let wh=labels(&white);let fa=labels(&faces);
        let mut primal=Vec::new();let mut dsu=Dsu::new(black.len()+white.len());
        for e in 0..q{if dsu.join(bl[e],black.len()+wh[e]){primal.push(e);}}
        if primal.len()!=black.len()+white.len()-1{return Err("disconnected quotient dessin".into());}
        let gi=inverse(&gamma);let mut dual=Vec::new();let mut adj=vec![Vec::new();faces.len()];let mut dsu=Dsu::new(faces.len());
        for e in 0..q{if primal.contains(&e){continue;}let(a,b)=(fa[e],fa[gi[e]as usize]);
            if dsu.join(a,b){dual.push(e);adj[a].push((b,e));adj[b].push((a,e));}
        }
        if dual.len()!=faces.len()-1{return Err("dual complement disconnected".into());}
        let handles:Vec<_>=(0..q).filter(|e|!primal.contains(e)&&!dual.contains(e)).collect();
        let handles:[usize;2]=handles.try_into().map_err(|_|"elliptic handle count")?;
        let mut queue=vec![0];let mut parent=vec![None;faces.len()];parent[0]=Some((0,usize::MAX));let mut at=0;
        while at<queue.len(){let v=queue[at];at+=1;for &(u,e)in &adj[v]{if parent[u].is_none(){parent[u]=Some((v,e));queue.push(u);}}}
        let leaves=queue[1..].iter().rev().map(|&v|(v,parent[v].unwrap().1)).collect();
        let words=faces.iter().map(|face|face.iter().flat_map(|&j|[Atom::Edge(j,false),Atom::Edge(gamma[j]as usize,true),Atom::Local(delta[j]as usize)]).collect()).collect();
        Ok(Self{beta,gamma,delta,primal,dual,handles,faces,leaves,words,s3:S3::new()})
    }
    fn value(&self,word:&[Atom],z:&[u8;MAX_Q],local:&[u8])->u8{
        let mut result=0;
        for atom in word{let value=match *atom{Atom::Local(j)=>local[j],Atom::Edge(j,inv)=>{
            assert!(z[j]<6,"unassigned intrinsic edge");if inv{self.s3.inv[z[j]as usize]}else{z[j]}
        }};result=self.s3.mul[value as usize][result as usize];}result
    }
    pub fn solve(&self,local:&[u8],handle_values:[u8;2])->Option<Lift>{
        let q=self.beta.len();assert_eq!(local.len(),q);assert!(local.iter().all(|&v|v<6)&&handle_values.iter().all(|&v|v<6));
        let mut z=[u8::MAX;MAX_Q];for &j in &self.primal{z[j]=0;}for(j,v)in self.handles.into_iter().zip(handle_values){z[j]=v;}
        for &(face,edge)in &self.leaves{
            let word=&self.words[face];let positions:Vec<_>=word.iter().enumerate().filter_map(|(j,a)|if matches!(a,Atom::Edge(i,_)if *i==edge){Some(j)}else{None}).collect();
            assert_eq!(positions.len(),1,"dual leaf edge must occur once");let cut=positions[0];
            let before=self.value(&word[..cut],&z,local);let after=self.value(&word[cut+1..],&z,local);
            let mut value=self.s3.mul[self.s3.inv[after as usize]as usize][self.s3.inv[before as usize]as usize];
            if matches!(word[cut],Atom::Edge(_,true)){value=self.s3.inv[value as usize];}z[edge]=value;
            assert_eq!(self.value(word,&z,local),0,"intrinsic leaf equation");
        }
        if self.value(&self.words[0],&z,local)!=0{return None;}
        let mut out:[[u8;MAX_N];4]=std::array::from_fn(|_|std::array::from_fn(|x|x as u8));
        for j in 0..q{let v=self.s3.mul[self.s3.inv[z[self.gamma[j]as usize]as usize]as usize][z[j]as usize];
            for x in 0..3{out[0][3*j+x]=(3*j+PERMS[local[j]as usize][x]as usize)as u8;
                out[1][3*j+x]=3*self.beta[j]+x as u8;
                out[2][3*j+x]=3*self.gamma[j]+PERMS[v as usize][x];
            }
        }
        for x in 0..3*q{let image=out[0][out[1][out[2][x]as usize]as usize]as usize;out[3][image]=x as u8;}
        Some(out)
    }
    pub fn image_mask(p:&[u8],mask:u32)->u32{
        let mut bits=mask;let mut out=0;while bits!=0{let x=bits.trailing_zeros()as usize;bits&=bits-1;out|=1<<p[x];}out
    }
    pub fn support_orbit(&self,start:u32)->(Vec<u32>,bool){
        let q=self.beta.len();assert_eq!(start.count_ones(),4);assert!(start<(1<<q));
        let masks:[Vec<u32>;3]=[cycles(&self.beta),cycles(&self.gamma),cycles(&self.delta)].map(|rows|rows.into_iter().map(|c|c.into_iter().fold(0,|m,x|m|1<<x)).collect());
        let mut seen=std::collections::HashSet::from([start]);let mut queue=VecDeque::from([start]);let mut orbit=Vec::new();let mut good=true;
        while let Some(mask)=queue.pop_front(){
            good&=masks.iter().flatten().all(|&cycle|(cycle&mask).count_ones()<=2);orbit.push(mask);
            for p in [&self.beta,&self.gamma]{let next=Self::image_mask(p,mask);if seen.insert(next){queue.push_back(next);}}
        }orbit.sort_unstable();(orbit,good)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn profile(p:&[u8])->Vec<usize>{let mut v:Vec<_>=cycles(p).into_iter().map(|c|c.len()).collect();v.sort_unstable();v}
    #[test]fn positive_genus_two_control_has_nine_exact_handle_frames(){
        let f=IntrinsicFrame::new(vec![1,2,0],vec![1,2,0]).unwrap();let mut count=0;
        for u in 0..6{for v in 0..6{if let Some(t)=f.solve(&[3,0,0],[u,v]){count+=1;
            assert_eq!(profile(&t[1][..9]),vec![3;3]);assert_eq!(profile(&t[2][..9]),vec![3;3]);assert_eq!(profile(&t[3][..9]),vec![3;3]);
        }}}assert_eq!(count,9);
    }
    #[test]fn nontrivial_dual_leaf_is_pure_and_tree_cotree_disjoint(){
        let f=IntrinsicFrame::new(vec![1,2,0,4,5,3],vec![1,3,4,0,5,2]).unwrap();
        assert_eq!(f.primal,vec![0,2,3]);assert_eq!(f.dual,vec![1]);assert_eq!(f.handles,[4,5]);let mut count=0;
        for u in 0..6{for v in 0..6{if let Some(t)=f.solve(&[1,2,0,1,2,0],[u,v]){count+=1;
            for p in &t[1..]{assert_eq!(profile(&p[..18]),vec![3;6]);}
        }}}assert_eq!(count,9);
    }
    #[test]fn disconnected_or_wrong_surface_never_constructs(){
        assert!(IntrinsicFrame::new(vec![1,0,2],vec![0,2,1]).is_err());
        assert!(IntrinsicFrame::new(vec![1,2,0,4,5,3],vec![1,2,0,4,5,3]).is_err());
    }
}
