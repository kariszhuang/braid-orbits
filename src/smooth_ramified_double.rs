// A ramified double lift of an all-smooth base. The base's localized entry
// lifts unramified; each other entry has at most one ramified cycle. Every
// additive incidence-forest component then has at most one branched vertex.
struct DoubleFrame { base:[Perm;4],tree:QuotientForest,ramified:usize }
struct RamifiedDoubleFrames { frames:Vec<DoubleFrame>,cursor:usize,initial:usize,r:usize }

fn rd_pure(q:Perm,r:usize)->Perm {
    let mut p=Perm::identity(2*r);
    for x in 0..r{for u in 0..2{p.img[2*x+u]=(2*q.img[x] as usize+u) as u8;}}p
}
fn rd_lift(base:[Perm;4],r:usize,tree:&QuotientForest,ramified:usize,rank:usize)->[Perm;4]{
    let free:Vec<_>=(0..tree.free_edges.len()).map(|i|{
        let mut p=Perm::identity(2);if rank&(1<<i)!=0{p.img.swap(0,1);}p
    }).collect();
    let (voltage,returned)=tree.voltages(&free,Some(ramified));
    assert_eq!(sm_index(returned,2),1);
    let a=rd_pure(base[0],r);let b=rd_pure(base[1],r);let mut c=Perm::identity(2*r);
    for x in 0..r{for u in 0..2{c.img[2*x+u]=(2*base[2].img[x] as usize+voltage[x].img[u] as usize) as u8;}}
    let t=[a,b,c,a.compose(b,2*r).compose(c,2*r).inverse(2*r)];
    debug_assert_eq!(t.map(|p|sm_index(p,2*r)),[2*sm_index(base[0],r),2*sm_index(base[1],r),2*sm_index(base[2],r)+1,2*sm_index(base[3],r)+1]);
    t
}
impl RamifiedDoubleFrames {
    fn new(path:&Path,n:usize,cursor:usize)->Result<Self,String>{
        if ![18,24].contains(&n){return Err("ramified double degree18/24".into());}
        let r=n/2;let mut frames=vec![];
        for line in fs::read_to_string(path).map_err(|e|e.to_string())?.lines(){
            let (degree,seed)=line.split_once('\t').ok_or("double base fields")?;
            if degree.parse::<usize>().map_err(|_|"base degree")?!=r{continue;}
            let base=sm_parse(seed,r)?;
            if ff_product_index(base[0],base[1],r)!=sm_index(base[0],r)+sm_index(base[1],r){return Err("base first pair geodesy".into());}
            let p=base[0].compose(base[1],r);let delta=p.compose(base[2],r);
            for root in 0..qf_cycles(delta,r).len(){for ramified in 0..qf_cycles(base[2],r).len(){
                let tree=QuotientForest::tree(p,base[2],r,2,root).ok_or("double voltage tree")?;
                frames.push(DoubleFrame{base,tree,ramified});
            }}
        }
        let total:usize=frames.iter().map(|f|1usize<<f.tree.free_edges.len()).sum();
        if frames.is_empty()||cursor>total{return Err("double base/cursor population".into());}
        Ok(Self{frames,cursor,initial:cursor,r})
    }
    fn total(&self)->usize{self.frames.iter().map(|f|1usize<<f.tree.free_edges.len()).sum()}
    fn next(&mut self)->Option<MechanismProposal>{
        if self.cursor==self.total(){return None;}
        let mut rank=self.cursor;
        for (source,f) in self.frames.iter().enumerate(){
            let size=1<<f.tree.free_edges.len();
            if rank<size {self.cursor+=1;return Some(MechanismProposal{tuple:rd_lift(f.base,self.r,&f.tree,f.ramified,rank),source,arm:5,points:vec![]});}
            rank-=size;
        }
        unreachable!()
    }
    fn save(&self,out:&Path)->Result<(),String>{
        fs::write(out.join("DOUBLE_LIFT.tsv"),format!("degree\tbase_degree\tframes\tstart\tnext\ttotal\tproducer_complete\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n",2*self.r,self.r,self.frames.len(),self.initial,self.cursor,self.total(),self.cursor==self.total())).map_err(|e|e.to_string())
    }
}

#[cfg(test)]
mod ramified_double_tests {
    use super::*;
    fn base()->[Perm;4]{sm_parse("[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]",9).unwrap()}
    #[test] fn solved_binary_population_equals_all512_direct_voltages(){
        let q=base();let p=q[0].compose(q[1],9);let mut solved=HashSet::new();
        for root in 0..3{for v in 0..3{let tree=QuotientForest::tree(p,q[2],9,2,root).unwrap();
            assert_eq!(tree.free_edges.len(),4);
            for rank in 0..16 {solved.insert(rd_lift(q,9,&tree,v,rank));}
        }}
        let mut direct=HashSet::new();let a=rd_pure(q[0],9);let b=rd_pure(q[1],9);
        for rank in 0..512 {let mut c=Perm::identity(18);
            for x in 0..9{for u in 0..2{c.img[2*x+u]=(2*q[2].img[x] as usize+(u^((rank>>x)&1))) as u8;}}
            let t=[a,b,c,a.compose(b,18).compose(c,18).inverse(18)];
            if sm_index(t[2],18)==13&&sm_index(t[3],18)==13{direct.insert(t);}
        }
        assert_eq!(solved,direct);assert_eq!(solved.len(),144);
    }
    #[test] fn every_double_lift_and_random_braid_translate_is_smooth(){
        let q=base();let mut rng=SplitMix64::new(2026101344);
        for root in 0..3{for v in 0..3{let tree=QuotientForest::tree(q[0].compose(q[1],9),q[2],9,2,root).unwrap();
            for rank in 0..16 {let mut t=rd_lift(q,9,&tree,v,rank);assert!(generators_transitive(18,&t));assert!(!mb_primitive(&t,18));
                for _ in 0..100 {let values=sm_values(&t,18);assert!(values.0);assert_eq!((values.1,values.2),(0,4));t=hurwitz_move_dynamic(18,t,rng.below(6));}
            }
        }}
    }
}
