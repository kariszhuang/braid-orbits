// Solve unramified voltage constraints for arbitrary connected quotient triples.
// Each block is an edge between its gamma-cycle and beta*gamma-cycle.
// Choose a spanning tree rooted at the latter's cycle containing block zero.
// Free non-tree voltages; eliminate leaves. No rejection of internal voltages.
struct QuotientForest {
    beta: Perm,
    gamma: Perm,
    r: usize,
    k: usize,
    orders: Vec<Vec<usize>>,
    root: usize,
    free_edges: Vec<usize>,
    leaves: Vec<(usize,usize)>,
    genus: usize,
}

fn qf_cycles(p: Perm, r: usize) -> Vec<Vec<usize>> {
    let mut seen=[false;MAX_D];let mut out=Vec::new();
    for x in 0..r {if !seen[x] {
        let mut cycle=Vec::new();let mut j=x;
        while !seen[j] {seen[j]=true;cycle.push(j);j=p.img[j] as usize;}
        out.push(cycle);
    }}
    out
}

impl QuotientForest {
    fn new(beta: Perm, gamma: Perm, r: usize, k: usize) -> Option<Self> {
        if r<3 || k<3 || r*k>MAX_D {return None;}
        Self::tree(beta,gamma,r,k,0)
    }

    fn tree(beta:Perm,gamma:Perm,r:usize,k:usize,root_cycle:usize)->Option<Self>{
        Self::tree_any(beta,gamma,r,k,root_cycle,false)
    }
    fn tree_any(beta:Perm,gamma:Perm,r:usize,k:usize,root_cycle:usize,rational:bool)->Option<Self>{
        if r<if rational{2}else{3}||k<2||r*k>MAX_D{return None;}
        let delta=beta.compose(gamma,r);
        if !rational && (beta.is_identity(r)||gamma.is_identity(r)||delta.is_identity(r)) {return None;}
        let mut orders=qf_cycles(gamma,r);let ng=orders.len();
        let ds=qf_cycles(delta,r);if root_cycle>=ds.len(){return None;}let root=ng+root_cycle;
        let total=sm_index(beta,r)+sm_index(gamma,r)+sm_index(delta,r);
        if total<if rational{2*r-2}else{2*r} || total%2!=0 {return None;}
        let genus=1+total/2-r;
        orders.extend(ds);
        let mut ends=vec![(0,0);r];
        for (v,cycle) in orders.iter().enumerate() {for &e in cycle {
            if v<ng {ends[e].0=v;}else{ends[e].1=v;}
        }}
        let mut parent=vec![usize::MAX;orders.len()];let mut edge=parent.clone();
        parent[root]=root;let mut queue=vec![root];let mut tree=vec![false;r];
        let mut cursor=0;
        while cursor<queue.len() {
            let v=queue[cursor];cursor+=1;
            for &e in &orders[v] {
                let (x,y)=ends[e];let u=if x==v {y}else{x};
                if parent[u]==usize::MAX {parent[u]=v;edge[u]=e;tree[e]=true;queue.push(u);}
            }
        }
        if queue.len()!=orders.len() {return None;}
        let free_edges=(0..r).filter(|&e|!tree[e]).collect();
        let leaves=queue[1..].iter().rev().map(|&v|(v,edge[v])).collect();
        Some(Self{beta,gamma,r,k,orders,root,free_edges,leaves,genus})
    }

    fn solve(&self, free: &[Perm]) -> [Perm;4] {
        let (voltage,return_map)=self.voltages(free,None);
        nc_lift(self.beta,self.gamma,self.k,&voltage,return_map.inverse(self.k))
    }

    fn voltages(&self,free:&[Perm],ramified_vertex:Option<usize>)->(Vec<Perm>,Perm){
        assert_eq!(free.len(),self.free_edges.len());
        let mut voltage=vec![Perm::identity(self.k);self.r];
        let mut assigned=vec![false;self.r];
        for (&e,&v) in self.free_edges.iter().zip(free) {voltage[e]=v;assigned[e]=true;}
        for &(v,e) in &self.leaves {
            let order=&self.orders[v];let cut=order.iter().position(|&x|x==e).unwrap();
            let mut before=Perm::identity(self.k);let mut after=before;
            for &j in &order[..cut] {assert!(assigned[j]);before=voltage[j].compose(before,self.k);}
            for &j in &order[cut+1..] {assert!(assigned[j]);after=voltage[j].compose(after,self.k);}
            // after * unknown * before = 1; composition acts rightmost first.
            let mut target=Perm::identity(self.k);
            if ramified_vertex==Some(v){assert_eq!(self.k,2);target.img.swap(0,1);}
            voltage[e]=after.inverse(self.k).compose(target,self.k).compose(before.inverse(self.k),self.k);
            assigned[e]=true;
        }
        assert!(assigned.iter().all(|&v|v));
        let mut return_map=Perm::identity(self.k);
        for &j in &self.orders[self.root] {return_map=voltage[j].compose(return_map,self.k);}
        (voltage,return_map)
    }
}

// More general than nc_hypothesis: quotient entries may have fixed blocks and
// unequal/multiple cycle lengths. Unramified means exactly k lifts per cycle.
fn qf_hypothesis(t: &[Perm;4], r: usize, k: usize) -> bool {
    let n=r*k;let mut local=0;
    for &p in t {
        let mut q=Perm::identity(r);
        for j in 0..r {
            q.img[j]=(p.img[j*k] as usize/k) as u8;
            if (0..k).any(|u|p.img[j*k+u] as usize/k!=q.img[j] as usize) {return false;}
        }
        if q.is_identity(r) && !p.is_identity(n) {
            let blocks:HashSet<_>=(0..n).filter(|&x|p.img[x] as usize!=x).map(|x|x/k).collect();
            if blocks.len()!=1 {return false;}
            local+=1;
        } else {
            let mut expected:Vec<_>=sp_type(q,r).into_iter().flat_map(|l|std::iter::repeat(l).take(k)).collect();
            expected.sort_unstable();let mut actual=sp_type(p,n);actual.sort_unstable();
            if expected!=actual {return false;}
        }
    }
    local==1
}

fn qf_shapes(n: usize) -> Vec<(usize,usize)> {
    (4..=n/3).filter(|r|n%r==0).map(|r|(r,n/r)).collect()
}

fn qf_tuple(n: usize, rng: &mut SplitMix64) -> Option<([Perm;4],usize,usize)> {
    let shapes=qf_shapes(n);if shapes.is_empty(){return None;}
    let (r,k)=shapes[rng.below(shapes.len())];
    for _ in 0..128 {
        let beta=if rng.below(2)==0 {dd_full_cycle(r,rng)}else{
            let mut p=Perm::identity(r);rng.shuffle(&mut p.img[..r]);p
        };
        let mut gamma=Perm::identity(r);rng.shuffle(&mut gamma.img[..r]);
        if [beta,gamma,beta.compose(gamma,r)].iter().all(|&p|sm_index(p,r)==r-1) {continue;}
        let Some(frame)=QuotientForest::new(beta,gamma,r,k) else{continue;};
        let free:Vec<_>=(0..frame.free_edges.len()).map(|_|{
            let mut p=Perm::identity(k);rng.shuffle(&mut p.img[..k]);p
        }).collect();
        let t=frame.solve(&free);
        debug_assert!(tuple_product_is_identity(n,&t));
        if !t[0].is_identity(n) {debug_assert!(qf_hypothesis(&t,r,k));}
        return Some((t,r,k));
    }
    None
}

#[cfg(test)]
mod quotient_forest_tests {
    use super::*;
    #[test]
    fn mixed_types_unramified_and_braid_stable() {
        let mut rng=SplitMix64::new(8761221);let mut good=0;let mut mixed=0;
        for n in [12,15,16,18,20,21,24,25] {for _ in 0..100 {
            let (t,r,k)=qf_tuple(n,&mut rng).unwrap();
            assert!(tuple_product_is_identity(n,&t));
            if t[0].is_identity(n)||!generators_transitive(n,&t){continue;}
            good+=1;assert!(qf_hypothesis(&t,r,k));
            let mut h=t;
            for _ in 0..100 {
                assert!(qf_hypothesis(&h,r,k));let v=sm_values(&h,n);
                assert!(v.0);assert_eq!(v.1,0);
                h=hurwitz_move_dynamic(n,h,rng.below(6));
            }
            if sp_type(t[1],n)!=sp_type(t[2],n){mixed+=1;}
        }}
        assert!(good>50 && mixed>25);
    }
    #[test]
    fn tree_solver_equals_exhaustive_internal_voltage_population() {
        // Even quotient degree, mixed (4),(3,1),(4): previously absent.
        let mut b=Perm::identity(4);b.img[..4].copy_from_slice(&[1,2,3,0]);
        let mut c=Perm::identity(4);c.img[..4].copy_from_slice(&[1,2,0,3]);
        let frame=QuotientForest::new(b,c,4,3).unwrap();assert_eq!(frame.genus,1);
        let local=nc_local_permutations(3);let mut solved=HashSet::new();
        for rank in 0usize..6usize.pow(frame.free_edges.len() as u32) {
            let mut v=rank;let free:Vec<_>=(0..frame.free_edges.len()).map(|_|{let p=local[v%6];v/=6;p}).collect();
            solved.insert(frame.solve(&free).map(|p|p.img[..12].to_vec()));
        }
        let mut direct=HashSet::new();
        for rank in 0usize..1296 {
            let mut v=rank;let volt:Vec<_>=(0..4).map(|_|{let p=local[v%6];v/=6;p}).collect();
            for &a in &local {
                let t=nc_lift(b,c,3,&volt,a);
                if sp_type(t[2],12)==vec![3,3,3] && sp_type(t[3],12)==vec![4,4,4] {
                    direct.insert(t.map(|p|p.img[..12].to_vec()));
                }
            }
        }
        assert_eq!(solved,direct);
    }
}
