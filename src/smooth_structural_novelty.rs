// New biased populations share the existing exact screens, keys and closure.
// Triangle budgets and prescribed returns are targets, never positive proofs.
fn sn_triangle_budget(t:&[Perm;4],n:usize)->bool {
    let mut ix=t.map(|p|sm_index(p,n));ix.sort_unstable();
    ix[0]+ix[3]>=n && ix[0]+ix[2]<n && ix[1]+ix[2]<n
}

fn sn_triangle_pair(n:usize,rng:&mut SplitMix64,min_low:usize,max_low:usize)->Option<[Perm;2]> {
    // Positive-genus triangle budgets force min index>=2, not balanced types.
    let high=(n-1)/2;let low=min_low.max(2);let top=max_low.min(high);
    if low>top{return None;}
    let ia=low+rng.below(top-low+1);let ib=ia+rng.below(high-ia+1);
    let mut colors=[0u8;MAX_D];colors[ia..ia+ib].fill(1);rng.shuffle(&mut colors[..ia+ib]);
    let mut parent=std::array::from_fn(|i|i);let mut pair=[Perm::identity(n);2];
    for &color in &colors[..ia+ib] {loop {
        let x=rng.below(n);let y=rng.below(n);
        if gd_join(&mut parent,x,y){pair[color as usize].img.swap(x,y);break;}
    }}
    debug_assert_eq!(sm_index(pair[0],n)+sm_index(pair[1],n),ff_product_index(pair[0],pair[1],n));
    Some(pair)
}

fn sn_triangle_passports(n:usize,min_r1:usize,max_r1:usize,min_genus:usize,target_genus:usize)->Vec<[usize;4]> {
    let g0=if target_genus>0{target_genus}else{min_genus.max(3)};
    let g1=if target_genus>0{target_genus}else{((n-1)/4).min(n.saturating_sub(1))};
    if g0>g1{return Vec::new();}
    let mut out=Vec::new();let r1max=max_r1.min((n-1)/2);
    for g in g0..=g1 {let total=2*n+2*g-2;
        for r1 in min_r1.max(2*g)..=r1max {
            for r2 in r1..=(n-1)/2 {for r3 in r2..=(n-1-r2) {
                let h=total-r1-r2-r3;
                // Reviewed C12/C13 high-index bound: in a positive-genus
                // TRIANGLE orbit, h >= floor(d/2) + 2g + 1. Apply it while
                // listing passports, before any permutation construction.
                let high_floor=n/2+2*g+1;
                if h>=r3&&h<n&&h>=(n+1)/2&&r1+h>=n&&h>=high_floor {out.push([r1,r2,r3,h]);}
            }}
        }
    }
    out
}

fn sn_is_involution(p:Perm,n:usize)->bool {(0..n).all(|x|p.img[p.img[x] as usize] as usize==x)}

fn sn_cycle_histogram(p:Perm,n:usize)->[usize;MAX_D+1] {
    let mut seen=[false;MAX_D];let mut hist=[0usize;MAX_D+1];
    for root in 0..n {if seen[root]{continue;}let mut x=root;let mut len=0;
        while !seen[x]{seen[x]=true;len+=1;x=p.img[x] as usize;}
        hist[len]+=1;
    }
    hist
}

fn sn_three_two_bulk_with_three_fixed(p:Perm,n:usize)->bool {
    if n%2!=0||n<8{return false;}let hist=sn_cycle_histogram(p,n);let half=n/2;
    hist[3]==1&&hist[2]==half-3&&hist[1]==3
}

// Independently reviewed C12 exclusion: a source-genus-three TRIANGLE
// component cannot contain a distinguished maximum-index entry of type
// 15 * 2^((d-15)/2). This cheap exact type check avoids futile closure work.
fn sn_fifteen_transposition_bulk(p:Perm,n:usize)->bool {
    if n<15||(n-15)%2!=0{return false;}
    let hist=sn_cycle_histogram(p,n);
    hist[15]==1&&hist[2]==(n-15)/2&&hist.iter().enumerate().all(|(len,&count)|
        len==2||len==15||count==0)
}

fn sn_centralizer_element(b:Perm,n:usize,rng:&mut SplitMix64)->Perm {
    let mut by_length:Vec<Vec<Vec<usize>>>=vec![Vec::new();n+1];let mut seen=[false;MAX_D];
    for root in 0..n {if seen[root]{continue;}let mut cycle=Vec::new();let mut x=root;
        while !seen[x]{seen[x]=true;cycle.push(x);x=b.img[x] as usize;}
        by_length[cycle.len()].push(cycle);
    }
    let mut c=Perm::identity(n);
    for cycles in by_length.iter_mut().skip(1) {
        let len=cycles.first().map_or(0,Vec::len);if len==0{continue;}
        let mut destination:Vec<usize>=(0..cycles.len()).collect();rng.shuffle(&mut destination);
        for source in 0..cycles.len(){let offset=rng.below(len);let target=destination[source];
            for pos in 0..len{c.img[cycles[source][pos]]=cycles[target][(pos+offset)%len] as u8;}
        }
    }
    debug_assert_eq!(c.compose(b,n),b.compose(c,n));c
}

fn sn_triangle_commuting_tuple(n:usize,rng:&mut SplitMix64,profiles:&[[usize;4]])->Option<MechanismProposal> {
    if profiles.is_empty(){return None;}
    let [r1,r2,r3,_]=profiles[rng.below(profiles.len())];
    let mut b=None;
    for _ in 0..128 {let candidate=ff_index_perm(n,r2,rng);
        if !sn_is_involution(candidate,n){b=Some(candidate);break;}}
    let b=b?;let mut c=None;
    for _ in 0..128 {let candidate=sn_centralizer_element(b,n,rng);
        if sm_index(candidate,n)==r3&&!sn_is_involution(candidate,n){c=Some(candidate);break;}}
    let c=c?;let frame=ForestFrame::new([b,b],n);let mut a=None;
    for _ in 0..16 {if let Some(candidate)=sn_joint_third(&frame,n,r1,rng,false){
        if !sn_is_involution(candidate,n){a=Some(candidate);break;}
    }}
    let a=a?;
    let d=a.compose(b,n).compose(c,n).inverse(n);let tuple=[a,b,c,d];
    debug_assert!(tuple_product_is_identity(n,&tuple));
    Some(MechanismProposal{tuple,source:r2,arm:5,points:vec![]})
}

// Exact d=16 construction for the remaining passport (7,7,7,15): build a
// geodesic triangle of three index-7 entries, then require their product to
// be the distinguished 16-cycle.
fn sn_triangle16_pairwise_low(rng:&mut SplitMix64)->Option<MechanismProposal> {
    let n=16;let mut b=Perm::identity(n);let mut c=Perm::identity(n);
    for i in 0..8{b.img[i]=((i+1)%8) as u8;c.img[8+i]=(8+(i+1)%8) as u8;}
    let bc=ForestFrame::new([b,c],n);let a=dd_joint_third(&bc,n,7,rng)?;
    if sm_index(a,n)!=7||sn_is_involution(a,n)
        ||ff_product_index(a,b,n)!=14||ff_product_index(a,c,n)!=14
        ||ff_product_index(b,c,n)!=14{return None;}
    let p=a.compose(b,n).compose(c,n);
    if sm_index(p,n)!=15{return None;}
    let d=p.inverse(n);let tuple=[a,b,c,d];
    debug_assert!(tuple_product_is_identity(n,&tuple));
    debug_assert_eq!(ff_product_index(b,c,n),14);
    Some(MechanismProposal{tuple,source:7,arm:6,points:vec![]})
}

// Build an additive degree-16 first pair from varied non-involution index-7
// cycle types. Each side has nine cycles, so the incidence graph has 18
// cycle-vertices and 16 point-edges; joining only distinct components
// constructs a two-component forest, hence ind(ab)=ind(a)+ind(b)=14. The
// types avoid both the reviewed low-involution obstruction and the excluded
// 3*2^5*1^3 special type.
fn sn_triangle16_forest_pair(rng:&mut SplitMix64)->Option<[Perm;2]> {
    const CYCLE_TYPES:[[usize;9];8]=[
        [4,2,2,2,2,1,1,1,1], [5,2,2,2,1,1,1,1,1],
        [8,1,1,1,1,1,1,1,1], [6,3,1,1,1,1,1,1,1],
        [4,4,2,1,1,1,1,1,1], [3,3,3,2,1,1,1,1,1],
        [7,2,1,1,1,1,1,1,1], [5,3,2,1,1,1,1,1,1],
    ];
    let n=16;
    for _ in 0..4096 {
        let mut left=CYCLE_TYPES[rng.below(CYCLE_TYPES.len())];
        let mut right=CYCLE_TYPES[rng.below(CYCLE_TYPES.len())];
        rng.shuffle(&mut left);rng.shuffle(&mut right);
        let mut remaining=[0usize;MAX_D];
        for j in 0..9 {remaining[j]=left[j];remaining[9+j]=right[j];}
        let mut parent=std::array::from_fn(|i|i);let mut edges=Vec::with_capacity(n);
        let mut failed=false;
        while edges.len()<n {
            let mut choices=[(0usize,0usize);81];let mut count=0usize;
            for a in 0..9 {if remaining[a]==0{continue;}
                for b in 9..18 {if remaining[b]>0 {
                    let ra=gd_find(&mut parent,a);let rb=gd_find(&mut parent,b);
                    if ra!=rb {choices[count]=(a,b);count+=1;}
                }}
            }
            if count==0 {failed=true;break;}
            let (a,b)=choices[rng.below(count)];
            if !gd_join(&mut parent,a,b) {failed=true;break;}
            remaining[a]-=1;remaining[b]-=1;edges.push((a,b));
        }
        if failed||remaining[..18].iter().any(|&x|x!=0){continue;}
        let mut incident:[Vec<usize>;18]=std::array::from_fn(|_|Vec::new());
        for (point,&(a,b)) in edges.iter().enumerate(){incident[a].push(point);incident[b].push(point);}
        let mut a=Perm::identity(n);let mut b=Perm::identity(n);
        for j in 0..9 {
            if incident[j].len()!=left[j]||incident[9+j].len()!=right[j]{failed=true;break;}
            rng.shuffle(&mut incident[j]);rng.shuffle(&mut incident[9+j]);
            for k in 0..incident[j].len(){a.img[incident[j][k]]=incident[j][(k+1)%incident[j].len()] as u8;}
            for k in 0..incident[9+j].len(){b.img[incident[9+j][k]]=incident[9+j][(k+1)%incident[9+j].len()] as u8;}
        }
        if failed||sm_index(a,n)!=7||sm_index(b,n)!=7||ff_product_index(a,b,n)!=14 {continue;}
        if sn_is_involution(a,n)||sn_is_involution(b,n){continue;}
        return Some([a,b]);
    }
    None
}

#[derive(Default)]
struct Triangle16RepairStats {
    pair_success:u64,pair_failure:u64,c_attempts:u64,c_generation_failure:u64,
    c_index_reject:u64,product_index_reject:u64,short_reject:[u64;5],
    deep_reject:u64,involution_reject:u64,special_pair_reject:u64,accepted:u64,
}

// Reuse the strong common-forest solver, but vary the third entry on one fixed
// low pair until the exact local short/deep braid screens pass. This directly
// targets the dominant sigma_1^2/sigma_4^2 witnesses seen in earlier runs.
fn sn_triangle16_middle_repair(rng:&mut SplitMix64,exact:Option<[usize;4]>,stats:&mut Triangle16RepairStats)->Option<MechanismProposal> {
    let n=16;if exact!=Some([7,7,7,15]){return None;}
    let Some(pair)=sn_triangle16_forest_pair(rng) else{stats.pair_failure+=1;return None;};
    stats.pair_success+=1;let frame=ForestFrame::new(pair,n);
    if frame.indices!=[7,7]{return None;}
    for _ in 0..64 {
        stats.c_attempts+=1;
        let c=match rng.below(3) {
            0=>dd_joint_third(&frame,n,7,rng),
            1=>sn_joint_third(&frame,n,7,rng,false),
            _=>Some(ff_index_perm(n,7,rng)),
        };
        let Some(c)=c else{stats.c_generation_failure+=1;continue;};
        let ic=sm_index(c,n);let id=ff_product_index(frame.p,c,n);
        if ic!=7 {stats.c_index_reject+=1;continue;}
        if id!=15 {stats.product_index_reject+=1;continue;}
        if let Some(word)=frame.short_word(c,ic,id,n) {stats.short_reject[word]+=1;continue;}
        if frame.deep_failure(c,ic,id,n,3).is_some(){stats.deep_reject+=1;continue;}
        let tuple=frame.tuple(c,n);
        let special=(0..3).filter(|&j|sn_three_two_bulk_with_three_fixed(tuple[j],n)).count();
        if (0..3).any(|j|sn_is_involution(tuple[j],n)){stats.involution_reject+=1;continue;}
        if special>=2{stats.special_pair_reject+=1;continue;}
        stats.accepted+=1;
        return Some(MechanismProposal{tuple,source:14,arm:7,points:vec![]});
    }
    None
}

fn sn_joint_third(frame:&ForestFrame,n:usize,target:usize,rng:&mut SplitMix64,strong:bool)->Option<Perm> {
    if strong{return dd_joint_third(frame,n,target,rng);}
    let mut c=Perm::identity(n);let mut products=frame.pair;
    for _ in 0..target {
        let labels=[dd_cycle_labels(c,n),dd_cycle_labels(products[0],n),dd_cycle_labels(products[1],n)];
        let mut choices=[(0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
        for x in 0..n {for y in x+1..n {
            if labels.iter().all(|v|v[x]!=v[y]){choices[count]=(x as u8,y as u8);count+=1;}
        }}
        if count==0{return None;}
        let (x,y)=choices[rng.below(count)];let (x,y)=(x as usize,y as usize);c.img.swap(x,y);
        for p in &mut products{p.img.swap(x,y);}
    }
    Some(c)
}

fn sn_triangle_tuple(n:usize,rng:&mut SplitMix64)->Option<MechanismProposal> {
    sn_triangle_targeted_tuple(n,rng,2,MAX_D,8)
}

fn sn_triangle_targeted_tuple(n:usize,rng:&mut SplitMix64,min_low:usize,max_low:usize,
                              strong_divisor:usize)->Option<MechanismProposal> {
    let mut stats=Triangle16RepairStats::default();
    sn_triangle_targeted_tuple_with_method(n,rng,min_low,max_low,strong_divisor,"mixed",&[],None,&mut stats)
}

fn sn_triangle_targeted_tuple_with_method(n:usize,rng:&mut SplitMix64,min_low:usize,max_low:usize,
                                          strong_divisor:usize,method:&str,profiles:&[[usize;4]],
                                          exact_profile:Option<[usize;4]>,stats:&mut Triangle16RepairStats)->Option<MechanismProposal> {
    if strong_divisor==0{return None;}
    if method=="middle_repair"{return sn_triangle16_middle_repair(rng,exact_profile,stats);}
    if method=="pairwise_low" {
        if n!=16||exact_profile!=Some([7,7,7,15]){return None;}
        return sn_triangle16_pairwise_low(rng);
    }
    if method=="centralizer"{return sn_triangle_commuting_tuple(n,rng,profiles);}
    if !["random","joint","strong","mixed"].contains(&method){return None;}
    let pair=sn_triangle_pair(n,rng,min_low,max_low)?;let frame=ForestFrame::new(pair,n);
    let [ia,ib]=frame.indices;let low=ia.max(ib).max(n+1-ia-ib);let high=n-1-ia.max(ib);
    if low>high{return None;}
    let target=if let Some(profile)=exact_profile {
        if ia!=profile[0]||ib!=profile[1]||profile[2]<low||profile[2]>high{return None;}
        profile[2]
    }else{low+rng.below(high-low+1)};
    let (c,method_id)=match method {
        "random"=>(ff_index_perm(n,target,rng),0),
        "joint"=>(sn_joint_third(&frame,n,target,rng,false)?,1),
        "strong"=>(dd_joint_third(&frame,n,target,rng)?,2),
        _=>{
            if rng.below(strong_divisor)==0 {
                if let Some(c)=dd_joint_third(&frame,n,target,rng){(c,2)}
                else if rng.below(2)==0 {(ff_index_perm(n,target,rng),0)}
                else {(sn_joint_third(&frame,n,target,rng,false)?,1)}
            }else if rng.below(2)==0 {(ff_index_perm(n,target,rng),0)}
            else {(sn_joint_third(&frame,n,target,rng,false)?,1)}
        }
    };
    let tuple=frame.tuple(c,n);
    if !sn_triangle_budget(&tuple,n){return None;}
    Some(MechanismProposal{tuple,source:ia+ib,arm:method_id,points:vec![]})
}

struct ReturnForestFrame { tree:QuotientForest,local:Vec<Perm> }
impl ReturnForestFrame {
    fn new(beta:Perm,gamma:Perm,r:usize,k:usize,local:Vec<Perm>)->Option<Self> {
        if local.len()!=r{return None;}
        let tree=QuotientForest::tree_any(beta,gamma,r,k,0,true)?;
        let mut a=Perm::identity(r*k);let mut b=a;
        for j in 0..r {for x in 0..k {
            a.img[j*k+x]=(j*k+local[j].img[x] as usize) as u8;
            b.img[j*k+x]=(beta.img[j] as usize*k+x) as u8;
        }}
        if ff_product_index(a,b,r*k)!=sm_index(a,r*k)+sm_index(b,r*k){return None;}
        Some(Self{tree,local})
    }
    fn solve(&self,free:&[Perm],targets:&[Perm])->[Perm;4] {
        let f=&self.tree;assert_eq!(targets.len(),f.orders.len());assert_eq!(free.len(),f.free_edges.len());
        let ng=qf_cycles(f.gamma,f.r).len();let delta=f.beta.compose(f.gamma,f.r);
        let mut voltage=vec![Perm::identity(f.k);f.r];let mut assigned=vec![false;f.r];
        for (&e,&v) in f.free_edges.iter().zip(free){voltage[e]=v;assigned[e]=true;}
        for &(v,e) in &f.leaves {
            let order=&f.orders[v];let cut=order.iter().position(|&x|x==e).unwrap();
            let factor=|j:usize|if v<ng{voltage[j]}else{self.local[delta.img[j] as usize].compose(voltage[j],f.k)};
            let mut before=Perm::identity(f.k);let mut after=before;
            for &j in &order[..cut]{assert!(assigned[j]);before=factor(j).compose(before,f.k);}
            for &j in &order[cut+1..]{assert!(assigned[j]);after=factor(j).compose(after,f.k);}
            let coefficient=if v<ng{Perm::identity(f.k)}else{self.local[delta.img[e] as usize]};
            voltage[e]=coefficient.inverse(f.k).compose(after.inverse(f.k),f.k)
                .compose(targets[v],f.k).compose(before.inverse(f.k),f.k);
            assigned[e]=true;
        }
        assert!(assigned.iter().all(|&v|v));
        let n=f.r*f.k;let mut t=[Perm::identity(n);4];
        for j in 0..f.r {for x in 0..f.k {
            t[0].img[j*f.k+x]=(j*f.k+self.local[j].img[x] as usize) as u8;
            t[1].img[j*f.k+x]=(f.beta.img[j] as usize*f.k+x) as u8;
            t[2].img[j*f.k+x]=(f.gamma.img[j] as usize*f.k+voltage[j].img[x] as usize) as u8;
        }}
        t[3]=t[0].compose(t[1],n).compose(t[2],n).inverse(n);
        t
    }
}

fn sn_return_tuple(n:usize,rng:&mut SplitMix64)->Option<MechanismProposal> {
    let mut shapes=[(0usize,0usize);MAX_D];let mut shape_count=0;
    for r in 2..=n/3 {if n%r==0 {shapes[shape_count]=(r,n/r);shape_count+=1;}}
    if shape_count==0{return None;}
    let (r,k)=shapes[rng.below(shape_count)];let beta=if rng.below(2)==0{dd_full_cycle(r,rng)}else{gd_shuffle(r,rng)};
    if beta.is_identity(r){return None;}
    let gamma=gd_shuffle(r,rng);let mut local=vec![Perm::identity(k);r];
    let mut parent=std::array::from_fn(|i|i);
    for j in 0..r {for x in 0..k {gd_join(&mut parent,j*k+x,beta.img[j] as usize*k+x);}}
    let target=2+rng.below(k-2);let mut active=0u64;
    for _ in 0..target {
        // In a disjoint union of r blocks, r*C(k,2) <= C(r*k,2).
        // Find each DSU root once rather than twice for every point pair.
        let roots:[usize;MAX_D]=std::array::from_fn(|i|
            if i<n {gd_find(&mut parent,i)} else {i});
        let mut choices=[(0u8,0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
        for j in 0..r {for x in 0..k {for y in x+1..k {
            if roots[j*k+x]!=roots[j*k+y]{choices[count]=(j as u8,x as u8,y as u8);count+=1;}
        }}}
        if count==0{return None;}
        let (j,x,y)=choices[rng.below(count)];let (j,x,y)=(j as usize,x as usize,y as usize);
        local[j].img.swap(x,y);active|=1u64<<j;
        assert!(gd_join(&mut parent,j*k+x,j*k+y));
    }
    if active.count_ones()<2{return None;}
    let frame=ReturnForestFrame::new(beta,gamma,r,k,local)?;let f=&frame.tree;
    let free:Vec<_>=(0..f.free_edges.len()).map(|_|gd_shuffle(k,rng)).collect();
    let arm=rng.below(3);let targets:Vec<_>=(0..f.orders.len()).map(|v|{
        if v==f.root||arm==0||rng.below(2)==0 {Perm::identity(k)}
        else if arm==1 {let mut p=Perm::identity(k);let x=rng.below(k);let y=(x+1+rng.below(k-1))%k;p.img.swap(x,y);p}
        else {gd_shuffle(k,rng)}
    }).collect();
    let tuple=frame.solve(&free,&targets);
    let points=(0..n).filter(|&x|tuple[0].img[x] as usize!=x).collect();
    debug_assert_eq!(ff_product_index(tuple[0],tuple[1],n),sm_index(tuple[0],n)+sm_index(tuple[1],n));
    Some(MechanismProposal{tuple,source:r,arm,points})
}

#[cfg(test)]
mod structural_novelty_tests {
    use super::*;
    // Deliberately retains the old allocations and repeated DSU finds. This
    // oracle checks proposal-law/RNG preservation, not just smooth survivors.
    fn literal_return(n:usize,rng:&mut SplitMix64)->Option<MechanismProposal> {
        let shapes:Vec<_>=(2..=n/3).filter(|r|n%r==0).map(|r|(r,n/r)).collect();
        if shapes.is_empty(){return None;}
        let (r,k)=shapes[rng.below(shapes.len())];
        let beta=if rng.below(2)==0{dd_full_cycle(r,rng)}else{gd_shuffle(r,rng)};
        if beta.is_identity(r){return None;}
        let gamma=gd_shuffle(r,rng);let mut local=vec![Perm::identity(k);r];
        let mut parent=std::array::from_fn(|i|i);
        for j in 0..r {for x in 0..k {gd_join(&mut parent,j*k+x,beta.img[j] as usize*k+x);}}
        let target=2+rng.below(k-2);let mut active=HashSet::new();
        for _ in 0..target {
            let mut choices=Vec::new();
            for j in 0..r {for x in 0..k {for y in x+1..k {
                if gd_find(&mut parent,j*k+x)!=gd_find(&mut parent,j*k+y){choices.push((j,x,y));}
            }}}
            if choices.is_empty(){return None;}
            let (j,x,y)=choices[rng.below(choices.len())];local[j].img.swap(x,y);active.insert(j);
            assert!(gd_join(&mut parent,j*k+x,j*k+y));
        }
        if active.len()<2{return None;}
        let frame=ReturnForestFrame::new(beta,gamma,r,k,local)?;let f=&frame.tree;
        let free:Vec<_>=(0..f.free_edges.len()).map(|_|gd_shuffle(k,rng)).collect();
        let arm=rng.below(3);let targets:Vec<_>=(0..f.orders.len()).map(|v|{
            if v==f.root||arm==0||rng.below(2)==0 {Perm::identity(k)}
            else if arm==1 {let mut p=Perm::identity(k);let x=rng.below(k);let y=(x+1+rng.below(k-1))%k;p.img.swap(x,y);p}
            else {gd_shuffle(k,rng)}
        }).collect();
        let tuple=frame.solve(&free,&targets);
        let points=(0..n).filter(|&x|tuple[0].img[x] as usize!=x).collect();
        Some(MechanismProposal{tuple,source:r,arm,points})
    }
    #[test]fn compact_returns_preserve_every_draw_and_rng_through_degree_25() {
        for n in 9..=25 {
            let mut old=SplitMix64::new(7210+n as u64);let mut new=SplitMix64::new(7210+n as u64);
            for _ in 0..1000 {
                let a=literal_return(n,&mut old);let b=sn_return_tuple(n,&mut new);
                assert_eq!(a.as_ref().map(|p|(&p.tuple,p.source,p.arm,&p.points)),
                           b.as_ref().map(|p|(&p.tuple,p.source,p.arm,&p.points)));
                assert_eq!(old.state,new.state);
            }
        }
    }
    #[test]fn compact_joint_choices_preserve_every_draw_and_rng() {
        fn literal(f:&ForestFrame,n:usize,target:usize,rng:&mut SplitMix64,strong:bool)->Option<Perm> {
            let mut c=Perm::identity(n);
            let mut products=if strong {f.factors.to_vec()} else {f.pair.to_vec()};
            for _ in 0..target {
                let labels:Vec<_>=std::iter::once(c).chain(products.iter().copied())
                    .map(|p|dd_cycle_labels(p,n)).collect();
                let mut choices=Vec::new();
                for x in 0..n {for y in x+1..n {
                    if labels.iter().all(|v|v[x]!=v[y]) {choices.push((x,y));}
                }}
                if choices.is_empty(){return None;}
                let (x,y)=choices[rng.below(choices.len())];c.img.swap(x,y);
                for p in &mut products {p.img.swap(x,y);}
            }
            Some(c)
        }
        let mut pairs=SplitMix64::new(20261021);
        for n in 9..=25 {for strong in [false,true] {for target in 1..n {
            let f=ForestFrame::new(gd_forest(n,&mut pairs),n);
            let seed=pairs.next_u64();let mut old=SplitMix64::new(seed);let mut new=SplitMix64::new(seed);
            assert_eq!(literal(&f,n,target,&mut old,strong),sn_joint_third(&f,n,target,&mut new,strong));
            assert_eq!(old.state,new.state);
        }}}
    }
    #[test]fn compact_choice_arrays_bound_all_supported_shapes() {
        assert!(MAX_D<=u8::MAX as usize);
        for n in 2..=MAX_D {for r in 2..=n {
            if n%r==0 {let k=n/r;assert!(r*k*(k-1)/2<=DD_PAIR_CAPACITY);}
        }}
    }
    #[test]fn triangle_draws_have_the_declared_geodesies_and_budget(){
        let mut rng=SplitMix64::new(202610133);let mut found=0;
        for n in 11..=18 {for _ in 0..200 {
            if let Some(p)=sn_triangle_tuple(n,&mut rng) {
                found+=1;assert!(sn_triangle_budget(&p.tuple,n));
                for (i,j) in [(0,1),(0,2),(1,2)] {
                    assert_eq!(ff_product_index(p.tuple[i],p.tuple[j],n),sm_index(p.tuple[i],n)+sm_index(p.tuple[j],n));
                }
                assert!(tuple_product_is_identity(n,&p.tuple));
            }
        }}assert!(found>0);
    }
    #[test]fn triangle16_forest_pairs_are_noninvolutive_and_additive(){
        let mut rng=SplitMix64::new(2026103116);
        for _ in 0..64 {
            let [a,b]=sn_triangle16_forest_pair(&mut rng).expect("degree-16 cycle forest exists");
            assert_eq!(sm_index(a,16),7);assert_eq!(sm_index(b,16),7);
            assert_eq!(ff_product_index(a,b,16),14);
            assert!(!sn_is_involution(a,16));assert!(!sn_is_involution(b,16));
            assert!(!sn_three_two_bulk_with_three_fixed(a,16));
            assert!(!sn_three_two_bulk_with_three_fixed(b,16));
        }
    }
    #[test]fn prescribed_nonroot_return_products_hold_with_multiblock_local_action(){
        let mut rng=SplitMix64::new(28445);
        let beta=dd_full_cycle(3,&mut rng);let gamma=beta;
        let mut x=Perm::identity(4);x.img.swap(0,1);let mut y=Perm::identity(4);y.img.swap(1,2);
        let f=ReturnForestFrame::new(beta,gamma,3,4,vec![x,y,Perm::identity(4)]).unwrap();
        for _ in 0..100 {
            let free:Vec<_>=(0..f.tree.free_edges.len()).map(|_|gd_shuffle(4,&mut rng)).collect();
            let targets:Vec<_>=(0..f.tree.orders.len()).map(|_|gd_shuffle(4,&mut rng)).collect();
            let t=f.solve(&free,&targets);let abc=t[0].compose(t[1],12).compose(t[2],12);
            let ng=qf_cycles(gamma,3).len();
            for (v,cycle) in f.tree.orders.iter().enumerate() {
                if v==f.tree.root{continue;}
                let p=if v<ng{t[2]}else{abc};
                let mut actual=Perm::identity(4);let block=cycle[0];
                for a in 0..4 {let mut point=block*4+a;for _ in cycle{point=p.img[point] as usize;}assert_eq!(point/4,block);actual.img[a]=(point%4) as u8;}
                assert_eq!(actual,targets[v]);
            }
            assert_eq!(ff_product_index(t[0],t[1],12),sm_index(t[0],12)+sm_index(t[1],12));
        }
    }
    #[test]fn rational_quotient_is_supported_not_silently_filtered(){
        let mut beta=Perm::identity(2);beta.img.swap(0,1);
        assert_eq!(QuotientForest::tree_any(beta,Perm::identity(2),2,6,0,true).unwrap().genus,0);
    }
}
