// Discovery only: failures have literal witnesses; resource caps are undecided.
// Numerical genus and quotient hypotheses never substitute for braid closure.
#[path="smooth_joint_completion.rs"]
mod joint_completion;
include!("smooth_cyclic_screen.rs");
include!("smooth_saturated_passport.rs");
include!("genus_three_frame_retirement.rs");

// P[s][b] counts partitions of s with largest part <=b. Subtracting the
// minimum part and conjugating partitions counts all lexicographic suffixes.
// This replaces exponential materialization without changing type weights,
// lexicographic ranks, or the old sampler's random-number consumption.
fn g3_partition_counts()->&'static [[usize;MAX_D+1];MAX_D+1] {
    static TABLE:std::sync::OnceLock<[[usize;MAX_D+1];MAX_D+1]>=std::sync::OnceLock::new();
    TABLE.get_or_init(||{
        let mut p=[[0usize;MAX_D+1];MAX_D+1];
        for b in 0..=MAX_D{p[0][b]=1;}
        for b in 1..=MAX_D{for s in 1..=MAX_D{
            p[s][b]=p[s][b-1]+if s>=b{p[s-b][b]}else{0};
        }}p
    })
}
fn g3_type_rank(n:usize,index:usize,mut rank:usize)->Vec<usize> {
    assert!(n<=MAX_D&&index<n);let p=g3_partition_counts();let mut remaining=n;
    let mut parts=n-index;let mut minimum=1;let mut out=Vec::with_capacity(parts);
    assert!(rank<p[index][parts]);
    while parts>0{
        let mut selected=None;
        for first in minimum..=remaining/parts{
            let suffix=remaining-parts*first;
            let count=if parts==1{usize::from(suffix==0)}else{p[suffix][parts-1]};
            if rank<count{selected=Some(first);break;}rank-=count;
        }
        let first=selected.expect("partition rank exhausted");out.push(first);
        remaining-=first;parts-=1;minimum=first;
    }
    assert_eq!(remaining,0);assert_eq!(rank,0);out
}
// Cache small, actually used index buckets in RAM. Large buckets stay ranked
// and are never materialized wholesale. This keeps old hot draws allocation-
// free without restoring the exponential all-partition startup at large n.
struct G3Types { n:usize,cached:Vec<std::cell::OnceCell<Vec<Vec<usize>>>> }
impl G3Types {
    fn new(n:usize)->Self {assert!(n<=MAX_D);Self{n,cached:(0..n).map(|_|std::cell::OnceCell::new()).collect()}}
    fn permutation(&self,n:usize,index:usize,rng:&mut SplitMix64)->Perm {
        assert_eq!(n,self.n);assert!(index<n);
        let count=g3_partition_counts()[index][n-index];
        let rank=rng.below(count);
        if count<=4096{
            let bucket=self.cached[index].get_or_init(||(0..count).map(|r|g3_type_rank(n,index,r)).collect());
            g3_type_perm(n,&bucket[rank],rng)
        }else{g3_type_perm(n,&g3_type_rank(n,index,rank),rng)}
    }
}
fn g3_type_perm(n:usize,shape:&[usize],rng:&mut SplitMix64)->Perm {
    let mut labels:[usize;MAX_D]=std::array::from_fn(|i|i);rng.shuffle(&mut labels[..n]);
    let mut p=Perm::identity(n);let mut offset=0;
    for &len in shape {
        for j in 0..len-1 {p.img[labels[offset+j]]=labels[offset+j+1] as u8;}
        p.img[labels[offset+len-1]]=labels[offset] as u8;offset+=len;
    }
    assert_eq!(offset,n);p
}
#[cfg(test)]
mod g3_partition_tests {
    use super::*;
    #[test]fn counted_unranking_is_identical_to_old_materialized_type_order(){
        for n in 1..=34.min(MAX_D){
            let mut parts=Vec::new();sc_partitions(n,1,&mut Vec::new(),&mut parts);
            for index in 0..n{
                let old:Vec<_>=parts.iter().filter(|p|p.len()==n-index).collect();
                assert_eq!(g3_partition_counts()[index][n-index],old.len());
                for (rank,shape)in old.iter().enumerate(){assert_eq!(&g3_type_rank(n,index,rank),*shape);}
                let mut before=SplitMix64::new(99177+index as u64);let mut after=before.clone();
                for _ in 0..10{
                    let shape=old[before.below(old.len())];let p=g3_type_perm(n,shape,&mut before);
                    assert_eq!(G3Types::new(n).permutation(n,index,&mut after),p);
                    assert_eq!(before.next_u64(),after.next_u64());
                }
            }
        }
    }
    #[test]fn high_degree_type_sampler_preserves_degree_and_index(){
        let mut rng=SplitMix64::new(87163);
        for n in 1..=MAX_D{let types=G3Types::new(n);for index in 0..n{for _ in 0..10{
            let p=types.permutation(n,index,&mut rng);assert_eq!(sm_index(p,n),index);
            let mut labels=p.img[..n].to_vec();labels.sort_unstable();
            assert_eq!(labels,(0..n as u8).collect::<Vec<_>>());
        }}}
    }
}
fn g3_passports(n:usize)->Vec<[usize;4]> {
    let mut out=Vec::new();
    for a in 1..n {for b in a..n {for c in b..n {
        let total=2*n+4;if a+b+c>=total {continue;}
        let d=total-a-b-c;
        if c<=d && d<n && a+c<n && (a+d<n || b+c<n) {out.push([a,b,c,d]);}
    }}}
    out
}
fn g3_index_profile(n:usize,text:&str)->Result<[usize;4],String>{
    let v=text.split(',').map(|x|x.parse::<usize>().map_err(|_|"index profile integer".to_string())).collect::<Result<Vec<_>,_>>()?;
    let p:[usize;4]=v.try_into().map_err(|_|"four sorted indices required".to_string())?;
    if !g3_passports(n).contains(&p){return Err("index profile outside genus3 necessary budgets".into());}Ok(p)
}
fn g3_involution_triangle_word(t:&[Perm;4],n:usize)->Option<Vec<usize>>{
    if n<5||n%4!=1||sm_index(t[3],n)!=n-1{return None;}
    let mut fixed=[0;3];
    for j in 0..3 {let mut count=0;for x in 0..n{
        if t[j].img[t[j].img[x]as usize]as usize!=x{return None;}
        if t[j].img[x]as usize==x{fixed[j]=x;count+=1;}
    }if count!=1{return None;}}
    let h=t[0].compose(t[1],n);if sm_index(h,n)!=n-1{return None;}
    let mut k=0;let mut x=fixed[0];while x!=fixed[2]{x=h.img[x]as usize;k+=1;assert!(k<n);}
    let mut word=if k<=n-k{vec![0;2*k]}else{vec![3;2*(n-k)]};word.push(4);
    assert!(ff_replay_negative(*t,n,&word));Some(word)
}
// Keep h fixed while splitting its cycles. a <- a*(xy), b <- (xy)*b.
// Different a-cycles and the same b-cycle preserve absolute-order geodesy.
fn g3_factor(h:Perm,n:usize,index:usize,rng:&mut SplitMix64)->Option<[Perm;2]> {
    let mut a=Perm::identity(n);let mut b=h;
    for _ in 0..index {
        let la=dd_cycle_labels(a,n);let lb=dd_cycle_labels(b,n);
        let mut choices=[(0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
        for x in 0..n {for y in x+1..n {if la[x]!=la[y] && lb[x]==lb[y] {
            choices[count]=(x as u8,y as u8);count+=1;
        }}}
        if count==0 {return None;}
        let (x,y)=choices[rng.below(count)];let (x,y)=(x as usize,y as usize);
        a.img.swap(x,y);
        for v in &mut b.img[..n] {if *v as usize==x {*v=y as u8;}else if *v as usize==y {*v=x as u8;}}
    }
    debug_assert_eq!(a.compose(b,n),h);
    debug_assert_eq!(sm_index(a,n)+sm_index(b,n),sm_index(h,n));
    Some([a,b])
}
struct G3Collision {
    n:usize,types:G3Types,passports:Vec<[usize;4]>,frame:Option<ForestFrame>,
    passport:[usize;4],used:usize,balanced:bool,retirement:G3FrameRetirement,
}
impl G3Collision {
    fn new(n:usize,balanced:bool)->Self {
        Self{n,types:G3Types::new(n),passports:g3_passports(n),frame:None,
             passport:[0;4],used:128,balanced,retirement:G3FrameRetirement::default()}
    }
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]> {
        if self.used>=128 || self.frame.is_none() {
            if self.passports.is_empty(){return None;}
            let ix=self.passports[rng.below(self.passports.len())];self.passport=ix;
            let e=ix[0]+ix[1];let h=if self.balanced {
                self.types.permutation(self.n,e,rng)
            }else{ff_index_perm(self.n,e,rng)};
            let pair=g3_factor(h,self.n,ix[0],rng)?;
            if self.retirement.retire(pair,ix,self.n){self.frame=None;self.used=128;return None;}
            self.frame=Some(ForestFrame::new(pair,self.n));self.used=0;
        }
        self.used+=1;let frame=self.frame.as_ref().unwrap();let n=self.n;
        let c=if self.balanced {self.types.permutation(n,self.passport[2],rng)}
              else {ff_index_perm(n,self.passport[2],rng)};
        if ff_product_index(frame.p,c,n)!=self.passport[3] {return None;}
        Some(frame.tuple(c,n))
    }
}

// Construct the third permutation while prescribing the indices of all five
// original short-word products AND (a*b)*c. Each right transposition joins or
// splits one cycle in each product, so the remaining split quotas are exact.
// No blocks are assumed. Dead ends are failed draws, not exclusions.
fn g3_joint_third(frame:&ForestFrame,n:usize,target:usize,fourth:usize,rng:&mut SplitMix64)->Option<Perm> {
    let[ia,ib]=frame.indices;let mut products=[Perm::identity(n);7];
    products[1..6].copy_from_slice(&frame.factors);products[6]=frame.p;
    let left=(ia+target).min(ib+fourth);let right=(ib+target).min(ia+fourth);
    if left>=n||right>=n{return None;}
    let goals=[target,left,left,right,right,right,fourth];let mut splits=[0usize;7];
    for j in 0..7 {
        let total=sm_index(products[j],n)+target;
        if total<goals[j]||(total-goals[j])%2!=0{return None;}
        splits[j]=(total-goals[j])/2;if splits[j]>target{return None;}
    }
    for step in 0..target {
        let remaining=target-step;let labels=products.map(|p|dd_cycle_labels(p,n));
        let mut choices=[(0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
        for x in 0..n {for y in x+1..n {
            if (0..7).all(|j|{
                let split=labels[j][x]==labels[j][y];
                !(split&&splits[j]==0||!split&&splits[j]==remaining)
            }) {choices[count]=(x as u8,y as u8);count+=1;}
        }}
        if count==0{return None;}
        let(x,y)=choices[rng.below(count)];let(x,y)=(x as usize,y as usize);
        for j in 0..7 {if labels[j][x]==labels[j][y]{splits[j]-=1;}products[j].img.swap(x,y);}
    }
    assert_eq!(splits,[0;7]);debug_assert_eq!(products.map(|p|sm_index(p,n)),goals);
    Some(products[0])
}
struct G3Joint {base:G3Collision,csp:bool,budget:u64,nodes:u64,found:u64,exhausted:u64,capped:u64}
impl G3Joint {
    fn new(n:usize)->Self {
        Self{base:G3Collision::new(n,true),csp:false,budget:5000,nodes:0,found:0,exhausted:0,capped:0}
    }
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]> {
        let b=&mut self.base;
        if b.frame.is_none()||b.used>=32 {
            if b.passports.is_empty(){return None;}
            b.passport=b.passports[rng.below(b.passports.len())];let ix=b.passport;
            let h=b.types.permutation(b.n,ix[0]+ix[1],rng);
            let pair=g3_factor(h,b.n,ix[0],rng)?;
            if b.retirement.retire(pair,ix,b.n){b.frame=None;b.used=32;return None;}
            b.frame=Some(ForestFrame::new(pair,b.n));b.used=0;
        }
        b.used+=1;let f=b.frame.as_ref().unwrap();
        let c=if self.csp {
            let [ia,ib,ic,id]=b.passport;
            let left=(ia+ic).min(ib+id);let right=(ib+ic).min(ia+id);
            if left>=b.n||right>=b.n{return None;}
            let mut products=[Perm::identity(b.n);7];products[1..6].copy_from_slice(&f.factors);products[6]=f.p;
            let images=products.map(|p|std::array::from_fn(|i|if i<b.n {p.img[i]}else{i as u8}));
            let answer=joint_completion::solve(b.n,&images,[ic,left,left,right,right,right,id],rng.next_u64(),self.budget);
            self.nodes+=answer.nodes;
            if let Some(img)=answer.image{self.found+=1;let mut c=Perm::identity(b.n);c.img[..b.n].copy_from_slice(&img[..b.n]);c}
            else {if answer.exhausted{self.exhausted+=1;b.used=32;}else{self.capped+=1;}return None;}
        }else{g3_joint_third(f,b.n,b.passport[2],b.passport[3],rng)?};
        debug_assert!(f.short_word(c,b.passport[2],b.passport[3],b.n).is_none());
        Some(f.tuple(c,b.n))
    }
}

struct G3Elliptic {n:usize,quotient_genus:usize,frame:Option<ReturnForestFrame>,first_frame:Option<ForestFrame>,last_indices:[usize;4],used:usize,minimum_moving:usize,local_types:Vec<G3Types>}
impl G3Elliptic {
    fn new(n:usize)->Self {Self::with_genus(n,1)}
    fn with_genus(n:usize,quotient_genus:usize)->Self {
        assert!(quotient_genus<=1);
        Self{n,quotient_genus,frame:None,first_frame:None,last_indices:[0;4],used:128,minimum_moving:1,local_types:(0..=n/2).map(G3Types::new).collect()}
    }
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]> {
        let n=self.n;
        if self.used>=128 || self.frame.is_none() {
            let shapes:Vec<_>=(if self.quotient_genus==1{3}else{2}..=n/3).filter(|&r|n%r==0).map(|r|(r,n/r)).collect();
            if shapes.is_empty(){return None;}
            let (r,k)=shapes[rng.below(shapes.len())];
            // A full beta reduces tree cycle-rank to the two elliptic handles.
            let mut beta=Perm::identity(r);for j in 0..r {beta.img[j]=((j+1)%r) as u8;}
            let gamma=gd_shuffle(r,rng);let Some(tree)=QuotientForest::tree_any(beta,gamma,r,k,0,true) else{return None;};
            if tree.genus!=self.quotient_genus{return None;}
            let mut local=vec![Perm::identity(k);r];let mut parent=std::array::from_fn(|i|i);
            for j in 0..r {for x in 0..k {gd_join(&mut parent,j*k+x,beta.img[j] as usize*k+x);}}
            let maximum=if self.quotient_genus==1 {3.min(k-1)} else {k-1};
            if maximum<self.minimum_moving{return None;}
            let moving=self.minimum_moving+rng.below(maximum-self.minimum_moving+1);
            for _ in 0..moving {
                let mut choices=[(0u8,0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
                let roots:[usize;MAX_D]=std::array::from_fn(|i|if i<n {gd_find(&mut parent,i)} else {i});
                for j in 0..r {for x in 0..k {for y in x+1..k {if roots[j*k+x]!=roots[j*k+y] {
                    choices[count]=(j as u8,x as u8,y as u8);count+=1;
                }}}}
                if count==0{return None;}
                let (j,x,y)=choices[rng.below(count)];let (j,x,y)=(j as usize,x as usize,y as usize);
                local[j].img.swap(x,y);assert!(gd_join(&mut parent,j*k+x,j*k+y));
            }
            self.frame=ReturnForestFrame::new(beta,gamma,r,k,local);self.first_frame=None;self.used=0;
        }
        self.used+=1;let frame=self.frame.as_ref()?;let f=&frame.tree;
        let moving:usize=frame.local.iter().map(|&p|sm_index(p,f.k)).sum();
        // Distribute a small exact budget among non-root returns; the root is
        // solved, not guessed. Final excess=4 is checked on the literal tuple.
        let desired_excess=if self.quotient_genus==1 {4}else{4+2*f.k};
        let mut budget=rng.below(desired_excess+1-moving);let mut targets=vec![Perm::identity(f.k);f.orders.len()];
        while budget>0 {
            let v=rng.below(f.orders.len()-1);let v=if v>=f.root {v+1}else{v};
            if !targets[v].is_identity(f.k){break;}
            let index=1+rng.below(budget.min(f.k-1));
            targets[v]=self.local_types[f.k].permutation(f.k,index,rng);budget-=index;
        }
        let free:Vec<_>=(0..f.free_edges.len()).map(|_|{
            if rng.below(3)==0 {gd_shuffle(f.k,rng)} else {
                self.local_types[f.k].permutation(f.k,1+rng.below(3.min(f.k-1)),rng)
            }
        }).collect();
        let t=frame.solve(&free,&targets);
        self.last_indices=t.map(|p|sm_index(p,n));let source_index_sum:usize=self.last_indices.iter().sum();
        if source_index_sum!=2*n+4 || (self.quotient_genus==1 && moving>=4){return None;}
        if self.first_frame.is_none(){self.first_frame=Some(ForestFrame::new([t[0],t[1]],n));}
        Some(t)
    }
}

// Every point of one elliptic fibre is a relative branch value. This deliberately
// removes the former moving-index<=3 restriction: total moving index is four,
// distributed over three/four blocks, with zero ordinary relative excess.
struct G3FullFiber {
    n:usize,frame:Option<ReturnForestFrame>,first_frame:Option<ForestFrame>,
    used:usize,types:Vec<G3Types>,last_indices:[usize;4],solved:bool,partial:bool,skip_covered:bool,
}
// All solutions of v*X*v^-1=Y map cycles of X to equal-length cycles of Y,
// choosing an arbitrary matching and cyclic offset. No factorial guess/reject.
fn g3_conjugator(x:Perm,y:Perm,k:usize,rng:&mut SplitMix64)->Option<Perm> {
    fn cycles(p:Perm,k:usize)->([(usize,usize);MAX_D],usize) {
        let mut result=[(0,0);MAX_D];let mut used=0u64;let mut count=0;
        for root in 0..k {if used&(1<<root)==0 {
            let mut len=0;let mut at=root;loop{used|=1<<at;len+=1;at=p.img[at]as usize;if at==root{break;}}
            result[count]=(len,root);count+=1;
        }}result[..count].sort_by_key(|v|v.0);(result,count)
    }
    let(mut xc,nx)=cycles(x,k);let(mut yc,ny)=cycles(y,k);
    let xc=&mut xc[..nx];let yc=&mut yc[..ny];
    if nx!=ny||xc.iter().map(|c|c.0).ne(yc.iter().map(|c|c.0)){return None;}
    let mut start=0;
    while start<yc.len(){let mut end=start+1;while end<yc.len()&&yc[end].0==yc[start].0{end+=1;}
        rng.shuffle(&mut yc[start..end]);start=end;}
    let mut v=Perm::identity(k);
    for (&(len,mut a),&(_,mut b)) in xc.iter().zip(yc.iter()) {
        for _ in 0..rng.below(len){b=y.img[b]as usize;}
        for _ in 0..len{v.img[a]=b as u8;a=x.img[a]as usize;b=y.img[b]as usize;}
    }
    debug_assert_eq!(v.compose(x,k).compose(v.inverse(k),k),y);Some(v)
}
fn g3_fiber_conjugacy(frame:&ReturnForestFrame,u:Perm)->(Perm,Perm) {
    let a=&frame.local;let k=frame.tree.k;let ui=u.inverse(k);
    if frame.tree.r==3 {
        (a[2].compose(ui,k),a[0].compose(u,k).compose(a[1],k).inverse(k))
    }else if sm_index(frame.tree.gamma,4)==3 {
        (a[3].compose(ui,k),a[0].compose(u,k).compose(a[2],k).compose(a[1],k).inverse(k))
    }else{
        (a[2].compose(ui,k),a[0].compose(a[3],k).compose(u,k).compose(a[1],k).inverse(k))
    }
}
impl G3FullFiber {
    fn new(n:usize)->Self {
        Self{n,frame:None,first_frame:None,used:128,
             types:(0..=n/3).map(G3Types::new).collect(),last_indices:[0;4],solved:true,partial:false,skip_covered:false}
    }
    fn partial(n:usize)->Self {let mut f=Self::new(n);f.partial=true;f}
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]> {
        let n=self.n;
        if self.frame.is_none()||self.used>=128 {
            let shapes:Vec<_>=[3,4].into_iter().filter(|&r|n%r==0&&n/r>=5&&
                !(self.skip_covered&&n/r<=6&&(!self.partial||r==3))).map(|r|(r,n/r)).collect();
            if shapes.is_empty(){return None;}
            let (r,k)=shapes[rng.below(shapes.len())];let mut beta=Perm::identity(r);
            for j in 0..r{beta.img[j]=((j+1)%r)as u8;}
            let mut gamma=beta;
            // The other degree-four elliptic quotient has S4 rather than C4
            // monodromy. Put a full cycle first to retain only two free edges.
            if r==4&&rng.below(2)==1{gamma=Perm::identity(r);gamma.img[0]=1;gamma.img[1]=2;gamma.img[2]=0;}
            let indices=if self.partial {
                // Missing middle population: all excess four occurs over the
                // moving fibre, but only two/three of its points are active.
                let mut ix=match rng.below(if r==3{2}else{3}) {
                    0=>vec![3,1],1=>vec![2,2],_=>vec![2,1,1],
                };ix.resize(r,0);rng.shuffle(&mut ix);ix
            }else{
                let double=if r==3{rng.below(r)}else{usize::MAX};
                (0..r).map(|j|if j==double{2}else{1}).collect()
            };
            let local=indices.into_iter().map(|index|self.types[k].permutation(k,index,rng)).collect();
            self.frame=ReturnForestFrame::new(beta,gamma,r,k,local);self.first_frame=None;self.used=0;
        }
        self.used+=1;let frame=self.frame.as_ref()?;let f=&frame.tree;
        assert_eq!(f.genus,1);assert_eq!(f.free_edges.len(),2);
        let u=gd_shuffle(f.k,rng);
        let v=if self.solved{let(x,y)=g3_fiber_conjugacy(frame,u);g3_conjugator(x,y,f.k,rng)?}else{gd_shuffle(f.k,rng)};
        let free=[u,v];
        let t=frame.solve(&free,&vec![Perm::identity(f.k);f.orders.len()]);
        let delta=f.beta.compose(f.gamma,f.r);
        let ix=[4,f.k*sm_index(f.beta,f.r),f.k*sm_index(f.gamma,f.r),f.k*sm_index(delta,f.r)];
        // Every relative excess is nonnegative; the leaf equations fix all but
        // one return. Matching both ordinary indices forces the last return=1.
        if !self.solved && (sm_index(t[2],n)!=ix[2]||sm_index(t[3],n)!=ix[3]){return None;}
        assert_eq!(sm_index(t[2],n),ix[2]);assert_eq!(sm_index(t[3],n),ix[3]);
        self.last_indices=ix;
        if self.first_frame.is_none(){self.first_frame=Some(ForestFrame::new([t[0],t[1]],n));}
        debug_assert_eq!(t.map(|p|sm_index(p,n)),ix);
        Some(t)
    }
}

// Finite gauge-reduced degree-five/six inner actions. No random tuple enumeration:
// normalize one local conjugacy class and enumerate the two remaining voltage
// variables. A complete parameter interval is not independent coverage until
// separately audited; every qualifying early rejection is persisted below.
struct G3FiberExact {
    frames:Vec<ReturnForestFrame>,perms:Vec<Perm>,cursor:usize,start:usize,total:usize,
    first_frame:Option<ForestFrame>,current:usize,last_indices:[usize;4],kind:String,n:usize,k:usize,
    outside_transitive:usize,
}
fn g3_local_key(local:&[Perm],perms:&[Perm],k:usize)->Vec<u8> {
    perms.iter().map(|&h|local.iter().flat_map(|&p|p.conjugate_by(h,k).img[..k].to_vec()).collect::<Vec<_>>()).min().unwrap()
}
impl G3FiberExact {
    fn new(n:usize,kind:&str,cursor:usize)->Result<Self,String> {
        let r=match kind{"C3"|"C3P31"|"C3P22"=>3,"C4"|"S4"=>4,_=>return Err("fiber kind C3/C3P31/C3P22/C4/S4".into())};
        if n%r!=0{return Err("finite fiber degree is not divisible by block count".into());}
        let k=n/r;if !(5..=6).contains(&k){return Err("finite fiber audit supports inner degrees five and six".into());}
        let perms:Vec<_>=sc_perms(k).into_iter().map(|v|v.0).collect();
        let trans:Vec<_>=perms.iter().copied().filter(|&p|sm_index(p,k)==1).collect();
        let mut beta=Perm::identity(r);for j in 0..r{beta.img[j]=((j+1)%r)as u8;}
        let gamma=if kind=="S4"{let mut p=Perm::identity(r);p.img[0]=1;p.img[1]=2;p.img[2]=0;p}else{beta};
        let mut raw=Vec::new();
        if kind=="C3P31"||kind=="C3P22" {
            let(first_index,second_index)=if kind=="C3P31"{(3,1)}else{(2,2)};
            let mut types=HashSet::new();
            for &first in &perms {if sm_index(first,k)!=first_index||!types.insert(sp_type(first,k)){continue;}
                for &second in &perms {if sm_index(second,k)!=second_index{continue;}
                    // Both cyclic orientations are retained; no unproved
                    // reflection of the marked quotient removes parameters.
                    raw.push(vec![first,second,Perm::identity(k)]);
                    raw.push(vec![first,Perm::identity(k),second]);
                }
            }
        }else if r==3 {
            let mut a=Perm::identity(k);a.img[0]=1;a.img[1]=2;a.img[2]=0;
            let mut b=Perm::identity(k);b.img.swap(0,1);b.img.swap(2,3);
            for first in [a,b]{for &x in &trans{for &y in &trans{raw.push(vec![first,x,y]);}}}
        }else{
            let mut first=Perm::identity(k);first.img.swap(0,1);
            for &x in &trans{for &y in &trans{for &z in &trans{raw.push(vec![first,x,y,z]);}}}
        }
        let mut seen=HashSet::new();let mut frames=Vec::new();
        for local in raw {
            if let Some(f)=ReturnForestFrame::new(beta,gamma,r,k,local.clone()){
                let key=g3_local_key(&local,&perms,k);
                if seen.insert(key.clone()){frames.push((key,f));}
            }
        }
        frames.sort_by(|a,b|a.0.cmp(&b.0));let frames:Vec<_>=frames.into_iter().map(|v|v.1).collect();
        let total=frames.len()*perms.len()*perms.len();if cursor>total{return Err("fiber cursor beyond total".into());}
        Ok(Self{frames,perms,cursor,start:cursor,total,first_frame:None,current:usize::MAX,
                last_indices:[0;4],kind:kind.into(),n,k,outside_transitive:0})
    }
    fn draw(&mut self)->Option<[Perm;4]> {
        if self.cursor==self.total{return None;}
        let size=self.perms.len();let fi=self.cursor/(size*size);let rank=self.cursor%(size*size);self.cursor+=1;
        if fi!=self.current{self.current=fi;self.first_frame=None;}
        let f=&self.frames[fi];let tree=&f.tree;
        let free=[self.perms[rank/size],self.perms[rank%size]];
        let t=f.solve(&free,&vec![Perm::identity(self.k);tree.orders.len()]);
        let ix=[4,self.k*sm_index(tree.beta,tree.r),self.k*sm_index(tree.gamma,tree.r),self.k*sm_index(tree.beta.compose(tree.gamma,tree.r),tree.r)];
        if sm_index(t[2],self.n)!=ix[2]||sm_index(t[3],self.n)!=ix[3]{return None;}
        // At inner degree five the additive first product is a full cycle.
        // This no longer holds at six: reject disconnected tuples explicitly,
        // and independently reconstruct this scope count in the coverage audit.
        if !generators_transitive(self.n,&t){self.outside_transitive+=1;return None;}
        self.last_indices=ix;
        if self.first_frame.is_none(){self.first_frame=Some(ForestFrame::new([t[0],t[1]],self.n));}
        Some(t)
    }
    fn save(&self,out:&Path)->Result<(),String> {
        fs::write(out.join("FIBER_COVERAGE.json"),format!("{{\"degree\":{},\"inner_degree\":{},\"kind\":\"{}\",\"local_frames\":{},\"internal_permutations\":{},\"start\":{},\"next\":{},\"total\":{},\"producer_complete\":{},\"outside_transitive\":{},\"scope\":\"normalized marked full-fiber parameter family, not a degree census\"}}\n",
            self.n,self.k,self.kind,self.frames.len(),self.perms.len(),self.start,self.cursor,self.total,self.cursor==self.total,self.outside_transitive)).map_err(|e|e.to_string())
    }
}
fn g3_early_negative(leads:&mut BufWriter<File>,resolutions:&mut BufWriter<File>,t:&[Perm;4],n:usize,id:u64,word:&[usize])->Result<(),String> {
    assert!(ff_replay_negative(*t,n,word));let text=tuple_string_dynamic(n,t);
    writeln!(leads,"{id}\t{text}").map_err(|e|e.to_string())?;leads.flush().map_err(|e|e.to_string())?;
    writeln!(resolutions,"{id}\tnonsmooth\t{}\t\t{text}",word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
    resolutions.flush().map_err(|e|e.to_string())
}

// Passport-preserving neighbourhoods of independently known smooth components.
// Break every old block system; do not infer primitivity or closure from that.
struct G3MutationSeed {tuple:[Perm;4],partitions:Vec<Vec<usize>>}
struct G3Mutation {n:usize,seeds:Vec<G3MutationSeed>,current:usize,frame:Option<ForestFrame>,third:Perm,third_index:usize,fourth_index:usize,used:usize}
fn g3_preserves(p:Perm,labels:&[usize],n:usize)->bool {
    let mut image=[usize::MAX;MAX_D];
    for x in 0..n {let a=labels[x];let b=labels[p.img[x] as usize];
        if image[a]==usize::MAX{image[a]=b;}else if image[a]!=b{return false;}
    }true
}
impl G3Mutation {
    fn load(path:&Path,n:usize)->Result<Self,String> {
        let text=fs::read_to_string(path).map_err(|e|e.to_string())?;let mut seeds=Vec::new();
        for line in text.lines() {
            let parts:Vec<_>=line.splitn(3,'\t').collect();if parts.len()!=3{return Err("mutation input shape".into());}
            let degree:usize=parts[0].parse().map_err(|_|"mutation degree")?;if degree!=n{continue;}
            let tuple=sm_parse(parts[1],n)?;
            if tuple.iter().map(|p|sm_index(*p,n)).sum::<usize>()!=2*n+4{return Err("mutation genus".into());}
            let partitions=parts[2].split(';').filter(|s|!s.is_empty()).map(|s|s.split(',').map(|v|v.parse::<usize>().map_err(|_|"mutation block label")).collect::<Result<Vec<_>,_>>()).collect::<Result<Vec<_>,_>>()?;
            for labels in &partitions {
                if labels.len()!=n||labels.iter().any(|&v|v>=n)||tuple.iter().any(|&p|!g3_preserves(p,labels,n)) {return Err("mutation partition".into());}
                let r=labels.iter().copied().max().unwrap()+1;let mut sizes=vec![0;r];for &v in labels{sizes[v]+=1;}
                if r<=1||r>=n||sizes.iter().any(|&v|v!=sizes[0]){return Err("proper equal blocks required".into());}
            }
            seeds.push(G3MutationSeed{tuple,partitions});
        }
        if seeds.is_empty(){return Err("no mutation seeds at degree".into());}
        Ok(Self{n,seeds,current:0,frame:None,third:Perm::identity(n),third_index:0,fourth_index:0,used:128})
    }
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]> {
        let n=self.n;
        if self.used>=128||self.frame.is_none() {
            self.current=rng.below(self.seeds.len());let mut t=self.seeds[self.current].tuple;
            for _ in 0..rng.below(16) {t=hurwitz_move_dynamic(n,t,rng.below(6));}
            if sm_index(t[0],n)+sm_index(t[1],n)!=ff_product_index(t[0],t[1],n) {t=[t[2],t[3],t[0],t[1]];}
            assert_eq!(sm_index(t[0],n)+sm_index(t[1],n),ff_product_index(t[0],t[1],n));
            self.frame=Some(ForestFrame::new([t[0],t[1]],n));self.third=t[2];self.third_index=sm_index(t[2],n);self.fourth_index=sm_index(t[3],n);self.used=0;
        }
        self.used+=1;let mut c=self.third;
        for _ in 0..1+rng.below(4) {
            let x=rng.below(n);let mut y=rng.below(n-1);if y>=x{y+=1;}
            // Conjugate c by a transposition without a heap permutation.
            c.img.swap(x,y);for v in &mut c.img[..n]{if *v as usize==x{*v=y as u8;}else if *v as usize==y{*v=x as u8;}}
        }
        if c==self.third||self.seeds[self.current].partitions.iter().any(|labels|g3_preserves(c,labels,n)) {return None;}
        let f=self.frame.as_ref().unwrap();if ff_product_index(f.p,c,n)!=self.fourth_index{return None;}
        Some(f.tuple(c,n))
    }
}

// Condition on the positive-genus right pair of a certified smooth control,
// rather than independently drawing a third entry. Re-factor p=(c*e)^-1
// with a low-index involution a. Fresh disjoint transpositions split a*p
// and simultaneously join a*c and a*e. Every accepted step therefore keeps
// (a,b),(a,c),(a,e) additive. This is only a proposal law: mixed braid frames
// can still escape, and a capped scan is never a positive certificate.
struct G3Refactor {
    n:usize,seeds:Vec<[Perm;4]>,current:Option<[Perm;4]>,used:usize,
    first_frame:Option<ForestFrame>,last_indices:[usize;4],cross:bool,joined:bool,min_index:usize,max_index:usize,
    adaptive:Vec<[Perm;4]>,adaptive_keys:HashSet<SmoothKey>,adaptive_limit:usize,adaptive_draws:u64,
}
impl G3Refactor {
    fn load(path:&Path,n:usize)->Result<Self,String>{
        let known=G3Mutation::load(path,n)?;
        Ok(Self{n,seeds:known.seeds.into_iter().map(|s|s.tuple).collect(),current:None,used:128,
            first_frame:None,last_indices:[0;4],cross:true,joined:false,min_index:1,max_index:4,
            adaptive:Vec::new(),adaptive_keys:HashSet::new(),adaptive_limit:0,adaptive_draws:0})
    }
    // These are untrusted CONSTRUCTOR templates, never known-positive states
    // and never used to skip a test. Every admission has a saved negative root
    // and an exact smooth prefix; its future proposals still need full closure.
    fn remember(&mut self,t:[Perm;4],word:&[usize])->Option<(usize,[Perm;4])>{
        if self.adaptive.len()>=self.adaptive_limit||word.len()<2{return None;}
        let n=self.n;let mut last=t;let mut prefix=0;
        for (j,&m) in word.iter().enumerate(){
            let u=hurwitz_move_dynamic(n,last,m);
            if !ff_smooth(&u,u.map(|p|sm_index(p,n)),n){break;}
            last=u;prefix=j+1;
        }
        if prefix==0||prefix>=word.len(){return None;}
        let key=ff_key(&last,n);
        if !self.adaptive_keys.insert(key){return None;}
        self.adaptive.push(last);Some((prefix,last))
    }
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]>{
        let n=self.n;
        if self.used>=128||self.current.is_none(){
            let mut t=if !self.adaptive.is_empty()&&rng.below(2)==1 {
                self.adaptive_draws+=1;self.adaptive[rng.below(self.adaptive.len())]
            }else{self.seeds[rng.below(self.seeds.len())]};
            for _ in 0..rng.below(24){
                let next=hurwitz_move_dynamic(n,t,rng.below(6));
                // A recycled prefix is smooth HERE, not everywhere. Never
                // assume that a random next braid remains smooth.
                if !ff_smooth(&next,next.map(|p|sm_index(p,n)),n){break;}t=next;
            }
            if sm_index(t[0],n)+sm_index(t[1],n)!=ff_product_index(t[0],t[1],n){t=[t[2],t[3],t[0],t[1]];}
            assert_eq!(sm_index(t[0],n)+sm_index(t[1],n),ff_product_index(t[0],t[1],n));
            self.current=Some(t);self.used=0;
        }
        self.used+=1;let old=self.current.unwrap();let p=old[0].compose(old[1],n);
        let ip=sm_index(p,n);let max_rank=self.max_index.min(n/2).min(ip-1);
        if max_rank<self.min_index{return None;}
        let rank=self.min_index+rng.below(max_rank-self.min_index+1);
        let ix=[rank,ip-rank,sm_index(old[2],n),sm_index(old[3],n)];
        let mut sorted=ix;sorted.sort_unstable();
        if ix[1]==0||!(sorted[0]+sorted[3]<n||sorted[1]+sorted[2]<n){return None;}
        let mut a=Perm::identity(n);let mut b=p;let mut ac=old[2];let mut ae=old[3];let mut moved=0u64;
        for _ in 0..rank{
            let lb=dd_cycle_labels(b,n);let lc=dd_cycle_labels(ac,n);let le=dd_cycle_labels(ae,n);
            let la=dd_cycle_labels(a,n);let ai=a.inverse(n);
            let mut choices=[(0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
            for x in 0..n{if !self.joined&&moved&(1<<x)!=0{continue;}for y in x+1..n{
                let boundary=if self.joined{la[x]!=la[y]&&lb[ai.img[x]as usize]==lb[ai.img[y]as usize]}
                             else{moved&(1<<y)==0&&lb[x]==lb[y]};
                if boundary&&(!self.cross||lc[x]!=lc[y]&&le[x]!=le[y]){
                    choices[count]=(x as u8,y as u8);count+=1;
                }
            }}
            if count==0{return None;}let(x,y)=choices[rng.below(count)];let(x,y)=(x as usize,y as usize);
            let(bx,by)=if self.joined{(ai.img[x]as usize,ai.img[y]as usize)}else{(x,y)};
            if self.joined{for v in &mut a.img[..n]{if *v as usize==x{*v=y as u8;}else if *v as usize==y{*v=x as u8;}}}
            else{moved|=(1<<x)|(1<<y);a.img.swap(x,y);}
            for v in &mut b.img[..n]{if *v as usize==bx{*v=by as u8;}else if *v as usize==by{*v=bx as u8;}}
            for q in [&mut ac,&mut ae]{for v in &mut q.img[..n]{if *v as usize==x{*v=y as u8;}else if *v as usize==y{*v=x as u8;}}}
        }
        let t=[a,b,old[2],old[3]];
        assert_eq!(a.compose(b,n),p);debug_assert_eq!(t.map(|q|sm_index(q,n)),ix);
        if self.cross{assert_eq!(sm_index(ac,n),rank+ix[2]);assert_eq!(sm_index(ae,n),rank+ix[3]);}
        self.last_indices=ix;self.first_frame=Some(ForestFrame::new([a,b],n));Some(t)
    }
}

// A centralizer element of a transitive action is uniquely fixed by the image
// of point zero. Propagate equivariance; no factorial group enumeration.
fn g3_has_deck(t:&[Perm;4],n:usize)->bool {
    for target in 1..n {
        let mut map=[u8::MAX;MAX_D];let mut used=0u64;let mut queue=[0usize;MAX_D];let mut size=1;let mut head=0;
        map[0]=target as u8;used|=1u64<<target;let mut valid=true;
        while head<size && valid {
            let x=queue[head];head+=1;let y=map[x] as usize;
            for p in &t[..3] {
                let a=p.img[x] as usize;let b=p.img[y] as usize;
                if map[a]!=u8::MAX {if map[a] as usize!=b {valid=false;break;}}
                else {if used&(1u64<<b)!=0 {valid=false;break;}map[a]=b as u8;used|=1u64<<b;queue[size]=a;size+=1;}
            }
        }
        if valid && size==n {return true;}
    }
    false
}

// Positive sufficient isotriviality witness, not a failure of smoothness.
// A genus-three single-local elliptic cover has moving index exactly four.
// Force its entire support into one block and close the equivalence under G.
fn g3_local_isotrivial(t:&[Perm;4],ix:[usize;4],n:usize)->Option<(usize,[usize;MAX_D])> {
    for special in 0..4 {
        if ix[special]!=4{continue;}
        let moved:Vec<_>=(0..n).filter(|&x|t[special].img[x] as usize!=x).collect();
        let mut parent=std::array::from_fn(|i|i);
        for &x in &moved[1..]{gd_join(&mut parent,moved[0],x);}
        loop {
            let mut first=[usize::MAX;MAX_D];let mut changed=false;
            for x in 0..n {
                let root=gd_find(&mut parent,x);
                if first[root]==usize::MAX{first[root]=x;continue;}
                for p in &t[..3]{changed|=gd_join(&mut parent,p.img[first[root]] as usize,p.img[x] as usize);}
            }
            if !changed{break;}
        }
        let mut names=[usize::MAX;MAX_D];let mut labels=[0usize;MAX_D];let mut reps=[0usize;MAX_D];let mut r=0;
        for x in 0..n {let root=gd_find(&mut parent,x);if names[root]==usize::MAX{names[root]=r;reps[r]=x;r+=1;}labels[x]=names[root];}
        if r<=1||r>=n||n%r!=0{continue;}
        let m=n/r;let mut sizes=[0usize;MAX_D];for &b in &labels[..n]{sizes[b]+=1;}
        if sizes[..r].iter().any(|&s|s!=m){continue;}
        let mut qi=[0usize;4];let mut valid=true;
        for j in 0..4 {
            let mut q=Perm::identity(r);for b in 0..r{q.img[b]=labels[t[j].img[reps[b]] as usize] as u8;}
            if (0..n).any(|x|labels[t[j].img[x] as usize]!=q.img[labels[x]] as usize){valid=false;break;}
            qi[j]=sm_index(q,r);
            if (j==special&&qi[j]!=0)||(j!=special&&ix[j]!=m*qi[j]){valid=false;break;}
        }
        if valid&&qi.iter().sum::<usize>()==2*r{return Some((special,labels));}
    }None
}
fn g3_iso_row(t:&[Perm;4],n:usize,candidate:u64,special:usize,labels:&[usize;MAX_D])->String {
    format!("{n}\t{candidate}\t{}\t{}\t{}\n",special+1,labels[..n].iter().map(|v|v.to_string()).collect::<Vec<_>>().join(","),tuple_string_dynamic(n,t))
}

// Negative-only outer orbit test. A bounded prefix cannot certify closure.
fn g3_outer_failure(t:[Perm;4],n:usize,rotations:usize)->Option<Vec<usize>> {
    let h=t[0].compose(t[1],n);let hi=h.inverse(n);let mut u=t;
    for k in 0..rotations {
        for (left,right) in [(false,false),(true,false),(false,true),(true,true)] {
            let mut v=u;let mut word=vec![0;2*k];
            if left {v=hurwitz_move_dynamic(n,v,0);word.push(0);}
            if right {v=hurwitz_move_dynamic(n,v,2);word.push(2);}
            let mut ix=v.map(|p|sm_index(p,n));ix.swap(1,2);
            for m in [1,4] {
                let w=hurwitz_move_dynamic(n,v,m);
                if !ff_smooth(&w,ix,n) {word.push(m);return Some(word);}
            }
        }
        u[0]=h.compose(u[0],n).compose(hi,n);u[1]=h.compose(u[1],n).compose(hi,n);
        if u[0]==t[0] && u[1]==t[1] {break;}
    }
    None
}
fn g3_path(parent:&[(usize,Vec<usize>)],mut i:usize)->Vec<usize> {
    let mut pieces=Vec::new();while i>0 {let (p,w)=&parent[i];pieces.push(w);i=*p;}
    pieces.into_iter().rev().flat_map(|w|w.iter().copied()).collect()
}
// All outer states remain in the certificate, but each whole H=<sigma1,sigma3>
// orbit enters the expensive middle scan just once. Squares differ by an inner
// conjugation, so 4 parity cases times the sigma1^2 period cover H exactly.
fn g3_macro_scan(seed:&[Perm;4],n:usize,cap:usize,deadline:Instant)->SmoothScan {
    assert!(cap>0);let root=ff_key(seed,n);let mut states=vec![root];let mut seen=HashMap::from([(root,0usize)]);
    let mut parents=vec![(0,Vec::new())];let mut expanded=HashSet::new();let mut head=0;
    if !sm_values(seed,n).0 {return SmoothScan{status:"nonsmooth",states,word:vec![]};}
    while head<states.len() {
        if expanded.contains(&states[head]) {head+=1;continue;}
        let start=states[head];let mut current=sm_unpack(start,n);let mut previous=head;
        let mut k=0usize;let mut members=Vec::new();let mut local=HashSet::new();
        loop {
            if Instant::now()>=deadline {return SmoothScan{status:"timeout",states,word:vec![]};}
            let base_key=ff_key(&current,n);
            if k>0 && base_key==start {break;}
            let base=if let Some(&i)=seen.get(&base_key){i}else{
                if states.len()==cap {return SmoothScan{status:"capped",states,word:vec![]};}
                let i=states.len();seen.insert(base_key,i);states.push(base_key);parents.push((previous,vec![0,0]));i
            };
            for (left,right) in [(false,false),(true,false),(false,true),(true,true)] {
                let mut v=current;let mut word=Vec::new();
                if left {v=hurwitz_move_dynamic(n,v,0);word.push(0);}
                if right {v=hurwitz_move_dynamic(n,v,2);word.push(2);}
                let key=if !left&&!right {base_key}else{ff_key(&v,n)};
                let id=if let Some(&i)=seen.get(&key){i}else{
                    if states.len()==cap {return SmoothScan{status:"capped",states,word:vec![]};}
                    let i=states.len();seen.insert(key,i);states.push(key);parents.push((base,word));i
                };
                if local.insert(key) {members.push(id);}
            }
            previous=base;current=hurwitz_move_dynamic(n,hurwitz_move_dynamic(n,current,0),0);k+=1;
        }
        for i in members {
            expanded.insert(states[i]);let t=sm_unpack(states[i],n);let mut ix=t.map(|p|sm_index(p,n));ix.swap(1,2);
            for m in [1,4] {
                if Instant::now()>=deadline {return SmoothScan{status:"timeout",states,word:vec![]};}
                let u=hurwitz_move_dynamic(n,t,m);
                if !ff_smooth(&u,ix,n) {let mut word=g3_path(&parents,i);word.push(m);return SmoothScan{status:"nonsmooth",states,word};}
                let key=ff_key(&u,n);
                if !seen.contains_key(&key) {
                    if states.len()==cap {return SmoothScan{status:"capped",states,word:vec![]};}
                    seen.insert(key,states.len());states.push(key);parents.push((i,vec![m]));
                }
            }
        }
        head+=1;
    }
    SmoothScan{status:"closed_smooth",states,word:vec![]}
}

struct G3TrieNode {children:[Option<usize>;6]}
struct G3Words {nodes:Vec<G3TrieNode>}
impl G3Words {
    fn new(words:&[Vec<usize>])->Self {
        let mut nodes=vec![G3TrieNode{children:[None;6]}];
        for word in words.iter().take(64) {let mut at=0;for &m in word.iter().take(24) {
            assert!(m<6);let next=match nodes[at].children[m] {Some(i)=>i,None=>{
                let i=nodes.len();nodes.push(G3TrieNode{children:[None;6]});nodes[at].children[m]=Some(i);i
            }};at=next;
        }}Self{nodes}
    }
    fn failure_reference(&self,t:[Perm;4],n:usize)->Option<Vec<usize>> {
        fn visit(tree:&G3Words,at:usize,t:[Perm;4],n:usize,path:&mut Vec<usize>)->Option<Vec<usize>> {
            for m in 0..6 {if let Some(next)=tree.nodes[at].children[m] {
                let u=hurwitz_move_dynamic(n,t,m);path.push(m);
                if !ff_smooth(&u,u.map(|p|sm_index(p,n)),n) {return Some(path.clone());}
                if let Some(w)=visit(tree,next,u,n,path){return Some(w);}path.pop();
            }}None
        }
        visit(self,0,t,n,&mut Vec::new())
    }
    fn failure_indices(&self,t:[Perm;4],indices:[usize;4],n:usize)->Option<Vec<usize>> {
        fn visit(tree:&G3Words,at:usize,t:[Perm;4],ix:[usize;4],n:usize,path:&mut Vec<usize>)->Option<Vec<usize>> {
            for m in 0..6 {if let Some(next)=tree.nodes[at].children[m] {
                let u=hurwitz_move_dynamic(n,t,m);let mut indices=ix;let k=m%3;indices.swap(k,k+1);path.push(m);
                if !ff_product_matches(u[0],u[1],n,indices[0]+indices[1],indices[2]+indices[3]){return Some(path.clone());}
                if let Some(w)=visit(tree,next,u,indices,n,path){return Some(w);}path.pop();
            }}None
        }
        visit(self,0,t,indices,n,&mut Vec::new())
    }
    fn failure(&self,t:[Perm;4],n:usize)->Option<Vec<usize>> {self.failure_indices(t,t.map(|p|sm_index(p,n)),n)}
}

pub fn genus_three_main(args:&[String])->Result<(),String> {
    let arg=|name:&str,default:&str|args.windows(2).find(|w|w[0]==name).map(|w|w[1].clone()).unwrap_or(default.into());
    if args.iter().any(|a|a=="--classify-frame-recipes") {
        let input=fs::read_to_string(arg("--classify-frame-recipes","")).map_err(|e|e.to_string())?;
        let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh --out required".into());}
        fs::create_dir(&out).map_err(|e|e.to_string())?;
        let mut rows=BufWriter::new(File::create(out.join("FRAME_CLASSIFICATION.tsv")).map_err(|e|e.to_string())?);
        writeln!(rows,"candidate\tdegree\tk\treason\tprofile\tpartial_pair_identity_padded_NOT_A_CANDIDATE").map_err(|e|e.to_string())?;
        for(i,line)in input.lines().enumerate(){
            let(raw_n,text)=line.split_once('\t').ok_or("fixture degree and four cycle types")?;let n=raw_n.parse::<usize>().map_err(|_|"fixture degree")?;
            if n==0||n>SM_MAX_DEGREE{return Err("frame fixture degree outside representation".into());}
            let shapes:Vec<_>=text.split('|').collect();if shapes.len()!=4{return Err("four fixture cycle types required".into());}
            let mut t=[Perm::identity(n);4];
            for(j,raw)in shapes.iter().enumerate(){let shape=raw.split_whitespace().map(|x|x.parse::<usize>().map_err(|_|"fixture cycle length")).collect::<Result<Vec<_>,_>>()?;
                if shape.contains(&0)||shape.iter().sum::<usize>()!=n{return Err("fixture cycle type degree".into());}
                let mut at=0;for length in shape{for x in 0..length{t[j].img[at+x]=(at+(x+1)%length)as u8;}at+=length;}
            }
            let ix=t.map(|p|sm_index(p,n));let cut=G3FrameRetirement::decision([t[0],t[1]],ix,n);
            let(k,reason)=if let Some((reason,k))=cut{(k.to_string(),G3FrameRetirement::reason_name(reason))}else{(String::new(),"not-retired")};
            writeln!(rows,"{}\t{n}\t{k}\t{reason}\t{}\t{}",i+1,ix.iter().map(|x|x.to_string()).collect::<Vec<_>>().join(","),tuple_string_dynamic(n,&[t[0],t[1],Perm::identity(n),Perm::identity(n)])).map_err(|e|e.to_string())?;
        }
        fs::write(out.join("SCOPE.txt"),"Counterfactual incomplete-frame diagnostics only. The two padding identities are NOT candidate entries, and fixture representatives need not satisfy the first-pair forest. Managed production uses actual additive construction frames.\n").map_err(|e|e.to_string())?;return Ok(());
    }
    if args.iter().any(|a|a=="--dump-saturated-rules") {
        println!("degree\tk\tprime\tmodulus\tmethod");
        for n in 2..=SM_MAX_DEGREE {let policy=G3SaturatedPassport::new(n);
            for k in 2..n {for rule in &policy.rules[k]{println!("{n}\t{k}\t{}\t{}\t{}",rule.prime,rule.modulus,rule.method());}}
        }return Ok(());
    }
    if args.iter().any(|a|a=="--classify-saturated") {
        let input=fs::read_to_string(arg("--classify-saturated","")).map_err(|e|e.to_string())?;
        let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh --out required".into());}
        fs::create_dir(&out).map_err(|e|e.to_string())?;
        let mut rows=BufWriter::new(File::create(out.join("CLASSIFICATION.tsv")).map_err(|e|e.to_string())?);
        write!(rows,"{G3_SATURATED_HEADER}").map_err(|e|e.to_string())?;let mut count=0;let mut rejected=0;
        for(i,line)in input.lines().enumerate(){let(n,text)=line.split_once('\t').ok_or("fixture degree and four cycle types")?;let n=n.parse::<usize>().map_err(|_|"fixture degree")?;
            if n==0||n>SM_MAX_DEGREE{return Err("fixture degree outside checked representation".into());}
            let shapes:Vec<_>=text.split('|').collect();if shapes.len()!=4{return Err("four fixture cycle types required".into());}
            let mut t=[Perm::identity(n);4];
            for(j,raw)in shapes.iter().enumerate(){let shape=raw.split_whitespace().map(|x|x.parse::<usize>().map_err(|_|"fixture cycle length")).collect::<Result<Vec<_>,_>>()?;
                if shape.contains(&0)||shape.iter().sum::<usize>()!=n{return Err("fixture cycle type degree".into());}
                let mut at=0;for length in shape{for x in 0..length{t[j].img[at+x]=(at+(x+1)%length)as u8;}at+=length;}
            }
            let cut=G3SaturatedPassport::new(n).rejection(&t,t.map(|p|sm_index(p,n)));
            rejected+=usize::from(cut.is_some());write!(rows,"{}",g3_sp_row((i+1)as u64,cut,&tuple_string_dynamic(n,&t))).map_err(|e|e.to_string())?;count+=1;
        }
        fs::write(out.join("SCOPE.txt"),"Arithmetic passport classifier only; fixture tuples need not have product one. Production invocation is restricted to identity-free source-genus-three product-one proposals. Theorem exclusions are not literal negative words or isotriviality labels.\n").map_err(|e|e.to_string())?;
        println!("classified={count} theorem_excluded={rejected}");return Ok(());
    }
    if args.iter().any(|a|a=="--check-controls") {return forest_fast_main(args);}
    if args.iter().any(|a|a=="--isotrivial-controls") {
        let input=fs::read_to_string(arg("--isotrivial-controls","")).map_err(|e|e.to_string())?;
        let path=PathBuf::from(arg("--out",""));fs::create_dir(&path).map_err(|e|e.to_string())?;
        let mut f=BufWriter::new(File::create(path.join("ISOTRIVIAL.tsv")).map_err(|e|e.to_string())?);
        writeln!(f,"degree\tcandidate\tspecial\tlabels\tseed").map_err(|e|e.to_string())?;let mut count=0;
        for (i,line) in input.lines().enumerate() {
            let p:Vec<_>=line.splitn(3,'\t').collect();let n=p[0].parse::<usize>().map_err(|_|"control degree")?;
            let t=sm_parse(p[1],n)?;let Some((s,l))=g3_local_isotrivial(&t,t.map(|p|sm_index(p,n)),n)else{return Err(format!("missing isotrivial control {i}"));};
            write!(f,"{}",g3_iso_row(&t,n,(i+1)as u64,s,&l)).map_err(|e|e.to_string())?;count+=1;
        }println!("isotrivial_controls={count}");return Ok(());
    }
    if args.iter().any(|a|a=="--controls") {
        let text=fs::read_to_string(arg("--controls","")).map_err(|e|e.to_string())?;
        let verifier=arg("--verifier","macro");
        if !["macro","ordinary"].contains(&verifier.as_str()){return Err("control verifier".into());}
        let out=if args.iter().any(|a|a=="--out") {
            let p=PathBuf::from(arg("--out",""));fs::create_dir(&p).map_err(|e|e.to_string())?;Some(p)
        }else{None};
        let mut total=0;let mut elapsed=0.0;let mut rows=String::from("control\tdegree\tstates\tclosure_seconds\n");
        for (i,line) in text.lines().enumerate() {
            let p:Vec<_>=line.splitn(3,'\t').collect();if p.len()!=3{return Err("control degree,states,tuple".into());}
            let n:usize=p[0].parse().map_err(|_|"control degree")?;let expected:usize=p[1].parse().map_err(|_|"control states")?;
            let saturated=G3SaturatedPassport::new(n);
            let t=sm_parse(p[2],n)?;let began=Instant::now();let limit=began+Duration::from_secs(120);
            let scan=if verifier=="macro"{g3_macro_scan(&t,n,1_000_000,limit)}else{ff_scan(&t,n,1_000_000,limit)};
            let seconds=began.elapsed().as_secs_f64();elapsed+=seconds;
            if scan.status!="closed_smooth"||scan.states.len()!=expected{return Err(format!("control {i} mismatch"));}
            if args.iter().any(|a|a=="--check-cyclic") {
                let mut cyclic=G3CyclicCache::default();
                let mut forest_cuts=G3ForestCutCache::default();
                for &key in &scan.states {
                    let state=sm_unpack(key,n);
                    assert!(cyclic.failure(&state,state.map(|p|sm_index(p,n)),n).is_none());
                    assert!(forest_cuts.failure(&state,state.map(|p|sm_index(p,n)),n,16).is_none());
                }
            }
            if args.iter().any(|a|a=="--check-saturated-passports") {
                for &key in &scan.states {let state=sm_unpack(key,n);assert!(saturated.rejection(&state,state.map(|p|sm_index(p,n))).is_none());}
            }
            if args.iter().any(|a|a=="--compare") {
                let ordinary=ff_scan(&t,n,1_000_000,Instant::now()+Duration::from_secs(120));
                if ordinary.status!=scan.status||ordinary.states.iter().copied().collect::<HashSet<_>>() != scan.states.iter().copied().collect() {
                    return Err(format!("control {i} ordinary state-set mismatch"));
                }
            }
            if let Some(p)=&out {sm_save(&t,n,&scan,&p.join(format!("control_{i:03}")))?;}
            rows.push_str(&format!("{i}\t{n}\t{expected}\t{seconds}\n"));total+=expected;
        }
        if let Some(p)=&out{fs::write(p.join("CONTROL_TIMINGS.tsv"),rows).map_err(|e|e.to_string())?;}
        println!("controls={} states={total} closure_seconds={elapsed} verifier={verifier}",text.lines().count());return Ok(());
    }
    let n:usize=arg("--d","14").parse().map_err(|_|"degree")?;
    let seconds:u64=arg("--seconds","45").parse().map_err(|_|"seconds")?;
    let attempts:u64=arg("--attempts","0").parse().map_err(|_|"attempts")?;
    let seed:u64=arg("--seed","1").parse().map_err(|_|"seed")?;
    let cap:usize=arg("--cap","250000").parse().map_err(|_|"cap")?;
    let closure:u64=arg("--closure-seconds","20").parse().map_err(|_|"closure")?;
    let mode=arg("--lane","collision");let verifier=arg("--verifier","macro");
    if !(9..=SM_MAX_DEGREE).contains(&n)||seconds==0||cap==0||closure==0||!["collision","legacy-types","elliptic","rational","mutation","refactor","full-fiber","partial-fiber","general-fiber","four-quotient","joint-genus3","involution-boundary"].contains(&mode.as_str())
        ||!["macro","ordinary"].contains(&verifier.as_str()) {return Err("invalid bounded genus3 configuration".into());}
    if args.iter().any(|a|a=="--tuple") {
        let t=sm_parse(&arg("--tuple",""),n)?;
        if t.iter().any(|p|p.is_identity(n))||t.iter().map(|p|sm_index(*p,n)).sum::<usize>()!=2*n+4 {
            return Err("replay must be identity-free genus three".into());
        }
        let scan=if verifier=="macro"{g3_macro_scan(&t,n,cap,Instant::now()+Duration::from_secs(closure))}
                 else{ff_scan(&t,n,cap,Instant::now()+Duration::from_secs(closure))};
        let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh --out required".into());}
        sm_save(&t,n,&scan,&out)?;println!("REPLAY status={} states={}",scan.status,scan.states.len());return Ok(());
    }
    let words=if args.iter().any(|a|a=="--words") {
        fs::read_to_string(arg("--words","")).map_err(|e|e.to_string())?.lines().filter(|l|!l.trim().is_empty()).map(|line|
            line.split(',').map(|s|s.parse::<usize>().map_err(|_|"word move")).collect::<Result<Vec<_>,_>>()
        ).collect::<Result<Vec<_>,_>>()?
    }else{Vec::new()};
    if words.len()>64||words.iter().any(|w|w.is_empty()||w.len()>24||w.iter().any(|&m|m>=6)){return Err("word bank bounds".into());}
    let trie=G3Words::new(&words);let mut collision=G3Collision::new(n,mode!="legacy-types");
    let mut boundary=if mode=="involution-boundary" {
        Some(G3InvolutionBoundary::new(n,arg("--cursor","0").parse().map_err(|_|"boundary cursor")?)?)
    }else{None};
    let mut joint=if mode=="joint-genus3"{Some(G3Joint::new(n))}else{None};
    let retire_frames=args.iter().any(|s|s=="--retire-proven-frames");
    collision.retirement.enabled=retire_frames;
    if let Some(j)=joint.as_mut(){j.base.retirement.enabled=retire_frames;}
    let exclude_transpositions=args.iter().any(|s|s=="--exclude-transpositions");
    let use_saturated=args.iter().any(|s|s=="--saturated-passport-cuts");
    let saturated_stage=arg("--saturated-stage","after-short");
    if !["early","after-short","late"].contains(&saturated_stage.as_str()){return Err("saturated stage early/after-short/late".into());}
    let saturated_trace_limit:usize=arg("--trace-saturated","0").parse().map_err(|_|"saturated trace bound")?;
    if saturated_trace_limit>4096{return Err("saturated trace bound at most4096".into());}
    let saturated=if use_saturated||saturated_trace_limit>0{Some(G3SaturatedPassport::new(n))}else{None};
    let index_profile=if args.iter().any(|s|s=="--index-profile") {
        if !["collision","legacy-types","joint-genus3"].contains(&mode.as_str()){return Err("index profile requires a block-free constructor".into());}
        let profile=g3_index_profile(n,&arg("--index-profile",""))?;
        collision.passports=vec![profile];if let Some(j)=joint.as_mut(){j.base.passports=vec![profile];}Some(profile)
    }else{None};
    let before_profiles=collision.passports.len();
    if exclude_transpositions{
        collision.passports.retain(|p|p.iter().all(|&i|i>1));
        if let Some(j)=joint.as_mut(){j.base.passports.retain(|p|p.iter().all(|&i|i>1));}
        if index_profile.is_some_and(|p|p.contains(&1)){return Err("selected genus3 profile excluded by the admitted transposition theorem".into());}
    }
    let avoided_profiles=before_profiles-collision.passports.len();
    if args.iter().any(|s|s=="--joint-csp") {
        let j=joint.as_mut().ok_or("--joint-csp requires joint-genus3")?;j.csp=true;
        j.budget=arg("--joint-nodes","5000").parse().map_err(|_|"joint node budget")?;
        if j.budget==0||j.budget>1000000{return Err("joint node budget outside1..1000000".into());}
    }
    let mut general=if mode=="general-fiber"{Some(G3GeneralFiber::new(n))}else{None};
    let mut four=if mode=="four-quotient"{Some(G3FourQuotient::new(n))}else{None};
    if args.iter().any(|s|s=="--return-rejection"){
        if let Some(q)=four.as_mut(){q.solved=false;}else if let Some(g)=general.as_mut(){g.solved=false;}
        else{return Err("return rejection requires general-fiber/four-quotient".into());}
    }
    let mut elliptic=if mode=="elliptic"||mode=="rational"{Some(G3Elliptic::with_genus(n,if mode=="elliptic"{1}else{0}))}else{None};
    if exclude_transpositions{if let Some(e)=elliptic.as_mut(){e.minimum_moving=2;}}
    let mut mutation=if mode=="mutation" {Some(G3Mutation::load(Path::new(&arg("--mutation-seeds","")),n)?)}else{None};
    let mut refactor=if mode=="refactor"{Some(G3Refactor::load(Path::new(&arg("--mutation-seeds","")),n)?)}else{None};
    if let Some(r)=refactor.as_mut(){
        r.cross=!args.iter().any(|s|s=="--refactor-free");
        r.joined=args.iter().any(|s|s=="--refactor-joined");
        r.max_index=arg("--refactor-max-index","4").parse().map_err(|_|"refactor max index")?;
        if !(1..=4).contains(&r.max_index){return Err("refactor max index1..4".into());}
        r.min_index=if exclude_transpositions{2}else{1};
        if r.max_index<r.min_index{return Err("refactor interval empty after transposition theorem".into());}
        r.adaptive_limit=arg("--adaptive-refactor","0").parse().map_err(|_|"adaptive template limit")?;
        if r.adaptive_limit>1024{return Err("adaptive template bound at most1024".into());}
    }
    let mut exact=if ["full-fiber","partial-fiber"].contains(&mode.as_str())&&args.iter().any(|s|s=="--fiber-kind") {
        Some(G3FiberExact::new(n,&arg("--fiber-kind",""),arg("--cursor","0").parse().map_err(|_|"fiber cursor")?)?)
    }else{None};
    let mut fiber=if exact.is_none()&&["full-fiber","partial-fiber"].contains(&mode.as_str()){
        Some(if mode=="partial-fiber"{G3FullFiber::partial(n)}else{G3FullFiber::new(n)})
    }else{None};
    let finite=exact.is_some()||boundary.is_some();
    let boundary_start=boundary.as_ref().map(|b|b.cursor);
    if args.iter().any(|s|s=="--skip-covered-fiber"){
        fiber.as_mut().ok_or("covered-family skip requires sampled fibre lane")?.skip_covered=true;
    }
    // Diagnostic baseline only; the admitted runner never supplies this flag.
    if args.iter().any(|s|s=="--fiber-rejection"){fiber.as_mut().ok_or("fiber rejection needs sampled full-fiber lane")?.solved=false;}
    let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh --out required".into());}
    fs::create_dir(&out).map_err(|e|e.to_string())?;
    if let Some(profile)=index_profile{fs::write(out.join("INDEX_PROFILE.txt"),format!("{profile:?}\n")).map_err(|e|e.to_string())?;}
    let mut leads=BufWriter::new(File::create(out.join("LEADS.tsv")).map_err(|e|e.to_string())?);
    let mut resolutions=BufWriter::new(File::create(out.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?);
    writeln!(resolutions,"candidate\tstatus\tword\tfolder\tseed").map_err(|e|e.to_string())?;
    fs::write(out.join("KNOWN.tsv"),"").map_err(|e|e.to_string())?;
    let mut iso=BufWriter::new(File::create(out.join("ISOTRIVIAL.tsv")).map_err(|e|e.to_string())?);
    writeln!(iso,"degree\tcandidate\tspecial\tlabels\tseed").map_err(|e|e.to_string())?;
    let mut theorem_log=if exclude_transpositions{
        let mut file=BufWriter::new(File::create(out.join("TARGET_EXCLUSIONS.tsv")).map_err(|e|e.to_string())?);
        writeln!(file,"candidate\tseed").map_err(|e|e.to_string())?;Some(file)
    }else{None};
    let mut theorem_skipped=0u64;let mut theorem_sampled=0u64;
    let mut saturated_log=if use_saturated{let mut f=BufWriter::new(File::create(out.join("SATURATED_EXCLUSIONS.tsv")).map_err(|e|e.to_string())?);write!(f,"{G3_SATURATED_HEADER}").map_err(|e|e.to_string())?;Some(f)}else{None};
    let mut saturated_trace=if saturated_trace_limit>0{let mut f=BufWriter::new(File::create(out.join("SATURATED_TRACE.tsv")).map_err(|e|e.to_string())?);write!(f,"{G3_SATURATED_HEADER}").map_err(|e|e.to_string())?;Some(f)}else{None};
    let mut saturated_checked=0u64;let mut saturated_sampled=0u64;let mut saturated_traced=0usize;
    let mut saturated_counts=std::collections::BTreeMap::<(u8,u8,u8,u8,u8),(u64,u64)>::new();
    let mut saturated_screen=|t:&[Perm;4],ix:[usize;4],candidate:u64|->Result<bool,String>{
        let Some(policy)=&saturated else{return Ok(false);};
        saturated_checked+=1;let cut=policy.rejection(t,ix);
        if saturated_traced<saturated_trace_limit{write!(saturated_trace.as_mut().unwrap(),"{}",g3_sp_row(candidate,cut,&tuple_string_dynamic(n,t))).map_err(|e|e.to_string())?;saturated_traced+=1;}
        if use_saturated{if let Some(cut)=cut{
            let entry=saturated_counts.entry(cut.key()).or_default();entry.0+=1;
            if saturated_sampled<64{write!(saturated_log.as_mut().unwrap(),"{}",g3_sp_row(candidate,Some(cut),&tuple_string_dynamic(n,t))).map_err(|e|e.to_string())?;entry.1+=1;saturated_sampled+=1;}
            return Ok(true);
        }}Ok(false)
    };
    let mut adaptive_log=if refactor.as_ref().is_some_and(|r|r.adaptive_limit>0){
        let mut file=BufWriter::new(File::create(out.join("ADAPTIVE.tsv")).map_err(|e|e.to_string())?);
        writeln!(file,"candidate\tprefix\tseed").map_err(|e|e.to_string())?;Some(file)
    }else{None};
    let mut exploration=if args.iter().any(|s|s=="--proposal-diagnostics"){
        let mut file=BufWriter::new(File::create(out.join("EXPLORATION.tsv")).map_err(|e|e.to_string())?);
        writeln!(file,"candidate\tstates_before_exit\tword_length\tstatus").map_err(|e|e.to_string())?;Some(file)
    }else{None};
    let mut rng=SplitMix64::new(seed);let began=Instant::now();let deadline=began+Duration::from_secs(seconds);
    let literal=args.iter().any(|s|s=="--literal-kernels");
    let literal_trie=args.iter().any(|s|s=="--literal-trie");
    let keep_deck=args.iter().any(|s|s=="--keep-deck");
    let keep_iso=args.iter().any(|s|s=="--keep-isotrivial");
    let mut counts=[0u64;8];let mut stage=[0u64;5];let mut seen=HashSet::new();
    let trace_limit:usize=arg("--trace-early","0").parse().map_err(|_|"trace early bound")?;
    if trace_limit>10000{return Err("early trace limit at most 10000".into());}
    let mut traced=0usize;
    let trie_trace_limit:usize=arg("--trace-trie","0").parse().map_err(|_|"trie trace bound")?;
    if trie_trace_limit>10000{return Err("trie trace limit at most10000".into());}
    let mut trie_traced=0usize;
    let mut trie_hits=std::collections::BTreeMap::<Vec<usize>,u64>::new();
    let mut cyclic=G3CyclicCache::default();
    let mut forest_cuts=G3ForestCutCache::default();
    let use_forest_cuts=args.iter().any(|s|s=="--cyclic-forest-cuts");
    let forest_steps:usize=arg("--forest-cut-steps","16").parse().map_err(|_|"forest cut steps")?;
    if !(1..=64).contains(&forest_steps){return Err("forest cut steps must be1..64".into());}
    let use_cyclic=!args.iter().any(|s|s=="--no-cyclic-screen");
    let early_cyclic=args.iter().any(|s|s=="--cyclic-before-deep");
    let cyclic_trace:usize=arg("--trace-cyclic","0").parse().map_err(|_|"cyclic trace bound")?;
    if cyclic_trace>10000{return Err("cyclic trace limit at most10000".into());}
    let mut cyclic_traced=0usize;
    let mut related_reject=0u64;
    while attempts==0 || counts[0]<attempts {
        if boundary.as_ref().is_some_and(|b|b.cursor==G3_BOUNDARY_TOTAL){break;}
        if exact.as_ref().is_some_and(|e|e.cursor==e.total){break;}
        if (literal||counts[0]&1023==0)&&Instant::now()>=deadline{break;}
        counts[0]+=1;let proposal=if let Some(b)=boundary.as_mut(){b.draw()}else if let Some(q)=four.as_mut(){q.draw(&mut rng)}else if let Some(g)=general.as_mut(){g.draw(&mut rng)}else if let Some(j)=joint.as_mut(){j.draw(&mut rng)}else if let Some(e)=elliptic.as_mut(){e.draw(&mut rng)}else if let Some(r)=refactor.as_mut(){r.draw(&mut rng)}else if let Some(m)=mutation.as_mut(){m.draw(&mut rng)}else if let Some(e)=exact.as_mut(){e.draw()}else if let Some(f)=fiber.as_mut(){f.draw(&mut rng)}else{collision.draw(&mut rng)};
        let Some(t)=proposal else {continue;};
        let ix=if literal {t.map(|p|sm_index(p,n))}else if boundary.is_some(){[7,7,7,15]}else if let Some(q)=&four{q.last_indices}else if let Some(g)=&general{g.last_indices}else if let Some(j)=&joint{j.base.passport}else if let Some(e)=&exact{e.last_indices}else if let Some(f)=&fiber{f.last_indices}else if let Some(e)=&elliptic {e.last_indices}
            else if let Some(r)=&refactor{r.last_indices}
            else if let Some(m)=&mutation{let f=m.frame.as_ref().unwrap();[f.indices[0],f.indices[1],m.third_index,m.fourth_index]}
            else{collision.passport};
        debug_assert_eq!(ix,t.map(|p|sm_index(p,n)));
        if ix.iter().any(|&i|i==0)||ix.iter().sum::<usize>()!=2*n+4 {continue;}counts[1]+=1;
        if exclude_transpositions&&ix.contains(&1){
            theorem_skipped+=1;
            if theorem_sampled<16{let file=theorem_log.as_mut().unwrap();
                writeln!(file,"{}\t{}",counts[0],tuple_string_dynamic(n,&t)).map_err(|e|e.to_string())?;theorem_sampled+=1;
            }continue;
        }
        if saturated_stage=="early"&&saturated.is_some()&&saturated_screen(&t,ix,counts[0])?{continue;}
        let passport=if literal{gd_passport(&t,n)}else{let mut i=ix;i.sort_unstable();i[0]+i[3]<n||i[1]+i[2]<n};
        if !passport {assert!(!finite);continue;}counts[3]+=1;
        let local_frame=if literal{Some(ForestFrame::new([t[0],t[1]],n))}else{None};
        let frame=if let Some(f)=&local_frame{f}else if let Some(b)=&boundary{b.frame()}else if let Some(q)=&four{q.first_frame.as_ref().unwrap()}else if let Some(g)=&general{g.first_frame.as_ref().unwrap()}else if let Some(j)=&joint{j.base.frame.as_ref().unwrap()}else if let Some(e)=&exact{e.first_frame.as_ref().unwrap()}else if let Some(f)=&fiber{f.first_frame.as_ref().unwrap()}else if let Some(e)=&elliptic{e.first_frame.as_ref().unwrap()}
            else if let Some(r)=&refactor{r.first_frame.as_ref().unwrap()}
            else if let Some(m)=&mutation{m.frame.as_ref().unwrap()}else{collision.frame.as_ref().unwrap()};
        if index_profile==Some([6,6,6,12]) {
            if let Some(word)=g3_involution_triangle_word(&t,n){
                if traced<trace_limit{g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&word)?;traced+=1;}
                counts[4]+=1;continue;
            }
        }
        if let Some(j)=frame.short_word(t[2],ix[2],ix[3],n) {
            if finite||(traced<trace_limit&&frame.transitive(t[2],n)){
                g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],FF_WORDS[j])?;traced+=1;
            }counts[4]+=1;continue;
        }
        if saturated_stage=="after-short"&&saturated.is_some()&&saturated_screen(&t,ix,counts[0])?{continue;}
        if use_cyclic&&early_cyclic {if let Some(failure)=cyclic.failure(&t,ix,n) {
            if finite||((traced<trace_limit||cyclic_traced<cyclic_trace)&&frame.transitive(t[2],n)) {
                let word=failure.word(&t,n);
                g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&word)?;traced+=1;cyclic_traced+=1;
            }counts[4]+=1;continue;
        }}
        if let Some(fail)=frame.deep_failure(t[2],ix[2],ix[3],n,3) {
            if finite||(traced<trace_limit&&frame.transitive(t[2],n)){
                g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&frame.deep_word(fail,n))?;traced+=1;
            }counts[4]+=1;continue;
        }
        if !frame.transitive(t[2],n) {assert!(!finite);continue;}counts[2]+=1;
        let failure=if literal_trie{trie.failure_reference(t,n)}else{trie.failure_indices(t,ix,n)};
        if let Some(w)=failure {
            if trie_trace_limit>0{*trie_hits.entry(w.clone()).or_default()+=1;}
            if finite||trie_traced<trie_trace_limit||traced<trace_limit {
                g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&w)?;
                trie_traced+=1;traced+=1;
            }stage[0]+=1;continue;
        }
        if let Some(w)=g3_outer_failure(t,n,32) {
            if finite||traced<trace_limit {
                g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&w)?;traced+=1;
            }stage[1]+=1;continue;
        }
        if !finite&&!keep_deck && g3_has_deck(&t,n) {stage[2]+=1;continue;}
        let key=ff_key(&t,n);if !finite&&seen.contains(&key){counts[5]+=1;continue;}
        if seen.len()<1000000 {seen.insert(key);}
        if !keep_iso&&!finite {if let Some((s,l))=g3_local_isotrivial(&t,ix,n) {
            write!(iso,"{}",g3_iso_row(&t,n,counts[0],s,&l)).map_err(|e|e.to_string())?;stage[4]+=1;continue;
        }}
        if use_forest_cuts {if let Some(word)=forest_cuts.failure(&t,ix,n,forest_steps) {
            if finite||traced<trace_limit||cyclic_traced<cyclic_trace {
                g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&word)?;traced+=1;cyclic_traced+=1;
            }
            counts[4]+=1;continue;
        }}
        if use_cyclic&&!early_cyclic {
            let word=if let Some(w)=cyclic.failure(&t,ix,n) {Some(w.word(&t,n))}else {
                let w=g3_related_cyclic_word(&t,ix,n);related_reject+=u64::from(w.is_some());w
            };
            if let Some(word)=word {
                if finite||traced<trace_limit||cyclic_traced<cyclic_trace {g3_early_negative(&mut leads,&mut resolutions,&t,n,counts[0],&word)?;traced+=1;cyclic_traced+=1;}
                counts[4]+=1;continue;
            }
        }
        if saturated_stage=="late"&&saturated.is_some()&&saturated_screen(&t,ix,counts[0])?{continue;}
        let id=counts[0];let text=tuple_string_dynamic(n,&t);
        writeln!(leads,"{id}\t{text}").map_err(|e|e.to_string())?;leads.flush().map_err(|e|e.to_string())?;
        let limit=Instant::now()+Duration::from_secs(closure);
        let scan=if verifier=="macro" {g3_macro_scan(&t,n,cap,limit)}else{ff_scan(&t,n,cap,limit)};
        let mut folder=String::new();
        if scan.status!="nonsmooth" {
            folder=format!("component_{id:09}");sm_save(&t,n,&scan,&out.join(&folder))?;
            if scan.status=="closed_smooth" {counts[6]+=1;for &key in &scan.states {if seen.len()<1000000{seen.insert(key);}}}
            else{counts[7]+=1;}
        }else{assert!(ff_replay_negative(t,n,&scan.word));}
        if scan.status=="nonsmooth"{if let Some(r)=refactor.as_mut(){if let Some((prefix,template))=r.remember(t,&scan.word){
            let file=adaptive_log.as_mut().expect("adaptive provenance sink");
            writeln!(file,"{id}\t{prefix}\t{}",tuple_string_dynamic(n,&template)).map_err(|e|e.to_string())?;
            file.flush().map_err(|e|e.to_string())?;
        }}}
        if let Some(file)=exploration.as_mut(){writeln!(file,"{id}\t{}\t{}\t{}",scan.states.len(),scan.word.len(),scan.status).map_err(|e|e.to_string())?;}
        stage[3]+=1;
        writeln!(resolutions,"{id}\t{}\t{}\t{folder}\t{text}",scan.status,scan.word.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
        resolutions.flush().map_err(|e|e.to_string())?;
        if Instant::now()>=deadline{break;}
    }
    drop(saturated_screen);
    if retire_frames {
        if let Some(j)=&joint{j.base.retirement.save(&out,n,"joint-genus3")?;}
        else{collision.retirement.save(&out,n,if ["collision","legacy-types"].contains(&mode.as_str()){mode.as_str()}else{"not-applicable"})?;}
    }
    fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\t{mode}\t{seed}\t{}\t{}\n",counts.iter().map(|v|v.to_string()).collect::<Vec<_>>().join("\t"),began.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
    if let Some(e)=&exact{e.save(&out)?;}
    if let Some(b)=&boundary {
        fs::write(out.join("BOUNDARY.tsv"),format!("start\tnext\ttotal\tqualified\tunresolved\n{}\t{}\t{}\t{}\t{}\n",boundary_start.unwrap(),b.cursor,G3_BOUNDARY_TOTAL,counts[1],counts[7])).map_err(|e|e.to_string())?;
    }
    fs::write(out.join("CYCLIC.tsv"),format!("enabled\tstage\ttested\tframes\trejected\trelated_rejected\tpersisted\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n",use_cyclic,if early_cyclic{"early"}else{"late"},cyclic.tested,cyclic.builds,cyclic.rejected,related_reject,cyclic_traced)).map_err(|e|e.to_string())?;
    fs::write(out.join("FOREST_CUTS.tsv"),format!("enabled\tmax_steps\ttested\tframes\trejected\n{}\t{}\t{}\t{}\t{}\n",use_forest_cuts,forest_steps,forest_cuts.tested,forest_cuts.builds,forest_cuts.rejected)).map_err(|e|e.to_string())?;
    if trie_trace_limit>0 {
        let mut text=String::from("word\thits\n");
        for (word,hits) in trie_hits {text.push_str(&format!("{}\t{hits}\n",word.iter().map(|x|x.to_string()).collect::<Vec<_>>().join(",")));}
        fs::write(out.join("TRIE_HITS.tsv"),text).map_err(|e|e.to_string())?;
    }
    if let Some(j)=&joint {if j.csp {
        fs::write(out.join("JOINT.tsv"),format!("nodes\tfound\texhausted_frames\tbudget_stops\tlimit\n{}\t{}\t{}\t{}\t{}\n",j.nodes,j.found,j.exhausted,j.capped,j.budget)).map_err(|e|e.to_string())?;
    }}
    fs::write(out.join("STAGES.tsv"),format!("trie_reject\touter_reject\tdeck_outside_target\tclosure_calls\tisotrivial_outside_target\trng_final\n{}\t{rng}\n",stage.iter().map(|v|v.to_string()).collect::<Vec<_>>().join("\t"),rng=rng.state)).map_err(|e|e.to_string())?;
    if let Some(r)=refactor.as_ref(){
        fs::write(out.join("ADAPTIVE_COUNTS.tsv"),format!("capacity\tadmitted\tadaptive_frame_draws\n{}\t{}\t{}\n",r.adaptive_limit,r.adaptive.len(),r.adaptive_draws)).map_err(|e|e.to_string())?;
    }
    if exclude_transpositions{
        fs::write(out.join("TARGET_THEOREM.tsv"),format!("theorem\texcluded_occurrences\tsampled_occurrences\tavoided_sorted_profiles\nno-positive-genus-all-smooth-transposition\t{theorem_skipped}\t{theorem_sampled}\t{avoided_profiles}\n")).map_err(|e|e.to_string())?;
    }
    if saturated.is_some(){
        let excluded:u64=saturated_counts.values().map(|v|v.0).sum();
        fs::write(out.join("SATURATED_SUMMARY.tsv"),format!("enabled\tchecked_scope_occurrences\ttheorem_excluded_occurrences\tsampled_occurrences\ttrace_records\tstage\n{use_saturated}\t{saturated_checked}\t{excluded}\t{saturated_sampled}\t{saturated_traced}\t{saturated_stage}\n")).map_err(|e|e.to_string())?;
        let mut text=String::from("method\tk\tprime\tmodulus\treason\texcluded_occurrences\tsampled_occurrences\n");
        for((method,k,prime,modulus,reason),(excluded,sampled))in saturated_counts{
            let c=G3SaturatedCut{special:0,partner:1,k,rule:G3SaturatedRule{method,prime,modulus},reason,violator:u8::MAX};
            text.push_str(&format!("{}\t{k}\t{prime}\t{modulus}\t{}\t{excluded}\t{sampled}\n",c.method(),c.reason()));
        }fs::write(out.join("SATURATED_THEOREMS.tsv"),text).map_err(|e|e.to_string())?;
    }
    println!("FINISHED degree={n} lane={mode} constructed={} closed={} unresolved={}",counts[0],counts[6],counts[7]);Ok(())
}

#[cfg(test)]
mod genus_three_tests {
    use super::*;
    #[test]fn joined_refactor_has_all_rank_two_and_three_cycle_shapes(){
        let n=15;
        let t=sm_parse("[(1 3 4 2 5), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 12)(2 8 15)(3 10 14)(4 7 11)(5 6 13), (1 8 14)(2 10 13)(3 7 15)(4 6 12)(5 9 11)]",n).unwrap();
        let mut r=G3Refactor{n,seeds:vec![t],current:None,used:128,first_frame:None,last_indices:[0;4],cross:true,joined:true,min_index:2,max_index:3,
            adaptive:Vec::new(),adaptive_keys:HashSet::new(),adaptive_limit:0,adaptive_draws:0};
        let mut rng=SplitMix64::new(143612);let mut types=HashSet::new();
        for _ in 0..10000{if let Some(v)=r.draw(&mut rng){
            let ix=v.map(|p|sm_index(p,n));assert!([2,3].contains(&ix[0]));assert!(tuple_product_is_identity(n,&v));
            for j in[1,2,3]{assert_eq!(ff_product_index(v[0],v[j],n),ix[0]+ix[j]);}
            types.insert(sp_type(v[0],n).into_iter().filter(|&x|x>1).collect::<Vec<_>>());
        }}
        for shape in[vec![3],vec![2,2],vec![4],vec![2,3],vec![2,2,2]]{assert!(types.contains(&shape),"missing{shape:?},got{types:?}");}
    }
    #[test]fn recycled_prefix_is_never_assumed_everywhere_smooth(){
        let n=12;
        let t=sm_parse("[(8 9), (1 12 7 3 2)(4 5 6 11 10), (1 8 11 3 10)(2 4 9 12 5)(6 7), (1 5 2 11 7 9)(3 6 12 10 8 4)]",n).unwrap();
        let word=[0,2,4,0,0,0,0,0,0,4,4];
        let mut r=G3Refactor{n,seeds:vec![t],current:None,used:128,first_frame:None,last_indices:[0;4],cross:true,joined:false,min_index:1,max_index:3,
            adaptive:Vec::new(),adaptive_keys:HashSet::new(),adaptive_limit:1,adaptive_draws:0};
        let(prefix,u)=r.remember(t,&word).unwrap();assert_eq!(prefix,10);
        assert!(ff_smooth(&u,u.map(|p|sm_index(p,n)),n));assert!(r.remember(t,&word).is_none());
        let mut rng=SplitMix64::new(911341);let mut accepted=0;
        for _ in 0..4000{if let Some(v)=r.draw(&mut rng){accepted+=1;
            let ix=v.map(|p|sm_index(p,n));assert_eq!(ix.iter().sum::<usize>(),2*n+4);
            for j in [1,2,3]{assert_eq!(ff_product_index(v[0],v[j],n),ix[0]+ix[j]);}
        }}assert!(accepted>0);assert!(r.adaptive_draws>0);
    }
    #[test]fn conditioned_refactor_preserves_genus_and_three_forests(){
        let n=15;
        let t=sm_parse("[(1 3 4 2 5), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 12)(2 8 15)(3 10 14)(4 7 11)(5 6 13), (1 8 14)(2 10 13)(3 7 15)(4 6 12)(5 9 11)]",n).unwrap();
        for cross in [false,true]{
            let mut r=G3Refactor{n,seeds:vec![t],current:None,used:128,first_frame:None,last_indices:[0;4],cross,joined:false,min_index:1,max_index:4,
                adaptive:Vec::new(),adaptive_keys:HashSet::new(),adaptive_limit:0,adaptive_draws:0};
            let mut rng=SplitMix64::new(329123);let mut ranks=[0usize;5];
            for _ in 0..4000{if let Some(u)=r.draw(&mut rng){
                let ix=u.map(|p|sm_index(p,n));ranks[ix[0]]+=1;
                assert!(tuple_product_is_identity(n,&u));assert_eq!(ix.iter().sum::<usize>(),2*n+4);
                assert_eq!(ff_product_index(u[0],u[1],n),ix[0]+ix[1]);
                if cross{for j in [2,3]{assert_eq!(ff_product_index(u[0],u[j],n),ix[0]+ix[j]);}}
                assert!(u[0].compose(u[0],n).is_identity(n));
            }}
            assert!(ranks[1..].iter().all(|&v|v>0),"cross={cross} ranks={ranks:?}");
        }
    }
    #[test]fn index_profile_scope_and_involution_obstruction_are_literal(){
        assert_eq!(g3_index_profile(13,"6,6,6,12").unwrap(),[6,6,6,12]);
        for raw in ["0,6,6,18","6,6,12,6","1,2,3","1,1,1,27"]{assert!(g3_index_profile(13,raw).is_err());}
        let n=13;let a=Perm{img:std::array::from_fn(|i|if i<n{((n-i)%n)as u8}else{i as u8})};
        let b=Perm{img:std::array::from_fn(|i|if i<n{((n-i-1)%n)as u8}else{i as u8})};let p=a.compose(b,n);
        let mut rng=SplitMix64::new(20813);let mut found=0;
        for _ in 0..5000{let c=g3_type_perm(n,&[1,2,2,2,2,2,2],&mut rng);let d=p.compose(c,n).inverse(n);
            let t=[a,b,c,d];let word=g3_involution_triangle_word(&t,n);
            if sm_index(d,n)==n-1{let w=word.unwrap();assert!(w.len()<=13&&ff_replay_negative(t,n,&w));found+=1;}else{assert!(word.is_none());}
        }assert!(found>100);
    }
    #[test]#[ignore]fn benchmark_cached_trie_indices(){
        let data=fs::read_to_string(std::env::var("G3_TRIE_FIXTURE").unwrap()).unwrap();
        let words=fs::read_to_string(std::env::var("G3_TRIE_WORDS").unwrap()).unwrap();
        let words:Vec<Vec<usize>>=words.lines().filter(|l|!l.is_empty()).map(|l|l.split(',').map(|v|v.parse().unwrap()).collect()).collect();
        let tree=G3Words::new(&words);let cases:Vec<_>=data.lines().map(|l|{let(n,text)=l.split_once('\t').unwrap();let n=n.parse().unwrap();let t=sm_parse(text,n).unwrap();(n,t,t.map(|p|sm_index(p,n)))}).collect();
        assert!(cases.len()>=10);
        for &(n,t,ix)in &cases{assert_eq!(tree.failure_reference(t,n),tree.failure_indices(t,ix,n));}
        for rep in 0..7 {let mut times=[0.0;2];for mode in [rep%2,1-rep%2]{let began=Instant::now();
            for _ in 0..20{for &(n,t,ix)in &cases{
                let word=if mode==0{tree.failure_reference(std::hint::black_box(t),n)}else{tree.failure_indices(std::hint::black_box(t),ix,n)};
                std::hint::black_box(word);
            }}times[mode]=began.elapsed().as_secs_f64();
        }println!("trie_pair {rep} {} {} cases={}",times[0],times[1],cases.len());}
    }
    #[test]fn cached_trie_indices_preserve_first_literal_witness(){
        let words=vec![vec![1,1,4,0,1],vec![4,0,4],vec![0,1,3,1],vec![2,4,1,0,5]];
        let tree=G3Words::new(&words);let mut rng=SplitMix64::new(2026102209);
        for n in [4,5,9,12,18,25,MAX_D.min(63)]{for _ in 0..500{
            let a=gd_shuffle(n,&mut rng);let b=gd_shuffle(n,&mut rng);let c=gd_shuffle(n,&mut rng);
            let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];let ix=t.map(|p|sm_index(p,n));
            let old=tree.failure_reference(t,n);assert_eq!(old,tree.failure_indices(t,ix,n));
            if let Some(word)=old{assert!(ff_replay_negative(t,n,&word));}
        }}
    }
    #[test]fn joint_quota_constructor_satisfies_all_literal_short_words() {
        let mut rng=SplitMix64::new(592276);let mut accepted=0;
        for n in [12,15,18,21,24,25] {
            let mut source=G3Joint::new(n);let mut local=0;
            for _ in 0..3000 {if let Some(t)=source.draw(&mut rng) {
                let ix=t.map(|p|sm_index(p,n));local+=1;
                assert!(tuple_product_is_identity(n,&t));assert_eq!(ix.iter().sum::<usize>(),2*n+4);
                assert!(ix.iter().all(|&i|i>0));
                assert_eq!(ff_product_index(t[0],t[1],n),ix[0]+ix[1]);
                for word in FF_WORDS {
                    let mut u=t;for &mv in word {u=hurwitz_move_dynamic(n,u,mv);}
                    let indices=u.map(|p|sm_index(p,n));let ip=sm_index(u[0].compose(u[1],n),n);
                    assert!(ip==indices[0]+indices[1]||ip==indices[2]+indices[3]);
                }
            }}assert!(local>0,"no joint-quota proposals in degree {n}");accepted+=local;
        }assert!(accepted>100);
    }
    #[test]fn joint_quota_failure_is_not_a_nonexistence_claim() {
        let mut rng=SplitMix64::new(15119);let p=Perm::identity(9);
        let frame=ForestFrame::new([p,p],9);
        assert!(g3_joint_third(&frame,9,1,8,&mut rng).is_none());
        // An impossible endpoint budget returns None; it neither journals a
        // braid witness nor claims that another initial pair cannot succeed.
    }
    #[test]fn finite_s6_scope_and_transitivity_are_explicit() {
        for(n,kind,expected_frames)in[(18,"C3",20),(24,"C4",44),(24,"S4",44)] {
            let mut exact=G3FiberExact::new(n,kind,0).unwrap();
            assert_eq!(exact.k,6);assert_eq!(exact.perms.len(),720);
            assert_eq!(exact.frames.len(),expected_frames);
            assert_eq!(exact.total,expected_frames*720*720);
            for _ in 0..5000 {if let Some(t)=exact.draw() {
                assert!(generators_transitive(n,&t));assert!(tuple_product_is_identity(n,&t));
                assert_eq!(t.map(|p|sm_index(p,n)),exact.last_indices);
                assert_eq!(exact.last_indices.iter().sum::<usize>(),2*n+4);
            }}
            let mut end=G3FiberExact::new(n,kind,exact.total).unwrap();assert!(end.draw().is_none());
            assert!(G3FiberExact::new(n,kind,exact.total+1).is_err());
        }
        assert!(G3FiberExact::new(21,"C3",0).is_err());
        assert!(G3FiberExact::new(18,"C4",0).is_err());
    }
    #[test]fn stack_conjugators_match_literal_cycles_and_rng() {
        fn literal(x:Perm,y:Perm,k:usize,rng:&mut SplitMix64)->Option<Perm>{
            let mut xc=qf_cycles(x,k);let mut yc=qf_cycles(y,k);
            xc.sort_by_key(|c|c.len());yc.sort_by_key(|c|c.len());
            if xc.iter().map(|c|c.len()).ne(yc.iter().map(|c|c.len())){return None;}
            let mut start=0;while start<yc.len(){let mut end=start+1;while end<yc.len()&&yc[end].len()==yc[start].len(){end+=1;}rng.shuffle(&mut yc[start..end]);start=end;}
            let mut v=Perm::identity(k);for(a,b)in xc.iter().zip(&yc){let offset=rng.below(a.len());for j in 0..a.len(){v.img[a[j]]=b[(j+offset)%b.len()]as u8;}}Some(v)
        }
        for(x,_)in sc_perms(5){for(y,_)in sc_perms(5){
            let mut a=SplitMix64::new(782);let mut b=SplitMix64::new(782);
            assert_eq!(literal(x,y,5,&mut a),g3_conjugator(x,y,5,&mut b));assert_eq!(a.state,b.state);
        }}
        let mut data=SplitMix64::new(1439);for k in 6..=8{for _ in 0..10000{
            let x=gd_shuffle(k,&mut data);let y=gd_shuffle(k,&mut data);let seed=data.next_u64();
            let mut a=SplitMix64::new(seed);let mut b=SplitMix64::new(seed);
            assert_eq!(literal(x,y,k,&mut a),g3_conjugator(x,y,k,&mut b));assert_eq!(a.state,b.state);
        }}
    }
    #[test]fn partitions_and_fixed_product_factorizations_are_exact() {
        let mut rng=SplitMix64::new(20261022);
        for n in 9..=25 {
            let types=G3Types::new(n);
            for e in 1..n {for a in 0..=e {
                let h=types.permutation(n,e,&mut rng);assert_eq!(sm_index(h,n),e);
                let pair=g3_factor(h,n,a,&mut rng).unwrap();assert_eq!(pair[0].compose(pair[1],n),h);
                assert_eq!([sm_index(pair[0],n),sm_index(pair[1],n)],[a,e-a]);
            }}
            for ix in g3_passports(n) {assert_eq!(ix.iter().sum::<usize>(),2*n+4);assert!(ix.windows(2).all(|p|p[0]<=p[1]));}
        }
    }
    #[test]fn outer_identity_and_literal_negative_filters() {
        let mut rng=SplitMix64::new(913);
        for n in 4..=25 {for _ in 0..100 {
            let a=gd_shuffle(n,&mut rng);let b=gd_shuffle(n,&mut rng);let c=gd_shuffle(n,&mut rng);
            let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];
            if !generators_transitive(n,&t) {continue;}
            let mut u=t;for m in [0,0,5,5] {u=hurwitz_move_dynamic(n,u,m);}
            let h=a.compose(b,n);assert_eq!(u,t.map(|p|p.conjugate_by(h,n)));
            if let Some(w)=g3_outer_failure(t,n,8) {assert!(ff_replay_negative(t,n,&w));}
            let trie=G3Words::new(&[vec![0,1,4],vec![0,1,2],vec![4,4,1]]);
            if let Some(w)=trie.failure(t,n) {assert!(ff_replay_negative(t,n,&w));}
        }}
    }
    #[test]fn deck_equivariance_matches_literal_s4_centralizers() {
        let all=sc_perms(4);
        for &(a,_) in &all {for &(b,_) in &all {
            let t=[a,b,a,b];if !generators_transitive(4,&t){continue;}
            let literal=all.iter().any(|(h,_)|!h.is_identity(4)&&[a,b].iter().all(|&p|p.compose(*h,4)==h.compose(p,4)));
            assert_eq!(g3_has_deck(&t,4),literal);
        }}
    }
    #[test]fn macro_equals_literal_complete_s4_orbits_and_negative_witnesses() {
        let all=sc_perms(4);let mut checked=0;
        for &(a,_) in &all {for &(b,_) in &all {for &(c,_) in &all {
            let t=[a,b,c,a.compose(b,4).compose(c,4).inverse(4)];
            if !generators_transitive(4,&t)||!sm_values(&t,4).0 {continue;}
            let deadline=Instant::now()+Duration::from_secs(10);
            let a=ff_scan(&t,4,10000,deadline);let b=g3_macro_scan(&t,4,10000,deadline);
            assert_eq!(a.status,b.status);
            if b.status=="closed_smooth" {assert_eq!(a.states.iter().copied().collect::<HashSet<_>>(),b.states.iter().copied().collect());}
            else {assert!(ff_replay_negative(t,4,&b.word));}checked+=1;
        }}}assert!(checked>100);
    }
    #[test]fn macro_positive_caps_and_targeted_constructor_scope() {
        let t=sm_parse("[(5 6 7), (1 2 5)(3 7 9)(4 6 8), (1 3 7)(2 6 9)(4 5 8), (1 4 5)(2 7 8)(3 6 9)]",9).unwrap();
        let d=Instant::now()+Duration::from_secs(10);assert_eq!(g3_macro_scan(&t,9,1000,d).states.len(),12);
        assert_eq!(g3_macro_scan(&t,9,1,d).status,"capped");
        assert_eq!(g3_macro_scan(&t,9,100,Instant::now()).status,"timeout");
        let mut rng=SplitMix64::new(77223);let mut found=0;
        for n in [12,15,18,20,21,24,25] {
            let mut c=G3Collision::new(n,true);let mut e=G3Elliptic::new(n);let mut r=G3Elliptic::with_genus(n,0);
            for _ in 0..2000 {for t in [c.draw(&mut rng),e.draw(&mut rng),r.draw(&mut rng)].into_iter().flatten() {
                found+=1;assert!(tuple_product_is_identity(n,&t));
                assert_eq!(t.iter().map(|p|sm_index(*p,n)).sum::<usize>(),2*n+4);
                assert_eq!(ff_product_index(t[0],t[1],n),sm_index(t[0],n)+sm_index(t[1],n));
            }}
        }assert!(found>100);
    }
    #[test]fn mutations_preserve_genus_geodesy_and_break_the_old_blocks() {
        let n=15;let t=sm_parse("[(1 3 4 2 5), (1 6 11)(2 7 12)(3 8 13)(4 9 14)(5 10 15), (1 9 12)(2 8 15)(3 10 14)(4 7 11)(5 6 13), (1 8 14)(2 10 13)(3 7 15)(4 6 12)(5 9 11)]",n).unwrap();
        let labels:Vec<_>=(0..n).map(|i|i/5).collect();assert!(t.iter().all(|&p|g3_preserves(p,&labels,n)));
        let (special,found)=g3_local_isotrivial(&t,t.map(|p|sm_index(p,n)),n).unwrap();
        assert_eq!(special,0);assert_eq!(&found[..n],labels.as_slice());
        let mut m=G3Mutation{n,seeds:vec![G3MutationSeed{tuple:t,partitions:vec![labels.clone()]}],current:0,
            frame:None,third:t[2],third_index:sm_index(t[2],n),fourth_index:sm_index(t[3],n),used:128};
        let mut rng=SplitMix64::new(3321);let mut found=0;
        for _ in 0..10000 {if let Some(u)=m.draw(&mut rng) {
            found+=1;assert!(tuple_product_is_identity(n,&u));
            assert_eq!(u.iter().map(|&p|sm_index(p,n)).sum::<usize>(),2*n+4);
            assert_eq!(sm_index(u[0],n)+sm_index(u[1],n),ff_product_index(u[0],u[1],n));
            assert!(!g3_preserves(u[2],&labels,n));
            let f=m.frame.as_ref().unwrap();assert_eq!(u.map(|p|sm_index(p,n)),[f.indices[0],f.indices[1],m.third_index,m.fourth_index]);
        }}assert!(found>10);
    }
    #[test]fn full_fibre_constructor_has_exact_genus_and_every_block_active() {
        let mut rng=SplitMix64::new(333718);let mut total=0;
        for n in [15,18,20,21,24] {
            let mut f=G3FullFiber::new(n);let mut count=0;
            for _ in 0..100000 {if let Some(t)=f.draw(&mut rng) {
                count+=1;let frame=f.frame.as_ref().unwrap();let(r,k)=(frame.tree.r,frame.tree.k);
                assert!(tuple_product_is_identity(n,&t));assert_eq!(t.map(|p|sm_index(p,n)),f.last_indices);
                assert_eq!(f.last_indices.iter().sum::<usize>(),2*n+4);
                assert_eq!(ff_product_index(t[0],t[1],n),sm_index(t[0],n)+sm_index(t[1],n));
                for j in 0..r{assert!((0..k).any(|x|t[0].img[j*k+x] as usize!=j*k+x));}
                for entry in 1..4{
                    let mut q=Perm::identity(r);for j in 0..r{q.img[j]=(t[entry].img[j*k]as usize/k)as u8;}
                    assert_eq!(sm_index(t[entry],n),k*sm_index(q,r));
                }
                assert!(g3_local_isotrivial(&t,f.last_indices,n).is_none());
            }}assert!(count>0,"no exact fiber proposals at degree {n}");total+=count;
        }assert!(total>100);
    }
    #[test]fn full_fibre_conjugacy_reduction_equals_literal_s5_voltage_pairs() {
        let mut rng=SplitMix64::new(847131);let mut count=0;
        for (n,kind) in [(15,"C3"),(20,"C4"),(20,"S4")] {
            let e=G3FiberExact::new(n,kind,0).unwrap();assert!(!e.frames.is_empty());
            for frame in &e.frames {
                let q=&frame.tree;assert_eq!(q.free_edges,if kind=="C4"{vec![2,3]}else{vec![1,2]});
                let targets=vec![Perm::identity(5);q.orders.len()];
                let target=5*sm_index(q.beta.compose(q.gamma,q.r),q.r);
                for &u in &e.perms {
                    let(x,y)=g3_fiber_conjugacy(frame,u);
                    let v=g3_conjugator(x,y,5,&mut rng);
                    let mut hits=0;
                    for &v in &e.perms {
                        let t=frame.solve(&[u,v],&targets);
                        let right=sm_index(t[3],n)==target;
                        assert_eq!(v.compose(x,5).compose(v.inverse(5),5)==y,right);
                        hits+=usize::from(right);count+=1;
                    }
                    assert_eq!(v.is_some(),hits>0);
                    if let Some(v)=v {assert_eq!(sm_index(frame.solve(&[u,v],&targets)[3],n),target);}
                }
            }
        }assert!(count>100000);
    }
}
