// Exact symbolic elimination of return equations. Constants and free voltages
// form words in S_k; no commutativity or sampled identity is assumed.
#[derive(Clone,Copy)]
enum ReturnAtom {Constant(Perm),Variable(usize,bool)}
fn return_inverse(word:&[ReturnAtom],k:usize)->Vec<ReturnAtom>{
    word.iter().rev().map(|a|match *a {
        ReturnAtom::Constant(p)=>ReturnAtom::Constant(p.inverse(k)),
        ReturnAtom::Variable(j,inverse)=>ReturnAtom::Variable(j,!inverse),
    }).collect()
}
fn return_value(word:&[ReturnAtom],free:&[Perm],k:usize)->Perm {
    word.iter().fold(Perm::identity(k),|p,a|p.compose(match *a {
        ReturnAtom::Constant(c)=>c,
        ReturnAtom::Variable(j,false)=>free[j],
        ReturnAtom::Variable(j,true)=>free[j].inverse(k),
    },k))
}
struct ReturnConjugacy {word:Vec<ReturnAtom>,variable:usize,first:usize,second:usize,inverted:bool,targets:Vec<Perm>}
impl ReturnConjugacy {
    fn new(frame:&ReturnForestFrame)->Option<Self> {
        Self::with_targets(frame,&vec![Perm::identity(frame.tree.k);frame.tree.orders.len()])
    }
    fn with_targets(frame:&ReturnForestFrame,targets:&[Perm])->Option<Self>{
        let f=&frame.tree;let k=f.k;let ng=qf_cycles(f.gamma,f.r).len();let delta=f.beta.compose(f.gamma,f.r);
        if targets.len()!=f.orders.len(){return None;}
        let mut voltage:Vec<Option<Vec<ReturnAtom>>>=vec![None;f.r];
        for(j,&e)in f.free_edges.iter().enumerate(){voltage[e]=Some(vec![ReturnAtom::Variable(j,false)]);}
        let factor=|v:usize,e:usize,voltage:&[Option<Vec<ReturnAtom>>]|{
            let mut word=Vec::new();
            if v>=ng{word.push(ReturnAtom::Constant(frame.local[delta.img[e]as usize]));}
            word.extend_from_slice(voltage[e].as_ref().unwrap());word
        };
        for &(v,e)in &f.leaves {
            let order=&f.orders[v];let cut=order.iter().position(|&j|j==e).unwrap();
            let mut before=Vec::new();let mut after=Vec::new();
            for &j in &order[..cut]{let mut next=factor(v,j,&voltage);next.extend(before);before=next;}
            for &j in &order[cut+1..]{let mut next=factor(v,j,&voltage);next.extend(after);after=next;}
            let coefficient=if v<ng{Perm::identity(k)}else{frame.local[delta.img[e]as usize]};
            let mut word=vec![ReturnAtom::Constant(coefficient.inverse(k))];
            word.extend(return_inverse(&after,k));word.push(ReturnAtom::Constant(targets[v]));
            word.extend(return_inverse(&before,k));voltage[e]=Some(word);
        }
        let mut word=Vec::new();
        for &e in &f.orders[f.root]{let mut next=factor(f.root,e,&voltage);next.extend(word);word=next;}
        word.insert(0,ReturnAtom::Constant(targets[f.root].inverse(k)));
        // In an orientable return graph a suitable non-tree edge occurs twice,
        // with opposite signs. Refuse unsupported words rather than guessing.
        for variable in 0..f.free_edges.len(){
            let positions:Vec<_>=word.iter().enumerate().filter_map(|(i,a)|match a {
                ReturnAtom::Variable(j,sign)if *j==variable=>Some((i,*sign)),_=>None,
            }).collect();
            if positions.len()==2&&positions[0].1!=positions[1].1 {
                return Some(Self{word,variable,first:positions[0].0,second:positions[1].0,inverted:positions[0].1,targets:targets.to_vec()});
            }
        }None
    }
    fn solve(&self,frame:&ReturnForestFrame,free:&mut[Perm],rng:&mut SplitMix64)->Option<[Perm;4]>{
        let k=frame.tree.k;
        let left=return_value(&self.word[..self.first],free,k);
        let x=return_value(&self.word[self.first+1..self.second],free,k);
        let right=return_value(&self.word[self.second+1..],free,k);
        // L*v^s*X*v^-s*R=1 iff v^s*X*v^-s=L^-1*R^-1.
        let y=left.inverse(k).compose(right.inverse(k),k);
        let v=g3_conjugator(x,y,k,rng)?;
        free[self.variable]=if self.inverted{v.inverse(k)}else{v};
        debug_assert!(return_value(&self.word,free,k).is_identity(k));
        Some(frame.solve(free,&self.targets))
    }
}

// A rational quotient with FOUR nonidentity entries, not the old identity-entry
// lift family. The first quotient pair is a forest and the third is a full
// cycle. Prescribed relative returns force total excess2k+4, hence source genus3.
struct G3FourQuotient {
    n:usize,pair:Option<[Perm;2]>,frame:Option<ReturnForestFrame>,solver:Option<ReturnConjugacy>,
    first_frame:Option<ForestFrame>,used:usize,last_indices:[usize;4],types:Vec<G3Types>,
    profiles:Vec<(usize,usize,usize,usize,usize)>,solved:bool,
}
fn return_excess_rows(mut pair:[Perm;2],d:Perm,cycles:&[Vec<usize>],r:usize)->Vec<u64>{
    let mut seen=HashSet::new();let mut rows=HashSet::new();
    while seen.insert(pair){
        for a in pair{
            let mut parent=std::array::from_fn(|i|i);
            for p in[a,d]{for x in 0..r{gd_join(&mut parent,x,p.img[x]as usize);}}
            let mut masks=[0u64;MAX_D];
            for(i,c)in cycles.iter().enumerate(){masks[gd_find(&mut parent,c[0])]|=1<<i;}
            for mask in masks{if mask!=0{rows.insert(mask);}}
        }
        pair=[pair[0].compose(pair[1],r).compose(pair[0].inverse(r),r),pair[0]];
    }
    // With nonnegative excess, a subset row is implied by any containing row
    // with the same cap. Removing it also sharpens the summed fractional bound.
    let all:Vec<_>=rows.into_iter().collect();
    let mut out:Vec<_>=all.iter().copied().filter(|&m|!all.iter().any(|&other|other!=m&&m&other==m)).collect();
    out.sort_unstable();out
}
fn return_excess_budget(rows:&[u64],cells:usize,budget:usize,cap:usize,rng:&mut SplitMix64)->Option<Vec<usize>>{
    // Summing exact inequalities gives this cheap rigorous upper bound.
    let cover=(0..cells).map(|i|rows.iter().filter(|&&m|m&(1<<i)!=0).count()).min()?;
    if cover>0&&budget*cover>rows.len()*cap{return None;}
    let mut values=vec![0;cells];let mut used=vec![0;rows.len()];
    for _ in 0..budget{
        let choices:Vec<_>=(0..cells).filter(|&i|values[i]<cap&&rows.iter().enumerate()
            .all(|(j,&m)|m&(1<<i)==0||used[j]<cap)).collect();
        if choices.is_empty(){return None;}let i=choices[rng.below(choices.len())];values[i]+=1;
        for(j,&m)in rows.iter().enumerate(){if m&(1<<i)!=0{used[j]+=1;}}
    }Some(values)
}
impl G3FourQuotient {
    fn new(n:usize)->Self{
        let mut profiles=Vec::new();
        for r in 4..=n/3{if n%r!=0{continue;}let k=n/r;
            for e in 2..r-1{for a in 1..e{for ec in 0..k{
                let ed=2*k+4-ec;
                // Exact decorated four-block quotient audit: all three pair
                // sums of the D-cycle excesses are <=k-1. Hence 2ed<=3(k-1).
                // This is only the marked (1,1,3,1), pure-A/B construction.
                if r==4&&2*ed>3*(k-1){continue;}
                // Every additive quotient pair orbit can absorb at most k-1
                // excess. The two low-entry/D pairings have a+1 and e-a+1
                // orbits; this is stronger than the source index passport.
                if ed>(a.min(e-a)+1)*(k-1){continue;}
                let mut ix=[k*a,k*(e-a),k*(r-1)+ec,k*(r-1-e)+ed];ix.sort_unstable();
                if ix[3]<n&&(ix[0]+ix[3]<n||ix[1]+ix[2]<n){profiles.push((r,k,e,a,ec));}
            }}}
        }
        Self{n,pair:None,frame:None,solver:None,first_frame:None,used:128,last_indices:[0;4],types:(0..=n/4).map(G3Types::new).collect(),profiles,solved:true}
    }
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]>{
        let n=self.n;
        if self.frame.is_none()||self.used>=128{
            if self.profiles.is_empty(){return None;}
            let(r,k,e,a,ec)=self.profiles[rng.below(self.profiles.len())];
            let gamma=dd_full_cycle(r,rng);
            let cut=g3_factor(gamma,r,e,rng)?;let product=cut[0].inverse(r);
            let pair=g3_factor(product,r,a,rng)?;
            let frame=ReturnForestFrame::new(product,gamma,r,k,vec![Perm::identity(k);r])?;
            if frame.tree.genus!=0{return None;}
            let f=&frame.tree;let mut target=vec![Perm::identity(k);f.orders.len()];
            assert_eq!(qf_cycles(gamma,r).len(),1);target[0]=self.types[k].permutation(k,ec,rng);
            // Include the root among ramified delta returns. The symbolic
            // solver enforces its prescribed value, not merely its index.
            let delta=product.compose(gamma,r);let dc=qf_cycles(delta,r);
            let rows=if r==4{vec![3,5,6]}else{return_excess_rows(pair,delta.inverse(r),&dc,r)};
            let budget=return_excess_budget(&rows,dc.len(),2*k+4-ec,k-1,rng)?;
            for(j,&ix)in budget.iter().enumerate(){target[1+j]=self.types[k].permutation(k,ix,rng);}
            let ix=[k*a,k*(e-a),k*sm_index(gamma,r)+ec,k*sm_index(product.compose(gamma,r),r)+2*k+4-ec];
            let solver=ReturnConjugacy::with_targets(&frame,&target)?;
            let mut lifted=[Perm::identity(n);2];
            for j in 0..2{for block in 0..r{for x in 0..k{lifted[j].img[block*k+x]=(pair[j].img[block]as usize*k+x)as u8;}}}
            self.pair=Some(lifted);self.last_indices=ix;self.frame=Some(frame);self.solver=Some(solver);self.first_frame=None;self.used=0;
        }
        self.used+=1;let frame=self.frame.as_ref().unwrap();let k=frame.tree.k;
        let mut free:Vec<_>=(0..frame.tree.free_edges.len()).map(|_|gd_shuffle(k,rng)).collect();
        let mut t=if self.solved{self.solver.as_ref().unwrap().solve(frame,&mut free,rng)?}
            else{frame.solve(&free,&self.solver.as_ref().unwrap().targets)};
        let pair=self.pair.unwrap();
        t[0]=pair[0];t[1]=pair[1];assert!(tuple_product_is_identity(n,&t));
        if !self.solved&&t.map(|p|sm_index(p,n))!=self.last_indices{return None;}
        assert_eq!(t.map(|p|sm_index(p,n)),self.last_indices);
        if self.first_frame.is_none(){self.first_frame=Some(ForestFrame::new(pair,n));}Some(t)
    }
}

struct G3GeneralFiber {
    n:usize,frame:Option<ReturnForestFrame>,solver:Option<ReturnConjugacy>,
    first_frame:Option<ForestFrame>,used:usize,last_indices:[usize;4],solved:bool,
}
impl G3GeneralFiber {
    fn new(n:usize)->Self{Self{n,frame:None,solver:None,first_frame:None,used:128,last_indices:[0;4],solved:true}}
    fn draw(&mut self,rng:&mut SplitMix64)->Option<[Perm;4]>{
        let n=self.n;
        if self.frame.is_none()||self.used>=128 {
            let shapes:Vec<_>=(3..=n/3).filter(|&r|n%r==0).map(|r|(r,n/r)).collect();
            if shapes.is_empty(){return None;}let(r,k)=shapes[rng.below(shapes.len())];
            let beta=gd_shuffle(r,rng);
            // This lane deliberately removes the full-beta restriction. It
            // admits inner degree3/4, impossible with beta a single cycle.
            if sm_index(beta,r)==r-1||qf_cycles(beta,r).len()*(k-1)<4{return None;}
            let gamma=gd_shuffle(r,rng);
            let tree=QuotientForest::tree_any(beta,gamma,r,k,0,true)?;
            if tree.genus!=1{return None;}
            let mut parent=std::array::from_fn(|i|i);let mut local=vec![Perm::identity(k);r];
            for j in 0..r{for x in 0..k{gd_join(&mut parent,j*k+x,beta.img[j]as usize*k+x);}}
            for _ in 0..4 {
                let roots:[usize;MAX_D]=std::array::from_fn(|x|if x<n{gd_find(&mut parent,x)}else{x});
                let mut choices=[(0u8,0u8,0u8);DD_PAIR_CAPACITY];let mut count=0;
                for j in 0..r{for x in 0..k{for y in x+1..k{if roots[j*k+x]!=roots[j*k+y]{
                    choices[count]=(j as u8,x as u8,y as u8);count+=1;
                }}}}
                if count==0{return None;}let(j,x,y)=choices[rng.below(count)];let(j,x,y)=(j as usize,x as usize,y as usize);
                local[j].img.swap(x,y);assert!(gd_join(&mut parent,j*k+x,j*k+y));
            }
            // With all ordinary returns unramified, one active block is the
            // established isotrivial family. It cannot meet the sole target.
            if local.iter().filter(|p|!p.is_identity(k)).count()==1{return None;}
            let frame=ReturnForestFrame::new(beta,gamma,r,k,local)?;
            let solver=ReturnConjugacy::new(&frame)?;
            self.frame=Some(frame);self.solver=Some(solver);self.first_frame=None;self.used=0;
        }
        self.used+=1;let frame=self.frame.as_ref().unwrap();let f=&frame.tree;
        let mut free:Vec<_>=(0..f.free_edges.len()).map(|_|gd_shuffle(f.k,rng)).collect();
        let t=if self.solved{self.solver.as_ref().unwrap().solve(frame,&mut free,rng)?}
            else{frame.solve(&free,&vec![Perm::identity(f.k);f.orders.len()])};
        let ix=[4,f.k*sm_index(f.beta,f.r),f.k*sm_index(f.gamma,f.r),f.k*sm_index(f.beta.compose(f.gamma,f.r),f.r)];
        if !self.solved&&t.map(|p|sm_index(p,n))!=ix{return None;}
        assert_eq!(t.map(|p|sm_index(p,n)),ix);self.last_indices=ix;
        if self.first_frame.is_none(){self.first_frame=Some(ForestFrame::new([t[0],t[1]],n));}
        Some(t)
    }
}

#[cfg(test)]
mod return_solver_tests {
    use super::*;
    #[test]fn maximal_return_rows_retain_exact_four_block_budget_obstructions(){
        let c=Perm{img:std::array::from_fn(|i|if i<4{((i+1)%4)as u8}else{i as u8})};
        let trans:Vec<_>=sc_perms(4).into_iter().filter_map(|(p,ix)|if ix==1{Some(p)}else{None}).collect();
        let mut marked=0;let mut overlap=0;let mut disjoint=0;
        for &a in &trans{for &b in &trans{
            let ab=a.compose(b,4);if sm_index(ab,4)!=2{continue;}
            let d=ab.compose(c,4).inverse(4);if sm_index(d,4)!=1{continue;}
            marked+=1;let dc=qf_cycles(d,4);let rows=return_excess_rows([a,b],d,&dc,4);
            if a.compose(b,4)==b.compose(a,4){disjoint+=1;assert_eq!(rows.len(),2);}
            else{overlap+=1;assert_eq!(rows,vec![3,5,6]);
                let mut rng=SplitMix64::new(2);assert!(return_excess_budget(&rows,3,17,11,&mut rng).is_none());
                // Above degree52 this obstruction disappears; sampling failure
                // is still not a nonexistence result.
                let some=(0..100).any(|_|return_excess_budget(&rows,3,18,12,&mut rng).is_some());assert!(some);
            }
        }}assert_eq!((marked,overlap,disjoint),(16,12,4));
        for n in [12,16,20,24,28,32]{assert!(!G3FourQuotient::new(n).profiles.iter().any(|p|p.0==4));}
        if MAX_D>=52{assert!(G3FourQuotient::new(52).profiles.iter().any(|p|p.0==4));}
    }
    #[test]fn symbolic_equation_matches_literal_voltage_elimination() {
        let mut rng=SplitMix64::new(931122);let mut frames=0;let mut completions=0;
        for r in 3..=7{for _ in 0..1000 {
            let beta=gd_shuffle(r,&mut rng);let gamma=gd_shuffle(r,&mut rng);
            let Some(frame)=ReturnForestFrame::new(beta,gamma,r,3,vec![Perm::identity(3);r])else{continue;};
            if frame.tree.genus!=1{continue;}let Some(solver)=ReturnConjugacy::new(&frame)else{continue;};
            frames+=1;let mut free:Vec<_>=(0..frame.tree.free_edges.len()).map(|_|gd_shuffle(3,&mut rng)).collect();
            let x=return_value(&solver.word[solver.first+1..solver.second],&free,3);
            let left=return_value(&solver.word[..solver.first],&free,3);
            let right=return_value(&solver.word[solver.second+1..],&free,3);
            let y=left.inverse(3).compose(right.inverse(3),3);let mut matches=0;
            for(v,_)in sc_perms(3){free[solver.variable]=v;
                let literal=frame.solve(&free,&vec![Perm::identity(3);frame.tree.orders.len()]);
                let root_zero=sm_index(literal[2],r*3)==3*sm_index(gamma,r)&&
                    sm_index(literal[3],r*3)==3*sm_index(beta.compose(gamma,r),r);
                let vs=if solver.inverted{v.inverse(3)}else{v};
                assert_eq!(vs.compose(x,3).compose(vs.inverse(3),3)==y,root_zero);
                assert_eq!(return_value(&solver.word,&free,3).is_identity(3),root_zero);
                matches+=usize::from(root_zero);
            }
            let solved=solver.solve(&frame,&mut free,&mut rng);assert_eq!(solved.is_some(),matches>0);
            completions+=matches;if frames>=40{break;}
        }if frames>=40{break;}}assert_eq!(frames,40);assert!(completions>0);
    }
    #[test]fn general_fibre_keeps_exact_genus_and_reaches_small_inner_actions() {
        let mut rng=SplitMix64::new(20261023);let mut accepted=0;let mut inner_small=0;
        for n in [12,15,16,18,20,21,24]{let mut g=G3GeneralFiber::new(n);
            for _ in 0..30000{if let Some(t)=g.draw(&mut rng){
                accepted+=1;let f=&g.frame.as_ref().unwrap().tree;
                assert_eq!(f.genus,1);assert!(qf_cycles(f.beta,f.r).len()>1);
                assert!(tuple_product_is_identity(n,&t));assert_eq!(t.map(|p|sm_index(p,n)),g.last_indices);
                assert_eq!(g.last_indices.iter().sum::<usize>(),2*n+4);
                assert_eq!(ff_product_index(t[0],t[1],n),g.last_indices[0]+g.last_indices[1]);
                inner_small+=usize::from(f.k<=4);
            }}
        }assert!(accepted>100);assert!(inner_small>100);
    }
    #[test]fn prescribed_ramified_returns_match_literal_solver(){
        let mut rng=SplitMix64::new(99174);let mut checked=0;
        for _ in 0..1000{
            let gamma=dd_full_cycle(6,&mut rng);let cut=g3_factor(gamma,6,4,&mut rng).unwrap();let beta=cut[0].inverse(6);
            let frame=ReturnForestFrame::new(beta,gamma,6,3,vec![Perm::identity(3);6]).unwrap();
            let targets:Vec<_>=(0..frame.tree.orders.len()).map(|_|gd_shuffle(3,&mut rng)).collect();
            let Some(solver)=ReturnConjugacy::with_targets(&frame,&targets)else{continue;};
            let mut free:Vec<_>=(0..frame.tree.free_edges.len()).map(|_|gd_shuffle(3,&mut rng)).collect();
            let mut matches=0;
            for(v,_)in sc_perms(3){free[solver.variable]=v;
                let t=frame.solve(&free,&targets);let symbolic=return_value(&solver.word,&free,3).is_identity(3);
                // Directly follow the final return around the chosen delta
                // cycle; this checks a nonidentity ROOT target as well.
                let x=frame.tree.orders[frame.tree.root][0]*3;
                let mut returned=Perm::identity(3);
                for point in 0..3{let mut y=x+point;
                    for _ in 0..frame.tree.orders[frame.tree.root].len(){y=t[0].compose(t[1],18).compose(t[2],18).img[y]as usize;}
                    returned.img[point]=(y%3)as u8;
                }
                // The chosen cycle may start away from block zero, but the
                // solver's target uses that very same sorted cycle start.
                assert_eq!(symbolic,returned==targets[frame.tree.root]);matches+=usize::from(symbolic);
            }
            assert_eq!(solver.solve(&frame,&mut free,&mut rng).is_some(),matches>0);checked+=1;
            if checked==20{break;}
        }assert_eq!(checked,20);
    }
    #[cfg(smooth_wide_keys)]
    #[test]fn four_nonidentity_quotient_entries_have_genus_zero_and_source_three(){
        let mut rng=SplitMix64::new(66201);let mut accepted=0;
        for n in [35,42,52,56,60,63]{
            if n>MAX_D{continue;}let mut q=G3FourQuotient::new(n);
            for _ in 0..15000{if let Some(t)=q.draw(&mut rng){
                accepted+=1;assert!(tuple_product_is_identity(n,&t));
                assert_eq!(t.map(|p|sm_index(p,n)),q.last_indices);assert_eq!(q.last_indices.iter().sum::<usize>(),2*n+4);
                let f=&q.frame.as_ref().unwrap().tree;
                let mut qi=[0usize;4];for j in 0..4{
                    let mut b=Perm::identity(f.r);for x in 0..f.r{b.img[x]=(t[j].img[x*f.k]as usize/f.k)as u8;}
                    assert!(!b.is_identity(f.r));qi[j]=sm_index(b,f.r);
                }assert_eq!(qi.iter().sum::<usize>(),2*f.r-2);
                assert_eq!(ff_product_index(t[0],t[1],n),q.last_indices[0]+q.last_indices[1]);
            }}
        }assert!(accepted>100);
    }
}
