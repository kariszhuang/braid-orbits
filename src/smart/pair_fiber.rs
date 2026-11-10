//! Elementary fixed-product moves, with exact orbit/genus/parity checks.
//! This is candidate construction, never a completeness or failure test.
use super::*;

fn cycle_labels(d:usize,p:Perm)->[u8;MAX_D] {
    let mut labels=[u8::MAX;MAX_D];
    for i in 0..d {if labels[i]==u8::MAX {
        let mut x=i;while labels[x]==u8::MAX {labels[x]=i as u8;x=p.img[x] as usize;}
    }}
    labels
}

// The full mask at every point makes partition equality independent of roots.
pub(super) fn pair_partition(d:usize,a:Perm,b:Perm)->[u64;MAX_D] {
    let mut masks=[0;MAX_D];let mut used=0u64;
    for root in 0..d {if used&(1u64<<root)==0 {
        let mut queue=[0usize;MAX_D];queue[0]=root;let (mut head,mut tail)=(0,1);
        let mut mask=1u64<<root;
        while head<tail {let x=queue[head];head+=1;
            for p in [a,b] {let y=p.img[x] as usize;let bit=1u64<<y;
                if mask&bit==0 {mask|=bit;queue[tail]=y;tail+=1;}
            }
        }
        used|=mask;for &x in &queue[..tail] {masks[x]=mask;}
    }}
    masks
}

pub(super) fn cut(d:usize,a:Perm,b:Perm,i:usize,j:usize)->(Perm,Perm) {
    let inverse=a.inverse(d);let (u,v)=(inverse.img[i] as usize,inverse.img[j] as usize);
    let binverse=b.inverse(d);let mut aa=a;let mut bb=b;
    // a'=s*a; b'=a^-1*s*a*b. Swap two outputs, not two domain labels.
    aa.img.swap(u,v);bb.img.swap(binverse.img[u] as usize,binverse.img[v] as usize);
    (aa,bb)
}

/// All one-transposition, index-sum- and partition-preserving cuts, in i,j order.
/// Unlike mutate_pair this is a finite population, with no trial cap or RNG.
pub(super) fn one_cut_pairs(d:usize,a:Perm,b:Perm)->Vec<(Perm,Perm)> {
    let target=pair_partition(d,a,b);let ca=cycle_labels(d,a);let cb=cycle_labels(d,b);
    let inverse=a.inverse(d);let mut out=Vec::new();
    for i in 0..d {for j in i+1..d {
        if target[i]!=target[j] || (ca[i]==ca[j])==(cb[inverse.img[i] as usize]==cb[inverse.img[j] as usize]) {continue;}
        let pair=cut(d,a,b,i,j);
        if pair_partition(d,pair.0,pair.1)==target {out.push(pair);}
    }}
    out
}

pub(super) fn mutate_pair(d:usize,mut a:Perm,mut b:Perm,steps:usize,rng:&mut SplitMix64)->(Option<(Perm,Perm)>,usize,usize) {
    let target=pair_partition(d,a,b);let mut previous=None;
    let (mut trials,mut made)=(0,0);
    for _ in 0..steps {
        let ca=cycle_labels(d,a);let cb=cycle_labels(d,b);let inverse=a.inverse(d);
        let mut choices=[0u16;MAX_D*(MAX_D-1)/2];let mut len=0;
        for i in 0..d {for j in i+1..d {
            if target[i]!=target[j] {continue;}
            let same_a=ca[i]==ca[j];let same_b=cb[inverse.img[i] as usize]==cb[inverse.img[j] as usize];
            if same_a!=same_b {choices[len]=((i<<6)|j) as u16;len+=1;}
        }}
        let mut chosen=None;
        // Bounded biased exploration, not enumeration. Do not undo the last cut.
        for _ in 0..len.min(16) {
            let k=rng.below(len);let packed=choices[k] as usize;len-=1;choices[k]=choices[len];trials+=1;
            let candidate=cut(d,a,b,packed>>6,packed&63);
            if previous==Some(candidate) || pair_partition(d,candidate.0,candidate.1)!=target {continue;}
            chosen=Some(candidate);break;
        }
        let Some((aa,bb))=chosen else {return (None,trials,made);};
        previous=Some((a,b));a=aa;b=bb;made+=1;
    }
    (Some((a,b)),trials,made)
}

/// Extra radius is even, so target parity depends only on the base profile.
/// Returns the proposal and measured trial/cut counts, including failed starts.
pub(super) fn mutate_tuple(d:usize,t:[Perm;4],ambient:Ambient,extra:usize,rng:&mut SplitMix64)->(Option<[Perm;4]>,usize,usize) {
    debug_assert!(extra%2==0);
    let mut profiles=[(0usize,0usize);6];let mut len=0;
    for (left,right) in [(2,0),(0,2),(1,1),(2,2),(1,0),(0,1)] {
        let odd=t.map(|p|!p.is_even(d));let parity=[odd[0]^(left%2==1),odd[1]^(left%2==1),
            odd[2]^(right%2==1),odd[3]^(right%2==1)];
        if (ambient==Ambient::Alternating && parity.iter().any(|&p|p))
            || (ambient==Ambient::Symmetric && !parity.iter().any(|&p|p)) {continue;}
        profiles[len]=(left,right);len+=1;
    }
    if len==0 {return (None,0,0);}
    let (left,right)=profiles[rng.below(len)];
    let (a,n,m)=mutate_pair(d,t[0],t[1],left+if left>0 {extra} else {0},rng);
    let Some((aa,bb))=a else {return (None,n,m);};
    let (b,nn,mm)=mutate_pair(d,t[2],t[3],right+if right>0 {extra} else {0},rng);
    let Some((cc,dd))=b else {return (None,n+nn,m+mm);};
    let u=[aa,bb,cc,dd];
    (if u==t {None} else {Some(u)},n+nn,m+mm)
}

#[cfg(test)]
mod tests {
use super::*;

#[test]
fn finite_one_cut_population_matches_literal_exhaustion() {
    for d in 3..=8 {let mut rng=SplitMix64::new(20261004+d as u64);
        for _ in 0..32 {let a=random_perm(d,Ambient::Symmetric,&mut rng);let b=random_perm(d,Ambient::Symmetric,&mut rng);
            let partition=pair_partition(d,a,b);let mut literal=Vec::new();
            for i in 0..d {for j in i+1..d {let mut s=Perm::identity(d);s.img.swap(i,j);
                let aa=s.compose(a,d);let bb=a.inverse(d).compose(s,d).compose(a,d).compose(b,d);
                if index(aa,d)+index(bb,d)==index(a,d)+index(b,d) && pair_partition(d,aa,bb)==partition {literal.push((aa,bb));}
            }}
            assert_eq!(one_cut_pairs(d,a,b),literal);
            assert_eq!(literal.iter().copied().collect::<HashSet<_>>().len(),literal.len());
        }
    }
}

#[test]
#[ignore = "paired one-worker exact cut-list benchmark; not a discovery multiplier"]
fn benchmark_finite_cut_listing() {
    let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt")).unwrap();
    let mut states=case.states.unwrap().into_iter().collect::<Vec<_>>();states.sort_unstable();
    let pairs=states.into_iter().flat_map(|k|{let t=unpack(k,18);[(t[0],t[1]),(t[2],t[3])]}).collect::<Vec<_>>();
    let literal=|a:Perm,b:Perm| {let target=pair_partition(18,a,b);let mut out=Vec::new();
        for i in 0..18 {for j in i+1..18 {let mut s=Perm::identity(18);s.img.swap(i,j);
            let aa=s.compose(a,18);let bb=a.inverse(18).compose(s,18).compose(a,18).compose(b,18);
            if index(aa,18)+index(bb,18)==index(a,18)+index(b,18) && pair_partition(18,aa,bb)==target {out.push((aa,bb));}
        }}out};
    for &(a,b) in &pairs {assert_eq!(literal(a,b),one_cut_pairs(18,a,b));}
    for block in 0..8 {for filtered in if block%2==0 {[false,true]} else {[true,false]} {
        let start=Instant::now();let mut total=0;
        for _ in 0..4 {for &(a,b) in &pairs {
            total+=std::hint::black_box(if filtered {one_cut_pairs(18,a,b)} else {literal(a,b)}).len();
        }}
        eprintln!("CUT_BENCH block={block} filtered={filtered} cuts={total} seconds={:.9}",start.elapsed().as_secs_f64());
    }}
}

#[test]
fn literal_cut_join_identities_and_all_control_invariants() {
    for d in 3..=7 {let mut rng=SplitMix64::new(20261004+d as u64);
        for _ in 0..32 {let a=random_perm(d,Ambient::Symmetric,&mut rng);let b=random_perm(d,Ambient::Symmetric,&mut rng);
            let ca=cycle_labels(d,a);let cb=cycle_labels(d,b);let inverse=a.inverse(d);
            for i in 0..d {for j in i+1..d {
                let (aa,bb)=cut(d,a,b,i,j);let mut s=Perm::identity(d);s.img.swap(i,j);
                assert_eq!(aa,s.compose(a,d));assert_eq!(bb,inverse.compose(s,d).compose(a,d).compose(b,d));
                assert_eq!(aa.compose(bb,d),a.compose(b,d));
                let opposite=(ca[i]==ca[j])!=(cb[inverse.img[i] as usize]==cb[inverse.img[j] as usize]);
                assert_eq!(index(aa,d)+index(bb,d)==index(a,d)+index(b,d),opposite);
                assert_ne!(aa.is_even(d),a.is_even(d));assert_ne!(bb.is_even(d),b.is_even(d));
            }}
        }
    }
    let mut changed=0;let mut alternating18=0;
    for path in fs::read_dir("evidence/certificates").unwrap().map(|r|r.unwrap().path()) {
        if !path.file_name().unwrap().to_string_lossy().ends_with("weak.txt") {continue;}
        let case=parse_certificate_case(&path).unwrap();let d=case.d;let mut rng=SplitMix64::new(8128);
        for key in case.states.unwrap() {let t=unpack(key,d);
            for ambient in [Ambient::Alternating,Ambient::Symmetric] {for _ in 0..8 {
                let (proposal,_,_)=mutate_tuple(d,t,ambient,0,&mut rng);
                if let Some(u)=proposal {changed+=1;
                    assert!(tuple_product_is_identity(d,&u));assert_ne!(u,t);
                    assert_eq!(u[0].compose(u[1],d),t[0].compose(t[1],d));
                    assert_eq!(u[2].compose(u[3],d),t[2].compose(t[3],d));
                    assert_eq!(pair_partition(d,u[0],u[1]),pair_partition(d,t[0],t[1]));
                    assert_eq!(pair_partition(d,u[2],u[3]),pair_partition(d,t[2],t[3]));
                    assert_eq!(genus_of(&u,d),genus_of(&t,d));assert_eq!(defect_of(&u,d),defect_of(&t,d));
                    assert_eq!(u.iter().all(|p|p.is_even(d)),ambient==Ambient::Alternating);
                    if d==18 && ambient==Ambient::Alternating {alternating18+=1;}
                }
            }}
        }
    }
    assert!(changed>0 && alternating18>0);
}
}
