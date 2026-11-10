//! Construct braid symmetries, then certify weak failures with the shared engine.
use super::*;
use std::sync::{Arc, OnceLock};
mod local_actions;
mod double_coset;
mod escape_repair;
mod joint_pure;

type ActionModel = (usize, Perm, Perm);

struct Frame { h: Perm, cycles: Vec<Vec<usize>>, pair: Option<(Perm,Perm)>, left_genus: i64 }
pub(super) struct Factory {
    frames: Vec<Frame>, active: Vec<usize>, discarded_frames: usize, central_frames: usize, twists: Vec<PairTwists>, pure: bool,
    period: usize, envelope: bool, envelope_cache: Mutex<HashMap<(Perm,Perm,Perm),bool>>,
    envelope_counts: [AtomicU64;3], glue: Option<(Vec<Arc<Vec<usize>>>,i64)>, excluded: HashSet<Passport>,
    completed_skips: AtomicU64, repair_steps: usize, repair_moves: Vec<Vec<Perm>>,
    pattern_counts: [AtomicU64;4], alias_frames: usize, pub counts: [AtomicU64;6],
    fiber_sources: Option<Vec<[Perm;4]>>, fiber_radius: Option<usize>, fiber_counts: [AtomicU64;4], screen_steps: usize,
    skip_known: bool, known: HashSet<Key>, known_skips: AtomicU64,
    local_audit: Option<local_actions::Audit>, boundary: Option<double_coset::Cache>, coordinated: Option<escape_repair::Repair>,
    joint: Option<joint_pure::Repair>,proposal_gates:[AtomicU64;6],
    witness_limit:u64,witness_exports:AtomicU64,witness:Option<Arc<positive_witness::Cache>>,
    witness_error:Mutex<Option<String>>
}

fn excluded_passports(args:&[String],d:usize,ambient:Ambient)->Result<HashSet<Passport>,String> {
    let mut excluded=HashSet::new();
    if let Some(path)=optional_arg(args,"--excluded-passports") {
        for line in fs::read_to_string(path).map_err(|e|e.to_string())?.lines() {
            if line.starts_with('#') || line.is_empty() {continue;}
            let fields=line.split('\t').collect::<Vec<_>>();
            if fields.len()!=4 {return Err("excluded-passports needs group, degree, passport, provenance".into());}
            if fields[0]==ambient.short() && fields[1]==d.to_string() {
                let mut p=parse_passport(fields[2],d,ambient)?;p.sort();excluded.insert(p);
            }
        }
    }
    Ok(excluded)
}

// A map of a transitive two-coloured orbit is determined by one root image.
// Cache all equivariant orbit maps once. Sampling then uses stack-only copies,
// including permutations of isomorphic blocks (not merely within-block twists).
struct PairTwists { orbits: Vec<Vec<usize>>, maps: Vec<Vec<(usize,Perm)>> }
impl PairTwists {
    fn new(d:usize,a:Perm,b:Perm)->Self {
        let orbits=point_orbits(d,&[a,b]);let mut orbit_of=[0usize;MAX_D];
        for (i,points) in orbits.iter().enumerate() {for &x in points {orbit_of[x]=i;}}
        let mut maps=Vec::new();
        for points in &orbits {
            let root=points[0];let mut choices=Vec::new();
            for image_root in 0..d {
                let mut images=[u8::MAX;MAX_D];let mut inverse=[u8::MAX;MAX_D];
                let mut queue=[0usize;MAX_D];queue[0]=root;
                images[root]=image_root as u8;inverse[image_root]=root as u8;
                let (mut head,mut tail,mut valid)=(0,1,true);
                while head<tail && valid {
                    let x=queue[head];head+=1;
                    for p in [a,b] {
                        let xx=p.img[x] as usize;let yy=p.img[images[x] as usize] as usize;
                        if (images[xx]!=u8::MAX && images[xx] as usize!=yy)
                            || (inverse[yy]!=u8::MAX && inverse[yy] as usize!=xx) {valid=false;break;}
                        if images[xx]==u8::MAX {images[xx]=yy as u8;inverse[yy]=xx as u8;queue[tail]=xx;tail+=1;}
                    }
                }
                if valid {let mut map=Perm::identity(d);for &x in points {map.img[x]=images[x];}
                    choices.push((orbit_of[image_root],map));}
            }
            assert!(!choices.is_empty());maps.push(choices);
        }
        Self {orbits,maps}
    }
    fn sample(&self,d:usize,rng:&mut SplitMix64)->Perm {
        let mut out=Perm::identity(d);let mut used=[false;MAX_D];
        for (points,maps) in self.orbits.iter().zip(&self.maps) {
            let count=maps.iter().filter(|(target,_)|!used[*target]).count();
            let (target,map)=maps.iter().filter(|(target,_)|!used[*target]).nth(rng.below(count)).unwrap();
            used[*target]=true;for &x in points {out.img[x]=map.img[x];}
        }
        out
    }
    fn enumerate(&self,d:usize,cap:usize)->(Vec<Perm>,bool) {
        fn visit(cache:&PairTwists,d:usize,i:usize,used:u64,p:Perm,cap:usize,out:&mut Vec<Perm>) {
            if out.len()>cap {return;}
            if i==cache.orbits.len() {out.push(p);return;}
            for &(target,map) in &cache.maps[i] {if used&(1u64<<target)!=0 {continue;}
                let mut next=p;for &x in &cache.orbits[i] {next.img[x]=map.img[x];}
                visit(cache,d,i+1,used|(1u64<<target),next,cap,out);
                if out.len()>cap {return;}
            }
        }
        let mut out=Vec::new();visit(self,d,0,0,Perm::identity(d),cap,&mut out);
        let complete=out.len()<=cap;out.truncate(cap);(out,complete)
    }
    // Internal automorphisms and exchanges of isomorphic coloured orbits
    // generate C(a,b). Reuse these for exact conjugacy reduction of h.
    fn generators(&self,d:usize)->Vec<Perm> {
        let mut out=HashSet::new();
        for (i,points) in self.orbits.iter().enumerate() {for &(j,map) in &self.maps[i] {
            let mut g=Perm::identity(d);
            for &x in points {g.img[x]=map.img[x];}
            if i!=j {for &x in points {g.img[map.img[x] as usize]=x as u8;}}
            if !g.is_identity(d) {out.insert(g);}
        }}
        let mut out=out.into_iter().collect::<Vec<_>>();out.sort_unstable_by_key(|p|p.img);out
    }
}

fn commutes(d:usize,a:Perm,b:Perm)->bool {a.compose(b,d)==b.compose(a,d)}
fn central_relabeling(d:usize,a:Perm,b:Perm,h:Perm)->bool {commutes(d,a,h)&&commutes(d,b,h)}

// C(h) has point orbits consisting of all cycles of a given length. If even
// <a,b,C(h)> is intransitive, no third entry c in C(h) can connect the tuple.
fn centralizer_can_connect(d:usize,a:Perm,b:Perm,cycles:&[Vec<usize>])->bool {
    let mut parent=std::array::from_fn::<_,MAX_D,_>(|i|i);
    fn root(parent:&mut [usize;MAX_D],mut x:usize)->usize {while parent[x]!=x {parent[x]=parent[parent[x]];x=parent[x];}x}
    fn join(parent:&mut [usize;MAX_D],a:usize,b:usize) {let x=root(parent,a);let y=root(parent,b);parent[x]=y;}
    for x in 0..d {for p in [a,b] {join(&mut parent,x,p.img[x] as usize);}}
    let mut representatives=[None;MAX_D+1];
    for cycle in cycles {let rep=*representatives[cycle.len()].get_or_insert(cycle[0]);for &x in cycle {join(&mut parent,rep,x);}}
    let r=root(&mut parent,0);(1..d).all(|x|root(&mut parent,x)==r)
}

// If a single entry already connects each <a,b> point orbit, adding c makes
// <a,b,c> transitive iff that entry together with c is transitive. Since d is
// determined by abc, this first pair can never start a generating weak failure.
fn redundant_pair_connectivity(d:usize,a:Perm,b:Perm)->bool {
    let orbits=point_orbits(d,&[a,b]).len();
    d-index(a,d)==orbits || d-index(b,d)==orbits
}

// Rotations and adjacent swaps of equal-length cycles generate C_S(h).
// For A_d use Schreier generators for the kernel of sign, not merely the
// individually even generators (products of odd generators matter).
fn centralizer_generators(d:usize,cycles:&[Vec<usize>],ambient:Ambient)->Vec<Perm> {
    let mut generators=Vec::new();
    for (i,cycle) in cycles.iter().enumerate() {
        if cycle.len()>1 {let mut rotation=Perm::identity(d);
            for j in 0..cycle.len() {rotation.img[cycle[j]]=cycle[(j+1)%cycle.len()] as u8;}
            generators.push(rotation);
        }
        if i>0 && cycles[i-1].len()==cycle.len() {let mut swap=Perm::identity(d);
            for (&a,&b) in cycles[i-1].iter().zip(cycle) {swap.img[a]=b as u8;swap.img[b]=a as u8;}
            generators.push(swap);
        }
    }
    if ambient==Ambient::Symmetric {return generators;}
    let Some(odd)=generators.iter().copied().find(|g|!g.is_even(d)) else {return generators;};
    let inverse=odd.inverse(d);let mut kernel=Vec::new();
    for t in [Perm::identity(d),odd] {for &g in &generators {
        let product=t.compose(g,d);let even=if product.is_even(d) {product} else {product.compose(inverse,d)};
        if !even.is_identity(d) {kernel.push(even);}
    }}
    kernel
}

fn envelope_order(d:usize,ambient:Ambient,a:Perm,b:Perm,cycles:&[Vec<usize>])->u128 {
    let mut generators=centralizer_generators(d,cycles,ambient);generators.extend([a,b]);
    subgroup_order(d,&generators)
}

fn power(d:usize,mut p:Perm,mut n:usize)->Perm {
    let mut out=Perm::identity(d);
    while n>0 {if n&1!=0 {out=out.compose(p,d);}p=p.compose(p,d);n>>=1;}out
}

// Full support on C(h): independently permute equal-length cycles and choose
// their cyclic offsets. No enumeration of S_d is needed.
fn centralizer(d:usize,cycles:&[Vec<usize>],rng:&mut SplitMix64)->Perm {
    let mut out=Perm::identity(d);let mut begin=0;
    while begin<cycles.len() {
        let l=cycles[begin].len();let mut end=begin+1;
        while end<cycles.len() && cycles[end].len()==l {end+=1;}
        let mut destinations=[0usize;MAX_D];
        for i in begin..end {destinations[i-begin]=i;}
        rng.shuffle(&mut destinations[..end-begin]);
        for i in begin..end {
            let target=&cycles[destinations[i-begin]];let shift=rng.below(l);
            for j in 0..l {out.img[cycles[i][j]]=target[(j+shift)%l] as u8;}
        }
        begin=end;
    }
    out
}

// Every square root joins two equal-length q-cycles, or takes the unique
// within-cycle square root of an odd cycle. Even cycles must all be paired.
fn square_root(d:usize,q:Perm,rng:&mut SplitMix64)->Option<Perm> {
    let mut cycles=permutation_cycles(d,q);let mut root=Perm::identity(d);let mut begin=0;
    while begin<cycles.len() {
        let l=cycles[begin].len();let mut end=begin+1;
        while end<cycles.len() && cycles[end].len()==l {end+=1;}
        if l%2==0 && (end-begin)%2!=0 {return None;}
        rng.shuffle(&mut cycles[begin..end]);let mut i=begin;
        while i<end {
            if i+1<end && (l%2==0 || rng.below(2)==0) {
                let shift=rng.below(l);
                for j in 0..l {
                    root.img[cycles[i][j]]=cycles[i+1][(j+shift)%l] as u8;
                    root.img[cycles[i+1][(j+shift)%l]]=cycles[i][(j+1)%l] as u8;
                }
                i+=2;
            } else {
                for j in 0..l {root.img[cycles[i][j]]=cycles[i][(j+(l+1)/2)%l] as u8;}
                i+=1;
            }
        }
        begin=end;
    }
    Some(root)
}

fn tuple_from_square(d:usize,h:Perm,x:Perm,c:Perm)->[Perm;4] {
    // a=xh, b=h^-1*x*h^2, ab=x^2*h^2. If [x^2,h]=[c,h]=1,
    // simultaneous conjugation by h is precisely sigma_1 on all four entries.
    let a=x.compose(h,d);let b=a.conjugate_by(h.inverse(d),d);
    [a,b,c,product_one_fourth(d,a,b,c)]
}

// Small triangle actions are built exactly on S4/H, rather than guessed by
// random point permutations. x^2=h^3=(xh)^4=1 in every action.
fn triangle_models()->Vec<(usize,Perm,Perm)> {
    let mut group=Vec::new();let mut p=Perm::identity(4);
    loop {group.push(p);if !next_permutation(&mut p.img[..4]) {break;}}
    let x=perm_of_type(&[2,1,1],4,None);
    let mut h=Perm::identity(4);h.img[..4].copy_from_slice(&[0,2,3,1]);
    let a=x.compose(h,4);let mut y=Perm::identity(4);y.img[..4].copy_from_slice(&[0,1,3,2]);
    let double=x.compose(y,4);let mut other=Perm::identity(4);other.img[..4].copy_from_slice(&[2,3,0,1]);
    let mut reflection=Perm::identity(4);reflection.img[..4].copy_from_slice(&[0,3,2,1]);
    let mut trans=Perm::identity(4);trans.img[..4].copy_from_slice(&[0,2,1,3]);
    let stabilizers=vec![vec![x,h],vec![h,double],vec![a,reflection],vec![h,trans],
        vec![a],vec![double,other],vec![x,y],vec![h],vec![x],vec![double],vec![]];
    let mut models=Vec::new();
    for generators in stabilizers {
        let mut subgroup=vec![Perm::identity(4)];let mut seen=HashSet::from([subgroup[0]]);let mut head=0;
        while head<subgroup.len() {let u=subgroup[head];head+=1;for &v in &generators {
            let w=u.compose(v,4);if seen.insert(w) {subgroup.push(w);}
        }}
        let key=|u:Perm|subgroup.iter().map(|&v|u.compose(v,4).img[..4].to_vec()).min().unwrap();
        let representatives=group.iter().map(|&u|key(u)).collect::<BTreeSet<_>>().into_iter().collect::<Vec<_>>();
        let n=representatives.len();let mut action=[Perm::identity(n);2];
        for (i,images) in representatives.iter().enumerate() {
            let mut u=Perm::identity(4);u.img[..4].copy_from_slice(images);
            for (j,g) in [x,h].iter().enumerate() {action[j].img[i]=representatives.binary_search(&key(g.compose(u,4))).unwrap() as u8;}
        }
        models.push((n,action[0],action[1]));
    }
    models
}

fn triangle_frames(d:usize,ambient:Ambient)->Vec<Frame> {
    fn extend(d:usize,ambient:Ambient,models:&[(usize,Perm,Perm)],begin:usize,
        used:usize,blocks:usize,x:Perm,h:Perm,out:&mut Vec<Frame>) {
        if used==d {
            let t=tuple_from_square(d,h,x,Perm::identity(d));
            if power_identity(d,t[0],2) || index(t[0],d)>d-3 || pair_transitive(d,t[0],t[1])
                || (ambient==Ambient::Alternating && !t[0].is_even(d)) {return;}
            let left_genus=defect_of(&t,d).1;
            if left_genus>d as i64-5 {return;}
            out.push(Frame {h,cycles:permutation_cycles(d,h),pair:Some((t[0],t[1])),left_genus});return;
        }
        // ponytail: at most eight point orbits; expand the representation
        // library if useful. This is not a complete degree census or time cap.
        if blocks==8 {return;}
        for i in begin..models.len() {
            let (n,xx,hh)=models[i];if used+n>d {continue;}
            let mut next_x=x;let mut next_h=h;
            for j in 0..n {next_x.img[used+j]=(used+xx.img[j] as usize) as u8;next_h.img[used+j]=(used+hh.img[j] as usize) as u8;}
            extend(d,ambient,models,i,used+n,blocks+1,next_x,next_h,out);
        }
    }
    let mut out=Vec::new();extend(d,ambient,&triangle_models(),0,0,0,Perm::identity(d),Perm::identity(d),&mut out);out
}

// Enumerate subgroups of the tiny abstract group, never subgroups/permutations
// of S_d. One representative per conjugacy class supplies every transitive
// G/H point action. u64 masks suffice for |A5|=60.
fn pure_models(model:&str)->Result<Vec<ActionModel>,String> {
    if model=="PSL27" {return Ok(fano_models());}
    let group=SmallPairGroup::new(model)?;
    Ok(group.actions(group.base))
}

// Work only inside a tiny abstract pair group. Enumerating generating pairs
// here is not enumerating permutations or subgroups of the ambient A_d/S_d.
struct SmallPairGroup {
    table: Vec<Vec<usize>>, subgroups: Vec<Vec<usize>>, base: [usize;2],
}

impl SmallPairGroup {
fn new(model:&str)->Result<Self,String> {
    let (n,even)=match model {"S3"=>(3,false),"A4"=>(4,true),"S4"=>(4,false),"A5"=>(5,true),
        _=>if let Some(m)=model.strip_prefix('D').and_then(|s|s.parse::<usize>().ok()) {
            if !(3..=15).contains(&m) {return Err("dihedral Dm requires 3<=m<=15".into());}(m,false)
        } else {return Err("pair-model must be S3, A4, S4, A5, PSL27 or D3..D15".into());}};
    let mut x=Perm::identity(n);let mut h=Perm::identity(n);
    if model.starts_with('D') {
        for i in 0..n {x.img[i]=((i+1)%n) as u8;h.img[i]=((n-i)%n) as u8;}
    } else {
        x.img[0]=1;x.img[1]=0;if even {x.img[2]=3;x.img[3]=2;}
        if n==3 {h.img[0]=1;h.img[1]=2;h.img[2]=0;}
        else if n==5 {h.img[1]=2;h.img[2]=4;h.img[4]=1;} else {h.img[1]=2;h.img[2]=3;h.img[3]=1;}
    }
    let mut group=vec![Perm::identity(n)];let mut seen=HashSet::from([group[0]]);let mut head=0;
    while head<group.len() {let p=group[head];head+=1;for generator in [x,h] {
        let q=p.compose(generator,n);if seen.insert(q) {group.push(q);}
    }}
    group.sort_unstable_by_key(|p|p.img);assert!(group.len()<=64);
    let lookup=group.iter().enumerate().map(|(i,&p)|(p,i)).collect::<HashMap<_,_>>();let size=group.len();
    let table=group.iter().map(|&a|group.iter().map(|&b|lookup[&a.compose(b,n)]).collect::<Vec<_>>()).collect::<Vec<_>>();
    let inverse=group.iter().map(|p|lookup[&p.inverse(n)]).collect::<Vec<_>>();
    let mut subgroups=vec![1u64];let mut seen=HashSet::from([1u64]);let mut head=0;
    while head<subgroups.len() {
        let subgroup=subgroups[head];head+=1;
        let members=(0..size).filter(|i|subgroup&(1u64<<i)!=0).collect::<Vec<_>>();
        for g in 0..size {if subgroup&(1u64<<g)!=0 {continue;}
            let mut mask=subgroup;let mut queue=members.clone();let mut pos=0;
            while pos<queue.len() {let u=queue[pos];pos+=1;
                for &v in members.iter().chain(std::iter::once(&g)) {let w=table[u][v];if mask&(1u64<<w)==0 {mask|=1u64<<w;queue.push(w);}}
            }
            if seen.insert(mask) {subgroups.push(mask);}
        }
    }
    let base=[lookup[&x],lookup[&h]];let mut representatives=Vec::new();
    for subgroup in subgroups {
        let members=(0..size).filter(|i|subgroup&(1u64<<i)!=0).collect::<Vec<_>>();
        if size/members.len()>MAX_D {continue;}
        if (0..size).any(|g|members.iter().fold(0u64,|mask,&u|mask|(1u64<<table[table[g][u]][inverse[g]]))<subgroup) {continue;}
        representatives.push(members);
    }
    Ok(Self {table,subgroups:representatives,base})
}

fn actions(&self,generators:[usize;2])->Vec<ActionModel> {
    let mut models=Vec::new();let size=self.table.len();
    for members in &self.subgroups {
        let key=|u:usize|members.iter().map(|&v|self.table[u][v]).min().unwrap();
        let reps=(0..size).map(key).collect::<BTreeSet<_>>().into_iter().collect::<Vec<_>>();let d=reps.len();
        let mut actions=[Perm::identity(d);2];
        for (i,&u) in reps.iter().enumerate() {for j in 0..2 {
            actions[j].img[i]=reps.binary_search(&key(self.table[generators[j]][u])).unwrap() as u8;
        }}
        models.push((d,actions[0],actions[1]));
    }
    models.sort_by_key(|m|m.0);models
}

// BFS rooted at 1, with ordered colours a,b, canonically labels the Cayley
// graph. Equality is precisely an isomorphism of marked abstract groups.
// Thus deduplication includes outer automorphisms, not just inner conjugacy.
fn marking_key(&self,generators:[usize;2])->Option<Vec<u8>> {
    let size=self.table.len();let mut labels=[usize::MAX;64];labels[0]=0;
    let mut queue=vec![0usize];let mut head=0;
    while head<queue.len() {let u=queue[head];head+=1;
        for g in generators {let v=self.table[u][g];if labels[v]==usize::MAX {labels[v]=queue.len();queue.push(v);}}
    }
    if queue.len()!=size {return None;}
    Some(queue.into_iter().flat_map(|u|generators.map(|g|labels[self.table[u][g]] as u8)).collect())
}

fn markings(&self,novel:bool)->Vec<[usize;2]> {
    let [x,y]=self.base;let xy=self.table[x][y];let mut seen=HashSet::new();
    if novel {for pair in [[x,y],[xy,y],[x,xy]] {seen.insert(self.marking_key(pair).unwrap());}}
    let mut result=Vec::new();
    for a in 1..self.table.len() {for b in 1..self.table.len() {
        if let Some(key)=self.marking_key([a,b]) {if seen.insert(key) {result.push([a,b]);}}
    }}
    result
}
}

fn marked_models(model:&str,novel:bool)->Result<Vec<Vec<ActionModel>>,String> {
    type Cache = HashMap<(String,bool),Vec<Vec<ActionModel>>>;
    static CACHE:OnceLock<Mutex<Cache>>=OnceLock::new();
    let cache=CACHE.get_or_init(||Mutex::new(HashMap::new()));let key=(model.to_string(),novel);
    if let Some(value)=cache.lock().unwrap().get(&key) {return Ok(value.clone());}
    let group=SmallPairGroup::new(model)?;
    let models=group.markings(novel).into_iter().map(|pair|group.actions(pair)).collect::<Vec<_>>();
    cache.lock().unwrap().insert(key,models.clone());Ok(models)
}

// Three exact (2,3,7) actions: nonzero vectors, dual covectors over F2,
// and P1(F7). The presentation also has [x,y]^4=1. No 168-subgroup census
// is asserted: this explicit point/line/projective library is the population.
fn fano_models()->Vec<(usize,Perm,Perm)> {
    let mut x=Perm::identity(7);let mut y=Perm::identity(7);
    for v in 1..=7 {x.img[v-1]=(v^((v&2)>>1)) as u8-1;y.img[v-1]=(((v<<1)&7)|(v>>2)) as u8-1;}
    let dot=|a:usize,b:usize|(a&b).count_ones()%2;
    let dual=|p:Perm| {let mut out=Perm::identity(7);
        for w in 1..=7 {out.img[w-1]=((1..=7).find(|&z|(1..=7).all(|v|dot(z,p.img[v-1] as usize+1)==dot(w,v))).unwrap()-1) as u8;}out};
    let mut xx=Perm::identity(8);let mut yy=Perm::identity(8);
    for v in 0..8 {for (p,offset) in [(&mut xx,0),(&mut yy,1)] {
        p.img[v]=if v==7 {0} else {let denominator=(v+offset)%7;
            if denominator==0 {7} else {let inverse=(1..7).find(|&u|u*denominator%7==1).unwrap();((7-inverse)%7) as u8}};
    }}
    vec![(1,Perm::identity(1),Perm::identity(1)),(7,x,y),(7,dual(x),dual(y)),(8,xx,yy)]
}

fn pure_frames(d:usize,ambient:Ambient,model:&str)->Result<Vec<Frame>,String> {
    pure_frames_filtered(d,ambient,model,false)
}

fn pure_frames_filtered(d:usize,ambient:Ambient,model:&str,post_cut:bool)->Result<Vec<Frame>,String> {
    let limit=if model=="S3" && d<=11 {d} else {8};
    Ok(action_frames(d,ambient,&pure_models(model)?,limit,false,post_cut))
}

fn marked_frames(args:&[String],d:usize,ambient:Ambient,model:&str)->Result<Vec<Frame>,String> {
    let mode=get_arg(args,"--marked-pairs",Some("novel"))?;
    if !matches!(mode.as_str(),"novel"|"all") {return Err("marked-pairs must be novel or all".into());}
    let limit:usize=parse_num(get_arg(args,"--point-orbit-max",Some(if d<=11 {"34"} else {"8"}))?,"--point-orbit-max")?;
    if limit==0 || limit>MAX_D {return Err("point-orbit-max must be 1..34".into());}
    let mut frames=Vec::new();let mut seen=HashSet::new();
    for models in marked_models(model,mode=="novel")? {
        for frame in action_frames(d,ambient,&models,limit,true,has_flag(args,"--post-cut-base-frames")) {
            if seen.insert(frame.pair.unwrap()) {frames.push(frame);}
        }
    }
    Ok(frames)
}

fn action_frames(d:usize,ambient:Ambient,models:&[ActionModel],limit:usize,marked:bool,post_cut:bool)->Vec<Frame> {
    fn extend(d:usize,ambient:Ambient,models:&[(usize,Perm,Perm)],block_limit:usize,begin:usize,used:usize,blocks:usize,
        x:Perm,y:Perm,marked:bool,post_cut:bool,out:&mut Vec<Frame>) {
        if used==d {
            let xy=x.compose(y,d);
            let pairs=[(x,y),(xy,y),(x,xy)];
            for &(a,b) in &pairs[..if marked {1} else {3}] {
                if pair_transitive(d,a,b) || (!post_cut && (index(a,d)==0 || index(b,d)==0 || index(a,d)>d-3 || index(b,d)>d-3
                    || (ambient==Ambient::Alternating && (!a.is_even(d)||!b.is_even(d))))) {continue;}
                let h=a.compose(b,d);let cycles=permutation_cycles(d,h);
                // Pure braid: sigma_1^2 acts by conjugation with p=ab on the
                // first pair. For any k in C(a,b), h=pk acts identically there;
                // c in C(h) makes the entire tuple fixed modulo conjugacy.
                let t=[a,b,Perm::identity(d),h.inverse(d)];let left_genus=defect_of(&t,d).1;
                if left_genus<=d as i64-5 {out.push(Frame {h,cycles,pair:Some((a,b)),left_genus});}
            }
            return;
        }
        // ponytail: eight point orbits, three generator-pair choices; expand
        // only if these constructions yield useful leads. Not a degree census.
        if blocks==block_limit {return;}
        for i in begin..models.len() {let (n,xx,yy)=models[i];if used+n>d {continue;}
            let mut a=x;let mut b=y;
            for j in 0..n {a.img[used+j]=(used+xx.img[j] as usize) as u8;b.img[used+j]=(used+yy.img[j] as usize) as u8;}
            extend(d,ambient,models,block_limit,i,used+n,blocks+1,a,b,marked,post_cut,out);
        }
    }
    let mut out=Vec::new();
    extend(d,ambient,models,limit,0,0,0,Perm::identity(d),Perm::identity(d),marked,post_cut,&mut out);out
}

fn dessin_frames(d:usize,ambient:Ambient)->Result<Vec<Frame>,String> {
    pair_frames(&[],d,ambient,0)
}

fn pair_frames(args:&[String],d:usize,ambient:Ambient,max_genus:i64)->Result<Vec<Frame>,String> {
    let post_cut=has_flag(args,"--post-cut-base-frames");
    let mut frames=triangle_frames(d,ambient);let mut seen=HashSet::new();
    for model in ["S3","A4","S4","A5","PSL27"] {frames.extend(pure_frames_filtered(d,ambient,model,post_cut)?);}
    if has_flag(args,"--expand-pair-actions") {
        for model in ["A4","S4","A5"] {frames.extend(marked_frames(args,d,ambient,model)?);}
        let mut all=args.to_vec();
        if let Some(i)=all.iter().position(|s|s=="--marked-pairs") {all[i+1]="all".into();}
        else {all.extend(["--marked-pairs".into(),"all".into()]);}
        for model in ["D5","D7","D9","D11","D13"] {frames.extend(marked_frames(&all,d,ambient,model)?);}
    }
    if let Some(raw)=optional_arg(args,"--known-certificates") {
        let mut paths=Vec::new();for path in raw.split(',') {collect_certificate_paths(Path::new(path),&mut paths)?;}
        for path in paths {let case=parse_certificate_case(&path)?;
            if case.d!=d || (!post_cut && ambient==Ambient::Alternating && case.ambient!=ambient) {continue;}
            verify(&path)?;
            for key in case.states.ok_or("gluing source requires closed states")? {
                let t=unpack(key,d);let h=t[0].compose(t[1],d);
                frames.push(Frame {h,cycles:permutation_cycles(d,h),pair:Some((t[0],t[1])),left_genus:defect_of(&t,d).1});
            }
        }
    }
    frames.retain(|f|f.left_genus<=max_genus && seen.insert(f.pair.unwrap()));
    // Normalize every boundary to p=ab; triangle frames originally stored a
    // one-braid relabeling, which is not the gluing boundary permutation.
    for f in &mut frames {let (a,b)=f.pair.unwrap();f.h=a.compose(b,d);f.cycles=permutation_cycles(d,f.h);}
    Ok(frames)
}

fn glue_tuple(d:usize,left:&Frame,right:&Frame,rng:&mut SplitMix64)->[Perm;4] {
    let (a,b)=left.pair.unwrap();let (_,y)=right.pair.unwrap();
    let mut q=Perm::identity(d);let mut begin=0;
    while begin<right.cycles.len() {let l=right.cycles[begin].len();let mut end=begin+1;
        while end<right.cycles.len() && right.cycles[end].len()==l {end+=1;}
        let mut destinations=[0usize;MAX_D];for i in begin..end {destinations[i-begin]=i;}
        rng.shuffle(&mut destinations[..end-begin]);
        for i in begin..end {let shift=rng.below(l);let target=&left.cycles[destinations[i-begin]];
            for j in 0..l {q.img[right.cycles[i][j]]=target[(j+shift)%l] as u8;}
        }begin=end;
    }
    let c=y.inverse(d).conjugate_by(q,d);[a,b,c,product_one_fourth(d,a,b,c)]
}

// These are explicit positive witnesses, never heuristic impossibility tests.
// Include the shared escape sigma2, sigma1^-1, sigma2 of the thirteen C107
// leads, both signs, and the intermediate middle-move states. No allocation,
// group-order computation or canonicalization is needed for this score.
fn bridge_score(d:usize,t:[Perm;4])->usize {
    let mut score=if generators_transitive(d,&t[..3]) {0} else {32};
    for i in 0..4 {for j in i+1..4 {score+=usize::from(pair_transitive(d,t[i],t[j]));}}
    for first in [1,4] {
        let u=hurwitz_move_dynamic(d,t,first);
        score+=usize::from(pair_transitive(d,u[0],u[1]));
        for inner in [0,3] {let v=hurwitz_move_dynamic(d,u,inner);
            for last in [1,4] {let w=hurwitz_move_dynamic(d,v,last);
                score+=usize::from(pair_transitive(d,w[0],w[1]));
            }
        }
    }
    score
}

fn repair_gluing(d:usize,mut t:[Perm;4],moves:&[Perm],steps:usize,rng:&mut SplitMix64)->(Option<[Perm;4]>,bool,usize) {
    let mut score=bridge_score(d,t);if score==0 {return (Some(t),false,0);}
    if moves.is_empty() {return (None,false,0);}
    for step in 0..steps {
        let mut s=moves[rng.below(moves.len())];
        if rng.below(2)==0 {s=s.inverse(d);}
        let mut u=t;u[2]=t[2].conjugate_by(s,d);u[3]=t[3].conjugate_by(s,d);
        let next=bridge_score(d,u);
        if next==0 {return (Some(u),true,step+1);}
        // Plateau moves explore offsets/block exchanges; occasional uphill
        // moves escape local minima. This is biased search, not a census.
        if next<=score || rng.below(16)==0 {t=u;score=next;}
    }
    (None,false,steps)
}

impl Factory {
    pub fn new(args:&[String],d:usize,ambient:Ambient,family:&str)->Result<Self,String> {
        let post_cut=has_flag(args,"--post-cut-base-frames");
        if post_cut && !matches!(family,"cut-join-glue"|"coordinated-glue") {return Err("post-cut-base-frames requires cut-join-glue or coordinated-glue".into());}
        let period:usize=parse_num(get_arg(args,"--pure-power",Some("1"))?,"--pure-power")?;
        if !(1..=7).contains(&period) {return Err("pure-power must be 1..7 (period divides 2r)".into());}
        let mut frames=Vec::new();let mut fiber_sources=None;let mut local_audit=None;
        if matches!(family,"pair-fiber"|"coordinated-fiber") {
            let path=PathBuf::from(get_arg(args,"--certificate",None)?);verify(&path)?;
            let case=parse_certificate_case(&path)?;
            if case.d!=d {return Err("fiber source degree mismatch".into());}
            // Source and target groups may differ: S18 -> A18, A25 -> S25.
            fiber_sources=Some(case.states.ok_or("fiber source needs closed states")?.into_iter().map(|key|unpack(key,d)).collect());
        } else if matches!(family,"independent-symmetry"|"independent-glue"|"joint-pure") {
            if !has_flag(args,"--joint-source-only") {
                let (library,audit)=local_actions::frames(args,d,ambient,post_cut)?;frames=library;local_audit=Some(audit);
            } else if family!="joint-pure" || optional_arg(args,"--certificate").is_none() {
                return Err("joint-source-only requires joint-pure and a certificate".into());
            }
            if family=="joint-pure" {
                if let Some(path)=optional_arg(args,"--certificate") {
                    verify(Path::new(&path))?;
                    let case=parse_certificate_case(Path::new(&path))?;
                    if case.d!=d {return Err("joint source degree mismatch".into());}
                    let mut seen=HashSet::new();
                    for key in case.states.ok_or("joint source needs closed states")? {
                        let t=unpack(key,d);
                        if ambient==Ambient::Alternating && (!t[0].is_even(d)||!t[1].is_even(d)) {continue;}
                        let moved=hurwitz_move_dynamic(d,hurwitz_move_dynamic(d,t,0),0);
                        if canonical(&moved,d)!=key {continue;}
                        let h=tuple_conjugator(d,t,moved)?;
                        if seen.insert(colored_representation_key(d,[t[0],t[1],h])) {
                            frames.push(Frame {h,cycles:permutation_cycles(d,h),pair:Some((t[0],t[1])),left_genus:defect_of(&t,d).1});
                        }
                    }
                }
            }
        } else if matches!(family,"pair-glue"|"resistant-glue"|"cut-join-glue"|"double-coset-glue"|"coordinated-glue") {
            let max_genus:i64=parse_num(get_arg(args,"--pair-genus-max",Some("2"))?,"--pair-genus-max")?;
            if !(0..=d as i64-5).contains(&max_genus) {return Err("pair-genus-max must be 0..d-5".into());}
            frames=pair_frames(args,d,ambient,max_genus)?;
            if matches!(family,"double-coset-glue"|"coordinated-glue") {
                let (extra,audit)=local_actions::frames(args,d,ambient,post_cut)?;frames.extend(extra.into_iter().filter(|f|f.left_genus<=max_genus));local_audit=Some(audit);
                let mut seen=HashSet::new();frames.retain(|f|seen.insert(colored_representation_key(d,[f.pair.unwrap().0,f.pair.unwrap().1,Perm::identity(d)])));
            }
        } else if family=="dessin-glue" {frames=dessin_frames(d,ambient)?;}
        else if family=="marked-symmetry" {frames=marked_frames(args,d,ambient,&get_arg(args,"--pair-model",Some("A5"))?)?;}
        else if family=="pure-symmetry" {frames=pure_frames(d,ambient,&get_arg(args,"--pair-model",Some("S4"))?)?;}
        else if family=="triangle-symmetry" {frames=triangle_frames(d,ambient);} else if matches!(family,"symmetry-centralizer"|"pure-centralizer") {
            let path=get_arg(args,"--certificate",None)?;verify(Path::new(&path))?;
            let case=parse_certificate_case(Path::new(&path))?;
            if case.d!=d {return Err("symmetry source degree mismatch".into());}
            let mut seen=HashSet::new();
            for key in case.states.ok_or("source requires closed states")? {
                let t=unpack(key,d);let moved=hurwitz_move_dynamic(d,t,0);
                if family=="symmetry-centralizer" && canonical(&moved,d)!=key {continue;}
                let h=if family=="pure-centralizer" {t[0].compose(t[1],d)} else {tuple_conjugator(d,t,moved)?};
                if seen.insert((t[0],t[1],h)) {
                    frames.push(Frame {h,cycles:permutation_cycles(d,h),pair:Some((t[0],t[1])),left_genus:defect_of(&t,d).1});
                }
            }
        } else {
            let mut profiles=BTreeSet::new();
            if let Some(raw)=optional_arg(args,"--symmetry-type") {
                let mut lengths=raw.split('.').map(|s|parse_num::<usize>(s.into(),"cycle length")).collect::<Result<Vec<_>,_>>()?;
                if lengths.iter().any(|&l|l==0 || l>d) || lengths.iter().sum::<usize>()>d {return Err("invalid symmetry type".into());}
                lengths.extend(vec![1;d-lengths.iter().sum::<usize>()]);lengths.sort_unstable();profiles.insert(lengths);
            } else {
                for l in 2..=6 {for m in l..=6 {for a in 1..=d/l {for b in 0..=(d-a*l)/m {
                    let fixed=d-a*l-b*m;
                    if fixed>4 || a+b+fixed<3 {continue;}
                    let mut p=vec![1;fixed];p.extend(vec![l;a]);p.extend(vec![m;b]);p.sort_unstable();profiles.insert(p);
                }}}}
            }
            for p in profiles {
                let h=perm_of_type(&p,d,None);
                if !h.is_identity(d) {frames.push(Frame {h,cycles:permutation_cycles(d,h),pair:None,left_genus:0});}
            }
        }
        let original_frames=frames.len();
        if !post_cut {frames.retain(|f|f.pair.is_none_or(|(a,b)|!redundant_pair_connectivity(d,a,b)));}
        let discarded_frames=original_frames-frames.len();
        if original_frames==0 && !matches!(family,"marked-symmetry"|"resistant-glue"|"cut-join-glue"|"pair-fiber"|"independent-symmetry"|"independent-glue"|"coordinated-fiber"|"coordinated-glue"|"double-coset-glue"|"joint-pure") {return Err("no nontrivial symmetry frames".into());}
        let pure=matches!(family,"pure-symmetry"|"pure-centralizer"|"marked-symmetry"|"independent-symmetry"|"joint-pure");
        if family=="joint-pure" && period!=1 {return Err("joint-pure imposes squares, so pure-power must be 1".into());}
        let mut alias_frames=0;
        if pure && has_flag(args,"--novel-power") {
            let before=frames.len();
            frames.retain(|f|{let (a,b)=f.pair.unwrap();
                // p^(r-s) in C(a,b) gives the same coset p^r C(a,b).
                // Inverting h preserves C(h), so +/-s are both aliases.
                !(1..=3).any(|s|central_relabeling(d,a,b,power(d,f.h,period.abs_diff(s)))
                    || central_relabeling(d,a,b,power(d,f.h,period+s)))
            });alias_frames=before-frames.len();
        }
        if pure {for frame in &mut frames {frame.h=power(d,frame.h,period);frame.cycles=permutation_cycles(d,frame.h);}}
        // k commutes with a,b. If p^r does too, every h=p^r k commutes
        // with the whole completion. A_d/S_d (d>=4) has trivial S_d centralizer;
        // h!=1 cannot generate the ambient group, and h=1 is not this lane.
        let before=frames.len();
        if pure && d>=4 {frames.retain(|f|{let (a,b)=f.pair.unwrap();!central_relabeling(d,a,b,f.h)});}
        let central_frames=before-frames.len();
        let twists=if pure || matches!(family,"triangle-symmetry"|"symmetry-centralizer") {frames.iter().map(|f|{let (a,b)=f.pair.unwrap();PairTwists::new(d,a,b)}).collect()} else {Vec::new()};
        let glue=if matches!(family,"dessin-glue"|"pair-glue"|"resistant-glue"|"cut-join-glue"|"independent-glue"|"double-coset-glue"|"coordinated-glue") {
            let genus:i64=parse_num(get_arg(args,"--glue-genus",Some("1"))?,"--glue-genus")?;
            if !(0..=d as i64-5).contains(&genus) {return Err("glue-genus outside nonidentity weak-search range".into());}
            let signatures=frames.iter().map(|f|f.cycles.iter().map(Vec::len).collect::<Vec<_>>()).collect::<Vec<_>>();
            let mut boundaries=BTreeMap::<Vec<usize>,Vec<usize>>::new();
            for (i,s) in signatures.iter().enumerate() {boundaries.entry(s.clone()).or_default().push(i);}
            let orbit_counts=frames.iter().map(|f|{let (a,b)=f.pair.unwrap();point_orbits(d,&[a,b]).len()}).collect::<Vec<_>>();
            let odd=frames.iter().map(|f|{let (a,b)=f.pair.unwrap();!a.is_even(d)||!b.is_even(d)}).collect::<Vec<_>>();
            // Eligibility depends on decorations, not the individual frame.
            // Share each matching list rather than allocate a quadratic number
            // of identical right-frame indices in the expanded library.
            let mut cache=HashMap::new();let mut lists=Vec::with_capacity(frames.len());
            for (i,f) in frames.iter().enumerate() {
                let key=(signatures[i].clone(),f.left_genus,orbit_counts[i],odd[i]);
                let matches=cache.entry(key).or_insert_with(||Arc::new(boundaries[&signatures[i]].iter().copied()
                    .filter(|&j|{let defect=f.cycles.len() as i64-orbit_counts[i] as i64-orbit_counts[j] as i64+1;
                        defect>=0 && defect+f.left_genus+frames[j].left_genus==genus
                            && (post_cut || ambient!=Ambient::Symmetric || odd[i] || odd[j])
                    }).collect::<Vec<_>>()));
                lists.push(Arc::clone(matches));
            }
            Some((lists,genus))
        } else {None};
        let repair_steps:usize=parse_num(get_arg(args,"--glue-repair-steps",Some(if family=="resistant-glue" {"16"} else {"0"}))?,"--glue-repair-steps")?;
        if repair_steps>128 || (repair_steps>0 && glue.is_none()) {return Err("glue-repair-steps needs gluing and must be 0..128".into());}
        let repair_moves=if repair_steps>0 {frames.iter().map(|f|centralizer_generators(d,&f.cycles,Ambient::Symmetric)).collect()} else {Vec::new()};
        let fiber_radius=if matches!(family,"pair-fiber"|"cut-join-glue") {
            let radius:usize=parse_num(get_arg(args,"--fiber-radius",Some("0"))?,"--fiber-radius")?;
            if radius>8 || radius%2!=0 {return Err("fiber-radius must be 0,2,4,6,8 extra cuts per changed pair".into());}
            Some(radius)
        } else {None};
        let screen_steps:usize=parse_num(get_arg(args,"--screen-steps",Some("16"))?,"--screen-steps")?;
        if screen_steps>512 {return Err("screen-steps must be 0..512".into());}
        let active=(0..frames.len()).filter(|&i|glue.as_ref().is_none_or(|(matches,_)|!matches[i].is_empty())).collect();
        let boundary=if matches!(family,"independent-glue"|"double-coset-glue"|"coordinated-glue") {
            let retire=has_flag(args,"--retire-boundary-plans");
            if retire && (family=="coordinated-glue" || repair_steps>0 || fiber_radius.is_some()) {
                return Err("boundary retirement is only valid for unchanged, unrepaired gluing".into());
            }
            Some(double_coset::Cache::new(parse_num(get_arg(args,"--boundary-centralizer-cap",Some("4096"))?,"--boundary-centralizer-cap")?,retire)?)
        } else {None};
        let coordinated=if matches!(family,"coordinated-glue"|"coordinated-fiber") {Some(escape_repair::Repair::new(args)?)} else {None};
        let joint=if family=="joint-pure" {Some(joint_pure::Repair::new(args)?)} else {None};
        if joint.is_none() && ["--joint-admission-aware","--joint-pure-basis","--joint-hidden-fiber"].iter().any(|flag|has_flag(args,flag)) {
            return Err("joint objectives require family joint-pure".into());
        }
        let witness_limit:u64=parse_num(get_arg(args,"--construction-witness-limit",Some("0"))?,"--construction-witness-limit")?;
        if witness_limit>1024 {return Err("construction-witness-limit must be 0..1024".into());}
        if has_flag(args,"--allow-partition-change") && coordinated.is_none() {return Err("allow-partition-change requires coordinated repair".into());}
        Ok(Self {frames,active,discarded_frames,central_frames,twists,pure,period,envelope:!has_flag(args,"--no-envelope"),
            envelope_cache:Mutex::new(HashMap::new()),envelope_counts:std::array::from_fn(|_|AtomicU64::new(0)),
            glue,excluded:excluded_passports(args,d,ambient)?,completed_skips:AtomicU64::new(0),
            repair_steps,repair_moves,pattern_counts:std::array::from_fn(|_|AtomicU64::new(0)),alias_frames,
            fiber_sources,fiber_radius,fiber_counts:std::array::from_fn(|_|AtomicU64::new(0)),screen_steps,
            skip_known:has_flag(args,"--skip-known-starts"),known:HashSet::new(),known_skips:AtomicU64::new(0),
            counts:std::array::from_fn(|_|AtomicU64::new(0)),local_audit,boundary,coordinated,joint,
            proposal_gates:std::array::from_fn(|_|AtomicU64::new(0)),witness_limit,
            witness_exports:AtomicU64::new(0),witness:None,witness_error:Mutex::new(None)})
    }

    pub fn is_empty(&self)->bool {self.fiber_sources.as_ref().map_or(self.active.is_empty(),Vec::is_empty)}

    pub fn set_known(&mut self,known:&HashSet<Key>) {if self.skip_known {self.known=known.clone();}}

    pub fn set_positive(&mut self,cache:Option<Arc<positive_witness::Cache>>)->Result<(),String> {
        if self.witness_limit>0 && cache.is_none() {return Err("construction witness export requires positive caching".into());}
        self.witness=cache;Ok(())
    }

    fn save_escape(&self,d:usize,t:[Perm;4],word:&[u8]) {
        let Some(cache)=&self.witness else {return;};
        if self.witness_exports.try_update(Ordering::Relaxed,Ordering::Relaxed,
            |n|(n<self.witness_limit).then_some(n+1)).is_err() {return;}
        if let Err(e)=cache.insert(d,t,word) {*self.witness_error.lock().unwrap()=Some(e.to_string());}
    }

    fn fiber_proposal(&self,d:usize,t:[Perm;4],ambient:Ambient,rng:&mut SplitMix64)->Option<[Perm;4]> {
        let Some(radius)=self.fiber_radius else {return Some(t);};
        self.fiber_counts[0].fetch_add(1,Ordering::Relaxed);
        let (proposal,trials,cuts)=pair_fiber::mutate_tuple(d,t,ambient,radius,rng);
        self.fiber_counts[1].fetch_add(trials as u64,Ordering::Relaxed);
        self.fiber_counts[2].fetch_add(cuts as u64,Ordering::Relaxed);
        self.fiber_counts[3].fetch_add(u64::from(proposal.is_some()),Ordering::Relaxed);proposal
    }

    fn indices_admit(d:usize,ambient:Ambient,t:[Perm;4],left_genus:i64)->bool {
        let indices=t.map(|p|index(p,d));
        !(indices.iter().any(|&i|i==0 || i>d-3) || indices.iter().sum::<usize>()<2*d-2
            || 1-d as i64+indices.iter().sum::<usize>() as i64/2<left_genus
            || (ambient==Ambient::Alternating && indices.iter().any(|i|i%2!=0))
            || (ambient==Ambient::Symmetric && indices.iter().all(|i|i%2==0)))
    }

    fn repair_grade(&self,d:usize,ambient:Ambient,t:[Perm;4],left_genus:i64)->u8 {
        if !Self::indices_admit(d,ambient,t,left_genus) {return 0;}
        if !generators_transitive(d,&t[..3]) || has_transitive_pair(d,&t) {return 1;}
        2
    }

    fn admit(&self,d:usize,ambient:Ambient,t:[Perm;4],left_genus:i64,check_pair:bool,rng:&mut SplitMix64)->Option<[Perm;4]> {
        self.counts[1].fetch_add(1,Ordering::Relaxed);
        if !Self::indices_admit(d,ambient,t,left_genus) {return None;}
        self.counts[2].fetch_add(1,Ordering::Relaxed);
        if !self.excluded.is_empty() {let mut p=t.map(|g|cycle_type(g,d));p.sort();
            if self.excluded.contains(&p) {self.completed_skips.fetch_add(1,Ordering::Relaxed);return None;}
        }
        if check_pair && redundant_pair_connectivity(d,t[0],t[1]) {return None;}
        if !generators_transitive(d,&t[..3]) || has_transitive_pair(d,&t) {return None;}
        self.counts[3].fetch_add(1,Ordering::Relaxed);
        if bridge_score(d,t)>0 {
            self.pattern_counts[3].fetch_add(1,Ordering::Relaxed);
            if self.witness_exports.load(Ordering::Relaxed)<self.witness_limit {
                'found: for first in [1,4] {
                    let u=hurwitz_move_dynamic(d,t,first);
                    if pair_transitive(d,u[0],u[1]) {self.save_escape(d,t,&[first as u8]);break;}
                    for inner in [0,3] {for last in [1,4] {
                        let w=hurwitz_move_dynamic(d,hurwitz_move_dynamic(d,u,inner),last);
                        if pair_transitive(d,w[0],w[1]) {self.save_escape(d,t,&[first as u8,inner as u8,last as u8]);break 'found;}
                    }}
                }
            }return None;
        }
        if !self.known.is_empty() && self.known.contains(&canonical(&t,d)) {
            // Already certified components are not new discoveries. This is
            // an exact membership skip, not a positive transitivity witness.
            self.known_skips.fetch_add(1,Ordering::Relaxed);return None;
        }
        if self.witness_limit>0 && self.witness.as_ref().is_some_and(|cache|cache.contains(canonical(&t,d))) {
            return None; // Imported/recorded literal positive proof, never a heuristic miss.
        }
        // A weak witness is valid even before full generation. New policies
        // can spend more cheap walk steps here to avoid expensive order checks.
        let trace=self.witness.is_some() && self.witness_exports.load(Ordering::Relaxed)<self.witness_limit;
        let mut word=Vec::new();
        if probe_with_word(d,t,rng,self.screen_steps,true,if trace {Some(&mut word)} else {None}).0 {
            if trace {self.save_escape(d,t,&word);}return None;
        }
        self.counts[4].fetch_add(1,Ordering::Relaxed);
        if !tuple_generates_dynamic(d,ambient,&t) {return None;}
        self.counts[5].fetch_add(1,Ordering::Relaxed);Some(t)
    }

    fn envelope_admits(&self,d:usize,ambient:Ambient,a:Perm,b:Perm,h:Perm,cycles:&[Vec<usize>])->bool {
        if !self.envelope {return true;}
        let key=(a,b,h);
        if let Some(&full)=self.envelope_cache.lock().unwrap().get(&key) {
            self.envelope_counts[1].fetch_add(1,Ordering::Relaxed);return full;
        }
        let full=envelope_order(d,ambient,a,b,cycles)==ambient.target_order(d);
        self.envelope_counts[0].fetch_add(1,Ordering::Relaxed);
        if !full {self.envelope_counts[2].fetch_add(1,Ordering::Relaxed);}
        let mut cache=self.envelope_cache.lock().unwrap();
        // ponytail: bounded per-process memoization; clearing costs recomputation
        // only, never loses proposals. Add an LRU only if measurements warrant it.
        if cache.len()==4096 {cache.clear();}cache.insert(key,full);full
    }

    pub fn generate(&self,d:usize,ambient:Ambient,rng:&mut SplitMix64)->Option<[Perm;4]> {
        if self.is_empty() {return None;}
        if let Some(sources)=&self.fiber_sources {
            for _ in 0..16 {
                self.counts[0].fetch_add(1,Ordering::Relaxed);
                let t=sources[rng.below(sources.len())];
                let proposal=if let Some(repair)=&self.coordinated {repair.repair(d,t,ambient,true,rng)} else {self.fiber_proposal(d,t,ambient,rng)};
                if let Some(u)=proposal {
                    if let Some(u)=self.admit(d,ambient,u,0,true,rng) {return Some(u);}
                }
            }
            return None;
        }
        // Amortize cycle extraction across 32 c proposals. All equivariant
        // block maps were solved at construction, not repeatedly in this loop.
        let variation=if !self.twists.is_empty() {
            self.proposal_gates[0].fetch_add(1,Ordering::Relaxed);
            let i=self.active[rng.below(self.active.len())];let (a,b)=self.frames[i].pair.unwrap();
            let h=self.frames[i].h.compose(self.twists[i].sample(d,rng),d);
            let cycles=permutation_cycles(d,h);
            // Identity relabelings impose no constraint on c. Keep those tuples
            // in the independent exploration lane, not this symmetry lane.
            if self.pure && h.is_identity(d) {self.proposal_gates[1].fetch_add(1,Ordering::Relaxed);return None;}
            if d>=4 && central_relabeling(d,a,b,h) && !h.is_identity(d) {self.proposal_gates[2].fetch_add(1,Ordering::Relaxed);return None;}
            if !centralizer_can_connect(d,a,b,&cycles) {self.proposal_gates[3].fetch_add(1,Ordering::Relaxed);return None;}
            if !self.envelope_admits(d,ambient,a,b,h,&cycles) {self.proposal_gates[4].fetch_add(1,Ordering::Relaxed);return None;}
            Some((i,h,cycles))
        } else {None};
        let joint_moves=if self.joint.is_some() {variation.as_ref().map(|v|centralizer_generators(d,&v.2,Ambient::Symmetric)).unwrap_or_default()} else {Vec::new()};
        for _ in 0..if self.repair_steps>0 || self.joint.is_some() {4} else {32} {
            self.counts[0].fetch_add(1,Ordering::Relaxed);
            let frame_index=variation.as_ref().map(|v|v.0).unwrap_or_else(||self.active[rng.below(self.active.len())]);
            let frame=&self.frames[frame_index];
            let (h,cycles)=variation.as_ref().map(|v|(v.1,&v.2)).unwrap_or((frame.h,&frame.cycles));
            if variation.is_none() && self.glue.is_none() {if let Some((a,b))=frame.pair {
                if !self.envelope_admits(d,ambient,a,b,h,cycles) {continue;}
            }}
            let t=if let Some((matches,_))=&self.glue {
                let eligible=&matches[frame_index];let j=eligible[rng.below(eligible.len())];let right=&self.frames[j];
                let t=if let Some(cache)=&self.boundary {
                    let Some(t)=cache.sample(d,frame_index,j,frame,right,rng) else {continue;};t
                } else {glue_tuple(d,frame,right,rng)};
                if self.repair_steps>0 {
                    self.pattern_counts[0].fetch_add(1,Ordering::Relaxed);
                    let (repaired,rescued,steps)=repair_gluing(d,t,&self.repair_moves[frame_index],self.repair_steps,rng);
                    self.pattern_counts[1].fetch_add(steps as u64,Ordering::Relaxed);
                    self.pattern_counts[2].fetch_add(u64::from(rescued),Ordering::Relaxed);
                    let Some(t)=repaired else {continue;};t
                } else {t}
            } else {
            let mut c=centralizer(d,cycles,rng);
            // Half of draws focus on small entry orders; half retain full C(h).
            if rng.below(2)==0 {
                let order=cycle_type(c,d).into_iter().fold(1,|a,b|a/gcd_usize(a,b)*b);
                let target=[2,3,4,6][rng.below(4)];c=power(d,c,order/gcd_usize(order,target));
            }
            if let Some((a,b))=frame.pair {[a,b,c,product_one_fourth(d,a,b,c)]} else {
                let q=match rng.below(3) {0=>Perm::identity(d),1=>power(d,h,1+rng.below(12)),_=>centralizer(d,&frame.cycles,rng)};
                let Some(x)=square_root(d,q,rng) else {continue;};
                let t=tuple_from_square(d,h,x,c);
                if power_identity(d,t[0],2) {continue;} // old involution-pair strata have their own engine
                t
            }};
            self.proposal_gates[5].fetch_add(1,Ordering::Relaxed);
            let Some(t)=self.fiber_proposal(d,t,ambient,rng) else {continue;};
            let proposal=if let Some(repair)=&self.coordinated {repair.repair(d,t,ambient,false,rng)} else {Some(t)};
            let Some(t)=proposal else {continue;};
            if let Some(repair)=&self.joint {if repair.admission_aware {
                if let Some(u)=repair.repair_admitted(d,t,&joint_moves,rng,
                    |u|self.repair_grade(d,ambient,u,frame.left_genus),
                    |u,rng|self.admit(d,ambient,u,frame.left_genus,frame.pair.is_none(),rng)) {return Some(u);}
                continue;
            }}
            let t=if let Some(repair)=&self.joint {
                let Some(u)=repair.repair(d,t,&joint_moves,rng) else {continue;};u
            } else {t};
            if let Some(t)=self.admit(d,ambient,t,if self.coordinated.is_some() {0} else {frame.left_genus},frame.pair.is_none(),rng) {return Some(t);}
        }
        None
    }

    pub fn summary(&self,out:&mut File)->io::Result<()> {
        self.summary_with_coverage(out,"biased symmetry sampling, not a census")
    }

    pub fn save_outputs(&self,dir:&Path,d:usize)->Result<(),String> {
        if let Some(repair)=&self.coordinated {repair.save_words(&dir.join("ESCAPE_WORDS.txt"))?;}
        if let Some(cache)=&self.boundary {cache.save(&dir.join("BOUNDARY_CLASSES.tsv"),d)?;}
        if self.local_audit.is_some() {
            let mut out=OpenOptions::new().write(true).create_new(true).open(dir.join("FRAME_LIBRARY.tsv")).map_err(|e|e.to_string())?;
            writeln!(out,"frame\tpair_genus\tpair_orbit_sizes\tboundary_cycle_type\ta_images_zero_based\tb_images_zero_based\tmarked_action_key").map_err(|e|e.to_string())?;
            for (i,f) in self.frames.iter().enumerate() {let (a,b)=f.pair.unwrap();
                writeln!(out,"{i}\t{}\t{}\t{:?}\t{}\t{}\t{}",f.left_genus,point_orbit_partition(d,&[a,b]),cycle_type(a.compose(b,d),d),
                    a.img[..d].iter().map(u8::to_string).collect::<Vec<_>>().join(","),b.img[..d].iter().map(u8::to_string).collect::<Vec<_>>().join(","),
                    colored_representation_key(d,[a,b,Perm::identity(d)])).map_err(|e|e.to_string())?;
            }
            out.sync_all().map_err(|e|e.to_string())?;
        }Ok(())
    }

    fn summary_with_coverage(&self,out:&mut File,coverage:&str)->io::Result<()> {
        if let Some(e)=self.witness_error.lock().unwrap().as_ref() {return Err(io::Error::other(e.clone()));}
        writeln!(out,"construction escape exports attempted = {}\nconstruction escape export limit = {}",
            self.witness_exports.load(Ordering::Relaxed),self.witness_limit)?;
        let construction=if self.fiber_sources.is_some() {"fixed-product cut/join pair fibers of certified sources; target parity may differ; exact pair partitions and genera preserved; no fixed-braid constraint"}
            else if self.glue.is_some() {"left/right pair covers; match p boundary cycles and offsets; genus=g_left+g_right+defect; no braid-fixed constraint"}
            else if !self.pure {"[x^2,h]=[c,h]=1; a=xh; b=h^-1*a*h"}
            else {"p=ab; k in C(a,b); h=p^r*k; c in C(h); sigma_1^(2r) is conjugation by h"};
        writeln!(out,"symmetry frames = {}\nconstruction = {construction}\ncoverage = {coverage}",self.frames.len())?;
        if let Some(audit)=&self.local_audit {audit.write(out)?;}
        if let Some(cache)=&self.boundary {cache.write(out)?;}
        if let Some(repair)=&self.coordinated {repair.write(out)?;}
        if let Some(repair)=&self.joint {repair.write(out)?;}
        for (label,count) in ["relabeling variations attempted","identity relabeling variations rejected",
            "central relabeling variations rejected","disconnected centralizer variations rejected",
            "proper envelope variations rejected","literal tuple proposals constructed"].iter().zip(&self.proposal_gates) {
            writeln!(out,"{label} = {}",count.load(Ordering::Relaxed))?;
        }
        writeln!(out,"pure braid power r = {}",self.period)?;
        writeln!(out,"active construction frames = {}\nredundant pair-connectivity frames excluded = {}",self.active.len(),self.discarded_frames)?;
        writeln!(out,"central pure-relabeling frames excluded = {}",self.central_frames)?;
        writeln!(out,"inverse/central equivalent power frames excluded = {}",self.alias_frames)?;
        writeln!(out,"boundary repair steps per start = {}",self.repair_steps)?;
        writeln!(out,"pre-generation weak-walk steps = {}\nfiber extra radius = {}\nfiber source states = {}",self.screen_steps,
            self.fiber_radius.map_or("disabled".into(),|r|r.to_string()),self.fiber_sources.as_ref().map_or(0,Vec::len))?;
        writeln!(out,"early known closed-state skips = {}",self.known_skips.load(Ordering::Relaxed))?;
        for (label,n) in ["fiber mutation starts","fiber cut trials","fiber accepted cuts","fiber nontrivial proposals"].iter().zip(&self.fiber_counts) {
            writeln!(out,"{label} = {}",n.load(Ordering::Relaxed))?;
        }
        for (label,n) in ["boundary repair starts","boundary repair local steps","boundary repair rescued starts","explicit middle-braid witnesses"].iter().zip(&self.pattern_counts) {
            writeln!(out,"{label} = {}",n.load(Ordering::Relaxed))?;
        }
        writeln!(out,"verified completed passports excluded = {}\ncompleted-passport proposal skips = {}",self.excluded.len(),self.completed_skips.load(Ordering::Relaxed))?;
        if let Some((lists,g))=&self.glue {
            let mut unique=HashSet::new();let mut stored=0;
            for list in lists {if unique.insert(Arc::as_ptr(list)) {stored+=list.len();}}
            writeln!(out,"gluing target genus = {g}\nunique shared boundary matching lists = {}\nstored right-frame references = {stored}\nlogical right-frame references = {}",unique.len(),lists.iter().map(|v|v.len()).sum::<usize>())?;
        }
        writeln!(out,"cached equivariant orbit maps = {}",self.twists.iter().flat_map(|t|&t.maps).map(Vec::len).sum::<usize>())?;
        let mut types=BTreeMap::new();
        for frame in &self.frames {*types.entry(frame.cycles.iter().map(Vec::len).collect::<Vec<_>>()).or_insert(0usize)+=1;}
        writeln!(out,"{} = {types:?}",if !self.twists.is_empty() {"base relabeling cycle types"} else {"symmetry cycle types"})?;
        let mut genera=BTreeMap::new();for frame in &self.frames {*genera.entry(frame.left_genus).or_insert(0usize)+=1;}
        writeln!(out,"fixed-pair genus lower bounds = {genera:?}")?;
        writeln!(out,"vary first-pair relabeling = {}",!self.twists.is_empty())?;
        for (label,n) in ["envelope exact checks","envelope cache hits","proper envelopes rejected"].iter().zip(&self.envelope_counts) {
            writeln!(out,"{label} = {}",n.load(Ordering::Relaxed))?;
        }
        for (label,n) in ["symmetry proposals","constructed tuples","index parity survivors","all-pair intransitive transitive tuples","preprobe survivors","full-group survivors"].iter().zip(&self.counts) {
            writeln!(out,"{label} = {}",n.load(Ordering::Relaxed))?;
        }
        Ok(())
    }
}

// Identical construction starts and per-start RNG streams. This is a bounded
// technical comparison, not a discovery rate or independent orbit count.
pub(super) fn benchmark_joint(args:&[String],d:usize,ambient:Ambient)->Result<(),String> {
    let cert=PathBuf::from(get_arg(args,"--certificate",None)?);verify(&cert)?;
    let case=parse_certificate_case(&cert)?;
    if case.d!=d || case.ambient!=ambient {return Err("benchmark source/ambient mismatch".into());}
    let starts:usize=parse_num(get_arg(args,"--benchmark-starts",Some("32"))?,"--benchmark-starts")?;
    let repeats:usize=parse_num(get_arg(args,"--benchmark-repeats",Some("3"))?,"--benchmark-repeats")?;
    if starts==0 || starts>256 || repeats==0 || repeats>5 {return Err("benchmark bounds: starts 1..256, repeats 1..5".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    if out.exists() {return Err("benchmark output must be new".into());}
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let seed:u64=parse_num(get_arg(args,"--seed",Some("20261010"))?,"--seed")?;
    let states=case.states.ok_or("benchmark needs closed source states")?;
    let known:HashSet<_>=states.iter().copied().collect();
    let sources=states.iter().copied().map(|k|unpack(k,d)).filter(|&t|
        pure_stability::equivariance_score(d,t,pure_stability::apply(d,t,&[0,0]))==4*d).collect::<Vec<_>>();
    if sources.is_empty() {return Err("source has no end-square-fixed state".into());}
    let mut panel=Vec::new();let mut rng=SplitMix64::new(seed);
    let mut input=File::create(out.join("STARTS.tsv")).map_err(|e|e.to_string())?;
    writeln!(input,"start\tseed\trepair_seed").map_err(|e|e.to_string())?;
    for i in 0..starts {
        let t=sources[i%sources.len()];
        let h=tuple_conjugator(d,t,pure_stability::apply(d,t,&[0,0]))?;
        let cycles=permutation_cycles(d,h);let moves=centralizer_generators(d,&cycles,Ambient::Symmetric);
        let c=if i%4==0 {t[2]} else {t[2].compose(centralizer(d,&cycles,&mut rng),d)};
        let u=[t[0],t[1],c,product_one_fourth(d,t[0],t[1],c)];
        let stream=seed ^ (i as u64+1).wrapping_mul(0xD1B54A32D192ED03);
        writeln!(input,"{i}\t{}\t{stream}",tuple_string_dynamic(d,&u)).map_err(|e|e.to_string())?;
        panel.push((u,moves,stream));
    }input.sync_all().map_err(|e|e.to_string())?;
    let mut results=File::create(out.join("BENCHMARK.tsv")).map_err(|e|e.to_string())?;
    writeln!(results,"repeat\torder\tarm\tstarts\tadmitted\tunique_admitted\tseconds\tknown_skips\tadmission_attempts\tfull_group_survivors").map_err(|e|e.to_string())?;
    for rep in 0..repeats {
        let mut order=[false,true];rng.shuffle(&mut order);
        for (slot,aware) in order.into_iter().enumerate() {
            let factory_args=vec!["--joint-source-only".into(),"--certificate".into(),cert.to_string_lossy().into(),"--skip-known-starts".into(),"--screen-steps".into(),"64".into()];
            let mut factory=Factory::new(&factory_args,d,ambient,"joint-pure")?;factory.set_known(&known);
            factory.joint=None; // The benchmark's separate repair owns these counters.
            let repair_args=if aware {vec!["--joint-admission-aware".into()]} else {Vec::new()};
            let repair=joint_pure::Repair::new(&repair_args)?;
            let begun=Instant::now();let mut accepted=0;let mut unique=HashSet::new();
            for (t,moves,stream) in &panel {
                let mut stream=SplitMix64::new(*stream);
                let got=if aware {repair.repair_admitted(d,*t,moves,&mut stream,
                    |u|factory.repair_grade(d,ambient,u,0),|u,r|factory.admit(d,ambient,u,0,false,r))}
                    else {repair.repair(d,*t,moves,&mut stream).and_then(|u|factory.admit(d,ambient,u,0,false,&mut stream))};
                if let Some(u)=got {accepted+=1;unique.insert(canonical(&u,d));}
            }
            let elapsed=begun.elapsed().as_secs_f64();
            writeln!(results,"{rep}\t{slot}\t{}\t{starts}\t{accepted}\t{}\t{elapsed:.9}\t{}\t{}\t{}",if aware {"admission-aware"} else {"legacy"},unique.len(),factory.known_skips.load(Ordering::Relaxed),factory.counts[1].load(Ordering::Relaxed),factory.counts[5].load(Ordering::Relaxed)).map_err(|e|e.to_string())?;
            let mut summary=File::create(out.join(format!("{rep}_{slot}_SUMMARY.txt"))).map_err(|e|e.to_string())?;
            repair.write(&mut summary).map_err(|e|e.to_string())?;
            factory.summary(&mut summary).map_err(|e|e.to_string())?;summary.sync_all().map_err(|e|e.to_string())?;
        }
    }results.sync_all().map_err(|e|e.to_string())?;
    Ok(())
}

// Finite marked S3 actions: every multiset of 1/2/3/6-point transitive
// actions, faithful involution pair, r=1,2, all k in C(a,b), all c in C(p^r k).
// This is not all tuples of a degree. Proper envelopes and witnessed weak
// successes are exact exclusions; surviving generating seeds use shared closure.
/// A bounded shard of the entire (one cut, one cut) S18 -> A18 population.
/// Construction leads are durable before the shared closure engine starts.
pub(super) fn fiber_population(args:&[String],d:usize,ambient:Ambient,out:&Path)->Result<Vec<[Perm;4]>,String> {
    if d!=18 || ambient!=Ambient::Alternating {return Err("pair-fiber-exact supports the S18-to-A18 one-cut bridge only".into());}
    if optional_arg(args,"--fiber-radius").is_some_and(|r|r!="0") {return Err("exact one-cut bridge has radius zero only".into());}
    let source=PathBuf::from(get_arg(args,"--certificate",None)?);verify(&source)?;
    let case=parse_certificate_case(&source)?;
    if case.d!=d || case.ambient!=Ambient::Symmetric {return Err("pair-fiber-exact needs an S18 source".into());}
    let mut states=case.states.ok_or("exact fiber source needs closed states")?.into_iter().collect::<Vec<_>>();
    states.sort_unstable();
    if states.iter().any(|&k|unpack(k,d).iter().any(|p|p.is_even(d))) {return Err("exact bridge needs four odd source entries".into());}
    let first:usize=parse_num(get_arg(args,"--fiber-state-start",Some("0"))?,"--fiber-state-start")?;
    let count:usize=parse_num(get_arg(args,"--fiber-state-count",Some("8"))?,"--fiber-state-count")?;
    if count==0 || count>8 || first>=states.len() {return Err("exact fiber requires an in-range shard of 1..8 source states".into());}
    let end=first.saturating_add(count).min(states.len());
    let create=|name:&str|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|e.to_string());
    let mut durable=create("construction_survivors.txt")?;
    let mut table=create("FIBER_STATES.tsv")?;
    writeln!(table,"source_state\tleft_cuts\tright_cuts\tcut_products\tnew_generating_leads").map_err(|e|e.to_string())?;
    let factory=Factory::new(args,d,ambient,"pair-fiber")?;
    let mut rng=SplitMix64::new(parse_num(get_arg(args,"--seed",Some("20261004"))?,"--seed")?);
    let mut seen=HashSet::new();let mut replay=Vec::new();let mut products=0u64;let mut duplicates=0u64;
    let start=Instant::now();
    for (i,&key) in states.iter().enumerate().take(end).skip(first) {
        let t=unpack(key,d);let left=pair_fiber::one_cut_pairs(d,t[0],t[1]);let right=pair_fiber::one_cut_pairs(d,t[2],t[3]);
        let before=replay.len();
        for &(a,b) in &left {for &(c,e) in &right {
            products+=1;let u=[a,b,c,e];let indices=u.map(|p|index(p,d));
            if indices.iter().any(|&x|x==0 || x>d-3) || has_transitive_pair(d,&u) {continue;}
            // Canonicalize only after cheap necessary tests, before repeated
            // witness walks/group checks. Partition preservation gives global
            // transitivity, genus and parity from the verified source.
            if !seen.insert(canonical(&u,d)) {duplicates+=1;continue;}
            if let Some(u)=factory.admit(d,ambient,u,0,true,&mut rng) {
                writeln!(durable,"group=A_18 seed={}",tuple_string_dynamic(d,&u)).map_err(|e|e.to_string())?;
                durable.sync_data().map_err(|e|e.to_string())?;replay.push(u);
            }
        }}
        writeln!(table,"{i}\t{}\t{}\t{}\t{}",left.len(),right.len(),left.len()*right.len(),replay.len()-before).map_err(|e|e.to_string())?;
        table.sync_data().map_err(|e|e.to_string())?;
    }
    let mut summary=create("FIBER_SUMMARY.txt")?;
    writeln!(summary,"finite construction complete = true\nfinite construction = all partition/index-sum preserving one-cut choices in each pair, on the specified sorted S18 source-state shard\nfinite source states = {}\nfinite source state start = {first}\nfinite source state end = {end}\nfinite raw choice upper bound = {}\nfinite cut products = {products}\nfinite canonical duplicates = {duplicates}\nfinite distinct necessary survivors = {}\nfinite distinct generating leads = {}\nfinite construction seconds = {:.3}\nwhole degree complete = false\nquotient = simultaneous S_d conjugation\nsource certificate = {}",states.len(),(end-first)*153*153,seen.len(),replay.len(),start.elapsed().as_secs_f64(),source.display()).map_err(|e|e.to_string())?;
    factory.summary_with_coverage(&mut summary,"exact specified one-cut construction shard; closure resolution reported separately").map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())?;Ok(replay)
}

pub(super) fn s3_population(args:&[String],d:usize,ambient:Ambient,out:&Path)->Result<Vec<[Perm;4]>,String> {
    if !(9..=11).contains(&d) {return Err("s3-exact supports degrees 9..11".into());}
    let excluded=excluded_passports(args,d,ambient)?;
    let create=|name:&str|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|e.to_string());
    let mut durable=create("construction_survivors.txt")?;
    let mut table=BufWriter::new(create("S3_FRAMES.tsv")?);
    writeln!(table,"frame\tr\tall_pair_twists\th_representatives\tfull_envelopes\tc_enumerated").map_err(|e|e.to_string())?;
    let mut frames=pure_frames(d,ambient,"S3")?;
    frames.retain(|f|{let (a,b)=f.pair.unwrap();power_identity(d,a,2)&&power_identity(d,b,2)
        && index(a,d)+index(b,d)>=4 && subgroup_order(d,&[a,b])==6 && !redundant_pair_connectivity(d,a,b)});
    let frame_cap:usize=parse_num(get_arg(args,"--frame-cap",Some("18446744073709551615"))?,"--frame-cap")?;
    if frame_cap==0 {return Err("positive frame-cap required".into());}
    let mut counts=[0u64;9];let mut replay=Vec::new();let mut seen=HashSet::new();
    let mut rng=SplitMix64::new(parse_num(get_arg(args,"--seed",Some("20261002"))?,"--seed")?);
    let start=Instant::now();
    for (i,frame) in frames.iter().take(frame_cap).enumerate() {
        let (a,b)=frame.pair.unwrap();let cache=PairTwists::new(d,a,b);
        let (twists,complete)=cache.enumerate(d,1_000_000);
        if !complete {return Err("S3 twist population exceeds validated finite limit".into());}
        let conjugators=cache.generators(d);
        for r in 1..=2 {
            let base=power(d,a.compose(b,d),r);let mut covered=HashSet::new();
            let (mut representatives,mut full,mut cs)=(0,0,0usize);
            for &k in &twists {
                let h=base.compose(k,d);if covered.contains(&h) {continue;}
                // Exact conjugacy in C_S(a,b), valid for simultaneous S_d
                // quotient even when all tuple entries are even.
                let mut orbit=vec![h];covered.insert(h);let mut head=0;
                while head<orbit.len() {let u=orbit[head];head+=1;
                    for &q in &conjugators {let v=u.conjugate_by(q,d);if covered.insert(v) {orbit.push(v);}}
                }
                representatives+=1;counts[0]+=1;
                if h.is_identity(d) || central_relabeling(d,a,b,h) {continue;}
                let cycles=permutation_cycles(d,h);
                if !centralizer_can_connect(d,a,b,&cycles) || envelope_order(d,ambient,a,b,&cycles)!=ambient.target_order(d) {continue;}
                full+=1;counts[1]+=1;
                let (centralizer,complete)=PairTwists::new(d,h,Perm::identity(d)).enumerate(d,1_000_000);
                if !complete {return Err("S3 completion population exceeds validated finite limit".into());}
                cs+=centralizer.len();
                for c in centralizer {
                    counts[2]+=1;let t=[a,b,c,product_one_fourth(d,a,b,c)];let indices=t.map(|p|index(p,d));
                    if indices.iter().any(|&x|x==0 || x>d-3) || indices.iter().sum::<usize>()<2*d-2
                        || (ambient==Ambient::Alternating && indices.iter().any(|x|x%2!=0))
                        || (ambient==Ambient::Symmetric && indices.iter().all(|x|x%2==0)) {continue;}
                    counts[3]+=1;
                    let mut passport=t.map(|p|cycle_type(p,d));passport.sort();
                    if excluded.contains(&passport) {counts[4]+=1;continue;}
                    if !generators_transitive(d,&t[..3]) || has_transitive_pair(d,&t) {continue;}
                    counts[5]+=1;
                    if probe(d,t,&mut rng,16,true).0 {continue;}
                    counts[6]+=1;
                    if !tuple_generates_dynamic(d,ambient,&t) {continue;}
                    counts[7]+=1;let key=canonical(&t,d);if !seen.insert(key) {counts[8]+=1;continue;}
                    writeln!(durable,"group={}_{} seed={}",ambient.short(),d,tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
                    durable.sync_data().map_err(|e|e.to_string())?;replay.push(t);
                }
            }
            writeln!(table,"{i}\t{r}\t{}\t{representatives}\t{full}\t{cs}",twists.len()).map_err(|e|e.to_string())?;
            table.flush().map_err(|e|e.to_string())?;
            eprintln!("S3 finite frame {}/{} r={r}: h={representatives} envelopes={full} c={cs} distinct leads={} elapsed={:.1}s",i+1,frames.len(),replay.len(),start.elapsed().as_secs_f64());
        }
    }
    let mut summary=create("S3_SUMMARY.txt")?;
    writeln!(summary,"finite construction complete = {}\nfinite construction = faithful S3 involution-pair action mixtures of sizes 1,2,3,6; r=1,2; nonidentity h=(ab)^r k; k in C_S(a,b); c in C_G(h)\nfinite frames = {}\nfinite frames inspected = {}\nfinite excluded passports = {}\nfinite distinct generating leads = {}\nfinite construction seconds = {:.3}\nwhole degree complete = false\nquotient = simultaneous S_d conjugation\nIdentity h is unrestricted, excluded from this symmetry population; r=3 is algebraically impossible for nonidentity h. Excluded passports refer to separately verified completed censuses.",frames.len()<=frame_cap,frames.len(),frames.len().min(frame_cap),excluded.len(),replay.len(),start.elapsed().as_secs_f64()).map_err(|e|e.to_string())?;
    for (label,count) in ["h representatives","full construction envelopes","centralizer completions","index parity survivors","completed-passport skips","all-pair intransitive transitive tuples","preprobe survivors","full-group occurrences","repeated generating tuples"].iter().zip(counts) {
        writeln!(summary,"finite {label} = {count}").map_err(|e|e.to_string())?;
    }
    summary.sync_all().map_err(|e|e.to_string())?;Ok(replay)
}

pub(super) fn audit(args:&[String],d:usize,ambient:Ambient)->Result<(),String> {
    let family=get_arg(args,"--family",Some("pure-symmetry"))?;
    let factory=Factory::new(args,d,ambient,&family)?;
    if factory.glue.is_some() || factory.frames.iter().any(|f|f.pair.is_none()) {
        return Err("envelope audit needs fixed first-pair construction frames".into());
    }
    let frame_cap:usize=parse_num(get_arg(args,"--frame-cap",Some("64"))?,"--frame-cap")?;
    let twist_cap:usize=parse_num(get_arg(args,"--twist-cap",Some("256"))?,"--twist-cap")?;
    if frame_cap==0 || !(1..=1_000_000).contains(&twist_cap) {return Err("positive frame-cap and 1..1000000 twist-cap required".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let create=|name:&str|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|e.to_string());
    let mut rows=BufWriter::new(create("ENVELOPES.tsv")?);
    writeln!(rows,"frame\ttwist\tall_frame_twists_enumerated\tidentity_h\tconnectivity_possible\tenvelope_order\tfull_ambient_possible\tcentralizer_order\th_cycle_type\tfirst_pair").map_err(|e|e.to_string())?;
    let (mut tested,mut rejected,mut complete_frames)=(0,0,0);
    for (i,f) in factory.frames.iter().take(frame_cap).enumerate() {
        let (a,b)=f.pair.unwrap();
        let (twists,complete)=if factory.twists.is_empty() {(vec![Perm::identity(d)],true)} else {factory.twists[i].enumerate(d,twist_cap)};
        complete_frames+=usize::from(complete);
        for (j,k) in twists.into_iter().enumerate() {
            let h=f.h.compose(k,d);let cycles=permutation_cycles(d,h);
            let connect=centralizer_can_connect(d,a,b,&cycles);let order=envelope_order(d,ambient,a,b,&cycles);
            let full=order==ambient.target_order(d);tested+=1;rejected+=usize::from(!full);
            let centralizer_order=subgroup_order(d,&centralizer_generators(d,&cycles,ambient));
            writeln!(rows,"{i}\t{j}\t{complete}\t{}\t{connect}\t{order}\t{full}\t{centralizer_order}\t{:?}\t{}",h.is_identity(d),cycle_type(h,d),tuple_string_dynamic(d,&[a,b,Perm::identity(d),a.compose(b,d).inverse(d)])).map_err(|e|e.to_string())?;
        }
    }
    rows.flush().map_err(|e|e.to_string())?;
    let mut summary=create("SUMMARY.txt")?;
    writeln!(summary,"population = {family}; {}{d}; selected finite frames and centralizer twists\nframes in library = {}\nframes inspected = {}\nframes with full twist enumeration = {complete_frames}\nenvelopes checked = {tested}\nproper envelopes excluded = {rejected}\nall library frames inspected = {}\nwhole-degree coverage = false\nA proper envelope is an exact exclusion of all c in C_G(h). A full envelope is only necessary. Identity h rows are reported but excluded from symmetry sampling.",ambient.short(),factory.frames.len(),factory.frames.len().min(frame_cap),factory.frames.len()<=frame_cap).map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn marked_groups_and_new_dihedral_actions_are_exact() {
        for (name,order) in [("S3",6),("A4",12),("S4",24),("A5",60),("D5",10),("D7",14),("D9",18),("D11",22),("D13",26)] {
            let group=SmallPairGroup::new(name).unwrap();assert_eq!(group.table.len(),order);
            let all=group.markings(false);let novel=group.markings(true);
            let keys=all.iter().map(|&p|group.marking_key(p).unwrap()).collect::<HashSet<_>>();assert_eq!(keys.len(),all.len());
            let [x,y]=group.base;let xy=group.table[x][y];
            let old=[[x,y],[xy,y],[x,xy]].map(|p|group.marking_key(p).unwrap()).into_iter().collect::<HashSet<_>>();
            assert_eq!(novel.len()+old.len(),all.len());
            for pair in all {for (d,a,b) in group.actions(pair) {
                assert!(generators_transitive(d,&[a,b]));assert_eq!(order as u128%subgroup_order(d,&[a,b]),0);
                // Literal multiplication on cosets must respect the marked
                // abstract generators' orders and their product order.
                for (g,p) in [(pair[0],a),(pair[1],b),(group.table[pair[0]][pair[1]],a.compose(b,d))] {
                    let mut v=g;let mut n=1;while v!=0 {v=group.table[v][g];n+=1;}
                    assert!(power_identity(d,p,n));
                }
            }}
            println!("MARKED_MODEL model={name} order={order} classes={} novel={} actions={}",keys.len(),novel.len(),group.subgroups.len());
        }
        let args=["--pair-model","D9","--marked-pairs","all","--pure-power","4","--novel-power"].map(str::to_string);
        let factory=Factory::new(&args,20,Ambient::Symmetric,"marked-symmetry").unwrap();
        assert!(!factory.is_empty());assert!(factory.alias_frames>0);
        for f in &factory.frames {let (a,b)=f.pair.unwrap();
            assert!(!central_relabeling(20,a,b,f.h));
            for r in 1..=3 {let p=a.compose(b,20);
                assert!(!central_relabeling(20,a,b,power(20,p,4-r)));
                assert!(!central_relabeling(20,a,b,power(20,p,4+r)));
            }
        }
    }

    #[test]
    fn post_cut_frames_delay_only_mutable_entry_filters() {
        let invalid=["--fiber-radius","2"].map(str::to_string);
        assert!(fiber_population(&invalid,18,Ambient::Alternating,Path::new(".")).is_err());
        assert!(fiber_population(&[],18,Ambient::Symmetric,Path::new(".")).is_err());
        let old=pair_frames(&[],9,Ambient::Alternating,3).unwrap();
        let args=["--post-cut-base-frames"].map(str::to_string);
        let expanded=pair_frames(&args,9,Ambient::Alternating,3).unwrap();
        let keys=expanded.iter().map(|f|f.pair.unwrap()).collect::<HashSet<_>>();
        assert!(old.iter().all(|f|keys.contains(&f.pair.unwrap())));
        assert!(expanded.iter().any(|f|{let (a,b)=f.pair.unwrap();!a.is_even(9)||!b.is_even(9)}));
        assert!(expanded.iter().all(|f|{let (a,b)=f.pair.unwrap();!pair_transitive(9,a,b)&&f.left_genus<=3}));
        let args=["--post-cut-base-frames","--glue-genus","1","--pair-genus-max","3"].map(str::to_string);
        let factory=Factory::new(&args,10,Ambient::Symmetric,"cut-join-glue").unwrap();
        assert_eq!(factory.discarded_frames,0);
        assert!(factory.active.iter().any(|&i|{let (a,b)=factory.frames[i].pair.unwrap();
            a.is_even(10)&&b.is_even(10)&&factory.glue.as_ref().unwrap().0[i].iter().any(|&j|{
                let (c,e)=factory.frames[j].pair.unwrap();c.is_even(10)&&e.is_even(10)})}));
        let mut rng=SplitMix64::new(20261004);let identity=Perm::identity(10);
        assert!(factory.admit(10,Ambient::Symmetric,[identity;4],0,true,&mut rng).is_none());
        let (a,b)=factory.frames.iter().find_map(|f|f.pair.filter(|(a,b)|a.is_even(10)&&b.is_even(10))).unwrap();
        let t=[a,b,identity,a.compose(b,10).inverse(10)];
        assert!(factory.admit(10,Ambient::Symmetric,t,0,true,&mut rng).is_none());
        assert!(Factory::new(&args,10,Ambient::Symmetric,"pair-glue").is_err());
    }

    #[test]
    fn fiber_factories_target_parity_and_known_skips_are_exact() {
        let args=["--certificate","evidence/certificates/CERTIFIED_S18_genus1_weak.txt",
            "--fiber-radius","0","--screen-steps","64"].map(str::to_string);
        let factory=Factory::new(&args,18,Ambient::Alternating,"pair-fiber").unwrap();
        assert!(!factory.is_empty());assert_eq!(factory.fiber_sources.as_ref().unwrap().len(),40);
        let mut rng=SplitMix64::new(20261004);
        for _ in 0..8 {if let Some(t)=factory.generate(18,Ambient::Alternating,&mut rng) {
            assert!(t.iter().all(|p|p.is_even(18)));assert_eq!(genus_of(&t,18),1);
        }}
        assert!(factory.fiber_counts[3].load(Ordering::Relaxed)>0);
        let args=["--certificate","evidence/certificates/CERTIFIED_A25_XXYY_genus4_weak.txt",
            "--skip-known-starts","--screen-steps","64"].map(str::to_string);
        let mut factory=Factory::new(&args,25,Ambient::Alternating,"pair-fiber").unwrap();
        let states=parse_certificate_case(Path::new(&args[1])).unwrap().states.unwrap();
        factory.set_known(&states.iter().copied().collect());
        for key in states {assert!(factory.admit(25,Ambient::Alternating,unpack(key,25),0,true,&mut rng).is_none());}
        assert_eq!(factory.known_skips.load(Ordering::Relaxed),6);
        assert_eq!(factory.counts[4].load(Ordering::Relaxed),0);
    }

    #[test]
    #[ignore = "identical-input paired work-order benchmark; explicitly run on one worker"]
    fn benchmark_witness_before_generation() {
        for d in [11,18,25] {
            let ambient=if d==25 {Ambient::Alternating} else {Ambient::Symmetric};
            let inputs=if d==25 {
                parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_A25_XXYY_genus4_weak.txt"))
                    .unwrap().states.unwrap().into_iter().map(|k|unpack(k,d)).collect::<Vec<_>>()
            } else {
                let args=["--pair-model","S4","--marked-pairs","novel","--pure-power","1"].map(str::to_string);
                let factory=Factory::new(&args,d,ambient,"marked-symmetry").unwrap();
                let mut rng=SplitMix64::new(20261004+d as u64);let mut inputs=Vec::new();
                for _ in 0..100_000 {if let Some(t)=factory.generate(d,ambient,&mut rng) {inputs.push(t);if inputs.len()==128 {break;}}}
                assert_eq!(inputs.len(),128,"benchmark needs a fixed full-generation input panel");inputs
            };
            for t in inputs.iter().take(4) {let mut rng=SplitMix64::new(9001);
                std::hint::black_box(probe(d,*t,&mut rng,64,true).0);std::hint::black_box(tuple_generates_dynamic(d,ambient,t));}
            let mut reference=None;
            for repeat in 0..6 {for witness_first in if repeat%2==0 {[false,true]} else {[true,false]} {
                let start=Instant::now();let mut checks=0;let mut decisions=Vec::new();
                for (i,t) in inputs.iter().enumerate() {let mut rng=SplitMix64::new(9001+i as u64);
                    let accepted=if witness_first {
                        if probe(d,*t,&mut rng,64,true).0 {false} else {checks+=1;tuple_generates_dynamic(d,ambient,t)}
                    } else {checks+=1;tuple_generates_dynamic(d,ambient,t) && !probe(d,*t,&mut rng,64,true).0};
                    decisions.push(accepted);
                }
                if let Some(ref r)=reference {assert_eq!(&decisions,r);} else {reference=Some(decisions.clone());}
                println!("ORDER_BENCH degree={d} repeat={repeat} witness_first={witness_first} starts={} accepted={} group_checks={checks} seconds={:.9}",
                    inputs.len(),decisions.iter().filter(|&&x|x).count(),start.elapsed().as_secs_f64());
            }}
        }
    }

    #[test]
    fn bridge_bank_and_boundary_repair_preserve_all_controls() {
        let mut rng=SplitMix64::new(20261003);
        for name in ["CERTIFIED_A25_weak.txt","CERTIFIED_A25_order3_weak.txt","CERTIFIED_A25_order6_genus12_weak.txt",
            "CERTIFIED_A25_order12_genus12_weak.txt","CERTIFIED_A25_mixed_genus0_weak.txt","CERTIFIED_A25_XXYY_genus4_weak.txt",
            "CERTIFIED_S18_genus1_weak.txt","CERTIFIED_S18_XXYY_genus1_weak.txt","CERTIFIED_S18_XXXX_genus1_weak.txt"] {
            let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();let d=case.d;
            for key in case.states.unwrap() {let t=unpack(key,d);assert_eq!(bridge_score(d,t),0,"{name}");
                let p=t[0].compose(t[1],d);let moves=centralizer_generators(d,&permutation_cycles(d,p),Ambient::Symmetric);
                let passport=t.map(|g|cycle_type(g,d));let invariant=defect_of(&t,d);
                for s in &moves {let mut u=t;u[2]=u[2].conjugate_by(*s,d);u[3]=u[3].conjugate_by(*s,d);
                    assert!(tuple_product_is_identity(d,&u));assert_eq!(u.map(|g|cycle_type(g,d)),passport);
                    assert_eq!(defect_of(&u,d),invariant);
                }
                let (repaired,rescued,steps)=repair_gluing(d,t,&moves,16,&mut rng);
                assert_eq!(repaired,Some(t));assert!(!rescued);assert_eq!(steps,0);
            }
        }
        for line in fs::read_to_string("evidence/weak_patterns_20261003/SHORT_WITNESSES.tsv").unwrap().lines().skip(1) {
            let fields=line.split('\t').collect::<Vec<_>>();let d=fields[0][1..].parse::<usize>().unwrap();let t=parse_tuple(fields[2],d).unwrap();
            assert!(tuple_product_is_identity(d,&t));assert!(!has_transitive_pair(d,&t));assert!(bridge_score(d,t)>0);
            let mut u=t;for mv in [1,3,1] {u=hurwitz_move_dynamic(d,u,mv);}
            assert!(pair_transitive(d,u[0],u[1]) && pair_transitive(d,u[2],u[3]));
            for first in 0..6 {let u=hurwitz_move_dynamic(d,t,first);assert!(!pair_transitive(d,u[0],u[1]));
                for last in 0..6 {let v=hurwitz_move_dynamic(d,u,last);assert!(!pair_transitive(d,v[0],v[1]));}
            }
        }
        let args=["--glue-genus","1","--pair-genus-max","3"].map(str::to_string);
        let factory=Factory::new(&args,10,Ambient::Symmetric,"resistant-glue").unwrap();
        for &i in factory.active.iter().take(24) {let left=&factory.frames[i];
            let literal=factory.frames.iter().enumerate().filter_map(|(j,right)|{
                if left.cycles.iter().map(Vec::len).ne(right.cycles.iter().map(Vec::len)) {return None;}
                let (a,b)=left.pair.unwrap();let (x,y)=right.pair.unwrap();
                let defect=left.cycles.len() as i64-point_orbits(10,&[a,b]).len() as i64-point_orbits(10,&[x,y]).len() as i64+1;
                (defect>=0 && defect+left.left_genus+right.left_genus==1 && [a,b,x,y].iter().any(|p|!p.is_even(10))).then_some(j)
            }).collect::<Vec<_>>();
            assert_eq!(factory.glue.as_ref().unwrap().0[i].as_ref(),&literal);
            for _ in 0..8 {let right=&factory.frames[literal[rng.below(literal.len())]];let t=glue_tuple(10,left,right,&mut rng);
                let (candidate,_,_)=repair_gluing(10,t,&factory.repair_moves[i],16,&mut rng);
                if let Some(u)=candidate {assert_eq!(bridge_score(10,u),0);assert!(tuple_product_is_identity(10,&u));
                    assert_eq!(u.map(|g|cycle_type(g,10)),t.map(|g|cycle_type(g,10)));assert_eq!(defect_of(&u,10),defect_of(&t,10));}
            }
        }
    }

    #[test]
    fn construction_escape_exports_are_literal_bounded_and_optional() {
        let t=parse_tuple("[(2 3)(5 7)(6 8)(9 10), (3 4)(5 6)(7 9)(8 10), (1 6 10)(2 9)(3 8)(4 5), (1 10 7)(2 5)(3 9)(4 8)]",10).unwrap();
        assert!(bridge_score(10,t)>0 && !has_transitive_pair(10,&t));
        let args=["--construction-witness-limit","1"].map(str::to_string);
        let mut factory=Factory::new(&args,10,Ambient::Symmetric,"pure-symmetry").unwrap();
        assert!(factory.set_positive(None).is_err());
        let out=std::env::temp_dir().join(format!("hurwitz-construction-word-{}-{}",std::process::id(),
            std::time::SystemTime::now().duration_since(std::time::SystemTime::UNIX_EPOCH).unwrap().as_nanos()));
        fs::create_dir(&out).unwrap();
        let cache=Arc::new(positive_witness::Cache::new(&[],10,&out).unwrap());
        factory.set_positive(Some(cache.clone())).unwrap();
        assert!(factory.admit(10,Ambient::Symmetric,t,0,false,&mut SplitMix64::new(7)).is_none());
        assert_eq!(factory.witness_exports.load(Ordering::Relaxed),1);
        assert!(cache.contains(canonical(&t,10)));
        factory.save_escape(10,t,&[6]); // Exhausted quota cannot introduce a false proof.
        assert!(factory.witness_error.lock().unwrap().is_none());
        let payload=fs::read_to_string(out.join("POSITIVE_WITNESSES.tsv")).unwrap();
        assert_eq!(payload.lines().count(),2);
        let fields=payload.lines().nth(1).unwrap().split('\t').collect::<Vec<_>>();
        let mut u=parse_tuple(fields[2],10).unwrap();
        for mv in fields[1].split(',') {u=hurwitz_move_dynamic(10,u,mv.parse().unwrap());}
        assert!(pair_transitive(10,u[0],u[1]));
        drop(factory);drop(cache);
        fs::remove_file(out.join("POSITIVE_WITNESSES.tsv")).unwrap();fs::remove_dir(out).unwrap();
    }

    #[test]
    #[ignore = "seeded paired construction benchmark; explicitly run on one worker"]
    fn benchmark_boundary_resistance() {
        for (d,ambient,genus) in [(10,Ambient::Symmetric,1),(18,Ambient::Symmetric,1),(25,Ambient::Alternating,4)] {
            let args=["--glue-genus".into(),genus.to_string(),"--pair-genus-max".into(),"3".into(),
                "--known-certificates".into(),"evidence/certificates".into()];
            let factory=Factory::new(&args,d,ambient,"resistant-glue").unwrap();assert!(!factory.is_empty());
            for repeat in 0..6 {for repaired in if repeat%2==0 {[false,true]} else {[true,false]} {
                let mut rng=SplitMix64::new(831+repeat);let mut raw=0;let mut survivors=0;let mut rescued=0;let mut evals=0;
                let start=Instant::now();
                for _ in 0..1000 {let i=factory.active[rng.below(factory.active.len())];let left=&factory.frames[i];
                    let right_indices=&factory.glue.as_ref().unwrap().0[i];let right=&factory.frames[right_indices[rng.below(right_indices.len())]];
                    let t=glue_tuple(d,left,right,&mut rng);raw+=1;
                    // Split proposal and repair RNGs so both arms receive the
                    // identical thousand starts, not merely the same seed.
                    let mut repair_rng=SplitMix64::new(rng.next_u64());
                    let (candidate,changed,n)=if repaired {repair_gluing(d,t,&factory.repair_moves[i],16,&mut repair_rng)}
                        else {(if bridge_score(d,t)==0 {Some(t)} else {None},false,0)};
                    evals+=n;rescued+=usize::from(changed);survivors+=usize::from(candidate.is_some());
                    std::hint::black_box(candidate);
                }
                println!("BRIDGE_BENCH degree={d} group={} genus={genus} repeat={repeat} repaired={repaired} starts={raw} survivors={survivors} rescued={rescued} steps={evals} seconds={:.9}",ambient.short(),start.elapsed().as_secs_f64());
            }}
        }
    }
    #[test]
    fn redundant_pair_exclusion_and_capped_twist_enumeration() {
        let d=5;let mut rng=SplitMix64::new(20261005);
        for _ in 0..30 {
            let a=random_perm(d,Ambient::Symmetric,&mut rng);let b=random_perm(d,Ambient::Symmetric,&mut rng);
            let cache=PairTwists::new(d,a,b);let (all,complete)=cache.enumerate(d,120);assert!(complete);
            let mut literal=HashSet::new();let mut c=Perm::identity(d);
            loop {
                if [a,b].iter().all(|&p|p.compose(c,d)==c.compose(p,d)) {literal.insert(c);}
                if redundant_pair_connectivity(d,a,b) && generators_transitive(d,&[a,b,c]) {
                    assert!(pair_transitive(d,a,c)||pair_transitive(d,b,c));
                }
                if !next_permutation(&mut c.img[..d]) {break;}
            }
            assert_eq!(all.into_iter().collect::<HashSet<_>>(),literal);
            assert_eq!(subgroup_order(d,&cache.generators(d)),literal.len() as u128);
            let (prefix,complete)=cache.enumerate(d,1);assert_eq!(prefix.len(),1);assert_eq!(complete,literal.len()==1);
        }
        for name in ["CERTIFIED_S18_XXXX_genus1_weak.txt","CERTIFIED_A25_mixed_genus0_weak.txt"] {
            let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();
            for key in case.states.unwrap() {let t=unpack(key,case.d);assert!(!redundant_pair_connectivity(case.d,t[0],t[1]));}
        }
    }
    #[test]
    fn s3_actions_central_pruning_and_positive_pair_gluing() {
        let models=pure_models("S3").unwrap();assert_eq!(models.iter().map(|m|m.0).collect::<Vec<_>>(),vec![1,2,3,6]);
        for &(d,a,b) in &models {assert!(generators_transitive(d,&[a,b]));assert_eq!(6%subgroup_order(d,&[a,b]),0);}
        for d in 9..=11 {
            let args=["--pair-model","S3","--pure-power","3"].map(str::to_string);
            let factory=Factory::new(&args,d,Ambient::Symmetric,"pure-symmetry").unwrap();
            assert!(factory.central_frames>0);
            assert!(factory.frames.iter().all(|f|{let (a,b)=f.pair.unwrap();!central_relabeling(d,a,b,f.h)}));
            let frames=pure_frames(d,Ambient::Alternating,"S3").unwrap();
            let expected=if d==9 {162} else {Ambient::Alternating.target_order(d)};
            assert!(frames.iter().filter(|f|{let (a,b)=f.pair.unwrap();power_identity(d,a,2)&&power_identity(d,b,2)
                && point_orbits(d,&[a,b]).iter().filter(|p|p.len()>1).map(Vec::len).collect::<BTreeSet<_>>()==BTreeSet::from([3,6])})
                .any(|f|{let (a,b)=f.pair.unwrap();PairTwists::new(d,a,b).enumerate(d,100).0.into_iter().any(|k|{
                    let h=f.h.compose(k,d);let cycles=permutation_cycles(d,h);
                    cycle_type(h,d).iter().filter(|&&x|x==3).count()==2 && envelope_order(d,Ambient::Alternating,a,b,&cycles)==expected})}));
        }
        let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_A25_XXYY_genus4_weak.txt")).unwrap();
        let states=case.states.unwrap();let mut defects=BTreeMap::new();let mut rng=SplitMix64::new(20261002);
        for key in states {let t=unpack(key,25);let (defect,left_g,right_g)=defect_of(&t,25);
            *defects.entry(defect).or_insert(0)+=1;assert_eq!(genus_of(&t,25),4);
            let left=Frame {h:t[0].compose(t[1],25),cycles:permutation_cycles(25,t[0].compose(t[1],25)),pair:Some((t[0],t[1])),left_genus:left_g};
            let right=Frame {h:left.h,cycles:permutation_cycles(25,t[3].inverse(25).compose(t[2].inverse(25),25)),pair:Some((t[3].inverse(25),t[2].inverse(25))),left_genus:right_g};
            for _ in 0..8 {let glued=glue_tuple(25,&left,&right,&mut rng);assert!(tuple_product_is_identity(25,&glued));
                assert_eq!(defect_of(&glued,25),(defect,left_g,right_g));assert_eq!(genus_of(&glued,25),4);}
            assert!(!central_relabeling(25,t[0],t[1],t[0].compose(t[1],25)));
        }
        assert_eq!(defects,BTreeMap::from([(2,2),(4,4)]));
    }
    #[test]
    fn asymmetric_planar_gluing_matches_both_boundaries() {
        let mut rng=SplitMix64::new(20261004);
        for d in [9,10,11,12,14,15,16,17,18] {
            let frames=dessin_frames(d,Ambient::Symmetric).unwrap();assert!(!frames.is_empty());
            for left in frames.iter().step_by((frames.len()/24).max(1)) {
                let signature=left.cycles.iter().map(Vec::len).collect::<Vec<_>>();
                for right in frames.iter().filter(|f|f.cycles.iter().map(Vec::len).collect::<Vec<_>>()==signature).take(6) {
                    let t=glue_tuple(d,left,right,&mut rng);assert!(tuple_product_is_identity(d,&t));
                    let (_,right_y)=right.pair.unwrap();assert_eq!(cycle_type(t[2],d),cycle_type(right_y,d));
                    assert_eq!(defect_of(&t,d).1,0);assert_eq!(defect_of(&t,d).2,0);
                    assert_eq!(genus_of(&t,d),defect_of(&t,d).0);
                }
            }
        }
    }
    #[test]
    fn fano_actions_and_longer_periods_are_exact() {
        let models=fano_models();assert_eq!(models.iter().map(|m|m.0).collect::<Vec<_>>(),vec![1,7,7,8]);
        for &(d,x,y) in &models[1..] {
            assert!(generators_transitive(d,&[x,y]));assert_eq!(subgroup_order(d,&[x,y]),168);
            assert!(power_identity(d,x,2) && power_identity(d,y,3) && power_identity(d,x.compose(y,d),7));
            let comm=x.compose(y,d).compose(x.inverse(d),d).compose(y.inverse(d),d);
            assert!(power_identity(d,comm,4));
        }
        // The two seven-point actions are inequivalent as marked actions.
        let (d,x,y)=models[1];let (_,xx,yy)=models[2];let twists=PairTwists::new(d*2,
            {let mut p=Perm::identity(14);for i in 0..7 {p.img[i]=x.img[i];p.img[7+i]=7+xx.img[i];}p},
            {let mut p=Perm::identity(14);for i in 0..7 {p.img[i]=y.img[i];p.img[7+i]=7+yy.img[i];}p});
        assert!(twists.maps[0].iter().all(|m|m.0==0));
        let mut rng=SplitMix64::new(20261003);
        for d in [9,10,11,14,15,16,17] {
            let frames=pure_frames(d,Ambient::Alternating,"PSL27").unwrap();assert!(!frames.is_empty());
            for frame in frames.iter().take(16) {let (a,b)=frame.pair.unwrap();let cache=PairTwists::new(d,a,b);
                for r in 1..=3 {for _ in 0..4 {
                    let h=power(d,a.compose(b,d),r).compose(cache.sample(d,&mut rng),d);
                    let c=centralizer(d,&permutation_cycles(d,h),&mut rng);let t=[a,b,c,product_one_fourth(d,a,b,c)];
                    let mut moved=t;for _ in 0..2*r {moved=hurwitz_move_dynamic(d,moved,0);}
                    assert_eq!(t.map(|p|p.conjugate_by(h,d)),moved);assert!(tuple_product_is_identity(d,&t));
                }}
            }
        }
    }
    #[test]
    fn full_and_even_centralizers_have_exact_support() {
        let d=5;let mut h=Perm::identity(d);
        loop {
            for ambient in [Ambient::Symmetric,Ambient::Alternating] {
                let generators=centralizer_generators(d,&permutation_cycles(d,h),ambient);
                let mut closure=HashSet::from([Perm::identity(d)]);let mut queue=vec![Perm::identity(d)];let mut head=0;
                while head<queue.len() {let u=queue[head];head+=1;for &g in &generators {let v=u.compose(g,d);if closure.insert(v) {queue.push(v);}}}
                let mut literal=HashSet::new();let mut p=Perm::identity(d);
                loop {if p.compose(h,d)==h.compose(p,d) && (ambient==Ambient::Symmetric || p.is_even(d)) {literal.insert(p);}
                    if !next_permutation(&mut p.img[..d]) {break;}}
                assert_eq!(closure,literal);assert_eq!(subgroup_order(d,&generators),literal.len() as u128);
            }
            if !next_permutation(&mut h.img[..d]) {break;}
        }
        for name in ["CERTIFIED_S18_XXXX_genus1_weak.txt","CERTIFIED_A25_mixed_genus0_weak.txt"] {
            let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();
            for key in case.states.unwrap() {let t=unpack(key,case.d);let moved=hurwitz_move_dynamic(case.d,t,0);
                if canonical(&moved,case.d)==key {let h=tuple_conjugator(case.d,t,moved).unwrap();
                    assert_eq!(envelope_order(case.d,case.ambient,t[0],t[1],&permutation_cycles(case.d,h)),case.ambient.target_order(case.d));}
            }
        }
    }
    #[test]
    fn pure_symmetry_models_twists_and_identity() {
        let mut rng=SplitMix64::new(20261002);
        for model in ["A4","S4","A5"] {
            let order=match model {"A4"=>12,"S4"=>24,_=>60};
            for (d,a,b) in pure_models(model).unwrap() {
                assert!(generators_transitive(d,&[a,b]));
                assert_eq!(order%subgroup_order(d,&[a,b]),0);
                assert!(power_identity(d,a,2) && power_identity(d,b,3));
                assert!(power_identity(d,a.compose(b,d),match model {"A4"=>3,"S4"=>4,_=>5}));
            }
            for d in [9,10,11,12,18,25] {
                let frames=pure_frames(d,Ambient::Symmetric,model).unwrap();assert!(!frames.is_empty());
                for frame in frames.iter().step_by((frames.len()/32).max(1)) {
                    let (a,b)=frame.pair.unwrap();let twists=PairTwists::new(d,a,b);
                    for _ in 0..4 {
                        let k=twists.sample(d,&mut rng);for p in [a,b] {assert_eq!(p.compose(k,d),k.compose(p,d));}
                        let h=frame.h.compose(k,d);let c=centralizer(d,&permutation_cycles(d,h),&mut rng);
                        let t=[a,b,c,product_one_fourth(d,a,b,c)];
                        assert_eq!(t.map(|p|p.conjugate_by(h,d)),hurwitz_move_dynamic(d,hurwitz_move_dynamic(d,t,0),0));
                        assert!(tuple_product_is_identity(d,&t));
                        assert_eq!(frame.left_genus,defect_of(&t,d).1);
                        if generators_transitive(d,&t[..3]) {assert!(centralizer_can_connect(d,a,b,&permutation_cycles(d,h)));}
                    }
                }
            }
            for d in [9,10,11,18,25,32,34,12,13,14,15,16,17,20,21,24,26,27,28,30] {
                assert!(!pure_frames(d,Ambient::Alternating,model).unwrap().is_empty(),"empty {model} A{d} job");
                assert!(!pure_frames(d,Ambient::Symmetric,model).unwrap().is_empty(),"empty {model} S{d} job");
            }
        }
        // Exhaustively compare the cached support with the literal centralizer
        // in S5. Membership is determined separately on each source orbit.
        let d=5;
        for (a,b) in [(perm_of_type(&[2,2,1],d,None),perm_of_type(&[3,1,1],d,None)),
            (perm_of_type(&[2,2,1],d,None),Perm::identity(d))] {
            let twists=PairTwists::new(d,a,b);let mut k=Perm::identity(d);
            loop {
                let cached=twists.orbits.iter().zip(&twists.maps).all(|(points,maps)|maps.iter().any(|(_,p)|points.iter().all(|&x|p.img[x]==k.img[x])));
                assert_eq!(cached,[a,b].iter().all(|p|p.compose(k,d)==k.compose(*p,d)));
                if !next_permutation(&mut k.img[..d]) {break;}
            }
        }
        for name in ["CERTIFIED_S18_XXXX_genus1_weak.txt","CERTIFIED_A25_mixed_genus0_weak.txt"] {
            let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();let d=case.d;
            for key in case.states.unwrap() {let t=unpack(key,d);let moved=hurwitz_move_dynamic(d,hurwitz_move_dynamic(d,t,0),0);
                if canonical(&moved,d)!=key {continue;}
                let h=tuple_conjugator(d,t,moved).unwrap();let k=t[0].compose(t[1],d).inverse(d).compose(h,d);
                let twists=PairTwists::new(d,t[0],t[1]);
                assert!(twists.orbits.iter().zip(&twists.maps).all(|(points,maps)|maps.iter().any(|(_,p)|points.iter().all(|&x|p.img[x]==k.img[x]))));
                assert!(centralizer_can_connect(d,t[0],t[1],&permutation_cycles(d,h)));
            }
        }
    }

    #[test]
    #[ignore = "technical microbenchmark; run explicitly on one worker"]
    fn benchmark_cached_pair_twists() {
        let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_XXXX_genus1_weak.txt")).unwrap();let d=case.d;
        let t=case.states.unwrap().iter().map(|&key|unpack(key,d)).find(|&t|canonical(&hurwitz_move_dynamic(d,t,0),d)==canonical(&t,d)).unwrap();
        let base=tuple_conjugator(d,t,hurwitz_move_dynamic(d,t,0)).unwrap();let cache=PairTwists::new(d,t[0],t[1]);
        for repeat in 0..6 {
            for cached in if repeat%2==0 {[false,true]} else {[true,false]} {
                let mut rng=SplitMix64::new(701+repeat);let start=Instant::now();
                for _ in 0..2000 {
                    let h=if cached {base.compose(cache.sample(d,&mut rng),d)} else {pair_braid_conjugator(d,t[0],t[1],&mut rng).unwrap()};
                    assert_eq!(t[..2].iter().map(|p|p.conjugate_by(h,d)).collect::<Vec<_>>(),vec![t[1].conjugate_by(t[0],d),t[0]]);
                    std::hint::black_box(h);
                }
                println!("TWIST_BENCH repeat={repeat} cached={cached} draws=2000 seconds={:.9}",start.elapsed().as_secs_f64());
            }
        }
    }
    #[test]
    fn square_roots_and_general_braid_identity() {
        let d=5;let mut q=Perm::identity(d);let mut rng=SplitMix64::new(91018);
        loop {
            let cycles=cycle_type(q,d);
            let exists=(2..=d).filter(|l|l%2==0).all(|l|cycles.iter().filter(|&&x|x==l).count()%2==0);
            for _ in 0..4 {let x=square_root(d,q,&mut rng);assert_eq!(x.is_some(),exists);if let Some(x)=x {assert_eq!(x.compose(x,d),q);}}
            if !next_permutation(&mut q.img[..d]) {break;}
        }
        for d in 4..=25 {for _ in 0..64 {
            let h=random_perm(d,Ambient::Symmetric,&mut rng);let cycles=permutation_cycles(d,h);
            let c=centralizer(d,&cycles,&mut rng);assert_eq!(h.compose(c,d),c.compose(h,d));
            let q=centralizer(d,&cycles,&mut rng);
            if let Some(x)=square_root(d,q,&mut rng) {
                let t=tuple_from_square(d,h,x,c);assert!(tuple_product_is_identity(d,&t));
                assert_eq!(t.map(|p|p.conjugate_by(h,d)),hurwitz_move_dynamic(d,t,0));
            }
        }}
        for name in ["CERTIFIED_S18_XXYY_genus1_weak.txt","CERTIFIED_S18_XXXX_genus1_weak.txt","CERTIFIED_A25_weak.txt","CERTIFIED_A25_mixed_genus0_weak.txt"] {
            let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();
            if name.contains("XXXX") {
                let closed=explore_goal(18,case.seed,1000,Instant::now()+Duration::from_secs(2),true);
                assert_eq!(closed.status,"certified_weak_and_strong_failure");
                assert_eq!(closed.states.iter().copied().collect::<HashSet<_>>(),case.states.as_ref().unwrap().iter().copied().collect());
                assert_eq!(closed.states.len(),10);assert_eq!(genus_of(&case.seed,18),1);
            }
            if name.contains("mixed") {
                let closed=explore_goal(25,case.seed,1000,Instant::now()+Duration::from_secs(2),true);
                assert_eq!(closed.status,"certified_weak_and_strong_failure");
                assert_eq!(closed.states.iter().copied().collect::<HashSet<_>>(),case.states.as_ref().unwrap().iter().copied().collect());
                assert_eq!(closed.states.len(),12);assert_eq!(genus_of(&case.seed,25),0);
                let passports=closed.states.iter().map(|&key|passport_string(&unpack(key,25),25)).collect::<HashSet<_>>();
                assert_eq!(passports.len(),12); // one class per ordering of the 3,3,2,6 entry types
                let mut pairs=BTreeMap::new();
                for &key in &closed.states {
                    let t=unpack(key,25);let mut partition=point_orbits(25,&t[..2]).iter().map(Vec::len).collect::<Vec<_>>();partition.sort_unstable();
                    *pairs.entry((subgroup_order(25,&t[..2]),partition)).or_insert(0usize)+=1;
                }
                assert_eq!(pairs,BTreeMap::from([((18,vec![1,1,2,6,6,9]),4),((576,vec![1,1,1,3,3,16]),2),
                    ((27993600,vec![10,15]),2),((111476736000,vec![5,20]),4)]));
                for &key in &closed.states {for mv in 0..3 {
                    let t=unpack(key,25);let moved=hurwitz_move_dynamic(25,t,mv);
                    assert_eq!(canonical(&hurwitz_move_dynamic(25,moved,mv),25),key);
                    let mut expected=t.map(|p|cycle_type(p,25));expected.swap(mv,mv+1);
                    assert_eq!(moved.map(|p|cycle_type(p,25)),expected);
                }}
            }
            let d=case.d;let mut fixed=0;
            for key in case.states.unwrap() {let t=unpack(key,d);let moved=hurwitz_move_dynamic(d,t,0);
                if canonical(&moved,d)!=key {continue;}
                let h=tuple_conjugator(d,t,moved).unwrap();let x=t[0].compose(h.inverse(d),d);let q=x.compose(x,d);
                assert_eq!(q.compose(h,d),h.compose(q,d));assert_eq!(tuple_from_square(d,h,x,t[2]),t);fixed+=1;
            }
            assert!(fixed>0);
        }
    }

    #[test]
    fn triangle_actions_and_mixed_frames_are_exact() {
        let models=triangle_models();let mut sizes=models.iter().map(|m|m.0).collect::<Vec<_>>();sizes.sort_unstable();
        assert_eq!(sizes,vec![1,2,3,4,6,6,6,8,12,12,24]);
        assert!(models.iter().any(|&(n,x,h)|n==6 && cycle_type(x.compose(h,n),n)==vec![1,1,4]));
        for &(n,x,h) in &models {
            assert!(power_identity(n,x,2) && power_identity(n,h,3) && power_identity(n,x.compose(h,n),4));
            assert!(generators_transitive(n,&[x,h]));assert_eq!(24%subgroup_order(n,&[x,h]),0);
        }
        let mut rng=SplitMix64::new(180910);
        for d in [9,10,11,12,18,24,32,34] {
            let frames=triangle_frames(d,Ambient::Alternating);assert!(!frames.is_empty());
            for frame in frames.iter().step_by((frames.len()/100).max(1)) {
                let (a,b)=frame.pair.unwrap();assert_eq!(power(d,a,4),Perm::identity(d));
                let c=centralizer(d,&frame.cycles,&mut rng);let t=[a,b,c,product_one_fourth(d,a,b,c)];
                assert_eq!(t.map(|p|p.conjugate_by(frame.h,d)),hurwitz_move_dynamic(d,t,0));
                let twists=PairTwists::new(d,a,b);
                let h=frame.h.compose(twists.sample(d,&mut rng),d);
                let c=centralizer(d,&permutation_cycles(d,h),&mut rng);let varied=[a,b,c,product_one_fourth(d,a,b,c)];
                assert_eq!(varied.map(|p|p.conjugate_by(h,d)),hurwitz_move_dynamic(d,varied,0));
                assert_eq!(frame.left_genus,defect_of(&t,d).1);assert!(frame.left_genus>=0);
                assert!(a.is_even(d) && !pair_transitive(d,a,b));
            }
        }
        let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_XXXX_genus1_weak.txt")).unwrap();
        for key in case.states.unwrap() {
            let t=unpack(key,18);let moved=hurwitz_move_dynamic(18,t,0);
            if canonical(&moved,18)!=key {continue;}
            let h=tuple_conjugator(18,t,moved).unwrap();let h0=power(18,t[0].compose(t[1],18),2);
            let x0=t[0].compose(h0.inverse(18),18);let twist=h0.inverse(18).compose(h,18);
            assert!(power_identity(18,x0,2) && power_identity(18,h0,3));
            assert_eq!(subgroup_order(18,&[x0,h0]),24);
            assert_eq!(point_orbits(18,&[x0,h0]).iter().map(Vec::len).collect::<Vec<_>>(),vec![6,6,6]);
            assert!(!twist.is_identity(18));assert!(power_identity(18,twist,3));
            for p in &t[..2] {assert_eq!(twist.compose(*p,18),p.compose(twist,18));}
            assert_eq!(subgroup_order(18,&[t[0],t[1],h]),72);
        }
    }
}
