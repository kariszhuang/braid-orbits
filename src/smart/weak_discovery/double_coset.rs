//! Exact inequivalent boundary gluings; randomized fallback above the cap.
use super::*;

pub(super) struct Plan { pub representatives: Vec<Perm>, pub masses: Vec<usize>, pub centralizer_size: usize }
struct Cached { plan: Option<Arc<Plan>>, cursor: AtomicU64 }
pub(super) struct Cache { cap: usize, plans: Mutex<HashMap<(usize,usize),Arc<Cached>>>, counts: [AtomicU64;5], retire: bool }

fn matching(d:usize,left:&Frame,right:&Frame)->Perm {
    let mut q=Perm::identity(d);
    for (a,b) in left.cycles.iter().zip(&right.cycles) {
        assert_eq!(a.len(),b.len());for (&x,&y) in a.iter().zip(b) {q.img[y]=x as u8;}
    }q
}

fn root(parents:&mut [usize],mut i:usize)->usize {
    while parents[i]!=i {parents[i]=parents[parents[i]];i=parents[i];}i
}

pub(super) fn classify(d:usize,left:&Frame,right:&Frame,cap:usize)->Option<Plan> {
    let p=left.h; let q=matching(d,left,right); let id=Perm::identity(d);
    let (elements,complete)=PairTwists::new(d,p,id).enumerate(d,cap);
    if !complete {return None;}
    let lookup=elements.iter().enumerate().map(|(i,&g)|(g,i)).collect::<HashMap<_,_>>();
    let (a,b)=left.pair.unwrap(); let (x,y)=right.pair.unwrap();
    let l=PairTwists::new(d,a,b).generators(d);
    let r=PairTwists::new(d,x,y).generators(d).into_iter().map(|g|g.conjugate_by(q,d)).collect::<Vec<_>>();
    let mut parents=(0..elements.len()).collect::<Vec<_>>();
    for (i,&s) in elements.iter().enumerate() {
        for next in l.iter().map(|&g|g.compose(s,d)).chain(r.iter().map(|&g|s.compose(g,d))) {
            let j=lookup[&next]; let ii=root(&mut parents,i); let jj=root(&mut parents,j); parents[ii]=jj;
        }
    }
    let mut classes=BTreeMap::<usize,Vec<usize>>::new();
    for i in 0..elements.len() {classes.entry(root(&mut parents,i)).or_default().push(i);}
    let mut classes=classes.into_values().map(|v| {
        let rep=*v.iter().min_by_key(|&&i|elements[i].img).unwrap(); (elements[rep].compose(q,d),v.len())
    }).collect::<Vec<_>>(); classes.sort_unstable_by_key(|(g,_)|g.img);
    let (representatives,masses)=classes.into_iter().unzip();
    Some(Plan {representatives,masses,centralizer_size:elements.len()})
}

pub(super) fn tuple(d:usize,left:&Frame,right:&Frame,q:Perm)->[Perm;4] {
    let (a,b)=left.pair.unwrap(); let (x,y)=right.pair.unwrap();
    [a,b,y.inverse(d).conjugate_by(q,d),x.inverse(d).conjugate_by(q,d)]
}

impl Cache {
    pub(super) fn new(cap:usize,retire:bool)->Result<Self,String> {
        if !(1..=65536).contains(&cap) {return Err("boundary-centralizer-cap must be 1..65536".into());}
        Ok(Self {cap,plans:Mutex::new(HashMap::new()),counts:std::array::from_fn(|_|AtomicU64::new(0)),retire})
    }
    pub(super) fn sample(&self,d:usize,i:usize,j:usize,left:&Frame,right:&Frame,rng:&mut SplitMix64)->Option<[Perm;4]> {
        // Compute under the lock: workers never duplicate an expensive exact
        // classification. Cache only 128 pairs; uncached choices still survive.
        let entry={let mut cache=self.plans.lock().unwrap();
            if let Some(entry)=cache.get(&(i,j)) {Some(Arc::clone(entry))}
            else if cache.len()<128 {
                let plan=classify(d,left,right,self.cap).map(Arc::new);
                self.counts[usize::from(plan.is_none())].fetch_add(1,Ordering::Relaxed);
                let entry=Arc::new(Cached {plan,cursor:AtomicU64::new(0)});cache.insert((i,j),Arc::clone(&entry));Some(entry)
            } else {None}
        };
        if let Some(entry)=entry {if let Some(plan)=&entry.plan {
            let k=entry.cursor.fetch_add(1,Ordering::Relaxed) as usize;
            if self.retire && k>=plan.representatives.len() {
                self.counts[4].fetch_add(1,Ordering::Relaxed);return None;
            }
            self.counts[2].fetch_add(1,Ordering::Relaxed);
            if k<plan.representatives.len() {self.counts[3].fetch_add(1,Ordering::Relaxed);}
            return Some(tuple(d,left,right,plan.representatives[k%plan.representatives.len()]));
        }}
        Some(glue_tuple(d,left,right,rng))
    }
    pub(super) fn write(&self,out:&mut File)->io::Result<()> {
        writeln!(out,"boundary centralizer enumeration cap = {}\nboundary plan cache cap = 128\ndouble cosets use full S_d relabeling; randomized fallback is not an exhaustion",self.cap)?;
        writeln!(out,"retire unchanged exact boundary plans = {}",self.retire)?;
        for (label,count) in ["exact boundary-pair plans","capped boundary-pair fallback plans","double-coset representative occurrences","first-cycle distinct boundary classes visited","retired boundary-plan repeat draws avoided"].iter().zip(&self.counts) {
            writeln!(out,"{label} = {}",count.load(Ordering::Relaxed))?;
        }
        let cache=self.plans.lock().unwrap();
        let mass:usize=cache.values().filter_map(|x|x.plan.as_ref()).map(|p|p.masses.iter().sum::<usize>()).sum();
        assert!(cache.values().filter_map(|x|x.plan.as_ref()).all(|p|p.masses.iter().sum::<usize>()==p.centralizer_size));
        writeln!(out,"classified boundary matching mass = {mass}")
    }
    pub(super) fn save(&self,path:&Path,d:usize)->Result<(),String> {
        let mut out=OpenOptions::new().write(true).create_new(true).open(path).map_err(|e|e.to_string())?;
        writeln!(out,"left_frame\tright_frame\tstatus\tclass\tclass_mass\tcentralizer_size\tmatching_zero_based").map_err(|e|e.to_string())?;
        let cache=self.plans.lock().unwrap();let mut keys=cache.keys().copied().collect::<Vec<_>>();keys.sort_unstable();
        for (i,j) in keys {if let Some(plan)=&cache[&(i,j)].plan {
            for (k,(&q,&mass)) in plan.representatives.iter().zip(&plan.masses).enumerate() {
                writeln!(out,"{i}\t{j}\texact\t{k}\t{mass}\t{}\t{}",plan.centralizer_size,q.img[..d].iter().map(u8::to_string).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
            }
        } else {writeln!(out,"{i}\t{j}\tcapped_random_fallback\t\t\t\t").map_err(|e|e.to_string())?;}}
        out.sync_all().map_err(|e|e.to_string())
    }
}

#[cfg(test)] mod tests {
    use super::*;
    fn frame(d:usize,a:Perm,b:Perm)->Frame {let h=a.compose(b,d);Frame {h,cycles:permutation_cycles(d,h),pair:Some((a,b)),left_genus:0}}
    #[test] fn double_cosets_partition_the_whole_boundary_group() {
        for d in 3..=6 {let mut rng=SplitMix64::new(d as u64);
            for _ in 0..8 {let a=random_perm(d,Ambient::Symmetric,&mut rng);let b=random_perm(d,Ambient::Symmetric,&mut rng);
                let left=frame(d,a,b);let p=left.h;let x=random_perm(d,Ambient::Symmetric,&mut rng);let y=x.inverse(d).compose(p,d);let right=frame(d,x,y);
                let plan=classify(d,&left,&right,720).unwrap();
                assert_eq!(plan.masses.iter().sum::<usize>(),plan.centralizer_size);
                let (all,complete)=PairTwists::new(d,p,Perm::identity(d)).enumerate(d,720);assert!(complete);
                let q=matching(d,&left,&right);
                let l=PairTwists::new(d,a,b).enumerate(d,720).0;
                let r=PairTwists::new(d,x,y).enumerate(d,720).0;
                let mut covered=HashSet::new();
                for (&rep,&mass) in plan.representatives.iter().zip(&plan.masses) {
                    let s=rep.compose(q.inverse(d),d);
                    let orbit=l.iter().flat_map(|&ll|r.iter().map(move |&rr|ll.compose(s,d).compose(rr.conjugate_by(q,d),d))).collect::<HashSet<_>>();
                    assert_eq!(orbit.len(),mass);for g in orbit {assert!(covered.insert(g));}
                    let t=tuple(d,&left,&right,rep);assert!(t.into_iter().fold(Perm::identity(d),|p,g|p.compose(g,d)).is_identity(d));
                }
                assert_eq!(covered,all.into_iter().collect());
            }
        }
    }
    #[test] fn capped_classification_cannot_be_mistaken_for_complete() {
        let d=6;let f=frame(d,Perm::identity(d),Perm::identity(d));assert!(classify(d,&f,&f,10).is_none());
    }
    #[test] fn retirement_visits_every_exact_class_once_but_fallback_keeps_sampling() {
        let d=4;let a=perm_of_type(&[2,2],d,None);let b=Perm::identity(d);let f=frame(d,a,b);
        let mut rng=SplitMix64::new(3);let cache=Cache::new(24,true).unwrap();
        let count=classify(d,&f,&f,24).unwrap().representatives.len();
        let mut tuples=HashSet::new();for _ in 0..count {tuples.insert(cache.sample(d,0,0,&f,&f,&mut rng).unwrap());}
        assert_eq!(tuples.len(),count);assert!(cache.sample(d,0,0,&f,&f,&mut rng).is_none());
        let fallback=Cache::new(1,true).unwrap();for _ in 0..16 {assert!(fallback.sample(d,0,0,&f,&f,&mut rng).is_some());}
        let cyclic=Cache::new(24,false).unwrap();for _ in 0..count+3 {assert!(cyclic.sample(d,0,0,&f,&f,&mut rng).is_some());}
    }
}
