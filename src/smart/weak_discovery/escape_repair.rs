//! Coordinated fixed-product repair scored by explicit positive braid witnesses.
use super::*;

pub(super) struct Repair {
    rounds:usize, beam:usize, extra:usize, relaxed:bool,
    words:Mutex<Vec<Vec<usize>>>, counts:[AtomicU64;6],
}

fn word_witness(d:usize,mut t:[Perm;4],word:&[usize])->bool {
    for &mv in word {t=hurwitz_move_dynamic(d,t,mv);}
    has_transitive_pair(d,&t)
}

fn parity_ok(d:usize,t:[Perm;4],ambient:Ambient)->bool {
    let all_even=t.iter().all(|p|p.is_even(d));
    if ambient==Ambient::Alternating {all_even} else {!all_even}
}

// Both pair products stay fixed. Default cuts preserve the exact partitions
// and index sums; relaxed cuts are a separately labeled new population.
fn stencil(d:usize,t:[Perm;4],ambient:Ambient,extra:usize,relaxed:bool,rng:&mut SplitMix64)->(Option<[Perm;4]>,usize) {
    let mut profiles=Vec::new();
    for (l,r) in [(2,2),(1,1),(1,2),(2,1)] {
        let parity=[t[0].is_even(d)^(l%2==1),t[1].is_even(d)^(l%2==1),
            t[2].is_even(d)^(r%2==1),t[3].is_even(d)^(r%2==1)];
        if (ambient==Ambient::Alternating && !parity.into_iter().all(|p|p))
            || (ambient==Ambient::Symmetric && parity.into_iter().all(|p|p)) {continue;}
        profiles.push((l+extra,r+extra));
    }
    if profiles.is_empty() {return (None,0);}
    let (l,r)=profiles[rng.below(profiles.len())];
    let (mut u,mut trials)=(t,0);
    for (start,n) in [(0,l),(2,r)] {
        if relaxed {
            for _ in 0..n {let i=rng.below(d);let mut j=rng.below(d-1);if j>=i {j+=1;}
                let (a,b)=pair_fiber::cut(d,u[start],u[start+1],i,j);u[start]=a;u[start+1]=b;trials+=1;
            }
        } else {
            let (pair,count,_)=pair_fiber::mutate_pair(d,u[start],u[start+1],n,rng);trials+=count;
            let Some((a,b))=pair else {return (None,trials);};u[start]=a;u[start+1]=b;
        }
    }
    if u[..2]==t[..2] || u[2..]==t[2..] || !parity_ok(d,u,ambient) {return (None,trials);}
    (Some(u),trials)
}

fn walk_witness(d:usize,mut t:[Perm;4],rng:&mut SplitMix64)->Option<Vec<usize>> {
    let mut word=Vec::new();let mut previous=None;
    for _ in 0..8 {
        let mut mv=rng.below(6);if previous.is_some_and(|old|mv==(old+3)%6) {mv=(mv+1)%6;}
        t=hurwitz_move_dynamic(d,t,mv);word.push(mv);previous=Some(mv);
        if has_transitive_pair(d,&t) {return Some(word);}
    }None
}

impl Repair {
    pub(super) fn new(args:&[String])->Result<Self,String> {
        let rounds:usize=parse_num(get_arg(args,"--coordinated-rounds",Some("8"))?,"--coordinated-rounds")?;
        let beam:usize=parse_num(get_arg(args,"--repair-beam",Some("8"))?,"--repair-beam")?;
        let extra:usize=parse_num(get_arg(args,"--coordinated-extra-cuts",Some("2"))?,"--coordinated-extra-cuts")?;
        if !(1..=128).contains(&rounds) || !(1..=64).contains(&beam) || extra>8 || extra%2!=0 {
            return Err("coordinated repair: rounds 1..128, beam 1..64, extra cuts 0,2,4,6,8".into());
        }
        let mut words=Vec::new();
        for first in [1,4] {words.push(vec![first]);for inner in [0,3] {for last in [1,4] {words.push(vec![first,inner,last]);}}}
        Ok(Self {rounds,beam,extra,relaxed:has_flag(args,"--allow-partition-change"),words:Mutex::new(words),
            counts:std::array::from_fn(|_|AtomicU64::new(0))})
    }
    fn score(&self,d:usize,t:[Perm;4],words:&[Vec<usize>])->usize {
        bridge_score(d,t)+words.iter().filter(|word|word_witness(d,t,word)).count()
    }
    pub(super) fn repair(&self,d:usize,t:[Perm;4],ambient:Ambient,force_change:bool,rng:&mut SplitMix64)->Option<[Perm;4]> {
        self.counts[0].fetch_add(1,Ordering::Relaxed);
        let mut words=self.words.lock().unwrap().clone();
        let (mut current,mut score,mut changed)=(t,self.score(d,t,&words),false);
        if score==0 && !force_change && parity_ok(d,t,ambient) {return Some(t);}
        for _ in 0..self.rounds {
            let mut best=None;
            for _ in 0..self.beam {
                let (candidate,trials)=stencil(d,current,ambient,self.extra,self.relaxed,rng);
                self.counts[1].fetch_add(trials as u64,Ordering::Relaxed);
                let Some(u)=candidate else {continue;};self.counts[2].fetch_add(1,Ordering::Relaxed);
                let next=self.score(d,u,&words);
                if best.as_ref().is_none_or(|(_,old)|next<*old) {best=Some((u,next));}
            }
            let Some((u,next))=best else {continue;};
            if next<=score || !changed || rng.below(16)==0 {current=u;score=next;changed=true;}
            if score==0 && changed {
                if let Some(word)=walk_witness(d,current,rng) {
                    // Store only words, never all sampled tuples. Each stored
                    // word is an independently checkable positive witness.
                    let mut bank=self.words.lock().unwrap();
                    if bank.len()<32 && !bank.contains(&word) {bank.push(word.clone());self.counts[3].fetch_add(1,Ordering::Relaxed);}
                    if !words.contains(&word) {words.push(word);}score=self.score(d,current,&words);continue;
                }
                self.counts[4].fetch_add(1,Ordering::Relaxed);return Some(current);
            }
        }
        self.counts[5].fetch_add(1,Ordering::Relaxed);None
    }
    pub(super) fn write(&self,out:&mut File)->io::Result<()> {
        writeln!(out,"coordinated repair rounds = {}\ncoordinated beam width = {}\ncoordinated extra cuts per changed pair = {}\ncoordinated partition-changing lane = {}\nrepair is biased; zero score is not a certificate",self.rounds,self.beam,self.extra,self.relaxed)?;
        for (label,count) in ["coordinated starts","coordinated cut trials","both-pair changed proposals","learned explicit escape words","repaired score-zero proposals","unrepaired starts"].iter().zip(&self.counts) {
            writeln!(out,"{label} = {}",count.load(Ordering::Relaxed))?;
        }Ok(())
    }
    pub(super) fn save_words(&self,path:&Path)->Result<(),String> {
        let mut out=fs::OpenOptions::new().create_new(true).write(true).open(path).map_err(|e|e.to_string())?;
        writeln!(out,"# move indices 0,1,2 positive; 3,4,5 inverse; positive witnesses only").map_err(|e|e.to_string())?;
        for word in self.words.lock().unwrap().iter() {writeln!(out,"{}",word.iter().map(usize::to_string).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;}
        out.sync_all().map_err(|e|e.to_string())
    }
}

#[cfg(test)] mod tests {
    use super::*;
    #[test] fn coordinated_stencils_preserve_products_partitions_and_genus() {
        let d=18;let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt")).unwrap();
        let t=unpack(*case.states.as_ref().unwrap().iter().next().unwrap(),d);let mut rng=SplitMix64::new(20261004);
        let mut made=0;
        for ambient in [Ambient::Alternating,Ambient::Symmetric] {for _ in 0..32 {
            if let (Some(u),_)=stencil(d,t,ambient,2,false,&mut rng) {
                made+=1;assert_eq!(u[0].compose(u[1],d),t[0].compose(t[1],d));assert_eq!(u[2].compose(u[3],d),t[2].compose(t[3],d));
                for i in [0,2] {assert_eq!(pair_fiber::pair_partition(d,u[i],u[i+1]),pair_fiber::pair_partition(d,t[i],t[i+1]));}
                assert_eq!(genus_of(&u,d),genus_of(&t,d));assert_eq!(defect_of(&u,d),defect_of(&t,d));assert!(parity_ok(d,u,ambient));
                assert_ne!(u[..2],t[..2]);assert_ne!(u[2..],t[2..]);
            }
        }}assert!(made>0);
    }
    #[test] fn identity_word_and_actual_braid_witness_are_distinguished() {
        let d=5;let mut a=Perm::identity(d);for i in 0..d {a.img[i]=((i+1)%d) as u8;}
        let t=[a,Perm::identity(d),Perm::identity(d),a.inverse(d)];assert!(word_witness(d,t,&[]));
        let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt")).unwrap();
        let t=unpack(*case.states.as_ref().unwrap().iter().next().unwrap(),18);
        assert!(!word_witness(18,t,&[1,3,1]));
    }
}
