//! Coupled pure-braid constraints, with moves confined to C(h).
//! A successful repair is a construction condition, NEVER a weak certificate.
use super::*;

pub(super) struct Repair {
    rounds: usize,
    source_only: bool,
    pub(super) admission_aware: bool,
    basis: bool,
    hidden: bool,
    counts: [AtomicU64; 6],
    funnel: [AtomicU64; 9],
}

// Test equivariance without canonicalizing two tuples or allocating per root.
// On a transitive tuple the image of 0 determines the entire relabeling.
// A partial score guides repair; only the full score is an exact equality test.
fn score(d: usize, source: [Perm; 4], target: [Perm; 4]) -> usize {
    pure_stability::equivariance_score(d,source,target)
}

fn middle_score(d: usize, t: [Perm; 4]) -> usize {
    let target = hurwitz_move_dynamic(d, hurwitz_move_dynamic(d, t, 1), 1);
    score(d, t, target)
}

impl Repair {
    pub(super) fn new(args: &[String]) -> Result<Self, String> {
        let rounds = parse_num(get_arg(args, "--joint-repair-rounds", Some("8"))?, "--joint-repair-rounds")?;
        if rounds > 64 { return Err("joint-repair-rounds must be 0..64".into()); }
        let hidden=has_flag(args,"--joint-hidden-fiber");
        let basis=hidden || has_flag(args,"--joint-pure-basis");
        Ok(Self { rounds, source_only:has_flag(args,"--joint-source-only"),
            admission_aware:basis || has_flag(args,"--joint-admission-aware"),basis,hidden,
            counts: std::array::from_fn(|_| AtomicU64::new(0)),
            funnel: std::array::from_fn(|_| AtomicU64::new(0)) })
    }

    fn evaluate(&self, d: usize, t: [Perm; 4]) -> usize {
        self.counts[1].fetch_add(1, Ordering::Relaxed);
        middle_score(d, t)
    }

    pub(super) fn repair(&self, d: usize, mut t: [Perm; 4], moves: &[Perm], rng: &mut SplitMix64) -> Option<[Perm; 4]> {
        self.counts[0].fetch_add(1, Ordering::Relaxed);
        // Nontransitive completions can be made transitive by a subsequent
        // centralizer move: do not reject them before the repair loop.
        let mut value = self.evaluate(d, t);
        if value == 4*d { self.counts[2].fetch_add(1, Ordering::Relaxed); return Some(t); }
        for _ in 0..self.rounds {
            if moves.is_empty() { break; }
            let start = rng.below(moves.len());
            let mut best = None;
            // Bounded best-neighbor search, plus inverses. Every move commutes
            // with h, so c stays in C(h); the end-square constraint survives.
            for j in 0..moves.len().min(12) {
                let g = moves[(start+j)%moves.len()];
                for g in [g, g.inverse(d)] {
                    let c = t[2].compose(g, d);
                    let u = [t[0], t[1], c, product_one_fourth(d,t[0],t[1],c)];
                    let next = self.evaluate(d,u);
                    if next == 4*d { self.counts[3].fetch_add(1, Ordering::Relaxed); return Some(u); }
                    if best.as_ref().is_none_or(|(v,_)|next>*v) { best=Some((next,u)); }
                }
            }
            let (next,u) = best?;
            if next >= value || rng.below(4)==0 { t=u; value=next; }
            else { self.counts[4].fetch_add(1, Ordering::Relaxed); }
        }
        self.counts[5].fetch_add(1, Ordering::Relaxed); None
    }

    // Rejected terminal matches remain traversable intermediates. The bounded
    // per-start cache only avoids rechecking a literally equivalent rejected
    // state; it never excludes a frame, passport or neighborhood.
    pub(super) fn repair_admitted<G,A>(&self,d:usize,mut t:[Perm;4],moves:&[Perm],rng:&mut SplitMix64,
        grade:G,mut admit:A)->Option<[Perm;4]>
        where G:Fn([Perm;4])->u8,A:FnMut([Perm;4],&mut SplitMix64)->Option<[Perm;4]> {
        self.counts[0].fetch_add(1,Ordering::Relaxed);
        let mut rejected=HashSet::new();
        let mut try_match=|u:[Perm;4],rng:&mut SplitMix64| {
            self.funnel[0].fetch_add(1,Ordering::Relaxed);
            let key=canonical(&u,d);
            if !rejected.insert(key) {self.funnel[2].fetch_add(1,Ordering::Relaxed);return None;}
            if self.basis {
                self.funnel[3].fetch_add(1,Ordering::Relaxed);
                if !pure_stability::fixes_basis(d,u) {self.funnel[4].fetch_add(1,Ordering::Relaxed);return None;}
            }
            let Some(u)=admit(u,rng) else {self.funnel[1].fetch_add(1,Ordering::Relaxed);return None;};
            if self.hidden {
                self.funnel[5].fetch_add(1,Ordering::Relaxed);
                if !pure_stability::has_hidden_fiber(d,u) {self.funnel[6].fetch_add(1,Ordering::Relaxed);return None;}
            }
            self.funnel[8].fetch_add(1,Ordering::Relaxed);Some(u)
        };
        let mut value=(grade(t),self.evaluate(d,t));
        if value.1==4*d {
            self.counts[2].fetch_add(1,Ordering::Relaxed);
            if let Some(u)=try_match(t,rng) {return Some(u);}
        }
        for _ in 0..self.rounds {
            if moves.is_empty() {break;}
            self.funnel[7].fetch_add(1,Ordering::Relaxed);
            let start=rng.below(moves.len());let mut best=None;
            for j in 0..moves.len().min(12) {
                let g=moves[(start+j)%moves.len()];
                for g in [g,g.inverse(d)] {
                    let c=t[2].compose(g,d);
                    let u=[t[0],t[1],c,product_one_fourth(d,t[0],t[1],c)];
                    let next=(grade(u),self.evaluate(d,u));
                    if next.1==4*d {
                        self.counts[3].fetch_add(1,Ordering::Relaxed);
                        if let Some(u)=try_match(u,rng) {return Some(u);}
                    }
                    if best.as_ref().is_none_or(|(v,_)|next>*v) {best=Some((next,u));}
                }
            }
            let (next,u)=best?;
            if next>=value || rng.below(4)==0 {t=u;value=next;}
            else {self.counts[4].fetch_add(1,Ordering::Relaxed);}
        }
        self.counts[5].fetch_add(1,Ordering::Relaxed);None
    }

    pub(super) fn write(&self, out: &mut File) -> io::Result<()> {
        writeln!(out,"joint pure constraint = end squares by construction; middle square by exact equivariance\njoint repair rounds = {}\njoint coverage = biased bounded centralizer repair; not a census or proof of pure stability over the component",self.rounds)?;
        writeln!(out,"joint certified source frames only = {}",self.source_only)?;
        writeln!(out,"joint admission aware = {}\njoint complete pure basis required = {}\njoint hidden outer passport fiber required = {}\njoint fitness = cheap admission grade then middle equivariance; rejected matches remain traversable",self.admission_aware,self.basis,self.hidden)?;
        for (i,label) in ["joint proposals","joint equivariance evaluations","joint direct matches","joint repaired matches","joint rejected downhill moves","joint unrepaired proposals"].iter().enumerate() {
            writeln!(out,"{label} = {}",self.counts[i].load(Ordering::Relaxed))?;
        }
        for (label,count) in ["joint gated matches","joint admission rejected matches","joint cached rejected matches",
            "joint complete basis checks","joint complete basis rejections","joint hidden fiber checks",
            "joint hidden fiber rejections","joint continued repair rounds","joint admitted matches"].iter().zip(&self.funnel) {
            writeln!(out,"{label} = {}",count.load(Ordering::Relaxed))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_score_agrees_with_canonical_states_and_certified_controls() {
        let root=Path::new(env!("CARGO_MANIFEST_DIR")).join("evidence/certificates");
        let mut checked=0;
        for entry in fs::read_dir(root).unwrap() {
            let path=entry.unwrap().path();
            if !path.file_name().unwrap().to_string_lossy().ends_with("_weak.txt") {continue;}
            let case=parse_certificate_case(&path).unwrap();
            for key in case.states.unwrap() {
                let t=unpack(key,case.d);
                for i in 0..3 {
                    let moved=hurwitz_move_dynamic(case.d,hurwitz_move_dynamic(case.d,t,i),i);
                    assert_eq!(score(case.d,t,moved)==4*case.d,canonical(&t,case.d)==canonical(&moved,case.d));
                }
                if case.d==25 {assert_eq!(middle_score(case.d,t),4*case.d);}
            }
            checked+=1;
        }
        assert_eq!(checked,9);
    }

    #[test]
    fn disconnected_agreement_is_not_a_full_match() {
        let t=[Perm::identity(9);4];
        assert!(score(9,t,t)<36);
    }

    #[test]
    fn rejected_exact_match_does_not_end_admission_aware_repair() {
        let path=Path::new(env!("CARGO_MANIFEST_DIR")).join("evidence/certificates/CERTIFIED_A25_XXYY_genus4_weak.txt");
        let case=parse_certificate_case(&path).unwrap();let t=case.seed;
        let h=tuple_conjugator(25,t,hurwitz_move_dynamic(25,hurwitz_move_dynamic(25,t,0),0)).unwrap();
        let moves=centralizer_generators(25,&permutation_cycles(25,h),Ambient::Symmetric);
        let repair=Repair::new(&["--joint-admission-aware".into(),"--joint-repair-rounds".into(),"2".into()]).unwrap();
        assert!(repair.repair_admitted(25,t,&moves,&mut SplitMix64::new(33),|_|0,|_,_|None).is_none());
        assert_eq!(repair.funnel[7].load(Ordering::Relaxed),2);
        assert!(repair.counts[1].load(Ordering::Relaxed)>1);
        let accepted=repair.repair_admitted(25,t,&moves,&mut SplitMix64::new(33),|_|0,|u,_|Some(u));
        assert_eq!(accepted,Some(t));
    }

    #[test]
    fn hidden_fiber_flag_implies_complete_basis_and_admission() {
        let repair=Repair::new(&["--joint-hidden-fiber".into()]).unwrap();
        assert!(repair.hidden && repair.basis && repair.admission_aware);
    }

    #[test]
    fn repair_preserves_product_one_and_end_constraint() {
        let path=Path::new(env!("CARGO_MANIFEST_DIR")).join("evidence/certificates/CERTIFIED_A25_XXYY_genus4_weak.txt");
        let case=parse_certificate_case(&path).unwrap();
        let t=case.seed;
        let p=t[0].compose(t[1],25);
        let target=hurwitz_move_dynamic(25,hurwitz_move_dynamic(25,t,0),0);
        let h=tuple_conjugator(25,t,target).unwrap();
        assert_eq!(t[0].conjugate_by(h,25),t[0].conjugate_by(p,25));
        let cycles=permutation_cycles(25,h);
        let moves=centralizer_generators(25,&cycles,Ambient::Symmetric);
        let repair=Repair::new(&[]).unwrap();
        let mut rng=SplitMix64::new(20261010);
        for g in moves.iter().take(6) {
            let c=t[2].compose(*g,25);
            let u=[t[0],t[1],c,product_one_fourth(25,t[0],t[1],c)];
            let end=hurwitz_move_dynamic(25,hurwitz_move_dynamic(25,u,0),0);
            assert_eq!(u.map(|x|x.conjugate_by(h,25)),end);
            if let Some(v)=repair.repair(25,u,&moves,&mut rng) {
                assert!(v.iter().fold(Perm::identity(25),|a,&b|a.compose(b,25)).is_identity(25));
                assert_eq!(middle_score(25,v),100);
                let end=hurwitz_move_dynamic(25,hurwitz_move_dynamic(25,v,0),0);
                assert_eq!(v.map(|x|x.conjugate_by(h,25)),end);
            }
        }
    }
}
