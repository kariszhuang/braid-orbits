//! Exact pure-braid stabilizers on simultaneous-S_d states. No new verifier.
use super::*;

// Literal execution order. These are A_12,A_23,A_34,A_13,A_24,A_14
// (using the inverse conjugation convention where necessary). The six
// conjugated squares generate P_4, not just its normal closure in B_4.
pub(super) const WORDS: [&[usize]; 6] = [
    &[0,0], &[1,1], &[2,2], &[4,0,0,1], &[5,1,1,2], &[5,4,0,0,1,2],
];

pub(super) fn apply(d:usize,mut t:[Perm;4],word:&[usize])->[Perm;4] {
    for &mv in word {t=hurwitz_move_dynamic(d,t,mv);}t
}

// Allocation-free root-image equivariance. Only full connected agreement is
// equality; partial agreement is a heuristic, never a negative certificate.
pub(super) fn equivariance_score(d:usize,source:[Perm;4],target:[Perm;4])->usize {
    let mut best=0;
    for root in 0..d {
        let mut image=[u8::MAX;MAX_D];let mut inverse=[u8::MAX;MAX_D];
        let mut queue=[0usize;MAX_D];image[0]=root as u8;inverse[root]=0;
        let (mut head,mut tail,mut matched,mut valid)=(0,1,0,true);
        while head<tail && valid {
            let x=queue[head];head+=1;
            for j in 0..4 {
                let xx=source[j].img[x] as usize;
                let yy=target[j].img[image[x] as usize] as usize;
                if (image[xx]!=u8::MAX && image[xx] as usize!=yy)
                    || (inverse[yy]!=u8::MAX && inverse[yy] as usize!=xx) {
                    valid=false;break;
                }
                if image[xx]==u8::MAX {image[xx]=yy as u8;inverse[yy]=xx as u8;
                    queue[tail]=xx;tail+=1;}
                matched+=1;
            }
        }
        best=best.max(matched);
        if valid && tail==d {return 4*d;}
    }best
}

pub(super) fn mask(d:usize,t:[Perm;4])->u8 {
    WORDS.iter().enumerate().fold(0,|mask,(i,w)| {
        mask | if equivariance_score(d,t,apply(d,t,w))==4*d {1<<i} else {0}
    })
}

pub(super) fn fixes_basis(d:usize,t:[Perm;4])->bool {
    WORDS.iter().all(|w|equivariance_score(d,t,apply(d,t,w))==4*d)
}

pub(super) fn passport_fibers(d:usize,states:&[Key])->BTreeMap<Passport,usize> {
    let mut fibers=BTreeMap::new();
    for &key in states {*fibers.entry(unpack(key,d).map(|g|cycle_type(g,d))).or_insert(0)+=1;}
    fibers
}

// Reuse the exact exploration engine. Capped/timed-out results do not qualify.
// Normality of P_4 bounds this orbit by 24 once its generators fix the seed.
pub(super) fn has_hidden_fiber(d:usize,t:[Perm;4])->bool {
    if !tuple_product_is_identity(d,&t) || !generators_transitive(d,&t[..3])
        || !fixes_basis(d,t) {return false;}
    let orbit=explore_goal(d,t,25,Instant::now()+Duration::from_secs(1),true);
    orbit.status=="certified_weak_and_strong_failure" && 24%orbit.states.len()==0
        && passport_fibers(d,&orbit.states).values().any(|&n|n>1)
}

pub(super) fn write_profile(file:&mut File,d:usize,states:&[Key])->io::Result<()> {
    let mut masks=BTreeMap::new();
    for &key in states {*masks.entry(mask(d,unpack(key,d))).or_insert(0usize)+=1;}
    let fibers=passport_fibers(d,states);
    let full=masks.get(&63).copied().unwrap_or(0)==states.len();
    if full && (states.is_empty() || 24%states.len()!=0) {
        return Err(io::Error::other("pure-braid orbit bound failed"));
    }
    writeln!(file,"\n# pure braid quotient = simultaneous S_d conjugation\n# pure generator order = A12,A23,A34,A13,A24,A14\n# pure generator mask histogram = {masks:?}\n# whole component factors through S4 = {full}\n# ordered passport fibers = {}\n# minimum passport fiber size = {}\n# maximum passport fiber size = {}\n# hidden outer passport fiber = {}",
        fibers.len(),fibers.values().min().copied().unwrap_or(0),
        fibers.values().max().copied().unwrap_or(0),fibers.values().any(|&n|n>1))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_nine_controls_match_literal_canonical_pure_words() {
        let root=Path::new(env!("CARGO_MANIFEST_DIR")).join("evidence/certificates");
        let mut cases=0;
        for entry in fs::read_dir(root).unwrap() {
            let path=entry.unwrap().path();
            if !path.file_name().unwrap().to_string_lossy().ends_with("_weak.txt") {continue;}
            let case=parse_certificate_case(&path).unwrap();let states=case.states.unwrap();
            for &key in &states {
                let t=unpack(key,case.d);
                for w in WORDS {assert_eq!(equivariance_score(case.d,t,apply(case.d,t,w))==4*case.d,
                    canonical(&apply(case.d,t,w),case.d)==key);}
                assert_eq!(fixes_basis(case.d,t),case.d==25);
            }
            let fibers=passport_fibers(case.d,&states);
            assert!(fibers.values().all(|&n|n==if case.d==25 {1} else {10}));
            if case.d==25 {assert!(!has_hidden_fiber(case.d,case.seed));}
            cases+=1;
        }assert_eq!(cases,9);
    }
    #[test]
    fn disconnected_match_cannot_certify_pure_stability() {
        let t=[Perm::identity(9);4];assert!(!fixes_basis(9,t));
        assert!(!has_hidden_fiber(9,t));
    }
}
