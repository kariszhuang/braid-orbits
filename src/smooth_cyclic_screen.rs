// Negative-only cyclic incidence sieve. A returned word is literally replayed;
// a missing collision proves nothing about smoothness or isotriviality.
// Kept separate from closure and from geometric sufficient criteria.
struct G3CyclicFrame {
    pair: [Perm; 2],
    letters: [u8; MAX_D],
    components: Vec<(usize, usize)>,
    // At start+delta store positions of a same-cycle pair, or the sentinel.
    differences: [[(u8, u8); MAX_D]; 2],
    steps: [[u64; MAX_D]; 2],
}
struct G3CyclicFailure { external:usize, side:usize, k:usize, len:usize }
impl G3CyclicFailure {
    // Materialize/replay only persisted proofs, not every discarded proposal.
    // The allocation-free plan is analogous to ForestFrame::deep_failure.
    fn word(&self,t:&[Perm;4],n:usize)->Vec<usize> {
        let Self {external,side,k,len}=*self;
        let mut word=Vec::with_capacity(2*k.min(len-k)+3);
        if external==3 {word.push(5);}
        if k<=len-k {word.extend(std::iter::repeat_n(0,2*k));}
        else {word.extend(std::iter::repeat_n(3,2*(len-k)));}
        if side==1 {word.push(3);}
        word.push(4);
        assert!(ff_replay_negative(*t,n,&word));
        word
    }
}
impl G3CyclicFrame {
    fn new(pair: [Perm; 2], n: usize) -> Self {
        let product = pair[0].compose(pair[1], n);
        let labels = pair.map(|p| dd_cycle_labels(p, n));
        let mut result = Self { pair, letters: [0; MAX_D], components: Vec::new(),
            differences: [[(u8::MAX, u8::MAX); MAX_D]; 2], steps:[[0;MAX_D];2] };
        let mut used = 0u64;
        let mut start = 0;
        for root in 0..n {
            if used & (1u64 << root) != 0 { continue; }
            let mut x = root;
            let mut len = 0;
            loop {
                result.letters[start + len] = x as u8;
                used |= 1u64 << x;
                len += 1;
                x = product.img[x] as usize;
                if x == root { break; }
            }
            for side in 0..2 { for x in 0..len { for y in x + 1..len {
                if labels[side][result.letters[start+x] as usize]
                    == labels[side][result.letters[start+y] as usize] {
                    result.differences[side][start+y-x] = (x as u8, y as u8);
                    result.differences[side][start+len-(y-x)] = (y as u8, x as u8);
                    result.steps[side][start] |= (1u64<<(y-x)) | (1u64<<(len-(y-x)));
                }
            } } }
            result.components.push((start, len));
            start += len;
        }
        debug_assert_eq!(start, n);
        result
    }

    fn failure(&self, t: &[Perm; 4], ix: [usize; 4], n: usize) -> Option<G3CyclicFailure> {
        for external in 2..4 {
            let other = 5-external;
            // Otherwise the complementary side might be additive: a cycle
            // incidence collision alone would not certify nonsmoothness.
            let forced = [ix[1]+ix[other] >= n, ix[0]+ix[other] >= n];
            if !forced[0] && !forced[1] { continue; }
            let labels = dd_cycle_labels(t[external], n);
            for &(start, len) in &self.components {
                let steps=(if forced[0]{self.steps[0][start]}else{0})
                    | (if forced[1]{self.steps[1][start]}else{0});
                if steps==0 {continue;}
                let mut masks=[0u64;MAX_D];
                for u in 0..len {masks[labels[self.letters[start+u] as usize] as usize]|=1u64<<u;}
                let full=(1u64<<len)-1;
                for u in 0..len {
                    let mask=masks[labels[self.letters[start+u] as usize] as usize];
                    if mask.is_power_of_two(){continue;}
                    let rotated=((mask>>u)|(mask<<(len-u)))&full;
                    let common=rotated&steps;
                    if common==0 {continue;}
                    let delta=common.trailing_zeros() as usize;
                    let side=usize::from(!forced[0] || self.steps[0][start]&(1u64<<delta)==0);
                    let (x,_) = self.differences[side][start+delta];
                    assert_ne!(x,u8::MAX);
                    // Q1^(2k) conjugates by p^k. Other component periods do
                    // not matter: one parallel incidence edge suffices.
                    return Some(G3CyclicFailure{external,side,k:(u+len-x as usize)%len,len});
                }
            }
        }
        None
    }
}

#[derive(Default)]
struct G3CyclicCache {
    frame: Option<G3CyclicFrame>,
    tested: u64,
    builds: u64,
    rejected: u64,
}
impl G3CyclicCache {
    fn failure(&mut self, t: &[Perm; 4], ix: [usize; 4], n: usize) -> Option<G3CyclicFailure> {
        // Two strict complementary-index tests, for either external entry.
        if (0..2).all(|side| (2..4).all(|e| ix[1-side]+ix[5-e] < n)) { return None; }
        self.tested += 1;
        let pair = [t[0], t[1]];
        if self.frame.as_ref().is_none_or(|f| f.pair != pair) {
            self.frame = Some(G3CyclicFrame::new(pair, n));
            self.builds += 1;
        }
        let word = self.frame.as_ref().unwrap().failure(t, ix, n);
        self.rejected += u64::from(word.is_some());
        word
    }
}

// Selected involution pairs need not occupy positions1/2. Inverse Hurwitz
// moves bring the chosen entries forward unchanged. This late-only extension
// also tests the full-product pairs in the degree14 triangle; it emits a
// literal witness, never a proof-only passport retirement.
fn g3_related_cyclic_word(t: &[Perm;4], ix: [usize;4], n: usize) -> Option<Vec<usize>> {
    let involution=std::array::from_fn::<_,4,_>(|j| (0..n).all(|x|t[j].img[t[j].img[x] as usize] as usize==x));
    let triangle=if n==14 {let mut ranks=ix;ranks.sort_unstable();ranks==[6,6,7,13]}else{false};
    if !triangle&&involution.iter().filter(|&&b|b).count()<2 {return None;}
    for i in 0..4 {for j in i+1..4 {
        if i==0&&j==1 {continue;}
        if !(involution[i]&&involution[j]||triangle&&ix[i]+ix[j]==n-1) {continue;}
        let mut labels=[0,1,2,3];let mut u=*t;let mut indices=ix;let mut prefix=Vec::new();
        for (destination,selected) in [i,j].into_iter().enumerate() {
            let mut at=labels.iter().position(|&x|x==selected).unwrap();
            while at>destination {
                let m=at-1+3;u=hurwitz_move_dynamic(n,u,m);prefix.push(m);
                labels.swap(at-1,at);indices.swap(at-1,at);at-=1;
            }
        }
        if !ff_smooth(&u,indices,n) {assert!(ff_replay_negative(*t,n,&prefix));return Some(prefix);}
        if let Some(word)=G3CyclicFrame::new([u[0],u[1]],n).failure(&u,indices,n) {
            prefix.extend(word.word(&u,n));assert!(ff_replay_negative(*t,n,&prefix));return Some(prefix);
        }
    }}
    None
}

// Exact C214 whole-cycle-subset cuts, evaluated in linear time by an
// incidence DSU. A repeated root detects any incidence loop, including loops
// involving several active cycles (which singleton difference masks miss).
// The finite prefix is a rejection sieve only, NOT an exhaustion of all k.
struct G3ForestCutFrame {
    pair: [Perm; 2],
    rotations: Vec<[Vec<u64>; 2]>,
}
fn g3_incidence_loop(supports: &[u64], labels: &[u8; MAX_D]) -> bool {
    let mut parent: [usize; MAX_D] = std::array::from_fn(|i| i);
    fn root(parent: &mut [usize; MAX_D], mut x: usize) -> usize {
        while parent[x] != x { parent[x] = parent[parent[x]]; x = parent[x]; }
        x
    }
    for &support in supports {
        let mut points = support;
        let first = points.trailing_zeros() as usize;
        points &= points-1;
        let anchor = labels[first] as usize;
        while points != 0 {
            let x = points.trailing_zeros() as usize;
            points &= points-1;
            let a = root(&mut parent, anchor);
            let b = root(&mut parent, labels[x] as usize);
            if a == b { return true; }
            parent[b] = a;
        }
    }
    false
}
impl G3ForestCutFrame {
    fn new(pair: [Perm; 2], n: usize, max_steps: usize) -> Self {
        assert!(n<=63 && (1..=64).contains(&max_steps));
        let initial = pair.map(|p| {
            let labels = dd_cycle_labels(p,n);
            let mut masks = [0u64; MAX_D];
            for x in 0..n { masks[labels[x] as usize] |= 1u64<<x; }
            masks.into_iter().filter(|m|m.count_ones()>1).collect::<Vec<_>>()
        });
        let p=pair[0].compose(pair[1],n);
        let mut rotations=vec![initial.clone()];
        while rotations.len()<max_steps {
            let next=rotations.last().unwrap().each_ref().map(|cycles| cycles.iter().map(|&mask| {
                let mut points=mask;let mut image=0u64;
                while points!=0 {let x=points.trailing_zeros()as usize;points&=points-1;image|=1u64<<p.img[x];}
                image
            }).collect::<Vec<_>>());
            // All product cycles rotate synchronously. Comparing the whole
            // support collection avoids inventing independent rotations.
            if next==initial {break;}
            rotations.push(next);
        }
        Self{pair,rotations}
    }
    fn failure(&self,t:&[Perm;4],ix:[usize;4],n:usize)->Option<Vec<usize>> {
        for external in 2..4 {
            let forced=[ix[1]+ix[5-external]>=n,ix[0]+ix[5-external]>=n];
            if !forced[0]&&!forced[1] {continue;}
            let labels=dd_cycle_labels(t[external],n);
            for (k,cycles) in self.rotations.iter().enumerate() {for side in 0..2 {
                if !forced[side]||!g3_incidence_loop(&cycles[side],&labels) {continue;}
                let mut word=Vec::with_capacity(2*k+3);
                if external==3 {word.push(5);}
                word.extend(std::iter::repeat_n(0,2*k));
                if side==1 {word.push(3);}
                word.push(4);
                assert!(ff_replay_negative(*t,n,&word));
                return Some(word);
            }}
        }
        None
    }
}
#[derive(Default)]
struct G3ForestCutCache { frame:Option<G3ForestCutFrame>,tested:u64,builds:u64,rejected:u64 }
impl G3ForestCutCache {
    fn failure(&mut self,t:&[Perm;4],ix:[usize;4],n:usize,steps:usize)->Option<Vec<usize>> {
        if (0..2).all(|s|(2..4).all(|e|ix[1-s]+ix[5-e]<n)) {return None;}
        self.tested+=1;
        let pair=[t[0],t[1]];
        if self.frame.as_ref().is_none_or(|f|f.pair!=pair) {
            self.frame=Some(G3ForestCutFrame::new(pair,n,steps));self.builds+=1;
        }
        let word=self.frame.as_ref().unwrap().failure(t,ix,n);
        self.rejected+=u64::from(word.is_some());word
    }
}

// First open three-involution boundary in genus three. This is a finite
// constructor, not a claim that every degree-16 passport has been exhausted.
// Two colored path frames, third entry a transversal partial matching.
const G3_BOUNDARY_TOTAL: u64 = (64+36)*5040;
struct G3InvolutionBoundary {
    frames: [ForestFrame; 2],
    cursor: u64,
}
impl G3InvolutionBoundary {
    fn new(n: usize, cursor: u64) -> Result<Self, String> {
        if n != 16 || cursor > G3_BOUNDARY_TOTAL {
            return Err("involution-boundary requires degree16 and cursor0..504000".into());
        }
        let frames = [8,7].map(|split| {
            let mut pair = [Perm::identity(n); 2];
            for (path, (start, len)) in [(0,split),(split,n-split)].into_iter().enumerate() {
                for x in 0..len-1 {
                    let side = (x + usize::from(split == 8 && path == 1))%2;
                    pair[side].img.swap(start+x, start+x+1);
                }
            }
            assert_eq!(pair.map(|p| sm_index(p,n)), [7,7]);
            assert_eq!(sm_index(pair[0].compose(pair[1],n),n),14);
            ForestFrame::new(pair, n)
        });
        Ok(Self { frames, cursor })
    }
    fn draw(&mut self) -> Option<[Perm; 4]> {
        if self.cursor == G3_BOUNDARY_TOTAL { return None; }
        let rank = self.cursor;
        self.cursor += 1;
        let shape = usize::from(rank >= 64*5040);
        let local = if shape == 0 { rank } else { rank-64*5040 };
        let omissions = local as usize/5040;
        let mut left = [0usize;7];
        let mut right = [0usize;7];
        if shape == 0 {
            let a = omissions/8;
            let b = omissions%8;
            for (j,x) in (0..8).filter(|&x| x!=a).enumerate() { left[j]=x; }
            for (j,x) in (0..8).filter(|&x| x!=b).enumerate() { right[j]=8+x; }
        } else {
            left = std::array::from_fn(|j| j);
            let mut k = 0;
            let mut omit = (0,1);
            'pairs: for a in 0..9 { for b in a+1..9 {
                if k == omissions { omit=(a,b); break 'pairs; }
                k += 1;
            } }
            for (j,x) in (0..9).filter(|&x| x!=omit.0 && x!=omit.1).enumerate() { right[j]=7+x; }
        }
        // Factoradic unranking: no rejection, RNG, or 7! permutation cache.
        let mut r = local as usize%5040;
        let mut c = Perm::identity(16);
        for j in 0..7 {
            let factorial = [720,120,24,6,2,1,1][j];
            let k = j+r/factorial;
            r %= factorial;
            let y=right[k];
            right[j..=k].rotate_right(1);
            c.img[left[j]]=y as u8;
            c.img[y]=left[j] as u8;
        }
        let f=&self.frames[shape];
        if ff_product_index(f.p,c,16)!=15 { return None; }
        Some(f.tuple(c,16))
    }
    fn frame(&self) -> &ForestFrame {
        &self.frames[usize::from(self.cursor > 64*5040)]
    }
}

#[cfg(test)]
mod cyclic_screen_tests {
    use super::*;
    #[test]
    fn whole_subset_dsu_matches_all_small_literal_product_indices() {
        for n in 1..=4 {let all=sc_perms(n);for &(a,_) in &all {for &(c,_) in &all {
            let frame=G3ForestCutFrame::new([a,Perm::identity(n)],n,1);
            assert_eq!(g3_incidence_loop(&frame.rotations[0][0],&dd_cycle_labels(c,n)),
                sm_index(a.compose(c,n),n)!=sm_index(a,n)+sm_index(c,n));
        }}}
    }
    #[test]
    fn whole_subset_catches_longer_loops_without_singleton_collisions() {
        for n in [6,8] {
            let mut a=Perm::identity(n);let mut c=a;
            for (x,y) in if n==6 {vec![(0,1),(2,3)]}else{vec![(0,1),(2,3),(4,5)]} {a.img.swap(x,y);}
            for (x,y) in if n==6 {vec![(0,2),(1,3)]}else{vec![(1,2),(3,4),(5,0)]} {c.img.swap(x,y);}
            let f=G3ForestCutFrame::new([a,Perm::identity(n)],n,16);
            let labels=dd_cycle_labels(c,n);
            assert!(f.rotations[0][0].iter().all(|m|!g3_incidence_loop(&[*m],&labels)));
            assert!(g3_incidence_loop(&f.rotations[0][0],&labels));
        }
    }
    #[test]
    fn whole_subset_rotations_equal_literal_conjugation() {
        let mut rng=SplitMix64::new(215);
        for n in [4,9,12,25,MAX_D.min(63)] {for _ in 0..30 {
            let a=gd_shuffle(n,&mut rng);let b=gd_shuffle(n,&mut rng);
            let c=gd_shuffle(n,&mut rng);let p=a.compose(b,n);
            let f=G3ForestCutFrame::new([a,b],n,16);
            let mut sides=[a,b];let labels=dd_cycle_labels(c,n);
            for rotation in &f.rotations {for side in 0..2 {
                assert_eq!(g3_incidence_loop(&rotation[side],&labels),
                    sm_index(sides[side].compose(c,n),n)!=sm_index(sides[side],n)+sm_index(c,n));
            }for side in 0..2 {sides[side]=p.compose(sides[side],n).compose(p.inverse(n),n);}}
        }}
    }
    #[test]
    fn whole_subset_words_literally_replay_with_both_external_entries() {
        let mut rng=SplitMix64::new(2151);let mut found=0;
        for n in [4,9,12,25,MAX_D.min(63)] {for _ in 0..100 {
            let a=gd_shuffle(n,&mut rng);let b=gd_shuffle(n,&mut rng);let c=gd_shuffle(n,&mut rng);
            let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];
            let f=G3ForestCutFrame::new([a,b],n,16);
            if let Some(w)=f.failure(&t,t.map(|p|sm_index(p,n)),n) {
                assert!(ff_replay_negative(t,n,&w));found+=1;
            }
        }}assert!(found>100);
    }
    #[test]
    fn every_reported_collision_is_a_literal_braid_witness() {
        let mut rng=SplitMix64::new(213);
        let mut found=0;
        for n in [3,4,5,9,12,16,25,MAX_D.min(63)] {
            let mut cache=G3CyclicCache::default();
            for _ in 0..1000 {
                let a=gd_shuffle(n,&mut rng); let b=gd_shuffle(n,&mut rng);
                let c=gd_shuffle(n,&mut rng);
                let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];
                let ix=t.map(|p|sm_index(p,n));
                if let Some(failure)=cache.failure(&t,ix,n) {
                    let word=failure.word(&t,n);
                    assert!(ff_replay_negative(t,n,&word)); found+=1;
                    assert!(word.len()<=n+3);
                }
            }
        }
        assert!(found>1000);
    }
    #[test]
    fn all_small_genus_zero_full_cycle_controls_are_retained() {
        let n=4; let all=sc_perms(n); let mut tested=0;
        for &(a,_) in &all { for &(b,_) in &all {
            let f=G3CyclicFrame::new([a,b],n);
            for &(c,_) in &all {
                let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];
                let ix=t.map(|p|sm_index(p,n));
                if ix.iter().sum::<usize>()!=2*n-2 || !ix.contains(&(n-1))
                    || !generators_transitive(n,&t[..3]) { continue; }
                assert!(f.failure(&t,ix,n).is_none()); tested+=1;
            }
        } }
        assert!(tested>100);
    }
    #[test]
    fn dihedral_both_parities_and_multicycle_paths_have_witnesses() {
        for n in 5..=MAX_D.min(63) {
            let a=Perm { img:std::array::from_fn(|x| if x<n { ((n-x)%n) as u8 } else {x as u8}) };
            let b=Perm { img:std::array::from_fn(|x| if x<n { ((n-x-1)%n) as u8 } else {x as u8}) };
            let c=a.compose(b,n);
            let t=[a,b,c,a.compose(b,n).compose(c,n).inverse(n)];
            let ix=t.map(|p|sm_index(p,n));
            assert!(ix.iter().sum::<usize>()>2*n-2);
            let word=G3CyclicFrame::new([a,b],n).failure(&t,ix,n).unwrap().word(&t,n);
            assert!(ff_replay_negative(t,n,&word));
            let moved=hurwitz_move_dynamic(n,t,4);
            let word=g3_related_cyclic_word(&moved,moved.map(|p|sm_index(p,n)),n).unwrap();
            assert!(ff_replay_negative(moved,n,&word));
        }
        // In a multicycle product no common global rotation period is needed.
        let source=G3InvolutionBoundary::new(16,0).unwrap();
        let [a,b]=source.frames[0].pair;
        let mut c=Perm::identity(16);
        for x in 0..16 { c.img[x]=((x+1)%16) as u8; }
        let t=[a,b,c,a.compose(b,16).compose(c,16).inverse(16)];
        let word=G3CyclicFrame::new([a,b],16).failure(&t,t.map(|p|sm_index(p,16)),16).unwrap().word(&t,16);
        assert!(ff_replay_negative(t,16,&word));
    }
    #[test]
    fn boundary_cursor_reconstructs_frames_and_never_repeats_a_third_entry() {
        assert!(G3InvolutionBoundary::new(15,0).is_err());
        assert!(G3InvolutionBoundary::new(16,G3_BOUNDARY_TOTAL+1).is_err());
        let mut source=G3InvolutionBoundary::new(16,0).unwrap();
        let mut qualified=0;
        let mut thirds=[HashSet::new(),HashSet::new()];
        while source.cursor<G3_BOUNDARY_TOTAL {
            let before=source.cursor;
            let t=source.draw();
            assert_eq!(source.cursor,before+1);
            if let Some(t)=t {
                let ix=t.map(|p|sm_index(p,16));
                assert_eq!(ix,[7,7,7,15]);
                assert!(tuple_product_is_identity(16,&t));
                assert!(generators_transitive(16,&t[..3]));
                let shape=usize::from(before>=64*5040);
                assert!(thirds[shape].insert(t[2]));
                if qualified%1024==0 {
                    let mut resumed=G3InvolutionBoundary::new(16,before).unwrap();
                    assert_eq!(resumed.draw(),Some(t));
                    assert_eq!(resumed.frame().pair,[t[0],t[1]]);
                }
                qualified+=1;
            }
        }
        assert!(qualified>10000);
        assert!(source.draw().is_none());
    }
}
