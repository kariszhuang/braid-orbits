// C251: incomplete proposal-FRAME target exclusions, not tuple negatives.
// Enabled only by a separately proof-bound non-isotrivial-genus-three plan.
#[derive(Default)]
struct G3FrameRetirement {
    enabled:bool,
    checked:u64,
    counts:[u64;6],
    samples:Vec<([Perm;2],[usize;4],u8,usize)>,
}
impl G3FrameRetirement {
    fn reason_name(reason:usize)->&'static str {match reason {
        0=>"C243-degree17-single3-saturated-no-all-smooth",
        1=>"C247-degree19-single3-saturated-no-all-smooth",
        2=>"C244-degree23-single4-saturated-no-all-smooth",
        3=>"C241-degree14-single3-middle-profile-no-all-smooth",
        4=>"C239-k-above5-Hodge-negative-no-all-smooth",
        5=>"C239-k5-Hodge-zero-nonisotrivial-target-excluded",
        _=>unreachable!(),
    }}
    fn decision(pair:[Perm;2],ix:[usize;4],n:usize)->Option<(usize,usize)> {
        if ix.contains(&0)||ix.iter().sum::<usize>()!=2*n+4{return None;}
        if sm_index(pair[0],n)!=ix[0]||sm_index(pair[1],n)!=ix[1]{return None;}
        let k=ix[0]+1;
        if !(3..n).contains(&k){return None;}
        // index+1 moved letters means EXACTLY ONE active cycle, not 2^2/3x2.
        if (0..n).filter(|&x|pair[0].img[x]as usize!=x).count()!=k{return None;}
        if !(1..4).any(|j|ix[j]==n-k){return None;}
        if (1..4).any(|i|(i+1..4).any(|j|ix[i]+ix[j]<n)){return None;}
        let reason=if n==17&&k==3{0}
            else if n==19&&k==3{1}
            else if n==23&&k==4{2}
            else if n==14&&k==3 {let mut sorted=ix;sorted.sort_unstable();if sorted!=[2,9,10,11]{return None;}3}
            else if k>5{4}else if k==5{5}else{return None;};
        Some((reason,k))
    }
    fn retire(&mut self,pair:[Perm;2],ix:[usize;4],n:usize)->bool {
        if !self.enabled{return false;}
        self.checked+=1;
        let Some((reason,k))=Self::decision(pair,ix,n)else{return false;};
        self.counts[reason]+=1;
        if self.samples.len()<48{self.samples.push((pair,ix,reason as u8,k));}
        true
    }
    fn save(&self,out:&Path,n:usize,constructor:&str)->Result<(),String>{
        let excluded:u64=self.counts.iter().sum();
        let summary=format!("constructor\tenabled\tchecked_new_frames\ttarget_excluded_frames\tsampled_frames\n{constructor}\t{}\t{}\t{excluded}\t{}\n",self.enabled,self.checked,self.samples.len());
        fs::write(out.join("FRAME_RETIREMENT_SUMMARY.tsv"),summary).map_err(|e|e.to_string())?;
        let mut counts=String::from("reason\ttarget_excluded_frames\n");
        for(reason,count)in self.counts.iter().enumerate(){counts.push_str(&format!("{}\t{count}\n",Self::reason_name(reason)));}
        fs::write(out.join("FRAME_RETIREMENT_REASONS.tsv"),counts).map_err(|e|e.to_string())?;
        let mut rows=String::from("sample\tk\tprofile\treason\tpartial_pair_identity_padded_NOT_A_CANDIDATE\n");
        for(j,(pair,ix,reason,k))in self.samples.iter().enumerate(){
            rows.push_str(&format!("{j}\t{k}\t{}\t{}\t{}\n",ix.iter().map(|x|x.to_string()).collect::<Vec<_>>().join(","),Self::reason_name(*reason as usize),tuple_string_dynamic(n,&[pair[0],pair[1],Perm::identity(n),Perm::identity(n)])));
        }
        fs::write(out.join("FRAME_RETIREMENT_SAMPLES.tsv"),rows).map_err(|e|e.to_string())?;
        fs::write(out.join("FRAME_RETIREMENT_SCOPE.txt"),"Counterfactual incomplete proposal FRAME exclusions for the NON-ISOTRIVIAL source-genus-three target. The remaining two entries were NOT constructed. Counts are frames, not transitive tuples, literal negative witnesses, components, or isotrivial labels. Identity padding is a serialization aid, NOT a candidate quadruple. All caps/raw unresolved work remain separate.\n").map_err(|e|e.to_string())
    }
}

#[cfg(test)]
mod frame_retirement_tests {
    use super::*;
    fn cycle(n:usize,length:usize)->Perm {let mut a=Perm::identity(n);for x in 0..length{a.img[x]=((x+1)%length)as u8;}a}
    fn pair(n:usize,k:usize,second:usize)->[Perm;2]{[cycle(n,k),cycle(n,second+1)]}
    #[test]fn finite_proven_single_cycles_retire_only_exact_scope(){
        for(n,k,ix,reason)in[(17,3,[2,9,13,14],0),(19,3,[2,12,12,16],1),(23,4,[3,12,16,19],2),(14,3,[2,9,10,11],3)]{
            assert_eq!(G3FrameRetirement::decision(pair(n,k,ix[1]),ix,n),Some((reason,k)));
        }
        assert!(G3FrameRetirement::decision(pair(14,3,10),[2,10,10,10],14).is_none());
        assert!(G3FrameRetirement::decision(pair(17,3,10),[2,10,13,13],17).is_none());
    }
    #[test]fn double_transpositions_and_other_multicycles_are_preserved(){
        let mut a=Perm::identity(17);a.img.swap(0,1);a.img.swap(2,3);
        assert!(G3FrameRetirement::decision([a,cycle(17,10)],[2,9,13,14],17).is_none());
        let mut b=Perm::identity(23);b.img[0]=1;b.img[1]=2;b.img[2]=0;b.img.swap(3,4);
        assert!(G3FrameRetirement::decision([b,cycle(23,13)],[3,12,16,19],23).is_none());
    }
    #[test]fn hodge_zero_is_target_exclusion_not_nonsmooth_classification(){
        assert_eq!(G3FrameRetirement::decision(pair(15,5,10),[4,10,10,10],15),Some((5,5)));
        // k6 saturated sourceg3, with all ordinary-pair sums beyond d-1.
        assert_eq!(G3FrameRetirement::decision(pair(20,6,12),[5,12,13,14],20),Some((4,6)));
        assert!(G3FrameRetirement::decision(pair(20,6,13),[5,13,13,13],20).is_none());
    }
    #[test]fn default_OFF_never_changes_sampling_or_counts(){
        let mut policy=G3FrameRetirement::default();let ix=[2,9,13,14];
        assert!(!policy.retire(pair(17,3,9),ix,17));assert_eq!(policy.checked,0);assert_eq!(policy.counts,[0;6]);
        policy.enabled=true;assert!(policy.retire(pair(17,3,9),ix,17));assert_eq!(policy.checked,1);assert_eq!(policy.counts[0],1);
    }
}
