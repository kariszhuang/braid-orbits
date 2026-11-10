// Target/theorem exclusions, NOT literal nonsmooth words or isotrivial labels.
// This passport predicate is scoped to identity-free source-genus-three draws.
// A hypothetical transitive all-smooth orbit permits reordering the special
// and designated classes; star smoothness then forces their product full.
#[derive(Clone,Copy,Debug,PartialEq,Eq,PartialOrd,Ord)]
struct G3SaturatedRule { method:u8, prime:u8, modulus:u8 }
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
struct G3SaturatedCut { special:u8, partner:u8, k:u8, rule:G3SaturatedRule, reason:u8, violator:u8 }
impl G3SaturatedRule {
    fn method(&self)->&'static str {match self.method{0=>"prime-three-tiling",1=>"consecutive-root-run",2=>"norm-one-four",3=>"prime-tile-valuation",4=>"four-tile-prime-periodicity",_=>unreachable!()}}
}
impl G3SaturatedCut {
    fn method(&self)->&'static str {self.rule.method()}
    fn reason(&self)->&'static str {match self.reason{0=>"length-bound",1=>"ordinary-cycle-congruence",2=>"ordinary-cycle-divisibility",_=>unreachable!()}}
    fn key(&self)->(u8,u8,u8,u8,u8){(self.rule.method,self.k,self.rule.prime,self.rule.modulus,self.reason)}
}
const G3_SATURATED_HEADER:&str="candidate\tstatus\tspecial\tpartner\tk\tprime\tmodulus\tmethod\treason\tviolator\tseed\n";
fn g3_sp_row(candidate:u64,cut:Option<G3SaturatedCut>,seed:&str)->String {
    if let Some(c)=cut {format!("{candidate}\ttheorem-excluded\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{seed}\n",
        c.special,c.partner,c.k,c.rule.prime,c.rule.modulus,c.method(),c.reason(),if c.violator==u8::MAX{String::new()}else{c.violator.to_string()})}
    else{format!("{candidate}\tnot-rejected\t\t\t\t\t\t\t\t\t{seed}\n")}
}
fn g3_sp_prime(n:usize)->bool {n>=2&&(2..).take_while(|p|p*p<=n).all(|p|n%p!=0)}
fn g3_sp_order(base:usize,prime:usize)->usize {
    let mut value=1;for order in 1..prime{value=value*base%prime;if value==1{return order;}}0
}
fn g3_sp_root_run(n:usize,k:usize,prime:usize)->bool {
    let mut orbit=1u64;let mut value=1;
    loop {orbit|=1u64<<value;value=value*prime%n;if value==1{break;}}
    let mut run=0;
    for i in 0..n+k-1 {if orbit&(1u64<<(i%n))!=0 {run+=1;if run>=k{return true;}}else{run=0;}}
    false
}
struct G3SaturatedPassport { n:usize,rules:Vec<Vec<G3SaturatedRule>>,active:bool }
impl G3SaturatedPassport {
    fn new(n:usize)->Self {
        let mut rules=vec![Vec::new();n];
        if g3_sp_prime(n)&&n<=63 {
            for k in 2..n {
                if k==3 {rules[k].push(G3SaturatedRule{method:0,prime:3,modulus:3});}
                else if g3_sp_prime(k){rules[k].push(G3SaturatedRule{method:3,prime:k as u8,modulus:k as u8});}
                if k==4 {rules[k].push(G3SaturatedRule{method:4,prime:2,modulus:4});}
                let mut remaining=k;
                for ell in 2..=k {if remaining%ell!=0{continue;}
                    let mut modulus=1;while remaining%ell==0{remaining/=ell;modulus*=ell;}
                    if g3_sp_root_run(n,k,ell)&&!rules[k].iter().any(|r|r.method!=3&&r.modulus as usize==modulus){
                        rules[k].push(G3SaturatedRule{method:1,prime:ell as u8,modulus:modulus as u8});
                    }
                }
            }
        }
        let active=rules.iter().any(|row|!row.is_empty());Self{n,rules,active}
    }
    fn uniform_modulo(p:Perm,n:usize,modulus:usize)->bool {
        let mut seen=0u64;let mut expected=None;
        for root in 0..n {if seen&(1u64<<root)!=0{continue;}
            let mut x=root;let mut length=0;
            loop {seen|=1u64<<x;length+=1;x=p.img[x]as usize;if x==root{break;}}
            let residue=length%modulus;
            if expected.is_some_and(|old|old!=residue){return false;}expected=Some(residue);
        }true
    }
    fn coprime_lengths(p:Perm,n:usize,prime:usize)->bool {
        let mut seen=0u64;
        for root in 0..n {if seen&(1u64<<root)!=0{continue;}
            let mut x=root;let mut length=0;
            loop {seen|=1u64<<x;length+=1;x=p.img[x]as usize;if x==root{break;}}
            if length%prime==0{return false;}
        }true
    }
    fn rejection(&self,t:&[Perm;4],indices:[usize;4])->Option<G3SaturatedCut> {
        let n=self.n;
        if !self.active||indices.contains(&0)||indices.iter().sum::<usize>()!=2*n+4{return None;}
        for special in 0..4 {
            let k=indices[special]+1;if k>=n||self.rules[k].is_empty(){continue;}
            if !(0..4).any(|partner|partner!=special&&indices[partner]==n-k){continue;}
            let mut ordinary=[0usize;3];let mut used=0;for j in 0..4{if j!=special{ordinary[used]=j;used+=1;}}
            if (0..3).any(|i|(i+1..3).any(|j|indices[ordinary[i]]+indices[ordinary[j]]<n)){continue;}
            if (0..n).filter(|&x|t[special].img[x]as usize!=x).count()!=k{continue;}
            for partner in 0..4 {if partner==special||indices[partner]!=n-k{continue;}
                for &rule in &self.rules[k] {
                    if rule.method!=3&&n/k<=rule.modulus as usize{return Some(G3SaturatedCut{special:special as u8,partner:partner as u8,k:k as u8,rule,reason:0,violator:u8::MAX});}
                    for violator in 0..4 {if violator==special||violator==partner{continue;}
                        let valid=if rule.method==3{Self::coprime_lengths(t[violator],n,rule.prime as usize)}else{Self::uniform_modulo(t[violator],n,rule.modulus as usize)};
                        if !valid{
                            return Some(G3SaturatedCut{special:special as u8,partner:partner as u8,k:k as u8,rule,reason:if rule.method==3{2}else{1},violator:violator as u8});
                        }
                    }
                }
            }
        }None
    }
}

#[cfg(test)]
mod saturated_passport_tests {
    use super::*;
    fn permutation(n:usize,shape:&[usize])->Perm {
        assert_eq!(shape.iter().sum::<usize>(),n);let mut p=Perm::identity(n);let mut start=0;
        for &length in shape{for x in 0..length{p.img[start+x]=(start+(x+1)%length)as u8;}start+=length;}p
    }
    fn types(n:usize,shapes:[&[usize];4])->[Perm;4]{shapes.map(|shape|permutation(n,shape))}
    #[test]fn arithmetic_rules_retain_prime_scope_and_distinct_strengths(){
        let a=G3SaturatedPassport::new(17);assert_eq!(g3_sp_order(2,17),8);
        assert!(a.rules[4].iter().any(|r|r.method==4&&r.modulus==4));
        let b=G3SaturatedPassport::new(23);assert_eq!(g3_sp_order(2,23),11);
        assert!(g3_sp_root_run(23,4,2));assert!(b.rules[4].iter().any(|r|r.method==4&&r.modulus==4));
        assert!(!G3SaturatedPassport::new(9).active);
        assert!(G3SaturatedPassport::new(31).rules[3].iter().any(|r|r.method==0));
        assert!(!g3_sp_root_run(31,4,2));assert_eq!(g3_sp_order(2,31),5);
        assert!(G3SaturatedPassport::new(31).rules[4].iter().any(|r|r.method==4));
    }
    #[test]fn designated_partner_is_never_subject_to_the_ordinary_congruence(){
        let n=23;let t=types(n,[&[4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],&[5,5,5,8],&[5,5,5,5,1,1,1],&[5,5,5,1,1,1,1,1,1,1,1]]);
        let ix=t.map(|p|sm_index(p,n));assert_eq!(ix,[3,19,16,12]);
        assert!(!G3SaturatedPassport::uniform_modulo(t[1],n,4));
        assert!(G3SaturatedPassport::new(n).rejection(&t,ix).is_none());
        let mut bad=t;bad[2]=permutation(n,&[2,2,2,2,5,5,5]);
        let cut=G3SaturatedPassport::new(n).rejection(&bad,bad.map(|p|sm_index(p,n))).unwrap();
        assert_eq!((cut.special,cut.partner,cut.violator,cut.reason),(0,1,2,1));
    }
    #[test]fn length_bound_star_and_genus_hypotheses_are_literal(){
        let n=17;let t=types(n,[&[4,1,1,1,1,1,1,1,1,1,1,1,1,1],&[4,4,4,5],&[1,1,1,1,1,1,11],&[1,1,1,1,13]]);
        let ix=t.map(|p|sm_index(p,n));assert_eq!(ix,[3,13,10,12]);
        let cut=G3SaturatedPassport::new(n).rejection(&t,ix).unwrap();assert_eq!(cut.reason,0);
        let mut wrong=ix;wrong[2]-=1;assert!(G3SaturatedPassport::new(n).rejection(&t,wrong).is_none());
        let multi=types(n,[&[2,2,2,1,1,1,1,1,1,1,1,1,1,1],&[4,4,4,5],&[1,1,1,1,1,1,11],&[1,1,1,1,13]]);
        assert!(G3SaturatedPassport::new(n).rejection(&multi,multi.map(|p|sm_index(p,n))).is_none());
        let outside=types(n,[&[9,1,1,1,1,1,1,1,1],&[9,1,1,1,1,1,1,1,1],&[17],&[7,1,1,1,1,1,1,1,1,1,1]]);
        assert_eq!(outside.map(|p|sm_index(p,n)),[8,8,16,6]);
        assert!(G3SaturatedPassport::new(n).rejection(&outside,outside.map(|p|sm_index(p,n))).is_none());
    }
    #[test]fn prime_five_valuation_does_not_invent_equal_unit_residues(){
        let n=31;let mut a=vec![5];a.extend([1;26]);let shape=[2,2,2,2,2,2,2,2,3,3,3,3,3];
        let t=types(n,[&a,&[6,6,6,6,7],&shape,&shape]);let ix=t.map(|p|sm_index(p,n));
        assert_eq!(ix,[4,26,18,18]);let policy=G3SaturatedPassport::new(n);
        assert_eq!(policy.rules[5],vec![G3SaturatedRule{method:3,prime:5,modulus:5}]);
        assert!(policy.rejection(&t,ix).is_none());
        let mut bad=t;bad[2]=permutation(n,&[1,1,1,1,1,1,1,1,2,3,4,5,9]);
        let cut=policy.rejection(&bad,bad.map(|p|sm_index(p,n))).unwrap();
        assert_eq!((cut.rule.method,cut.reason,cut.violator),(3,2,2));
    }
}
