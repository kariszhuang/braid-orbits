// Histograms stay in RAM. All output is buffered stdout, once per cell/block.
fn cp_cycles(p: Perm, d: usize) -> usize {
    let mut seen=0u64; let mut count=0;
    for i in 0..d {if seen&(1u64<<i)==0 {
        count+=1;let mut x=i;
        while seen&(1u64<<x)==0 {seen|=1u64<<x;x=p.img[x] as usize;}
    }}count
}
fn cp_stats(t:&[Perm;4],d:usize)->(usize,usize,usize,usize,usize) {
    let total:usize=t.iter().map(|p|d-cp_cycles(*p,d)).sum();
    assert_eq!(total%2,0);
    let g=total/2+1-d;
    let cycles=cp_cycles(t[0].compose(t[1],d),d);
    let l=point_orbit_count_perms(d,t[0],t[1]);
    let r=point_orbit_count_perms(d,t[2],t[3]);
    assert!(cycles+1>=l+r);let delta=cycles+1-l-r;assert!(delta<=g);
    if l==1&&r==1 {assert_eq!(delta,cycles-1);}
    (g,delta,cycles,l,r)
}
fn cp_tuple(raw:&str,d:usize)->Result<[Perm;4],String> {
    let values:Vec<usize>=raw.split(',').map(|v|v.parse().map_err(|_|"bad point".to_string())).collect::<Result<_,_>>()?;
    if values.len()!=4*d {return Err("expected4*d points".into());}
    let mut t=[Perm::identity(d);4];
    for j in 0..4 {let mut seen=0u64;for x in 0..d {
        let v=values[j*d+x];if v>=d || seen&(1u64<<v)!=0 {return Err("not a permutation".into());}
        seen|=1u64<<v;t[j].img[x]=v as u8;
    }}
    if !tuple_product_is_identity(d,&t)||!generators_transitive(d,&t[..3]) {return Err("product/transitivity".into());}
    Ok(t)
}
fn cp_row(out:&mut impl Write,kind:&str,id:&str,d:usize,a:&str,rep:usize,block:usize,g:usize,h:&[u64],both:u64)->Result<(),String> {
    let count:u64=h.iter().sum(); if count==0 {return Ok(());}
    for (k,&n) in h.iter().enumerate() {if n!=0 {
        writeln!(out,"{kind}\t{id}\t{d}\t{a}\t{rep}\t{block}\t{g}\t{k}\t{n}\t{count}\t{both}").map_err(|e|e.to_string())?;
    }}Ok(())
}
pub fn coin_probe_main(args:&[String])->Result<(),String> {
    let arg=|name:&str,default:&str|args.windows(2).find(|w|w[0]==name).map(|w|w[1].clone()).unwrap_or(default.into());
    let trials:usize=arg("--trials","100000").parse().map_err(|_|"trials")?;
    let seed:u64=arg("--seed","2026101202").parse().map_err(|_|"seed")?;
    let input=arg("--seeds","");let mut out=BufWriter::with_capacity(1<<20,io::stdout().lock());
    writeln!(out,"kind\tid\tdegree\tambient\treplicate\tblock\tgenus\tdefect\tcount\ttotal\tboth_pairs_transitive").map_err(|e|e.to_string())?;
    for d in [8usize,9,12,18,25] {for ambient in [Ambient::Symmetric,Ambient::Alternating] {for rep in 0..3 {
        let mut rng=SplitMix64::new(seed+(d as u64)*100+(rep as u64)*13+if ambient==Ambient::Alternating {10000}else{0});
        let mut hist=[[0u64;MAX_D+1];MAX_D+1];let mut both=[0u64;MAX_D+1];
        for _ in 0..trials {
            let a=random_perm(d,ambient,&mut rng);let b=random_perm(d,ambient,&mut rng);let c=random_perm(d,ambient,&mut rng);
            let t=[a,b,c,product_one_fourth(d,a,b,c)];
            if t.iter().any(|p|p.is_identity(d))||!generators_transitive(d,&t[..3]) {continue;}
            let(g,delta,_,l,r)=cp_stats(&t,d);if g==0 {continue;}
            hist[g][delta]+=1;if l==1&&r==1 {both[g]+=1;}
        }
        for g in 1..=d-1 {cp_row(&mut out,"iid","uniform",d,ambient.short(),rep,0,g,&hist[g],both[g])?;}
    }}}
    if !input.is_empty() {
        let text=fs::read_to_string(input).map_err(|e|e.to_string())?;
        for line in text.lines().skip(1) {
            let f:Vec<_>=line.split('\t').collect();if f.len()!=3 {return Err("seed TSV columns".into());}
            let d:usize=f[1].parse().map_err(|_|"degree")?;if !(3..=MAX_D).contains(&d) {return Err("degree range".into());}
            let t=cp_tuple(f[2],d)?;let genus=cp_stats(&t,d).0;
            for rep in 0..4 {
                let mut u=t;let mut rng=SplitMix64::new(seed+(rep as u64)*31+f[0].bytes().fold(0u64,|a,b|a.wrapping_mul(37).wrapping_add(b as u64)));
                for _ in 0..4000 {u=hurwitz_move_dynamic(d,u,rng.below(6));}
                for block in 0..8 {
                    let mut hist=[0u64;MAX_D+1];let mut both=0;
                    for _ in 0..4000 {
                        for _ in 0..8 {u=hurwitz_move_dynamic(d,u,rng.below(6));}
                        let(g,delta,_,l,r)=cp_stats(&u,d);assert_eq!(g,genus);
                        hist[delta]+=1;if l==1&&r==1 {both+=1;}
                    }
                    cp_row(&mut out,"walk",f[0],d,"S",rep,block,genus,&hist,both)?;
                }
            }
        }
    }
    out.flush().map_err(|e|e.to_string())
}
#[cfg(test)]
mod coin_probe_tests {
    use super::*;
    #[test] fn every_s3_transitive_product_one_tuple_checks_defect() {
        let d=3;let mut p=Perm::identity(d);let mut ps=vec![p];while next_permutation(&mut p.img[..d]) {ps.push(p);}
        for &a in &ps {for &b in &ps {for &c in &ps {
            let t=[a,b,c,product_one_fourth(d,a,b,c)];if !generators_transitive(d,&t[..3]) {continue;}
            let(g,_,_,_,_)=cp_stats(&t,d);
            for m in 0..6 {let u=hurwitz_move_dynamic(d,t,m);assert!(tuple_product_is_identity(d,&u));assert_eq!(cp_stats(&u,d).0,g);}
        }}}
    }
    #[test] fn full_cycles_force_defect_parity() {
        for d in [5usize,6] {
            let mut a=Perm::identity(d);for x in 0..d {a.img[x]=((x+1)%d) as u8;}
            let t=[a,a.inverse(d),a,a.inverse(d)];let(g,delta,_,l,r)=cp_stats(&t,d);
            assert_eq!((g,l,r),(d-1,1,1));assert_eq!(delta%2,(d-1)%2);
        }
    }
}
