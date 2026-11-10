// Isolated <=128-point statistic kernel: production Perm is capped at34.
// No group certification, canonicalization or discovery engine is duplicated.
use std::{collections::BTreeMap, env, fs, io::{self, Write, BufWriter}};
const N:usize=128;
type P=[u8;N];
type T=[P;4];
struct R(u64);
impl R {
    fn next(&mut self)->u64 {self.0=self.0.wrapping_add(0x9e3779b97f4a7c15);let mut z=self.0;z=(z^(z>>30)).wrapping_mul(0xbf58476d1ce4e5b9);z=(z^(z>>27)).wrapping_mul(0x94d049bb133111eb);z^(z>>31)}
    fn below(&mut self,n:usize)->usize {let n=n as u64;let limit=n.wrapping_neg()%n;loop{let v=self.next();if v>=limit{return(v%n)as usize}}}
}
fn id()->P {std::array::from_fn(|i|i as u8)}
fn mul(a:&P,b:&P,d:usize)->P {let mut p=id();for i in 0..d{p[i]=a[b[i]as usize]}p}
fn inv(a:&P,d:usize)->P {let mut p=id();for i in 0..d{p[a[i]as usize]=i as u8}p}
fn cyc(p:&P,d:usize)->usize {let mut seen=0u128;let mut c=0;for i in 0..d{if seen&(1u128<<i)==0{c+=1;let mut x=i;while seen&(1u128<<x)==0{seen|=1u128<<x;x=p[x]as usize}}}c}
fn solve(t:&mut T,d:usize){t[3]=inv(&mul(&mul(&t[0],&t[1],d),&t[2],d),d)}
fn genus(t:&T,d:usize)->i32 {let s:usize=t.iter().map(|p|d-cyc(p,d)).sum();assert_eq!(s%2,0);1-d as i32+(s/2)as i32}
fn root(pa:&mut P,mut x:usize)->usize {while pa[x]as usize!=x{pa[x]=pa[pa[x]as usize];x=pa[x]as usize}x}
fn orbits(ps:&[P],d:usize)->usize {let mut pa=id();for p in ps{for(x,&y)in p[..d].iter().enumerate(){let a=root(&mut pa,x);let b=root(&mut pa,y as usize);pa[a]=b as u8}}(0..d).filter(|&x|root(&mut pa,x)==x).count()}
fn valid(t:&T,d:usize)->bool {t.iter().all(|p|p[..d]!=id()[..d])&&orbits(&t[..3],d)==1}
fn stat(t:&T,d:usize)->(usize,usize,usize){let c=cyc(&mul(&t[0],&t[1],d),d);let l=orbits(&t[..2],d);let r=orbits(&t[2..],d);assert!(c+1>=l+r);assert!((c+1-l-r)as i32<=genus(t,d));(c,l,r)}
fn braid(t:&T,d:usize,m:usize)->T {let mut u=*t;let i=m/2;let a=t[i];let b=t[i+1];if m%2==0{u[i]=b;u[i+1]=mul(&mul(&inv(&b,d),&a,d),&b,d)}else{u[i]=mul(&mul(&a,&b,d),&inv(&a,d),d);u[i+1]=a}u}
fn random(d:usize,a:bool,r:&mut R)->P {let mut p=id();for i in(1..d).rev(){p.swap(i,r.below(i+1))}if a&&(d-cyc(&p,d))%2==1{p.swap(0,1)}p}
fn seed(d:usize,g:usize,r:&mut R)->T {
    assert!(g>0&&g<d);let labels=random(d,false,r);let mut a=id();for i in 0..d{a[labels[i]as usize]=labels[(i+1)%d]}
    let mut c=id();let mut k=0;while k<g{let x=r.below(d);let mut y=r.below(d-1);if y>=x{y+=1}let mut q=c;q.swap(x,y);if cyc(&q,d)+1==cyc(&c,d){c=q;k+=1}}
    let mut t=[a,inv(&a,d),c,inv(&c,d)];for _ in 0..100{t=braid(&t,d,r.below(6))}assert!(valid(&t,d));assert_eq!(genus(&t,d),g as i32);t
}
fn step(t:&mut T,d:usize,g:usize,r:&mut R)->bool {
    if r.below(4)==0 {*t=braid(t,d,r.below(6));return true}
    let mut u=*t;let j=r.below(3);let x=r.below(d);let mut y=r.below(d-1);if y>=x{y+=1}u[j].swap(x,y);solve(&mut u,d);
    if genus(&u,d)==g as i32&&valid(&u,d){*t=u;true}else{false}
}
fn smooth_control(d:usize)->T {
    assert_eq!(d%3,0);let k=d/3;let mut x=id();for(i,v)in x[..k].iter_mut().enumerate(){*v=((i+1)%k)as u8}let mut y=id();y.swap(0,1);
    let z=inv(&mul(&y,&x,k),k);let tau=mul(&mul(&mul(&inv(&x,k),&y,k),&x,k),&inv(&y,k),k);
    let mut t=[id();4];t[0][..k].copy_from_slice(&tau[..k]);for b in 0..3{for i in 0..k{t[1][b*k+i]=(((b+1)%3)*k+i)as u8;t[2][b*k+i]=(((b+1)%3)*k+[x,y,z][b][i]as usize)as u8}}solve(&mut t,d);assert_eq!(genus(&t,d),2);assert!(valid(&t,d));t
}
fn save_tuple(w:&mut impl Write,job:usize,phase:&str,t:&T,d:usize)->io::Result<()> {write!(w,"{job}\t{phase}\t")?;for(j,p)in t.iter().enumerate(){for(i,v)in p[..d].iter().enumerate(){if j!=0||i!=0{write!(w,",")?}write!(w,"{v}")?}}writeln!(w)}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=env::args().collect();if args.len()!=4{return Err("usage: defect_degree_probe JOBS.tsv TUPLES.tsv METRICS.tsv".into())}
    let jobs=fs::read_to_string(&args[1])?;let mut out=BufWriter::with_capacity(1<<20,io::stdout().lock());let mut tuples=BufWriter::new(fs::File::create(&args[2])?);let mut metrics=BufWriter::new(fs::File::create(&args[3])?);
    writeln!(out,"job\tkind\tdegree\tambient\treplicate\tblock\tgenus\tcycles\tleft\tright\tcount")?;writeln!(tuples,"job\tphase\tpermutations")?;writeln!(metrics,"job\taccepted\tproposals\tretained\tseconds")?;
    for(row,line)in jobs.lines().skip(1).enumerate(){let f:Vec<_>=line.split('\t').collect();if f.len()!=9{return Err("nine job columns required".into())}let kind=f[0];let d:usize=f[1].parse()?;let g:usize=f[2].parse()?;let a=f[3];let rep:usize=f[4].parse()?;let mut r=R(f[5].parse()?);let samples:usize=f[6].parse()?;let burn:usize=f[7].parse()?;let stride:usize=f[8].parse()?;assert!((3..=N).contains(&d));assert!(samples>0&&stride>0);let timer=std::time::Instant::now();let mut accepted=0;let mut proposals=0;let mut retained=0;
        if kind=="iid"{let mut hist=BTreeMap::new();for _ in 0..samples{let mut t=[random(d,a=="A",&mut r),random(d,a=="A",&mut r),random(d,a=="A",&mut r),id()];solve(&mut t,d);if !valid(&t,d){continue}let gg=genus(&t,d);assert!(gg>=0);let(c,l,rr)=stat(&t,d);*hist.entry((gg,c,l,rr)).or_insert(0u64)+=1;retained+=1;if retained==1{save_tuple(&mut tuples,row,"initial",&t,d)?}}
            for((gg,c,l,rr),n)in hist{writeln!(out,"{row}\t{kind}\t{d}\t{a}\t{rep}\t0\t{gg}\t{c}\t{l}\t{rr}\t{n}")?}proposals=samples;
        }else{
            let mut t=if kind=="smooth_control"{smooth_control(d)}else if kind=="parity_control"{let mut p=id();for(i,v)in p[..d].iter_mut().enumerate(){*v=((i+1)%d)as u8}[p,inv(&p,d),p,inv(&p,d)]}else{assert_eq!(kind,"conditional");seed(d,g,&mut r)};assert_eq!(genus(&t,d),g as i32);save_tuple(&mut tuples,row,"initial",&t,d)?;
            for _ in 0..burn{proposals+=1;if kind=="conditional"{accepted+=usize::from(step(&mut t,d,g,&mut r))}else{t=braid(&t,d,r.below(6));accepted+=1}}
            for block in 0..8{let mut hist=BTreeMap::new();for _ in 0..samples{for _ in 0..stride{proposals+=1;if kind=="conditional"{accepted+=usize::from(step(&mut t,d,g,&mut r))}else{t=braid(&t,d,r.below(6));accepted+=1}}assert_eq!(genus(&t,d),g as i32);let(c,l,rr)=stat(&t,d);if kind=="smooth_control"{assert_eq!(c+1-l-rr,0)}if kind=="parity_control"{assert_eq!((c+1-l-rr)%2,(d-1)%2)}*hist.entry((c,l,rr)).or_insert(0u64)+=1;retained+=1}
                for((c,l,rr),n)in hist{writeln!(out,"{row}\t{kind}\t{d}\t{a}\t{rep}\t{block}\t{g}\t{c}\t{l}\t{rr}\t{n}")?}
            }save_tuple(&mut tuples,row,"final",&t,d)?;
        }
        writeln!(metrics,"{row}\t{accepted}\t{proposals}\t{retained}\t{:.6}",timer.elapsed().as_secs_f64())?;out.flush()?;metrics.flush()?;tuples.flush()?;eprintln!("job={row} kind={kind} d={d} g={g} retained={retained} seconds={:.3}",timer.elapsed().as_secs_f64());
    }Ok(())
}
#[cfg(test)]mod tests{use super::*;
    #[test]fn kernels_and_controls(){let mut r=R(13);for d in[3,4,9,60,61,120,128]{for _ in 0..100{let p=random(d,false,&mut r);assert_eq!(mul(&p,&inv(&p,d),d),id());let t=seed(d,(d-1).min(4),&mut r);for m in 0..6{let u=braid(&t,d,m);assert_eq!(braid(&u,d,m^1),t);assert_eq!(genus(&u,d),genus(&t,d));assert!(valid(&u,d));let _=stat(&u,d);}}}for d in[9,60,90,120]{let t=smooth_control(d);assert_eq!(genus(&t,d),2);for m in 0..6{let(c,l,r)=stat(&braid(&t,d,m),d);assert_eq!(c+1-l-r,0)}}}
    #[test]fn symmetric_proposal(){let mut r=R(17);for d in[4,9,60]{let mut t=seed(d,2,&mut r);for _ in 0..1000{let old=t;let _=step(&mut t,d,2,&mut r);assert_eq!(genus(&t,d),2);assert!(valid(&t,d));for j in 0..3{for x in 0..d.min(6){let y=(x+1)%d;let mut u=old;u[j].swap(x,y);solve(&mut u,d);u[j].swap(x,y);solve(&mut u,d);assert_eq!(u,old)}}}}}
}
