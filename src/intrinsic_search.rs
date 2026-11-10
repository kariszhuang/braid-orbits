// Standalone finite, explicitly marked S3-cover experiment. Coordinate algebra
// is separately reproduced in Python; only literal words/full closures certify.
fn intrinsic_pair_key(a:&[u8],b:&[u8])->Vec<u8>{
    let n=a.len();let mut best=vec![u8::MAX;2*n];
    for root in 0..n{let mut labels=vec![u8::MAX;n];labels[root]=0;let mut queue=vec![root];let mut head=0;
        while head<queue.len(){let x=queue[head];head+=1;for p in [a,b]{let y=p[x]as usize;if labels[y]==u8::MAX{labels[y]=queue.len()as u8;queue.push(y);}}}
        if queue.len()!=n{continue;}let mut key=Vec::with_capacity(2*n);
        for p in [a,b]{for &x in &queue{key.push(labels[p[x]as usize]);}}
        if key<best{best=key;}
    }best
}
fn intrinsic_next_permutation(p:&mut[u8])->bool{
    let Some(i)=(0..p.len()-1).rev().find(|&i|p[i]<p[i+1])else{return false;};
    let j=(i+1..p.len()).rev().find(|&j|p[i]<p[j]).unwrap();p.swap(i,j);p[i+1..].reverse();true
}
fn intrinsic_profile(p:&[u8])->Vec<usize>{let mut out:Vec<_>=crate::intrinsic_s3::cycles(p).iter().map(|c|c.len()).collect();out.sort_unstable();out}
fn intrinsic_cycle_stats(p:&[u8])->(usize,usize){
    let mut seen=0u32;let(mut count,mut longest)=(0,0);
    for start in 0..p.len(){if seen&(1<<start)!=0{continue;}count+=1;let mut x=start;let mut size=0;
        while seen&(1<<x)==0{seen|=1<<x;size+=1;x=p[x]as usize;}longest=longest.max(size);
    }(count,longest)
}
fn intrinsic_elliptic_betas(q:usize)->Vec<Vec<u8>>{
    fn shapes(left:usize,count:usize,minimum:usize,maximum:usize,row:&mut Vec<usize>,out:&mut Vec<Vec<usize>>){
        if count==0{if left==0{out.push(row.clone());}return;}
        for k in minimum..=maximum.min(left/count){if left-k>(count-1)*maximum{continue;}row.push(k);shapes(left-k,count-1,k,maximum,row,out);row.pop();}
    }
    let mut profiles=Vec::new();for c in 1..=q/3{shapes(q,c,1,q/2,&mut Vec::new(),&mut profiles);}
    profiles.into_iter().map(|parts|{let mut beta=vec![0;q];let mut at=0;for k in parts{for j in 0..k{beta[at+j]=(at+(j+1)%k)as u8;}at+=k;}beta}).collect()
}
fn intrinsic_frames(q:usize,all:bool,deadline:Instant)->Result<Vec<crate::intrinsic_s3::IntrinsicFrame>,String>{
    use crate::intrinsic_s3::{IntrinsicFrame,compose};
    if all{
        let mut out=Vec::new();let mut seen=HashSet::new();
        for beta in intrinsic_elliptic_betas(q){let cb=intrinsic_cycle_stats(&beta).0;let mut gamma:Vec<u8>=(0..q as u8).collect();let mut tested=0u64;
            loop{let(cg,lg)=intrinsic_cycle_stats(&gamma);
                tested+=1;if tested&1023==0&&Instant::now()>=deadline{return Err("absolute deadline during quotient enumeration".into());}
                if cg>=cb&&cg+2*cb<=q&&lg<=q/2{
                    let delta=compose(&beta,&gamma);let(cd,ld)=intrinsic_cycle_stats(&delta);
                    if cd>=cb&&cb+cg+cd==q&&ld<=q/2{if let Ok(f)=IntrinsicFrame::new(beta.clone(),gamma.clone()){
                        if seen.insert(intrinsic_pair_key(&beta,&gamma)){out.push(f);}
                    }}
                }
                if !intrinsic_next_permutation(&mut gamma){break;}
            }
        }return Ok(out);
    }
    let(mut beta,shape,delta_shape):(Vec<u8>,Vec<usize>,Vec<usize>)=if q==8{((0..q).map(|j|(4*(j/4)+(j+1)%4)as u8).collect(),vec![2;4],vec![4;2])}
        else{((0..q).map(|j|(3*(j/3)+(j+1)%3)as u8).collect(),vec![3;q/3],vec![3;q/3])};
    let mut gamma:Vec<u8>=(0..q as u8).collect();let mut seen=HashSet::new();let mut out=Vec::new();
    loop{if intrinsic_profile(&gamma)==shape&&intrinsic_profile(&compose(&beta,&gamma))==delta_shape{
        if let Ok(frame)=IntrinsicFrame::new(beta.clone(),gamma.clone()){
            if seen.insert(intrinsic_pair_key(&beta,&gamma)){out.push(frame);}
        }
    }if !intrinsic_next_permutation(&mut gamma){break;}}
    // Keep the control's exact original marking, not an accidental inversion.
    if q==3{beta=vec![1,2,0];assert_eq!(out.len(),1);assert_eq!(out[0].beta,beta);}
    Ok(out)
}
fn intrinsic_list<T:std::fmt::Display>(v:&[T])->String{v.iter().map(|x|x.to_string()).collect::<Vec<_>>().join(",")}
pub fn intrinsic_main(args:&[String])->Result<(),String>{
    use crate::intrinsic_s3::{TRANS,PERMS};
    let arg=|k:&str,d:&str|args.windows(2).find(|w|w[0]==k).map(|w|w[1].clone()).unwrap_or(d.into());
    let q:usize=arg("--q","9").parse().map_err(|_|"quotient degree")?;
    let all=args.iter().any(|s|s=="--all-elliptic");
    if if all{!(6..=11).contains(&q)}else{![3,6,8,9].contains(&q)}{return Err("bounded uniform q3,6,8,9 or all-elliptic q6..11".into());}
    let n=3*q;let out=PathBuf::from(arg("--out",""));if out.as_os_str().is_empty(){return Err("fresh --out required".into());}
    let wall: f64=arg("--deadline-epoch","1792762200").parse().map_err(|_|"deadline")?;
    let now=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_err(|e|e.to_string())?.as_secs_f64();
    let remaining=wall-now;if !remaining.is_finite()||remaining<=225.0{return Err("absolute research deadline/drain".into());}
    let began=Instant::now();let deadline=began+Duration::from_secs_f64(remaining-225.0);
    let cap:usize=arg("--cap","1000000").parse().map_err(|_|"cap")?;if cap==0{return Err("positive cap".into());}
    let seconds:u64=arg("--closure-seconds","15").parse().map_err(|_|"closure seconds")?;if seconds==0{return Err("positive closure seconds".into());}
    let words=if args.iter().any(|v|v=="--words"){
        fs::read_to_string(arg("--words","")).map_err(|e|e.to_string())?.lines().filter(|s|!s.is_empty()).map(|line|line.split(',').map(|x|x.parse::<usize>().map_err(|_|"literal move")).collect::<Result<Vec<_>,_>>()).collect::<Result<Vec<_>,_>>()?
    }else{Vec::new()};
    if words.len()>64||words.iter().any(|w|w.is_empty()||w.len()>24||w.iter().any(|&m|m>=6)){return Err("bounded literal word bank".into());}
    let trie=G3Words::new(&words);let frames=intrinsic_frames(q,all,deadline)?;
    fs::create_dir(&out).map_err(|e|e.to_string())?;
    let mut quotients=BufWriter::new(File::create(out.join("QUOTIENTS.tsv")).map_err(|e|e.to_string())?);
    let mut cases=BufWriter::new(File::create(out.join("CASES.tsv")).map_err(|e|e.to_string())?);
    let mut coords=BufWriter::new(File::create(out.join("COORDINATES.tsv")).map_err(|e|e.to_string())?);
    let mut leads=BufWriter::new(File::create(out.join("LEADS.tsv")).map_err(|e|e.to_string())?);
    let mut resolutions=BufWriter::new(File::create(out.join("RESOLUTIONS.tsv")).map_err(|e|e.to_string())?);
    writeln!(quotients,"quotient\tbeta\tgamma\tprimal\tdual\thandles").map_err(|e|e.to_string())?;
    writeln!(cases,"quotient\tsupport\torbit\teligible\tparameters\tsolved\ttransitive\tsmooth").map_err(|e|e.to_string())?;
    writeln!(coords,"candidate\tquotient\tsupport\tlocal\thandles\ttransitive\tsmooth\tseed").map_err(|e|e.to_string())?;
    writeln!(resolutions,"candidate\tstatus\tword\tfolder\tseed").map_err(|e|e.to_string())?;
    let(mut attempted,mut solved,mut transitive,mut smooth,mut negative,mut closed,mut unresolved)=(0u64,0u64,0u64,0u64,0u64,0u64,0u64);
    let mut complete=true;
    'quotient:for(fi,f)in frames.iter().enumerate(){
        writeln!(quotients,"{}\t{}\t{}\t{}\t{}\t{}",fi,intrinsic_list(&f.beta),intrinsic_list(&f.gamma),intrinsic_list(&f.primal),intrinsic_list(&f.dual),intrinsic_list(&f.handles)).map_err(|e|e.to_string())?;
        let mut visited=HashSet::new();let masks:Vec<u32>=if q==3{vec![1]}else{(0..1u32<<q).filter(|m|m.count_ones()==4).collect()};
        for mask in masks{if !visited.insert(mask){continue;}
            let(orbit,eligible)=if q==3{(vec![1],true)}else{f.support_orbit(mask)};
            visited.extend(orbit.iter().copied());
            let support:Vec<_>=(0..q).filter(|&j|mask&(1<<j)!=0).collect();let limit=if q==3{1}else{81};
            let(mut cparam,mut csolved,mut ctrans,mut csmooth)=(0,0,0,0);
            if eligible{for assignment in 0..limit{let mut local=vec![0;q];let mut code=assignment;
                for &j in &support{local[j]=if q==3{3}else{let value=TRANS[code%3];code/=3;value};}
                for u in 0..6{for v in 0..6{
                    if Instant::now()>=deadline{complete=false;break 'quotient;}
                    attempted+=1;cparam+=1;
                    let Some(lift)=f.solve(&local,[u,v])else{continue;};solved+=1;csolved+=1;
                    let t:[Perm;4]=std::array::from_fn(|j|{let mut p=Perm::identity(n);p.img[..n].copy_from_slice(&lift[j][..n]);p});
                    assert!(tuple_product_is_identity(n,&t));assert_eq!(t.iter().map(|p|sm_index(*p,n)).sum::<usize>(),2*n+if q==3{2}else{4});
                    let tr=generators_transitive(n,&t[..3]);let sm=sm_values(&t,n).0;let text=tuple_string_dynamic(n,&t);
                    writeln!(coords,"{attempted}\t{fi}\t{mask}\t{}\t{u},{v}\t{tr}\t{sm}\t{text}",intrinsic_list(&local)).map_err(|e|e.to_string())?;
                    if !tr{continue;}transitive+=1;ctrans+=1;if !sm{continue;}smooth+=1;csmooth+=1;
                    writeln!(leads,"{attempted}\t{text}").map_err(|e|e.to_string())?;leads.flush().map_err(|e|e.to_string())?;
                    let frame=ForestFrame::new([t[0],t[1]],n);let ix=t.map(|p|sm_index(p,n));
                    let word=frame.short_word(t[2],ix[2],ix[3],n).map(|j|FF_WORDS[j].to_vec())
                        .or_else(||trie.failure_reference(t,n)).or_else(||g3_outer_failure(t,n,32));
                    let scan=if let Some(word)=word{assert!(ff_replay_negative(t,n,&word));SmoothScan{status:"nonsmooth",states:Vec::new(),word}}
                        else{g3_macro_scan(&t,n,cap,deadline.min(Instant::now()+Duration::from_secs(seconds)))};
                    let mut folder=String::new();if scan.status=="nonsmooth"{assert!(ff_replay_negative(t,n,&scan.word));negative+=1;}
                    else{folder=format!("component_{attempted:09}");sm_save(&t,n,&scan,&out.join(&folder))?;
                        if scan.status=="closed_smooth"{closed+=1;}else{unresolved+=1;}}
                    writeln!(resolutions,"{attempted}\t{}\t{}\t{folder}\t{text}",scan.status,intrinsic_list(&scan.word)).map_err(|e|e.to_string())?;resolutions.flush().map_err(|e|e.to_string())?;
                    if q!=3&&scan.status=="closed_smooth"{complete=false;break 'quotient;}
                }}
            }}
            writeln!(cases,"{fi}\t{mask}\t{}\t{eligible}\t{cparam}\t{csolved}\t{ctrans}\t{csmooth}",intrinsic_list(&orbit)).map_err(|e|e.to_string())?;
        }
    }
    fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\tintrinsic\t0\t{attempted}\t{solved}\t{transitive}\t{solved}\t{negative}\t0\t{closed}\t{unresolved}\t{}\n",began.elapsed().as_secs_f64())).map_err(|e|e.to_string())?;
    fs::write(out.join("POPULATION.tsv"),format!("quotient_degree\tall_elliptic\tquotient_frames\tcomplete\tparameters\tsolved\ttransitive\tsmooth\tnegative\tclosed\tunresolved\tdeadline_epoch\n{q}\t{all}\t{}\t{complete}\t{attempted}\t{solved}\t{transitive}\t{smooth}\t{negative}\t{closed}\t{unresolved}\t{wall}\n",frames.len())).map_err(|e|e.to_string())?;
    println!("intrinsic q={q} frames={} parameters={attempted} solved={solved} smooth={smooth} negative={negative} closed={closed} unresolved={unresolved} complete={complete}",frames.len());
    let _=PERMS;Ok(())
}
