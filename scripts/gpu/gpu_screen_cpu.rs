//! Native C ABI shared by discovery and benchmarks; exact research primitives.
#[allow(dead_code)]
mod production {
    include!("../src/main.rs");

    type Session=smart::hybrid::Session;
    #[no_mangle]
    pub unsafe extern "C" fn hybrid_open(out:*const std::ffi::c_char,known:*const std::ffi::c_char)->*mut Session {
        if out.is_null() || known.is_null() {return std::ptr::null_mut();}
        let result=(|| {
            let out=std::ffi::CStr::from_ptr(out).to_str().map_err(|e|e.to_string())?;
            let known=std::ffi::CStr::from_ptr(known).to_str().map_err(|e|e.to_string())?;
            smart::hybrid::open(Path::new(out),Path::new(known))
        })();
        match result {Ok(s)=>Box::into_raw(Box::new(s)),Err(e)=>{eprintln!("hybrid init: {e}");std::ptr::null_mut()}}
    }
    #[no_mangle]
    pub unsafe extern "C" fn hybrid_close(s:*mut Session) {if !s.is_null() {drop(Box::from_raw(s));}}

    #[no_mangle]
    pub unsafe extern "C" fn hybrid_fill(s:*const Session,d:u32,n:u32,family:u32,seed:u64,workers:u32,out:*mut u8)->u32 {
        if s.is_null() || out.is_null() || workers==0 || workers>10 || n==0 || n>1048576
           || !matches!((d,family),(9..=17,0..=1)|(18,2..=3)|(25,4)) {return 1;}
        let s=&*s;let d=d as usize;let n=n as usize;let address=out as usize;let next=AtomicU64::new(0);
        thread::scope(|scope| {
            for _ in 0..workers {let next=&next;scope.spawn(move || {
                loop {
                    let shard=next.fetch_add(1,Ordering::Relaxed) as usize;let start=shard*4096;if start>=n {break;}
                    let count=(n-start).min(4096);
                    // Atomic shard ownership gives each worker disjoint mutable bytes.
                    let bytes=std::slice::from_raw_parts_mut((address+start*3*d) as *mut u8,count*3*d);
                    smart::hybrid::generate(s,d,family,seed^(shard as u64).wrapping_mul(0xD1B54A32D192ED03),bytes);
                }
            });}
        });0
    }
    #[no_mangle]
    pub unsafe extern "C" fn hybrid_cpu_screen(d:u32,n:u32,workers:u32,input:*const u8,out:*mut u32)->u32 {
        parallel_screen(d,n,workers,1,input,out)
    }
    #[no_mangle]
    pub unsafe extern "C" fn hybrid_cpu_probe(d:u32,n:u32,workers:u32,input:*const u8,out:*mut u32)->u32 {
        parallel_screen(d,n,workers,2,input,out)
    }
    unsafe fn parallel_screen(d:u32,n:u32,workers:u32,mode:u32,input:*const u8,out:*mut u32)->u32 {
        if !(4..=34).contains(&d) || workers==0 || workers>10 || n==0 || n>1048576 || input.is_null() || out.is_null() {return 1;}
        let bytes=std::slice::from_raw_parts(input,n as usize*3*d as usize);
        let output=std::slice::from_raw_parts_mut(out,n as usize);let stride=(n as usize+workers as usize-1)/workers as usize;
        thread::scope(|scope| {
            for (rows,flags) in bytes.chunks(stride*3*d as usize).zip(output.chunks_mut(stride)) {
                scope.spawn(move || {assert_eq!(bench_screen(d,flags.len() as u32,mode,rows.as_ptr(),flags.as_mut_ptr()),0);});
            }
        });0
    }
    #[no_mangle]
    pub unsafe extern "C" fn hybrid_validate(s:*const Session,d:u32,n:u32,workers:u32,cap:u32,millis:u64,input:*const u8,flags:*const u32,stats:*mut u64)->u32 {
        if s.is_null() || input.is_null() || flags.is_null() || stats.is_null() || !(9..=25).contains(&d)
           || n==0 || n>1048576 || workers==0 || workers>10 || cap==0 || millis==0 {return 1;}
        let bytes=std::slice::from_raw_parts(input,n as usize*3*d as usize);
        let flags=std::slice::from_raw_parts(flags,n as usize);
        match smart::hybrid::validate(&*s,d as usize,bytes,flags,workers as usize,cap as usize,millis) {
            Ok(values)=>{std::slice::from_raw_parts_mut(stats,12).copy_from_slice(&values);0},
            Err(e)=>{eprintln!("hybrid validation: {e}");2}
        }
    }

    fn indices(t: &[Perm;4], d:usize)->[u32;4] {
        t.map(|p| {
            let mut seen=[false;MAX_D];let mut cycles=0;
            for x in 0..d {
                if seen[x] {continue;}cycles+=1;let mut y=x;
                while !seen[y] {seen[y]=true;y=p.img[y] as usize;}
            }
            (d-cycles) as u32
        })
    }

    fn braid_probe(d:usize,t:&mut[Perm;4],stop:bool)->bool {
        let mut state=0xcbf29ce484222325u64;
        for p in &t[..3] {for &x in &p.img[..d] {state=(state^x as u64).wrapping_mul(0x100000001b3);}}
        let mut rng=SplitMix64::new(state);
        for _ in 0..64 {
            let mv=rng.below(6);*t=hurwitz_move_dynamic(d,*t,mv);
            if stop && mv%3==1 && generators_transitive(d,&t[..2]) {return true;}
        }
        false
    }

    fn candidate(d:usize,family:u32,rng:&mut SplitMix64,i:usize,base:Option<[Perm;4]>)->[Perm;4] {
        let mut t=[Perm::identity(d);4];
        if let Some(seed)=base {
            t=seed;
            // Include exact failures as controls, then independent local mutations.
            if i%16!=0 {
                let entry=rng.below(3);let mut h=Perm::identity(d);
                for _ in 0..[1,2,4,8][i%4] {let a=rng.below(d);let b=rng.below(d);h.img.swap(a,b);}
                t[entry]=t[entry].conjugate_by(h,d);
            }
        } else {
            for p in &mut t[..3] {
                if family==0 {*p=random_perm(d,Ambient::Symmetric,rng);}
                else {
                    let mut labels:Vec<_>=(0..d).collect();rng.shuffle(&mut labels);
                    let minimum=(d+2)/3;let k=minimum+rng.below(d/2-minimum+1);
                    for pair in labels[..2*k].chunks_exact(2) {p.img[pair[0]]=pair[1] as u8;p.img[pair[1]]=pair[0] as u8;}
                }
            }
        }
        t[3]=product_one_fourth(d,t[0],t[1],t[2]);t
    }

    #[no_mangle]
    pub unsafe extern "C" fn bench_generate(d:u32,n:u32,family:u32,seed:u64,out:*mut u8)->u32 {
        if !(4..=34).contains(&d) || n==0 || out.is_null() || family>2 {return 1;}
        let d=d as usize;let mut rng=SplitMix64::new(seed);
        let base=if family==2 {
            let text=match d {
                18=>include_str!("../evidence/certificates/CERTIFIED_S18_XXYY_genus1_weak.txt"),
                25=>include_str!("../evidence/certificates/CERTIFIED_A25_weak.txt"),
                _=>return 2,
            };
            let seed=text.lines().find_map(|s|s.strip_prefix("seed=")).unwrap();
            let mut t=[Perm::identity(d);4];
            for (j,part) in seed.trim_matches(['[',']']).split(',').enumerate() {
                let mut used=[false;MAX_D];
                for cycle in part.split('(').skip(1) {
                    let points:Vec<usize>=cycle.split(')').next().unwrap().split_whitespace().map(|x|x.parse::<usize>().unwrap()-1).collect();
                    for &x in &points {assert!(x<d && !used[x]);used[x]=true;}
                    for i in 0..points.len() {t[j].img[points[i]]=points[(i+1)%points.len()] as u8;}
                }
            }
            assert!(tuple_product_is_identity(d,&t));Some(t)
        } else {None};
        let bytes=std::slice::from_raw_parts_mut(out,n as usize*3*d);
        for (i,row) in bytes.chunks_exact_mut(3*d).enumerate() {
            let t=candidate(d,family,&mut rng,i,base);
            for j in 0..3 {row[j*d..(j+1)*d].copy_from_slice(&t[j].img[..d]);}
        }
        0
    }

    #[no_mangle]
    pub unsafe extern "C" fn bench_screen(d:u32,n:u32,mode:u32,input:*const u8,out:*mut u32)->u32 {
        if !(4..=34).contains(&d) || n==0 || mode>3 || (mode==3 && d!=4) || input.is_null() || out.is_null() {return 1;}
        let d=d as usize;let bytes=std::slice::from_raw_parts(input,n as usize*3*d);
        let output=std::slice::from_raw_parts_mut(out,n as usize);
        for (row,result) in bytes.chunks_exact(3*d).zip(output.iter_mut()) {
            let mut t=[Perm::identity(d);4];
            for j in 0..3 {t[j].img[..d].copy_from_slice(&row[j*d..(j+1)*d]);}
            t[3]=product_one_fourth(d,t[0],t[1],t[2]);
            if mode==3 {
                braid_probe(d,&mut t,false);let mut packed=0;
                for j in 0..4 {let mut rank=0;for x in 0..d {let digit=(x+1..d).filter(|&y|t[j].img[y]<t[j].img[x]).count();rank=rank*(d-x)+digit;}packed|=(rank as u32)<<(5*j);}
                *result=packed;continue;
            }
            let ix=indices(&t,d);
            if mode>=1 {
                *result=0;
                if ix.iter().any(|&i|i>d as u32-3) || ix.iter().sum::<u32>()<2*d as u32-2 {continue;}
                if mode==2 && !generators_transitive(d,&t[..3]) {continue;}
                let transitive=(0..4).any(|i|(i+1..4).any(|j|generators_transitive(d,&[t[i],t[j]])));
                *result=u32::from(!transitive && (mode!=2 || !braid_probe(d,&mut t,true)));continue;
            }
            let mut mask=0;let mut bit=0;
            for i in 0..4 {for j in i+1..4 {
                if generators_transitive(d,&[t[i],t[j]]) {mask|=1<<bit;}bit+=1;
            }}
            *result=ix[0]|ix[1]<<6|ix[2]<<12|ix[3]<<18|mask<<24;
        }
        0
    }

    #[no_mangle]
    pub unsafe extern "C" fn bench_compact(n:u32,input:*const u32,out:*mut u32)->u32 {
        let input=std::slice::from_raw_parts(input,n as usize);let out=std::slice::from_raw_parts_mut(out,n as usize);
        let mut used=0;
        for (i,&value) in input.iter().enumerate() {if value==1 {out[used]=i as u32;used+=1;}}
        used as u32
    }
}
