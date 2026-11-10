//! Benchmark adapter only; certificate admission uses the independent auditor.
use std::slice;
#[no_mangle]
pub unsafe extern "C" fn finite_prefix(frames:u32,count:u32,workers:u32,f:*const u64,c:*const u64,out:*mut u16)->i32{
    if frames==0||count==0||workers==0||workers>10||f.is_null()||c.is_null()||out.is_null(){return 1;}
    let(nc,nf)=(count as usize,frames as usize);
    let factors=slice::from_raw_parts(f,nf*5);let cs=slice::from_raw_parts(c,nc);
    let output=slice::from_raw_parts_mut(out,nf*nc);let chunk=nf.div_ceil(workers as usize)*nc;
    std::thread::scope(|scope|{for(part,block)in output.chunks_mut(chunk).enumerate(){
        scope.spawn(move||{for(local,row)in block.chunks_mut(nc).enumerate(){let frame=part*(chunk/nc)+local;
            for(&c,result)in cs.iter().zip(row){*result=0;
                'factors:for j in 0..5{let f=factors[5*frame+j];let mut x=0;
                    for len in 1..=15{let y=(c>>(4*x))&15;x=(f>>(4*y))&15;
                        if x==0{if len!=15{*result=j as u16+1;break 'factors;}break;}
                    }
                }
            }
        }});
    }});0
}
