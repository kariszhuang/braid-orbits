#include <metal_stdlib>
using namespace metal;
kernel void prefix(device const ulong* cs [[buffer(0)]],
                   device const ulong* factors [[buffer(1)]],
                   device ushort* output [[buffer(2)]],
                   constant uint* params [[buffer(3)]],uint tid [[thread_position_in_grid]]){
    uint count=params[0],states=params[1];if(tid>=states)return;
    ulong c=cs[tid%count];uint frame=tid/count;
    for(uint j=0;j<5;++j){ulong f=factors[frame*5+j];uint x=0;
        for(uint len=1;len<=15;++len){uint y=uint((c>>(4*x))&15);x=uint((f>>(4*y))&15);
            if(x==0){if(len!=15){output[tid]=ushort(j+1);return;}break;}
        }
    }output[tid]=0; // unresolved by these five words, NEVER a rejection
}
