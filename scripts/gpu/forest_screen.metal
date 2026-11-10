#include <metal_stdlib>
using namespace metal;

// Separate smoothness predicate, NOT the historical weak-search screen.
inline uint composed_index(device const uchar *a, device const uchar *b, uint d) {
    ulong seen=0;uint cycles=0;
    for(uint root=0;root<d;++root){
        if(seen&(ulong(1)<<root))continue;
        ++cycles;uint x=root;
        while(!(seen&(ulong(1)<<x))){seen|=ulong(1)<<x;x=a[b[x]];}
    }
    return d-cycles;
}
kernel void forest_screen(device const uchar *input [[buffer(0)]],
                          device uchar *output [[buffer(1)]],
                          constant uint *params [[buffer(2)]],
                          uint tid [[thread_position_in_grid]]) {
    uint d=params[0],count=params[1];if(tid>=count)return;
    device const uchar *row=input+tid*(6*d+5);
    device const uchar *c=row+5*d;
    uint ia=row[6*d],ib=row[6*d+1],ic=row[6*d+2],id=row[6*d+3];
    for(uint j=0;j<5;++j){
        uint left=j<2?ia+ic:ib+ic,right=j<2?ib+id:ia+id;
        if(left>=d&&right>=d){output[tid]=uchar(j+1);return;}
        uint ip=composed_index(row+j*d,c,d);
        if(ip!=left&&ip!=right){output[tid]=uchar(j+1);return;}
    }
    output[tid]=0;
}
