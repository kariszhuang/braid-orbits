#include <metal_stdlib>
using namespace metal;

// Rows: a,b,b^-1,c,c^-1,ia,ib,ic,id,expected. The seed has already
// passed the five original smoothness filters. A geodesic first pair gives
// the unique positive-genus target ia+ib for the two middle squares.
inline bool square_matches(device const uchar *r,uint d,bool forward) {
    device const uchar *a=r,*b=r+d,*bi=r+2*d,*c=r+3*d,*ci=r+4*d;
    uint expected=d-r[5*d]-r[5*d+1],cycles=0,remaining=(1u<<d)-1;
    while(remaining) {
        if(++cycles>expected)return false;
        uint x=ctz(remaining);
        while(remaining&(1u<<x)) {
            remaining&=~(1u<<x);
            uint y=forward?b[c[b[ci[bi[x]]]]]:ci[b[c[x]]];
            x=a[y];
        }
        if(popcount(remaining)<expected-cycles)return false;
    }
    return cycles==expected;
}
kernel void forest_middle_screen(device const uchar *input [[buffer(0)]],
                                 device uchar *output [[buffer(1)]],
                                 constant uint *params [[buffer(2)]],
                                 uint tid [[thread_position_in_grid]]) {
    uint d=params[0],count=params[1];if(tid>=count)return;
    device const uchar *r=input+tid*(5*d+5);
    output[tid]=!square_matches(r,d,true)?1:!square_matches(r,d,false)?2:0;
}
