#include <metal_stdlib>
using namespace metal;

inline uint permutation_index(thread const uchar *p, uint d) {
    ulong seen=0;uint cycles=0;
    for (uint x=0;x<d;++x) {
        if (seen&(ulong(1)<<x)) continue;
        ++cycles;uint y=x;
        while (!(seen&(ulong(1)<<y))) {seen|=ulong(1)<<y;y=p[y];}
    }
    return d-cycles;
}

inline bool pair_transitive(thread const uchar *a,thread const uchar *b,uint d) {
    const ulong full=(ulong(1)<<d)-1;ulong seen=1,todo=1;
    while (todo) {
        uint x=uint(ctz(todo));todo&=todo-1;
        ulong bit=ulong(1)<<a[x];if (!(seen&bit)) {seen|=bit;todo|=bit;}
        bit=ulong(1)<<b[x];if (!(seen&bit)) {seen|=bit;todo|=bit;}
        if (seen==full) return true;
    }
    return false;
}

inline bool tuple_transitive(thread uchar p[4][34],uint d) {
    const ulong full=(ulong(1)<<d)-1;ulong seen=1,todo=1;
    while (todo) {
        uint x=uint(ctz(todo));todo&=todo-1;
        for (uint j=0;j<3;++j) {
            ulong bit=ulong(1)<<p[j][x];if (!(seen&bit)) {seen|=bit;todo|=bit;}
        }
        if (seen==full) return true;
    }
    return false;
}

inline ulong splitmix_next(thread ulong &state) {
    state+=0x9E3779B97F4A7C15ul;ulong z=state;
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ul;z=(z^(z>>27))*0x94D049BB133111EBul;
    return z^(z>>31);
}

inline bool braid_probe(thread uchar p[4][34],uint d,bool stop_on_witness) {
    ulong state=0xcbf29ce484222325ul;
    for (uint j=0;j<3;++j) for (uint x=0;x<d;++x) state=(state^ulong(p[j][x]))*0x100000001b3ul;
    uchar inv[34],changed[34];
    for (uint step=0;step<64;++step) {
        ulong draw;do {draw=splitmix_next(state);} while (draw>=0xfffffffffffffffcUL);
        uint mv=uint(draw%6),i=mv%3;
        if (mv<3) {
            for (uint x=0;x<d;++x) inv[p[i][x]]=uchar(x);
            for (uint x=0;x<d;++x) changed[x]=p[i][p[i+1][inv[x]]];
            for (uint x=0;x<d;++x) {p[i+1][x]=p[i][x];p[i][x]=changed[x];}
        } else {
            for (uint x=0;x<d;++x) inv[p[i+1][x]]=uchar(x);
            for (uint x=0;x<d;++x) changed[x]=inv[p[i][p[i+1][x]]];
            for (uint x=0;x<d;++x) {p[i][x]=p[i+1][x];p[i+1][x]=changed[x];}
        }
        // Only the middle move changes the first-pair subgroup/transitivity.
        if (stop_on_witness && i==1 && pair_transitive(p[0],p[1],d)) return true;
    }
    return false;
}

kernel void screen(device const uchar *input [[buffer(0)]],
                   device uint *output [[buffer(1)]],
                   constant uint *params [[buffer(2)]],
                   uint tid [[thread_position_in_grid]]) {
    uint d=params[0],n=params[1],mode=params[2];if (tid>=n) return;
    uchar p[4][34];uint ix[4];
    for (uint j=0;j<3;++j) for (uint x=0;x<d;++x) p[j][x]=input[tid*3*d+j*d+x];
    // Product-one fourth entry: (a*b*c)^-1, with the production convention.
    for (uint x=0;x<d;++x) p[3][p[0][p[1][p[2][x]]]]=uchar(x);
    if (mode==3) {
        // Exact, collision-free diagnostic for independent exhaustive S4 oracle:
        // all four Lehmer ranks after 64 moves, five bits per permutation.
        braid_probe(p,d,false);uint packed=0;
        for (uint j=0;j<4;++j) {
            uint rank=0;for (uint x=0;x<d;++x) {uint digit=0;for (uint y=x+1;y<d;++y) digit+=p[j][y]<p[j][x];rank=rank*(d-x)+digit;}
            packed|=rank<<(5*j);
        }
        output[tid]=packed;return;
    }
    uint sum=0;bool allowed=true;
    for (uint j=0;j<4;++j) {ix[j]=permutation_index(p[j],d);sum+=ix[j];allowed&=ix[j]<=d-3;}
    if (mode>=1 && (!allowed || sum<2*d-2)) {output[tid]=0;return;}
    if (mode==2 && !tuple_transitive(p,d)) {output[tid]=0;return;}
    uint mask=0,bit=0;
    for (uint i=0;i<4;++i) for (uint j=i+1;j<4;++j) {
        bool transitive=pair_transitive(p[i],p[j],d);
        if (mode>=1 && transitive) {output[tid]=0;return;}
        if (transitive) mask|=1<<bit;++bit;
    }
    if (mode==2 && braid_probe(p,d,true)) {output[tid]=0;return;}
    output[tid]=mode>=1 ? 1 : ix[0]|(ix[1]<<6)|(ix[2]<<12)|(ix[3]<<18)|(mask<<24);
}
