import json,sys
def load(f):
    L=[json.loads(l) for l in open(f)]
    s=[o for o in L if o.get('summary')]; o=[o for o in L if not o.get('summary')]
    return sorted((x['states'],x['clusters'],x['bad_clusters'],x['lambda12_LL'],x['comb_smooth_states'],tuple(sorted(x['types'])),tuple(x['delta_LL']),tuple(x['bad_cluster_sizes']),x['sigma1_orbits'],x['bad_sigma1_orbits'],x['defect_states'],x['primitive'],tuple(x['block_sizes'])) for x in o), s[0]
a,sa=load(sys.argv[1]); b,sb=load(sys.argv[2])
print(sys.argv[2].split('/')[-1], 'orbits', len(a), len(b), 'MATCH' if a==b else 'MISMATCH', 'time old %.1fs new %.1fs'%(sa['seconds'],sb['seconds']))
