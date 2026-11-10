"""C251 source-bound admission and independent incomplete-FRAME accounting.

The third/fourth entries do not exist at this stage. No tuple, negative
word, transitive start or isotrivial component is fabricated from a frame.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from hashlib import sha256
import json
from pathlib import Path

from verify_smooth_geometry import parse_tuple,cycles,compose

REASONS=(
    'C243-degree17-single3-saturated-no-all-smooth',
    'C247-degree19-single3-saturated-no-all-smooth',
    'C244-degree23-single4-saturated-no-all-smooth',
    'C241-degree14-single3-middle-profile-no-all-smooth',
    'C239-k-above5-Hodge-negative-no-all-smooth',
    'C239-k5-Hodge-zero-nonisotrivial-target-excluded',
)


def sha(path):
    h=sha256()
    with Path(path).open('rb')as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def admit(path,project):
    path=Path(path).resolve();project=Path(project).resolve();value=json.loads(path.read_text())
    if not(value.get('schema')==1 and value.get('claim')=='C251_proposal_admission'and value.get('passed')is True
           and value.get('theorem')=='retire-completed-or-objective-excluded-single-cycle-construction-frames'
           and value.get('reviewed_by_parent')is True and value.get('default_activation')is False
           and value.get('performance_validation_required')is True):raise ValueError('frame retirement admission scope')
    target=dict(source_genus=3,nonisotrivial=True,all_smooth=True,orbit='FULLsixmoves modulo simultaneousSdconjugacy',transitive=True,product_one=True)
    if value.get('production_target')!=target:raise ValueError('frame retirement target must be nonisotrivial genus3')
    expected=[dict(degree=14,single_cycle=3,sorted_indices=[2,9,10,11],dependency='C241'),
              dict(degree=17,single_cycle=3,ordinary_index='d-k',dependency='C243'),
              dict(degree=19,single_cycle=3,ordinary_index='d-k',dependency='C247'),
              dict(degree=23,single_cycle=4,ordinary_index='d-k',dependency='C244')]
    if value.get('absolute_negative_rules')!=expected:raise ValueError('frame retirement finite scopes differ')
    only=value.get('objective_only_rule',{})
    if not(only.get('single_cycle')=='k>=5'and only.get('ordinary_index')=='d-k'and only.get('star_required')is True
           and only.get('dependency')=='C239'and only.get('conclusion')=='no NON-ISOTRIVIAL all-smooth completion; not classification of this incomplete frame as an actual isotrivial component'):
        raise ValueError('frame Hodge rule must remain objective-only')
    bindings=value.get('bindings',[]);required={'C239','C241','C243','C244','C247'}
    if len(bindings)!=5 or {Path(item['path']).stem for item in bindings}!={claim+'_ADMISSION'for claim in required}:
        raise ValueError('frame retirement prerequisite admissions missing')
    checked={};expanded=set();admissions=[]
    def bound_file(leaf,digest):
        leaf=Path(leaf).resolve()
        if leaf in checked and checked[leaf]!=digest:raise ValueError('conflicting recursive proof bindings: '+str(leaf)+' '+checked[leaf]+' versus '+digest)
        if leaf not in checked:
            if sha(leaf)!=digest:raise ValueError('recursive frame proof binding changed: '+str(leaf))
            checked[leaf]=digest
        if leaf.suffix=='.json'and leaf not in expanded:
            expanded.add(leaf);walk(json.loads(leaf.read_text()),leaf.parent)

    def walk(node,context=project):
        if isinstance(node,list):
            for child in node:walk(child,context)
        elif isinstance(node,dict):
            if isinstance(node.get('path'),str)and isinstance(node.get('sha256'),str):
                leaf=Path(node['path']);leaf=(leaf if leaf.is_absolute()else project/leaf).resolve();digest=node['sha256']
                bound_file(leaf,digest)
            for key in('bindings','raw_bindings'):
                mapping=node.get(key)
                if isinstance(mapping,dict):
                    for name,digest in mapping.items():
                        if not isinstance(digest,str)or len(digest)!=64:continue
                        raw=Path(name);candidates=[raw]if raw.is_absolute()else[context/raw,project/raw]
                        if not raw.is_absolute():
                            if node.get('producer_root'):candidates.insert(0,Path(node['producer_root'])/raw)
                            # Audit INPUTS/FINISHED maps use producer-relative
                            # names, whereas audit PLAN maps use audit-relative
                            # names. Resolve this ambiguity by the EXPECTED
                            # content hash, never by existence alone.
                            for sibling in('PLAN.json','INPUTS.json'):
                                file=context/sibling
                                if file.is_file():
                                    metadata=json.loads(file.read_text())
                                    if metadata.get('producer_root'):candidates.insert(0,Path(metadata['producer_root'])/raw)
                        matched=False
                        for leaf in dict.fromkeys(p.resolve()for p in candidates):
                            if not leaf.is_file():continue
                            if leaf in checked:
                                if checked[leaf]!=digest:continue
                            elif sha(leaf)!=digest:continue
                            else:checked[leaf]=digest
                            bound_file(leaf,digest);matched=True;break
                        if not matched:raise ValueError('recursive named proof binding changed: '+str(name))
            for key,child in node.items():
                if key not in('path','sha256'):walk(child,context)
    for item in bindings:
        child=Path(item['path']);child=child if child.is_absolute()else project/child
        record=json.loads(child.read_text());admissions.append(record.get('claim'))
        if record.get('passed')is not True:raise ValueError('unfinished child frame admission')
    if set(admissions)!=required:raise ValueError('frame dependency claim mismatch')
    walk(bindings)
    return dict(path=str(path),sha256=sha(path),**value,
        recursively_checked_files=[dict(path=str(p),sha256=digest)for p,digest in sorted(checked.items())])


def decision(n,pair,profile):
    actual=[n-len(cycles(p))for p in pair]
    if any(x==0 for x in profile)or sum(profile)!=2*n+4 or actual!=profile[:2]:return None
    active=[len(c)for c in cycles(pair[0])if len(c)>1]
    if len(active)!=1:return None
    k=active[0]
    if not 3<=k<n or not any(i==n-k for i in profile[1:]):return None
    if any(profile[i]+profile[j]<=n-1 for i in range(1,4)for j in range(i+1,4)):return None
    if(n,k)==(17,3):reason=0
    elif(n,k)==(19,3):reason=1
    elif(n,k)==(23,4):reason=2
    elif(n,k)==(14,3)and sorted(profile)==[2,9,10,11]:reason=3
    elif k>5:reason=4
    elif k==5:reason=5
    else:return None
    return REASONS[reason],k


def verify(folder):
    folder=Path(folder)
    def rows(name):
        with(folder/name).open()as stream:return list(csv.DictReader(stream,delimiter='\t'))
    summary=rows('FRAME_RETIREMENT_SUMMARY.tsv');counts=rows('FRAME_RETIREMENT_REASONS.tsv');samples=rows('FRAME_RETIREMENT_SAMPLES.tsv')
    if len(summary)!=1 or {r['reason']for r in counts}!=set(REASONS)or len(counts)!=6:raise ValueError('frame summary shape')
    s=summary[0];enabled=s['enabled']=='true'
    if s['enabled']not in('true','false'):raise ValueError('frame enabled flag')
    mass=sum(int(r['target_excluded_frames'])for r in counts)
    if not 0<=mass<=int(s['checked_new_frames'])or mass!=int(s['target_excluded_frames']):raise ValueError('frame counter mass')
    if not len(samples)==int(s['sampled_frames'])<=min(48,mass):raise ValueError('frame sample mass')
    if not enabled and(mass or int(s['checked_new_frames'])):raise ValueError('OFF policy changed frame counts')
    with(folder/'COUNTS.tsv').open()as stream:n=int(next(csv.DictReader(stream,delimiter='\t'))['degree'])
    for j,row in enumerate(samples):
        t=parse_tuple(row['partial_pair_identity_padded_NOT_A_CANDIDATE'],n)
        if int(row['sample'])!=j or t[2:]!=(tuple(range(n)),)*2:raise ValueError('partial frame serialization')
        profile=list(map(int,row['profile'].split(',')))
        if n-len(cycles(compose(t[0],t[1])))!=sum(profile[:2]):raise ValueError('frame pair not an additive construction')
        if decision(n,t[:2],profile)!=(row['reason'],int(row['k'])):raise ValueError('unsupported frame target exclusion')
    return dict(passed=True,enabled=enabled,constructor=s['constructor'],checked_new_frames=int(s['checked_new_frames']),
        target_excluded_frames=mass,sampled_frames=len(samples),
        qualification='incomplete construction FRAME exclusions only; NOT tuples/literal negatives/transitive starts/isotrivial components',
        bindings={name:sha(folder/name)for name in('FRAME_RETIREMENT_SUMMARY.tsv','FRAME_RETIREMENT_REASONS.tsv','FRAME_RETIREMENT_SAMPLES.tsv','FRAME_RETIREMENT_SCOPE.txt')})
