"""Run bounded, exact fixed-first-entry A9 slices and independently check their mass."""
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile

PROJECT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('campaign',PROJECT/'scripts/weak/run_campaign.py')
runner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
(PROJECT/'runs/experiments/exact-a9-slice-checks').mkdir(parents=True,exist_ok=True)
root=Path(tempfile.mkdtemp(prefix='exact-a9-slices-',dir=PROJECT/'runs/experiments/exact-a9-slice-checks'))
binary=PROJECT/'target/smart/release/hurwitz_breaker'
known={p.name:p for p in (PROJECT/'evidence/certificates').glob('CERTIFIED_*.txt')}
headers=[runner.certificate_header(p) for p in known.values()]
rows=[]
for i,(passport,total) in enumerate(runner.exact_slices(headers)):
    out=root/f'p{i:02d}'
    cmd=[binary,'passport-exact','--group','A','--d','9','--passport',passport,
         '--start','0','--count',str(total),'--threads','1','--seconds','60',
         '--seed','951357','--out',out]
    result=subprocess.run(list(map(str,cmd)),text=True,capture_output=True,check=True)
    (out/'COMMAND.json').write_text(json.dumps(list(map(str,cmd))))
    (out/'RUN.log').write_text(result.stdout+result.stderr)
    summary=dict(line.split(' = ',1) for line in (out/'SUMMARY.txt').read_text().splitlines() if ' = ' in line)
    assert int(summary['scan total candidate pairs'])==total
    certificates=[]
    for certificate in out.glob('CERTIFIED_*.txt'):
        subprocess.run([binary,'verify','--certificate',certificate],check=True)
        certificates.append(runner.certificate_header(certificate))
    # Generating A9 has trivial centralizer in S9. Divide the labelled mass
    # equally among distinct passport orderings, then fix the first entry.
    if summary['scan range resolved']=='true':
        centralizer=362880//runner.class_size(passport.split('/')[0],9)
        types=passport.split('/')
        orderings=math.factorial(4)
        for t in set(types): orderings//=math.factorial(types.count(t))
        numerator=sum(int(h['states']) for h in certificates)*centralizer
        assert numerator%orderings==0
        mass=numerator//orderings
        assert mass==int(summary['accepted generating seeds'])-int(summary['solved seeds']),(mass,summary)
    rows.append({'summary':summary,'certificates':certificates,'directory':str(out)})
    print(passport,'pairs',total,'generating',summary['accepted generating seeds'],'solved',summary['solved seeds'],
          'bad_components',len(certificates),'resolved',summary['scan range resolved'],'seconds',summary['elapsed seconds'],flush=True)
    (root/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
print(root,flush=True)
