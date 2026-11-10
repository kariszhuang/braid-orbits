#!/usr/bin/env python3
"""Differential literal mathematics plus malformed-evidence rejection controls."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools, json, random, struct, subprocess, tempfile, unittest
from pathlib import Path
from odd_order_spin import hurwitz, inverse, cycle_lengths
from analyze_fixed_braid_constraints import product, point_orbits
from run_smooth_mission import values
from run_smooth_degree9 import rows

PROJECT=Path(__file__).resolve().parents[2]
BINARY=PROJECT/'target/audit-fast/release/smooth_negative_audit'

class Server:
    def __init__(self,binary,n,literal=False):
        self.p=subprocess.Popen([str(binary),'--literal-server' if literal else '--server',str(n)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if self.p.stdout.readline().strip()!='READY':raise ValueError('startup')
    def ask(self,fields):
        self.p.stdin.write('\t'.join(map(str,fields))+'\n');self.p.stdin.flush()
        line=self.p.stdout.readline()
        if not line:raise ValueError(self.p.stderr.read())
        return json.loads(line)
    def close(self):
        self.p.stdin.close();self.p.wait(timeout=10);self.p.stdout.close();self.p.stderr.close()

class Checks(unittest.TestCase):
    def test_resume_and_tampered_frame_report(self):
        from run_negative_audit import audit
        source=Path('data/smooth_frontier_20261012_123427_295724/validation/focused')
        with tempfile.TemporaryDirectory(dir='data',prefix='audit-resume-test-') as temp:
            root=Path(temp);first=audit(root,source,BINARY,1);second=audit(root,source,BINARY,1)
            self.assertEqual(first['coverage'],second['coverage'])
            report=root/'audit_frames/frame_00000.json';data=json.loads(report.read_text());data['count'][0]+=1;report.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'coverage mismatch'):audit(root,source,BINARY,1)

    def test_all_product_one_s4_and_larger_words(self):
        rng=random.Random(19);number=0
        for n in [4,5,9,10]:
            server=Server(BINARY,n,True)
            try:
                perms=list(itertools.permutations(range(n))) if n==4 else None
                starts=itertools.product(perms,repeat=3) if n==4 else (tuple(tuple(rng.sample(range(n),n)) for _ in range(3)) for _ in range(1000))
                for a,b,c in starts:
                    t=(a,b,c,inverse(product(a,b,c))); word=[rng.randrange(6) for _ in range(rng.randrange(8))];moved=t
                    for m in word:moved=hurwitz(moved,m)
                    got=server.ask([*(','.join(map(str,p)) for p in t),','.join(map(str,word))])
                    self.assertEqual(got['tuple'],list(map(list,moved)))
                    self.assertEqual(got['indices'],[n-len(cycle_lengths(p)) for p in moved])
                    self.assertEqual(got['pair_orbits'],[len(point_orbits(moved[i:i+2])) for i in range(3)])
                    self.assertEqual(got['transitive'],len(point_orbits(moved))==1)
                    self.assertEqual(got['smooth'],values(moved)['smooth']);number+=1
            finally:server.close()
        print('LITERAL_DIFFERENTIAL',number)

    def test_corrupt_and_unsupported_fail_closed(self):
        source=Path('data/smooth_frontier_20261012_123427_295724/degree_10')
        from odd_order_spin import parse_seed
        row=rows(source/'FRAMES.tsv')[0];a,b,_,_=parse_seed(f"[{row['a']}, {row['b']}, (), ()]",10)
        original=source/'frames/frame_00000'
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            with (original/'OUTCOMES.bin').open('rb') as stream:raw=stream.read(8*10000)
            cases=[(raw[:7],(original/'WITNESSES.tsv').read_text()),
                   (struct.pack('<II',0,0)+raw[8:],(original/'WITNESSES.tsv').read_text()),
                   (raw,'token\tstatus\tword_or_component\n0\tclosed_smooth\tcomponent_0000\n'),
                   (raw,'token\tstatus\tword_or_component\n0\tnonsmooth\t6\n'),
                   (raw,'token\tstatus\tword_or_component\n0\tnonsmooth\t0,3\n')]
            for data,words in cases:
                (folder/'OUTCOMES.bin').write_bytes(data);(folder/'WITNESSES.tsv').write_text(words)
                server=Server(BINARY,10)
                try:
                    with self.assertRaises(ValueError):server.ask([folder,','.join(map(str,a)),','.join(map(str,b)),1,10000])
                finally:server.close()
            # Full population mode must reject even one trailing byte, including
            # a frame with no eligible candidates (identity first entries).
            (folder/'OUTCOMES.bin').write_bytes(b'x');(folder/'WITNESSES.tsv').write_text('token\tstatus\tword_or_component\n')
            server=Server(BINARY,5)
            try:
                with self.assertRaisesRegex(ValueError,'extra outcome'):server.ask([folder,'0,1,2,3,4','0,1,2,3,4',1,0])
            finally:server.close()

if __name__=='__main__':unittest.main()
