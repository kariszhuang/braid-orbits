"""Exact small independent comparisons and evidence-corruption regressions."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import permutations
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from verify_genus_three_chords import admissible_partitions,audit,cycle_permutations,frames
from verify_genus_three_finite import mul,profile


class ChordTests(unittest.TestCase):
    def test_scalar_degree_eleven_is_complete(self):
        budgets=[]
        for a in range(1,11):
            for b in range(a,11):
                for c in range(b,11):
                    e=26-a-b-c
                    if c<=e<=10 and a+c<=10 and min(a+e,b+c)<=10:budgets.append((a,b,c,e))
        self.assertEqual(budgets,[(1,7,9,9),(1,8,8,9),(2,8,8,8)])

    def test_small_colouring_criterion_against_all_permutations(self):
        n=7
        for mode in ('transposition','double'):
            f=frames(n,mode);self.assertEqual(len(f),3 if mode=='transposition' else 10)
            for a in f:
                built=[]
                for blocks in admissible_partitions(a):built.extend(cycle_permutations(blocks,n))
                self.assertEqual(len(built),len(set(built)))
                brute=set()
                for b in permutations(range(n)):
                    if len(profile(b))!=3:continue
                    shifted=[]
                    for k in range(n):
                        u=[0]*n
                        for x in range(n):u[(x+k)%n]=(a[x]+k)%n
                        shifted.append(tuple(u))
                    rank=1 if mode=='transposition' else 2
                    if all(n-len(profile(mul(u,b)))==rank+n-3 for u in shifted):brute.add(b)
                self.assertEqual(set(built),brute)

    def test_single_cycle_balancing_small_literal_cases(self):
        for n,k in ((6,3),(7,3),(7,2)):
            a=tuple(list(range(1,k))+[0]+list(range(k,n)))
            count=0
            for b in permutations(range(n)):
                if len(profile(b))!=k:continue
                good=True
                for shift in range(n):
                    u=[0]*n
                    for x in range(n):u[(x+shift)%n]=(a[x]+shift)%n
                    if profile(mul(tuple(u),b))!=(n,):good=False;break
                if good:
                    count+=1;self.assertEqual(n%k,0);self.assertEqual(profile(b),(n//k,)*k)
            if n%k:self.assertEqual(count,0)

    @unittest.skipUnless(os.environ.get('G3_CHORD_BINARY')and os.environ.get('G3_CHORD_AUDITOR'),'explicit frozen native fixtures required')
    def test_complete_native_control_and_corruptions(self):
        with tempfile.TemporaryDirectory()as tmp:
            base=Path(tmp)/'base';binary=os.environ['G3_CHORD_BINARY'];native=Path(os.environ['G3_CHORD_AUDITOR'])
            subprocess.run([binary,'7',str(base),'double'],check=True,capture_output=True,text=True)
            self.assertEqual(audit(base,native)['status'],'complete_negative_slice')
            for kind in ('missing','false_word','frame_omission','count','duplicate'):
                folder=Path(tmp)/kind;shutil.copytree(base,folder)
                if kind in ('frame_omission','count'):
                    c=json.loads((folder/'COVERAGE.json').read_text())
                    if kind=='frame_omission':c['frames'].pop()
                    else:c['constructed']+=1
                    (folder/'COVERAGE.json').write_text(json.dumps(c))
                else:
                    rows=(folder/'RESOLUTIONS.tsv').read_text().splitlines();leads=(folder/'LEADS.tsv').read_text().splitlines()
                    if kind=='missing':rows.pop();leads.pop()
                    elif kind=='false_word':
                        fields=rows[1].split('\t');fields[2]='0';rows[1]='\t'.join(fields)
                    else:
                        fields=rows[-1].split('\t');first=rows[1].split('\t');fields[2]=first[2];fields[4]=first[4];rows[-1]='\t'.join(fields);leads[-1]=fields[0]+'\t'+fields[4]
                    (folder/'RESOLUTIONS.tsv').write_text('\n'.join(rows)+'\n');(folder/'LEADS.tsv').write_text('\n'.join(leads)+'\n')
                with self.subTest(kind=kind),self.assertRaises((AssertionError,subprocess.CalledProcessError,ValueError,KeyError)):
                    audit(folder,native)


if __name__=='__main__':unittest.main()
