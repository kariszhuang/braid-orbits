"""Actual producer/auditor controls; fresh disposable folders, no pool changes."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from audit_guided_native import audit_guided_native
from verify_smooth_geometry import parse_tuple, compose, cycles
from verify_genus_three_adaptive import verify_adaptive


@unittest.skipUnless(os.environ.get('G3_REFACTOR_BUILD'), 'requires validated frozen refactor build')
class RefactorTests(unittest.TestCase):
    def setUp(self):
        self.build=Path(os.environ['G3_REFACTOR_BUILD'])
        self.root=Path(os.environ['G3_REFACTOR_ROOT'])

    def test_actual_low_index_refactors_and_independent_witnesses(self):
        with tempfile.TemporaryDirectory()as tmp:
            for n in (18,21,24,25):
                out=Path(tmp)/str(n)
                subprocess.run([str(self.build/'bin/compact'),'--d',str(n),'--lane','refactor',
                    '--mutation-seeds',str(self.root/'MUTATION_SEEDS.tsv'),'--refactor-max-index','3',
                    '--seed','329123','--attempts','3000','--seconds','60','--closure-seconds','10',
                    '--cap','250000','--trace-early','32','--words',str(self.root/'banks/words_00000000.txt'),
                    '--out',str(out)],check=True,capture_output=True,text=True)
                audited=audit_guided_native(out,self.build/'bin/auditor',native_components=True)
                self.assertGreater(audited['independent_native']['records'],0)
                with (out/'RESOLUTIONS.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
                for row in rows:
                    t=parse_tuple(row['seed'],n)
                    ix=[n-len(cycles(p))for p in t]
                    self.assertIn(ix[0],(1,2,3));self.assertEqual(sum(ix),2*n+4)
                    self.assertEqual(compose(t[0],t[0]),tuple(range(n)))
                    for j in (1,2,3):self.assertEqual(n-len(cycles(compose(t[0],t[j]))),ix[0]+ix[j])

    def test_bounded_trie_trace_preserves_proposals_and_rng(self):
        with tempfile.TemporaryDirectory()as tmp:
            counts=[];stages=[];persisted=[]
            for bound in (0,64):
                out=Path(tmp)/str(bound)
                subprocess.run([str(self.build/'bin/wide'),'--d','60','--lane','four-quotient',
                    '--seed','998201','--attempts','250000','--seconds','60','--closure-seconds','10',
                    '--cap','250000','--literal-trie','--trace-trie',str(bound),
                    '--words',str(self.root/'banks/words_00000000.txt'),'--out',str(out)],
                    check=True,capture_output=True,text=True)
                with (out/'COUNTS.tsv').open()as stream:count=next(csv.DictReader(stream,delimiter='\t'))
                count.pop('seconds');counts.append(count)
                with (out/'STAGES.tsv').open()as stream:stages.append(next(csv.DictReader(stream,delimiter='\t')))
                checked=audit_guided_native(out,self.build/'bin/auditor',native_components=True)
                persisted.append(checked['independent_native']['records'])
            self.assertEqual(counts[0],counts[1]);self.assertEqual(stages[0],stages[1])
            self.assertGreater(persisted[1],persisted[0])

    def test_actual_adaptive_templates_have_literal_negative_prefix_provenance(self):
        if not os.environ.get('G3_ADAPTIVE_BUILD'):self.skipTest('requires adaptive build')
        with tempfile.TemporaryDirectory()as tmp:
            out=Path(tmp)/'adaptive'
            words=Path(tmp)/'words';words.write_text('')
            subprocess.run([str(self.build/'bin/compact'),'--d','18','--lane','refactor',
                '--mutation-seeds',str(self.root/'MUTATION_SEEDS.tsv'),'--refactor-max-index','3',
                '--adaptive-refactor','4','--seed','911341','--attempts','100000','--seconds','60',
                '--closure-seconds','10','--cap','250000','--words',str(words),
                '--out',str(out)],check=True,capture_output=True,text=True)
            audit_guided_native(out,self.build/'bin/auditor',native_components=True)
            checked=verify_adaptive(out)
            self.assertTrue(checked['passed']);self.assertGreater(checked['templates'],0)
            self.assertLessEqual(checked['templates'],4)
            self.assertGreater(checked['adaptive_frame_draws'],0)


if __name__=='__main__':unittest.main()
