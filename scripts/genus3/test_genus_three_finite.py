"""Independent count, kernel and evidence-corruption boundaries."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import verify_genus_three_finite as v


class FiniteTests(unittest.TestCase):
    def test_character_counts_match_literal_small_groups(self):
        for n in range(2, 8):
            target = tuple((i+1) % n for i in range(n))
            counts = {}
            profiles = set()
            for a in itertools.permutations(range(n)):
                mu = v.profile(a)
                nu = v.profile(v.mul(v.inv(a), target))
                profiles.add(mu)
                counts[mu, nu] = counts.get((mu, nu), 0) + 1
            for mu in profiles:
                for nu in profiles:
                    self.assertEqual(v.connection(mu, nu), counts.get((mu, nu), 0))

    def test_target_counts_and_tree_frame_counts(self):
        self.assertEqual(v.connection(v.A, v.A), 2800)
        self.assertEqual(v.connection(v.A, v.B), 63000)
        self.assertEqual(v.connection(v.B, v.B), 1428000)
        self.assertEqual(v.connection(v.B, (2,2,3,4,4)), 2425500)
        for p, frames in [(v.A, 2), (v.B, 40), ((2,2,3,4,4),60),
                          ((2,2,3,3,5),60), ((2,2,2,4,5),40), ((2,2,2,3,6),40)]:
            self.assertEqual(v.connection(v.SPECIAL, p), frames*15)

    def test_all_braid_inverses_and_product(self):
        a = tuple((i+1) % 15 for i in range(15))
        b = tuple(reversed(range(15)))
        t = (a, b, v.inv(a), v.inv(v.mul(v.mul(a,b),v.inv(a))))
        for m in range(6):
            u = v.braid(t,m)
            self.assertEqual(v.braid(u,(m+3)%6),t)
            self.assertEqual(v.mul(v.mul(u[0],u[1]),v.mul(u[2],u[3])),v.IDENTITY)

    def test_empty_frame_table_cannot_claim_complete(self):
        # A tiny malformed compact certificate must never be called complete.
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw)
            (path/'FRAMES.tsv').write_text('frame\ta\tb\n')
            summary = dict(schema=2,first_ordinary_profile=[3]*5,frames=0,completions=0)
            with self.assertRaises(AssertionError):
                v.audit_slice(path,summary)


@unittest.skipUnless(os.environ.get('GENUS_THREE_FINITE_BIN'), 'set GENUS_THREE_FINITE_BIN for artifact mutation controls')
class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tmp.name)
        cls.good = cls.root/'good'
        cls.binary = os.environ['GENUS_THREE_FINITE_BIN']
        subprocess.run([cls.binary, str(cls.good), '--slice', 'A'], check=True, capture_output=True)
        assert v.audit(cls.good)['status'] == 'verified_complete_negative_passport'

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def altered(self):
        path = Path(tempfile.mkdtemp(dir=self.root))/'certificate'
        shutil.copytree(self.good,path)
        return path

    def test_truncated_and_invalid_tokens_rejected(self):
        path = self.altered()
        data = (path/'WITNESSES.u16le').read_bytes()
        (path/'WITNESSES.u16le').write_bytes(data[:-1])
        with self.assertRaises(AssertionError):
            v.audit(path)
        (path/'WITNESSES.u16le').write_bytes(b'\xff\xff'+data[2:])
        with self.assertRaisesRegex(AssertionError,'dictionary token'):
            v.audit(path)

    def test_false_negative_word_rejected(self):
        path = self.altered()
        # q1 preserves the first-pair product; it cannot be a negative witness.
        lines = (path/'WORDS.tsv').read_text().splitlines()
        lines[1] = '1\t0'
        (path/'WORDS.tsv').write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(AssertionError,'false negative'):
            v.audit(path)

    def test_invalid_and_duplicate_cached_completions_rejected(self):
        path = self.altered()
        file = path/'COMPLETIONS.tsv'
        lines = file.read_text().splitlines()
        parts = lines[1].split('\t')
        parts[1] = '255,' + parts[1].split(',',1)[1]
        lines[1] = '\t'.join(parts)
        file.write_text('\n'.join(lines)+'\n')
        result = subprocess.run([self.binary,str(path/'invalid'),'--slice','A','--library',str(file)],capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn(b'library images',result.stderr)
        lines = (self.good/'COMPLETIONS.tsv').read_text().splitlines()
        lines[2] = '1\t'+lines[1].split('\t')[1]
        file.write_text('\n'.join(lines)+'\n')
        result = subprocess.run([self.binary,str(path/'duplicate'),'--slice','A','--library',str(file)],capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn(b'library coverage/duplicates',result.stderr)

    def test_cache_and_constructor_produce_identical_certificates(self):
        out = self.root/'cached'
        subprocess.run([self.binary,str(out),'--slice','A','--library',str(self.good/'COMPLETIONS.tsv')],check=True,capture_output=True)
        for name in ('FRAMES.tsv','COMPLETIONS.tsv','WORDS.tsv','WITNESSES.u16le'):
            self.assertEqual((out/name).read_bytes(),(self.good/name).read_bytes())
        self.assertEqual(v.audit(out)['status'],'verified_complete_negative_passport')

    def test_existing_directory_is_never_overwritten(self):
        result = subprocess.run([self.binary,str(self.good),'--slice','A'],capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(v.audit(self.good)['status'],'verified_complete_negative_passport')

    @unittest.skipUnless(os.environ.get('GENUS_THREE_FINITE_AUDIT_BIN'),'native auditor not selected')
    def test_native_matches_python_and_rejects_corrupt_witnesses(self):
        native=Path(os.environ['GENUS_THREE_FINITE_AUDIT_BIN'])
        python=v.audit(self.good);result=v.audit(self.good,native,4)
        self.assertEqual(python['negative'],result['negative'])
        self.assertEqual({tuple(x['word']):x['states'] for x in python['word_histogram']},
                         {tuple(x['word']):x['states'] for x in result['word_histogram']})
        for kind in ('false_word','bad_token','bad_images'):
            path=self.altered()
            if kind=='false_word':
                lines=(path/'WORDS.tsv').read_text().splitlines();lines[1]='1\t0'
                (path/'WORDS.tsv').write_text('\n'.join(lines)+'\n')
            elif kind=='bad_token':
                data=(path/'WITNESSES.u16le').read_bytes()
                (path/'WITNESSES.u16le').write_bytes(b'\xff\xff'+data[2:])
            else:
                text=(path/'COMPLETIONS.tsv').read_text()
                (path/'COMPLETIONS.tsv').write_text(text.replace('0\t','0\t255,',1))
            with self.assertRaises(AssertionError):v.audit(path,native,4)


if __name__ == '__main__':
    unittest.main()
