try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
import tempfile
import unittest
from build_public_smooth_catalog import release
from package_current_math_release import compact
import test_build_public_smooth_catalog as fixtures


class CompactReleaseTests(unittest.TestCase):
    def make_history(self,root):
        tests=fixtures.PublicationTests();source=tests.source()
        a=root/'v1';b=root/'v2'
        release(source,a,now='2026-10-14T17:00:00Z')
        source['components'][0]['tuple'][0]='(6 7 5)'
        release(source,b,a,now='2026-10-14T18:00:00Z')
        return a,b

    def test_retention_alias_and_incremental_reconstruction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);old,full=self.make_history(root);out=root/'compact'
            m=compact(full,out)
            self.assertEqual(m['minimum_incremental_base_version'],1)
            self.assertFalse(m['versions'][0]['snapshot']['retained'])
            self.assertFalse(m['versions'][0]['changes']['retained'])
            self.assertFalse((out/'versions/1.json').exists())
            self.assertFalse((out/'changes/1.json').exists())
            self.assertEqual((out/'all-smooth.json').read_bytes(),(full/'all-smooth.json').read_bytes())
            self.assertEqual((out/'all-smooth.json').stat().st_ino,(out/'versions/2.json').stat().st_ino)
            self.assertEqual((out/'changes/2.json').read_bytes(),(full/'changes/2.json').read_bytes())
            state={e['id']:e for e in json.loads((old/'all-smooth.json').read_text())['components']}
            change=json.loads((out/'changes/2.json').read_text())
            state.update({e['id']:e for e in change['upserted']})
            for e in change['deleted']:state.pop(e['id'],None)
            self.assertEqual(list(state.values()),json.loads((out/'all-smooth.json').read_text())['components'])
            with self.assertRaises(ValueError):compact(full,out)

    def test_corruption_is_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);_,full=self.make_history(root)
            with (full/'changes/2.json').open('ab') as stream:stream.write(b' ')
            with self.assertRaises(ValueError):compact(full,root/'bad')
            self.assertFalse((root/'bad').exists())

    def test_bootstrap_only_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);_,full=self.make_history(root)
            m=compact(full,root/'bootstrap',changes_from=3)
            self.assertEqual(m['minimum_incremental_base_version'],2)
            self.assertFalse((root/'bootstrap/changes').exists())


if __name__=='__main__':unittest.main()
