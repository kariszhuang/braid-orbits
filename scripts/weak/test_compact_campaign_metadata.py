try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
import struct
import tempfile
import unittest
from compact_campaign_metadata import metadata_plan
from compact_smooth_outputs import compact,restore,sha


class MetadataTests(unittest.TestCase):
    def fixture(self,root):
        p=root/'jobs/job/attempt_000';p.mkdir(parents=True)
        (p.parent/'DONE').write_text('done');(p/'VERIFIED').write_text('checked')
        (p/'SUMMARY.txt').write_text('scientific evidence')
        (p/'._SUMMARY.txt').write_bytes(struct.pack('>II',0x00051607,0x00020000)+bytes(18)+b'real metadata')
        return p.resolve()

    def test_roundtrip_preserves_scientific_file_and_resource_fork(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'root';p=self.fixture(root);data=sha(p/'SUMMARY.txt');metadata=sha(p/'._SUMMARY.txt')
            files,report=metadata_plan(root)
            self.assertIn(p/'._SUMMARY.txt',files)
            self.assertEqual(report['metadata_files'],len(files))
            out=Path(tmp)/'archive';compact(root,out,files,[])
            self.assertEqual(sha(p/'SUMMARY.txt'),data);self.assertFalse((p/'._SUMMARY.txt').exists())
            restore(out/'REPEATED_OUTPUTS.tar.gz',root);self.assertEqual(sha(p/'._SUMMARY.txt'),metadata)

    def test_protected_partial_foreign_and_symlink_files_stay(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=self.fixture(root)
            self.assertNotIn(p/'._SUMMARY.txt',metadata_plan(root,[p/'SUMMARY.txt'])[0])
            (p/'hard_candidates.txt').write_text('saved lead');self.assertEqual(metadata_plan(root)[0],[])
            (p/'hard_candidates.txt').write_text('');(p/'._SUMMARY.txt').write_bytes(b'not AppleDouble')
            self.assertNotIn(p/'._SUMMARY.txt',metadata_plan(root)[0])
            (p/'._SUMMARY.txt').unlink();(p/'._SUMMARY.txt').symlink_to(p/'SUMMARY.txt')
            self.assertNotIn(p/'._SUMMARY.txt',metadata_plan(root)[0])

    def test_unverified_or_unfinished_job_stays(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=self.fixture(root);(p/'VERIFIED').unlink();self.assertEqual(metadata_plan(root)[0],[])
            (p/'VERIFIED').write_text('checked');(p.parent/'DONE').unlink();self.assertEqual(metadata_plan(root)[0],[])

    def test_jobs_symlink_is_not_followed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'root';outside=Path(tmp)/'outside';outside.mkdir();root.mkdir()
            (root/'jobs').symlink_to(outside,target_is_directory=True)
            with self.assertRaises(ValueError):metadata_plan(root)


if __name__=='__main__':unittest.main()
