try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from export_genus_three_public import export, validate_folder


class GenusThreeExportTests(unittest.TestCase):
    def test_certificate_path_is_directly_in_jobs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / 'jobs/job_1/component_1'
            self.assertEqual(validate_folder(root, dict(raw_folder=str(folder))), folder.resolve())
            for path in [root/'component_1', root/'other/job_1/component_1', root/'jobs/job_1/../component_1']:
                with self.assertRaises(ValueError):
                    validate_folder(root, dict(raw_folder=str(path)))

    def test_active_campaign_is_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'STATE.json').write_text(json.dumps(dict(status='searching', active_jobs=[])))
            with self.assertRaisesRegex(ValueError, 'inactive'):
                export(root, root/'previous.json', root/'out.json')

    def test_inactive_empty_source_preserves_previous(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'STATE.json').write_text(json.dumps(dict(status='stopped', active_jobs=[])))
            (root/'NEW_COMPONENTS.json').write_text('[]')
            previous = dict(component_count=1, components=[dict(id='old', states=12, revision=3)])
            (root/'previous.json').write_text(json.dumps(previous))
            checks = export(root, root/'previous.json', root/'out.json')
            self.assertEqual(checks['added'], 0)
            self.assertEqual(json.loads((root/'out.json').read_text()), previous)
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                export(root, root/'previous.json', root/'out.json')


if __name__ == '__main__':
    unittest.main()
