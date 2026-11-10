"""Real launcher admission under a launchd-style minimal PATH; no service."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import fcntl
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_forest_fast as runner
from run_smooth_degree9 import sha


class RuntimeTests(unittest.TestCase):
    def test_only_verified_native_zero_pending_guided_history_can_skip_journals(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder=root/'guided_0_0000000_18_returns';folder.mkdir()
            checks=dict(construction=dict(unresolved='0'),independent_native=dict(unresolved=0))
            self.assertFalse(runner.resolved_native_history(folder,checks))
            (folder/'VERIFIED').write_text('independent audit completed')
            self.assertTrue(runner.resolved_native_history(folder,checks))
            for changed in [dict(construction=dict(unresolved='1'),independent_native=dict(unresolved=1)),
                            dict(construction=dict(unresolved='0'),independent_native=dict(unresolved=1)),
                            dict(construction=dict(unresolved='0'))]:
                self.assertFalse(runner.resolved_native_history(folder,changed))
            replay=root/'replay_000000';replay.mkdir();(replay/'VERIFIED').write_text('done')
            self.assertFalse(runner.resolved_native_history(replay,checks))

    def test_isolated_compilation_preserves_dirty_workspace_and_binds_every_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);clean=root/'clean';(clean/'src').mkdir(parents=True)
            (clean/'src/main.rs').write_text('clean compilation input')
            dirty=root/'main.rs';dirty.write_text('unrelated user edit')
            expected={'src/main.rs':sha(clean/'src/main.rs')}
            files=runner.compiled_source_files(clean,expected)
            self.assertEqual(files,{'src/main.rs':(clean/'src/main.rs').resolve()})
            self.assertEqual(dirty.read_text(),'unrelated user edit')
            (clean/'src/new.rs').write_text('unexpected input')
            with self.assertRaisesRegex(ValueError,'file set'):runner.compiled_source_files(clean,expected)
            (clean/'src/new.rs').unlink();(clean/'src/main.rs').write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed after'):runner.compiled_source_files(clean,expected)
            for bad in [None,{}, {'../src/main.rs':'fake'}, {'/src/main.rs':'fake'}]:
                with self.assertRaises(ValueError):runner.compiled_source_files(clean,bad)

    def fixture(self,root):
        prior=root/'prior';prior.mkdir()
        (prior/'FINISHED.json').write_text(json.dumps(dict(status='stopped')))
        (prior/'STATE.json').write_text(json.dumps(dict(status='stopped')))
        runtime=root/'uv';runtime.write_text('#!/bin/sh\nexit 0\n');runtime.chmod(0o700)
        (prior/'PLAN.json').write_text(json.dumps(dict(profile_uv=str(runtime))))
        for name in ['binary','auditor','known_states.tsv','mutation.tsv','profile.json','primitive_exclusions.json','double_bases.tsv']:
            (root/name).write_text('fixture')
        checks=dict(passed=True,design_passed=True,mechanism_passed=True,native_components_passed=True,
            binary_sha256=sha(root/'binary'),audit_binary_sha256=sha(root/'auditor'),controls=66,control_states=25192,
            known_bank_sha256=sha(root/'known_states.tsv'),mutation_seeds_sha256=sha(root/'mutation.tsv'),
            profile_queue_sha256=sha(root/'profile.json'),primitive_exclusions_sha256=sha(root/'primitive_exclusions.json'),
            double_bases_sha256=sha(root/'double_bases.tsv'))
        (root/'RESULTS.json').write_text(json.dumps(checks))
        return prior,runtime,['runner','--run','--launchd','--inherit-stopped','--predecessor',str(prior),
            '--binary',str(root/'binary'),'--audit-binary',str(root/'auditor'),'--checks',str(root/'RESULTS.json'),
            '--design-search','--native-components','--mechanism-search','--mutation-seeds',str(root/'mutation.tsv'),
            '--profile-input',str(root/'profile.json')]

    def test_minimal_path_reuses_saved_runtime_before_shared_admission(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve();prior,runtime,argv=self.fixture(root)
            control=root/'.run-control';control.mkdir()
            with (control/'.campaign.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                with patch.object(runner,'PROJECT',root),patch.object(sys,'argv',argv),patch.object(runner.shutil,'which',return_value=None):
                    # The valid runtime passes; shared admission then blocks us.
                    with self.assertRaises(BlockingIOError):runner.main()

    def test_missing_saved_runtime_fails_before_admission(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve();prior,runtime,argv=self.fixture(root);runtime.unlink()
            with patch.object(runner,'PROJECT',root),patch.object(sys,'argv',argv),patch.object(runner.shutil,'which',return_value=None):
                with self.assertRaisesRegex(ValueError,'isolated SymPy runtime required'):runner.main()
            self.assertFalse((root/'.run-control').exists())


if __name__=='__main__':unittest.main()
