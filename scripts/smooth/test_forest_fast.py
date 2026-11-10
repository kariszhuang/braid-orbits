try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
import fcntl
from unittest.mock import patch
from types import SimpleNamespace

PROJECT=Path(__file__).resolve().parents[2]

class FastRunnerTests(unittest.TestCase):
    def test_storage_guard_has_ten_gib_boundary_and_data_folder_requirement(self):
        from run_forest_fast import storage_available, MIN_FREE_BYTES
        self.assertEqual(MIN_FREE_BYTES,10*2**30)
        with patch('run_forest_fast.Path.is_dir',return_value=True):
            for free,expected in [(MIN_FREE_BYTES-1,False),(MIN_FREE_BYTES,True),(20*2**30,True)]:
                with patch('run_forest_fast.shutil.disk_usage',return_value=SimpleNamespace(free=free)):
                    self.assertEqual(storage_available(Path('data')),expected)
        with patch('run_forest_fast.Path.is_dir',return_value=False):
            self.assertFalse(storage_available(Path('data')))
    def test_stopped_handoff_requires_explicit_opt_in_and_final_state(self):
        from run_forest_fast import predecessor_finished
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'FINISHED.json').write_text(json.dumps(dict(status='stopped',epochs=[7,8])))
            (root/'STATE.json').write_text(json.dumps(dict(status='searching')))
            with self.assertRaises(RuntimeError):predecessor_finished(root,True)
            (root/'STATE.json').write_text(json.dumps(dict(status='stopped')))
            with self.assertRaises(RuntimeError):predecessor_finished(root)
            self.assertEqual(predecessor_finished(root,True)['epochs'],[7,8])
    def test_dead_queue_pointer_replacement_is_atomic(self):
        from run_forest_fast import replace_queue_guard
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'guard.json';path.write_text(json.dumps(dict(pid=1,output='old')))
            replace_queue_guard(path,dict(pid=2,output='new'))
            self.assertEqual(json.loads(path.read_text()),dict(pid=2,output='new'))
            self.assertFalse(path.with_name(path.name+'.tmp').exists())
    def test_plan_is_read_only_and_bounded(self):
        result=subprocess.run([sys.executable,str(PROJECT/'scripts/smooth/run_forest_fast.py')],capture_output=True,text=True,check=True)
        plan=json.loads(result.stdout)
        self.assertEqual(plan['workers'],6)
        self.assertFalse(plan['gpu'])
        self.assertIn('shared lock',plan['activation'])

    def test_queued_worker_does_not_start_and_incomplete_audit_blocks(self):
        with tempfile.TemporaryDirectory(prefix='forest-fast-lock-test-')as temp:
            root=Path(temp);predecessor=root/'predecessor';predecessor.mkdir()
            (predecessor/'FINISHED.json').write_text(json.dumps(dict(status='incomplete')))
            (root/'PLAN.json').write_text(json.dumps(dict(predecessor=str(predecessor))))
            guard=(root/'lock').open('a');fcntl.flock(guard,fcntl.LOCK_EX)
            waiter=(root/'lock').open('a')
            child=subprocess.Popen([sys.executable,str(PROJECT/'scripts/smooth/run_forest_fast.py'),'--worker',str(root),'--lock-fd',str(waiter.fileno())],pass_fds=(waiter.fileno(),),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            try:
                limit=time.monotonic()+10
                while not(root/'STATE.json').exists()and time.monotonic()<limit:time.sleep(.02)
                self.assertEqual(json.loads((root/'STATE.json').read_text())['status'],'queued')
                self.assertIsNone(child.poll())
                self.assertFalse((root/'INHERITED_LEADS.json').exists())
                fcntl.flock(guard,fcntl.LOCK_UN);guard.close()
                self.assertNotEqual(child.wait(timeout=10),0)
                self.assertEqual(json.loads((root/'STATE.json').read_text())['status'],'failed')
                self.assertIn('did not complete',(root/'FAILURE.json').read_text())
            finally:
                if child.poll()is None:child.terminate();child.wait(timeout=10)
                waiter.close();guard.close()

if __name__=='__main__':unittest.main()
