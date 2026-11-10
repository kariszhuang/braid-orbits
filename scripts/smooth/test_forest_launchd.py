"""Launch ownership, admission and failure cleanup controls (no real service)."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import os
from pathlib import Path
import plistlib
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from run_forest_fast import launchd_pid, spawn_launchd


class LaunchdTests(unittest.TestCase):
    @unittest.skipUnless(sys.platform=='darwin','macOS admission CLI')
    def test_managed_queue_cannot_bypass_parent_admission_lock(self):
        script=Path(__file__).resolve().parents[2]/'scripts/smooth/run_forest_fast.py'
        result=subprocess.run([sys.executable,str(script),'--queue','--launchd'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0);self.assertIn('requires --run admission',result.stderr)

    def test_service_must_have_a_live_pid(self):
        with patch('run_forest_fast.subprocess.run') as run:
            run.return_value=subprocess.CompletedProcess([],0,'{\n "PID" = 123;\n}\n','')
            self.assertEqual(launchd_pid('test.research'),123)
            run.return_value=subprocess.CompletedProcess([],0,'{\n "LastExitStatus" = 0;\n}\n','')
            with self.assertRaises(RuntimeError):launchd_pid('test.research')

    def test_submission_binds_service_to_actual_worker_state(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);calls=[]
            def run(command,**kwargs):
                calls.append(command)
                if command[1]=='bootstrap':
                    (root/'STATE.json').write_text(json.dumps(dict(pid=123,status='queued')))
                    return subprocess.CompletedProcess(command,0,'','')
                return subprocess.CompletedProcess(command,0,'{ "PID" = 123; }','')
            with patch('run_forest_fast.subprocess.run',side_effect=run):
                self.assertEqual(spawn_launchd(root,['python','worker'],'test.research'),123)
            self.assertEqual(calls[0][:3],['/bin/launchctl','bootstrap',f'gui/{os.getuid()}'])
            config=plistlib.loads(Path(calls[0][-1]).read_bytes())
            self.assertEqual((root/'LAUNCHD_CONFIG.plist').read_bytes(),Path(calls[0][-1]).read_bytes())
            self.assertEqual(config['ProgramArguments'],['python','worker'])
            self.assertTrue(config['RunAtLoad']);self.assertFalse(config['KeepAlive'])
            self.assertNotIn('StartInterval',config);self.assertNotIn('StartCalendarInterval',config)

    def test_missing_or_mismatched_registration_stops_only_new_service(self):
        for mismatched in [False,True]:
            with tempfile.TemporaryDirectory() as temp:
                root=Path(temp);calls=[]
                if mismatched:(root/'STATE.json').write_text(json.dumps(dict(pid=999,status='queued')))
                def run(command,**kwargs):
                    calls.append(command)
                    return subprocess.CompletedProcess(command,0,'{ "PID" = 123; }','')
                with patch('run_forest_fast.subprocess.run',side_effect=run):
                    with self.assertRaises(RuntimeError):spawn_launchd(root,['python','worker'],'test.research',timeout=.01)
                self.assertTrue((root/'STOP').exists())
                self.assertEqual(calls[-1],['/bin/launchctl','bootout',f'gui/{os.getuid()}/test.research'])


if __name__=='__main__':unittest.main()
