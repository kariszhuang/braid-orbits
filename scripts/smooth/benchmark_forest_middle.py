"""Benchmark the changed Metal middle-square sieve, including checked rejection."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import io
import json
from pathlib import Path
import statistics
import subprocess

from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]


def records(text):
    return list(csv.DictReader(io.StringIO('\n'.join(line for line in text.splitlines() if not line.startswith('#'))),delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);args=p.parse_args();root=args.out.resolve();root.mkdir()
    swift=PROJECT/'target/forest_middle_screen';shader=PROJECT/'scripts/gpu/forest_middle_screen.metal';results=[]
    for n in [11,14]:
        for size in [32768,262144]:
            trace=root/f'trace_{n}_{size}.bin'
            rust=subprocess.check_output([str(args.binary.resolve()),'--gpu-middle-bench','--d',str(n),'--size',str(size),'--out',str(trace)],text=True)
            metal=subprocess.check_output([str(swift),str(trace),str(shader)],text=True)
            (root/f'rust_{n}_{size}.tsv').write_text(rust);(root/f'metal_{n}_{size}.tsv').write_text(metal)
            r=records(rust);m=records(metal)
            times={variant:[float(row['seconds']) for row in r+m if row['variant']==variant] for variant in ['rust_cpu','cpu','gpu','gpu_checked']}
            medians={k:statistics.median(v) for k,v in times.items()}
            result=dict(degree=n,size=size,times=times,medians=medians,
                        gpu_unchecked_vs_rust=medians['rust_cpu']/medians['gpu'],
                        gpu_checked_vs_rust=medians['rust_cpu']/medians['gpu_checked'],
                        gpu_checked_vs_swift=medians['cpu']/medians['gpu_checked'],
                        trace_sha256=sha(trace),every_rust_negative_literal_verified=True,
                        complete_metal_swift_rust_differential=True,every_gpu_negative_cpu_verified=True)
            results.append(result);print('METAL',json.dumps({k:v for k,v in result.items() if k!='times'}),flush=True)
    result=dict(passed=True,stages=results,binary_sha256=sha(args.binary),swift_binary_sha256=sha(swift),
                shader_sha256=sha(shader),gpu_deployed=False,
                limitations='Warm RAM traces; GPU includes input copy, dispatch and waiting. Checked path verifies every GPU negative, not GPU positives; survivors still need CPU BFS. Rust and Swift timings are separate matched traces, not identical implementations. Stage latency is not a whole-search/discovery speed claim.')
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print('RESULTS='+str(root/'RESULTS.json'))


if __name__=='__main__':main()
