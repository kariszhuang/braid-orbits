"""Attach the checked static UI and pinned proof to a local release (no deploy)."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import hashlib
from pathlib import Path
from classify_smooth_isotrivial import PROOF_SHA256


def package(release):
    project=Path(__file__).resolve().parents[2]
    proof=project/'deploy/math/sources/isotrivial-proof.pdf'
    if hashlib.sha256(proof.read_bytes()).hexdigest()!=PROOF_SHA256:
        raise ValueError('proof bytes differ from the classification reference')
    sources=[(project/'deploy/math/ui'/name,'ui/'+name) for name in ['index.html','app.mjs','catalog.mjs','style.css']]
    sources.append((proof,'sources/isotrivial-proof.pdf'))
    for source,relative in sources:
        target=release/relative;data=source.read_bytes()
        if target.exists():
            if target.read_bytes()!=data:raise ValueError('refusing to overwrite packaged assets')
            continue
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as stream:stream.write(data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--release',type=Path,required=True)
    package(p.parse_args().release)
