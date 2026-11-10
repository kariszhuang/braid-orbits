"""Put every scripts/ topic folder on the import path, so scripts can import each other.

Each topic folder holds an identical copy of this file. A script imports it first
(`import _scripts_path`); Python finds the copy in the script's own folder.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOLDERS = sorted(d for d in HERE.parent.iterdir() if d.is_dir() and (d / '_scripts_path.py').exists())
for _d in FOLDERS:
    if str(_d) not in sys.path:
        sys.path.append(str(_d))


def find(name):
    """Path of a script (or other file) by name, whichever topic folder holds it."""
    for d in [HERE, *FOLDERS]:
        if (d / name).exists():
            return d / name
    return HERE / name
