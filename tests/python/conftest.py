import sys
import pathlib

# ensure build python module is discoverable during ctest runs
ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD_ROOT = ROOT / 'build'
if BUILD_ROOT.exists():
    sys.path.insert(0, str(BUILD_ROOT))
    sys.path.insert(0, str(BUILD_ROOT / 'python'))  # out-of-tree module dir
