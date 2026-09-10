"""
_path_helper.py
Utility imported at the top of every analysis script to locate
config/paths.py regardless of the script's working directory.

Usage in any script under scripts/*:
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..', 'config'))
    from paths import *

Or simply:
    from _path_helper import *
"""
import sys
import os

# Walk up from this file's location until we find the 'config' directory
_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = None
for _candidate in [_here,
                   os.path.join(_here, '..'),
                   os.path.join(_here, '..', '..')]:
    if os.path.isdir(os.path.join(_candidate, 'config')):
        _repo_root = os.path.abspath(_candidate)
        break

if _repo_root is None:
    raise RuntimeError(
        "Cannot find repository root (expected a 'config/' directory). "
        "Ensure you are running from within the nature_aging_code_release folder."
    )

_config_dir = os.path.join(_repo_root, 'config')
if _config_dir not in sys.path:
    sys.path.insert(0, _config_dir)

from paths import *  # noqa: F401, F403  expose all path constants
