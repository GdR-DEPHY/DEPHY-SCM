"""Root pytest configuration.

Ensures that the project's root directory (containing the
``dephycf`` package) is present in ``sys.path``, regardless of the
directory ``pytest`` is invoked from.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
