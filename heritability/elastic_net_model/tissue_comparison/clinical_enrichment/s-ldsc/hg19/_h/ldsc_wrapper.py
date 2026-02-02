#!/usr/bin/env python
"""
Wrapper for LDSC scripts to fix pandas/Python 3 compatibility issues.

LDSC was written for Python 2 and older pandas versions. This wrapper patches:
1. pd.set_option() - 'precision' -> 'display.precision', etc.
2. DataFrame.ix - removed in pandas 1.0, restored as alias to .loc/.iloc
3. np.matrix deprecation warnings - suppressed
4. gzip.open() - use text mode by default for Python 2 compatibility

Usage:
    python ldsc_wrapper.py <ldsc_dir> <script_name> [args...]

Examples:
    python ldsc_wrapper.py /path/to/ldsc munge_sumstats.py --sumstats ...
    python ldsc_wrapper.py /path/to/ldsc ldsc.py --h2 ...
"""
import sys
import warnings
import gzip

# Suppress numpy matrix deprecation warnings
warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
warnings.filterwarnings('ignore', message='.*matrix subclass.*')

# 0. Patch gzip.open to use text mode by default (Python 2 compatibility)
# In Python 2, gzip.open returned strings; in Python 3, it returns bytes by default
_original_gzip_open = gzip.open

def _patched_gzip_open(filename, mode='rb', *args, **kwargs):
    """Patch gzip.open to use text mode for compatibility with Python 2-era code."""
    # If mode is default 'rb' or just 'r', convert to text mode
    if mode in ('rb', 'r'):
        mode = 'rt'
    elif mode == 'wb':
        mode = 'wt'
    return _original_gzip_open(filename, mode, *args, **kwargs)

gzip.open = _patched_gzip_open

# Patch pandas before ldsc imports it
import pandas as pd
import numpy as np

# 1. Patch pd.set_option for renamed options
_original_set_option = pd.set_option

def _patched_set_option(key, *args, **kwargs):
    """Translate old pandas option names to new ones."""
    option_map = {
        'precision': 'display.precision',
        'max_rows': 'display.max_rows',
        'max_columns': 'display.max_columns',
        'max_colwidth': 'display.max_colwidth',
    }
    if key in option_map:
        key = option_map[key]
    return _original_set_option(key, *args, **kwargs)

pd.set_option = _patched_set_option


# 2. Restore deprecated .ix accessor (removed in pandas 1.0)
# .ix was a mixed integer/label indexer; we approximate with .loc
class _IXIndexer:
    """Compatibility shim for deprecated pandas .ix accessor."""
    def __init__(self, df):
        self._df = df

    def __getitem__(self, key):
        # Handle tuple keys (row, col)
        if isinstance(key, tuple):
            row_key, col_key = key
            # Use iloc for integer slices, loc otherwise
            if isinstance(row_key, slice) and isinstance(row_key.start, (int, type(None))):
                if isinstance(col_key, slice) and isinstance(col_key.start, (int, type(None))):
                    return self._df.iloc[row_key, col_key]
                elif isinstance(col_key, list) and all(isinstance(c, int) for c in col_key):
                    return self._df.iloc[row_key, col_key]
            return self._df.loc[key]
        return self._df.loc[key]

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self._df.loc[key] = value
        else:
            self._df.loc[key] = value


# Add .ix property to DataFrame if it doesn't exist
if not hasattr(pd.DataFrame, 'ix'):
    pd.DataFrame.ix = property(lambda self: _IXIndexer(self))

if not hasattr(pd.Series, 'ix'):
    pd.Series.ix = property(lambda self: self.loc)


# Now run the requested LDSC script
if __name__ == '__main__':
    import os

    if len(sys.argv) < 3:
        print("Usage: python ldsc_wrapper.py <ldsc_dir> <script.py> [args...]")
        print("       python ldsc_wrapper.py <ldsc_dir> /path/to/script.py [args...]")
        sys.exit(1)

    ldsc_dir = sys.argv[1]
    script_name = sys.argv[2]

    # If script_name is an absolute path, use it directly
    # Otherwise, combine with ldsc_dir
    if os.path.isabs(script_name):
        script_path = script_name
    else:
        script_path = f'{ldsc_dir}/{script_name}'

    # Update argv so the script sees correct arguments
    sys.argv = [script_path] + sys.argv[3:]

    # Add ldsc directory to path (for imports like ldscore, ldsc modules)
    sys.path.insert(0, ldsc_dir)

    # Execute the requested script
    exec(open(script_path).read())
