#!/usr/bin/env python
"""
Wrapper for LDSC scripts to fix remaining Python 3 compatibility issues.

The main LDSC code has been updated for Python 3 and pandas 2.x compatibility.
This wrapper handles remaining compatibility issues:
1. np.matrix deprecation warnings - suppressed
2. gzip.open() - use text mode by default for Python 2-era code

Usage:
    python ldsc_wrapper.py <ldsc_dir> <script_name> [args...]

Examples:
    python ldsc_wrapper.py /path/to/ldsc munge_sumstats.py --sumstats ...
    python ldsc_wrapper.py /path/to/ldsc ldsc.py --h2 ...
"""
import sys
import warnings
import gzip

# Suppress numpy matrix deprecation warnings (still used for string formatting in some places)
warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
warnings.filterwarnings('ignore', message='.*matrix subclass.*')

# Patch gzip.open to use text mode by default (Python 2 compatibility)
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
