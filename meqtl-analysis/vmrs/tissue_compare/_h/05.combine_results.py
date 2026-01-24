"""
Combine text files into one (tab-separated format).
This uses a streaming approach to reduce memory usage by
loading one file at a time.
"""
import gzip
import argparse
import session_info
from pathlib import Path

def combine_files(input_dir, output_dir, label):
    data_dir = Path(input_dir)
    input_files = sorted(data_dir.glob(f"{label}*.txt.gz"))
    output_file = Path(output_dir) / f"{label}_ancestry_combined.txt.gz"

    with gzip.open(output_file, 'wt') as outfile:
        for i, infile_path in enumerate(input_files):
            with gzip.open(infile_path, 'rt') as infile:
                for j, line in enumerate(infile):
                    if i > 0 and j == 0:
                        continue # Skip header in subsequent files
                    outfile.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=str, default="temp")
    parser.add_argument('--output-dir', type=str, default="local")
    parser.add_argument('--label', type=str, default="lfsr")
    args = parser.parse_args()

    combine_files(args.input_dir, args.output_dir, args.label)


if __name__ == "__main__":
    main()
