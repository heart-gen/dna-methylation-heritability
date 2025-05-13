import glob
import re
from collections import defaultdict
import os

# -----------------------
# Define pipeline stages
# -----------------------

pipeline_steps = {
    "extract_snp": {
        "region": re.compile(r'on ([0-9XY]+): ([0-9]+)-([0-9]+)'),
        "error": [
            re.compile(r'No variants remaining after main filters', re.IGNORECASE)
        ]
    },
    "greml": {
        "region": re.compile(r'Performing GREML analysis on ([0-9XY]+): ([0-9]+)-([0-9]+)', re.IGNORECASE),
        "error": [
            re.compile(r'analysis stopped because more than half of the variance components are constrained', re.IGNORECASE),
            re.compile(r'\.fam not found', re.IGNORECASE),
            re.compile(r'Log-likelihood not converged \(stop after 100 iteractions\)', re.IGNORECASE)
        ]
    }
}

# -----------------------
# Initialize results
# -----------------------

errors_by_stage_chr = defaultdict(set)

# Get all log files
log_files = glob.glob("./logs/*/*.log")
print(f"Found {len(log_files)} log files")

# -----------------------
# Process each log
# -----------------------

for filepath in log_files:
    filename = os.path.basename(filepath)

    # Get stage based on filename
    stage = next((step for step in pipeline_steps if step in filename), None)

    step_conf = pipeline_steps[stage]
    region_regex = step_conf["region"]
    error_regexes = step_conf["error"]

    regions = set()
    error_msgs = []

    with open(filepath) as f:
        for line in f:
            # Match region
            region_match = region_regex.search(line)
            if region_match:
                chr_val = f"chr{region_match.group(1)}"
                print(f"Matched region: {region_match.groups()} in file {filename}")
                if len(region_match.groups()) == 3:
                    start = region_match.group(2)
                    end = region_match.group(3)
                else:
                    start = end = ""
                regions.add((chr_val, start, end))

            # Match errors
            for err_re in error_regexes:
                if err_re.search(line):
                    error_msgs.append(line.strip())
                    print(f"Matched error: {line.strip()} in file {filename}")
                    break  # stop after first matched error per line

    # If there are error messages, save them
    if error_msgs:
        if regions:
            for chr_val, start, end in regions:
                for err in error_msgs:
                    errors_by_stage_chr[(stage, chr_val)].add((start, end, err))
        else:
            # Try to infer chromosome from filename
            chr_match = re.search(r'chr[_\-]?([0-9XY]+)', filepath, re.IGNORECASE)
            chr_val = f"chr{chr_match.group(1)}" if chr_match else "chrNA"
            for err in error_msgs:
                errors_by_stage_chr[(stage, chr_val)].add(("", "", err))

# -----------------------
# Output results
# -----------------------

with open("summary/pipeline_errors.tsv", "w") as out:
    out.write("stage\tchr\tstart\tend\terror\n")
    for (stage, chr_val), entries in sorted(errors_by_stage_chr.items()):
        for start, end, error_msg in sorted(entries):
            out.write(f"{stage}\t{chr_val}\t{start}\t{end}\t{error_msg}\n")

print("Summary of errors per stage and chromosome:\n")
print("Stage\tChromosome\t# Regions with Errors")

# Now include all stages in the summary print
for (stage, chr_val), entries in sorted(errors_by_stage_chr.items()):
    print(f"{stage}\t{chr_val}\t{len(entries)}")