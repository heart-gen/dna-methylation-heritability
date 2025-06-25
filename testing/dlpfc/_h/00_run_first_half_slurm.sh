#!/bin/bash
#SBATCH --output=run_first_half_slurm.log
#SBATCH --mem=25G

#  After running 'install_software.sh', this should point to the directory
#  where this repo was cloned, and not say "$PWD"
ORIG_DIR=$PWD

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

$ORIG_DIR/Software/bin/nextflow $ORIG_DIR/first_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "paired" \
    --reference "hg38" \
    --all_alignments \
    --input /projects/p32505/users/elisa/projects/dna-methylation-heritability/testing/dlpfc/_h/BiocMAP/benchmark/samples.manifest \
    --output "./results" \
    --trim_mode "adaptive" \
    -profile first_half_slurm
