#!/bin/bash
# Slice PGC3 hg19 summary statistics to the hg38 schizophrenia locus windows.
#
#   slice_gwas_sumstats.sh <sumstats.gz> <chain> <windows.bed> <outdir> \
#                          <liftover-bin> <bedtools-bin>
#
# Called by _h/01_define_scz_loci.R, which runs on the submit host. Reading the
# 7.6M-variant table, lifting it, and joining in memory was OOM-killed there, so
# the coordinates are streamed instead: awk emits them, liftOver maps them,
# bedtools selects the windows, and only surviving rows are ever read into R.
#
# The join key is the LINE NUMBER of the de-commented stream, so no coordinate
# or allele string round-trips through two formats. Outputs, in $OUTDIR:
#   sumstats.hg19.bed        every variant's hg19 coordinate + line number
#   sumstats.inwindow.bed    lifted hg38 coordinate + line number + locus_id
#   sumstats.selected.tsv    the original rows for those line numbers, headered

set -euo pipefail

SUMSTATS=${1:?sumstats.gz}
CHAIN=${2:?chain file}
WIN=${3:?windows bed}
OUTDIR=${4:?output dir}
LIFTOVER=${5:?liftOver binary}
BEDTOOLS=${6:?bedtools binary}

mkdir -p "$OUTDIR"
HEADERED="$OUTDIR/sumstats.headered.tsv"
trap 'rm -f "$HEADERED"' EXIT

zcat "$SUMSTATS" | grep -v '^##' > "$HEADERED"

# Column positions are read from the header rather than assumed. CHROM ID POS
# is the documented PGC layout, but a silently reordered file would otherwise
# yield coordinates taken from the wrong field -- which would liftOver cleanly
# and be wrong everywhere downstream.
awk -F'\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        if (!(("CHROM" in col) && ("POS" in col))) {
            print "ERROR: sumstats header lacks CHROM/POS" > "/dev/stderr"
            exit 1
        }
        next
    }
    {
        p = $col["POS"] + 0
        if (p > 0) printf "chr%s\t%d\t%d\t%d\n", $col["CHROM"], p - 1, p, NR
    }
' "$HEADERED" > "$OUTDIR/sumstats.hg19.bed"

"$LIFTOVER" "$OUTDIR/sumstats.hg19.bed" "$CHAIN" \
    "$OUTDIR/sumstats.hg38.bed" "$OUTDIR/sumstats.unmapped.bed"

sort -k1,1 -k2,2n "$OUTDIR/sumstats.hg38.bed" > "$OUTDIR/sumstats.hg38.sorted.bed"
sort -k1,1 -k2,2n "$WIN" > "$OUTDIR/locus-windows.sorted.bed"

"$BEDTOOLS" intersect -a "$OUTDIR/sumstats.hg38.sorted.bed" \
    -b "$OUTDIR/locus-windows.sorted.bed" -wa -wb -sorted \
    > "$OUTDIR/sumstats.inwindow.bed"

# Pull the surviving rows out of the original stream by line number, keeping
# the header. A variant in two overlapping windows appears once here; the locus
# assignment is carried by sumstats.inwindow.bed, not by this file.
#
# row_key is emitted as a COLUMN, not left implicit in the row order. The
# selected file is a subset, so its Nth data row is not line N+1 of the stream,
# and reconstructing the key from position would join every variant to the
# wrong locus -- silently, since most keys still exist.
awk -F'\t' -v OFS='\t' '
    NR == FNR { keep[$4] = 1; next }
    FNR == 1 { print "row_key", $0; next }
    (FNR in keep) { print FNR, $0 }
' "$OUTDIR/sumstats.inwindow.bed" "$HEADERED" > "$OUTDIR/sumstats.selected.tsv"

echo "sliced $(wc -l < "$OUTDIR/sumstats.inwindow.bed") window hits into $OUTDIR"
