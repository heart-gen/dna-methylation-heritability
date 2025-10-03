## Script to prepare inputs for meQTL analysis of VMRs
suppressMessages({
    library(here)
    library(dplyr)
    library(SummarizedExperiment)
    library(data.table)
})

#### Main

                                        # Create output dir
output_path <- here("meqtl-analysis", "vmrs", "hippocampus", "_m")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}


#### Phenotype data

                                        # Read in pheno data
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-AA.tsv")
pheno <- fread(pheno_file_path, header = TRUE)

sample_file_path <- here("heritability/hippocampus/_m/samples.txt")
samples <- fread(sample_file_path, header = FALSE)
colnames(samples) <- c("brnum", "fid")

                                        # Basic filtering
pheno <- pheno |>
        select(brnum, agedeath, sex, race, primarydx, region) |>
        filter(race == "AA", agedeath >= 17, region == "hippocampus") |>
        inner_join(samples, by = "brnum")

                                        # Save filtered phenotype file
out_pheno   <- file.path(output_path, "phenotypes.csv")
as.data.frame(pheno) |>
    write.csv(file = out_pheno)

#### Combine residualized methylation matrix

                                        # Read in .phen files
meth_file_path <- here("heritability/hippocampus/_m/vmr")
meth_files <- list.files(path = meth_file_path, pattern = "_meth\\.phen$", 
                         recursive = TRUE, full.names = TRUE)

feature_ids <- paste0("VMR", seq_along(meth_files))
meth_list <- vector("list", length(meth_files))
for (i in seq_along(meth_files)) {
  df <- fread(meth_files[i], select = c("V1", "V3"))
  colnames(df) <- c("brnum", "meth")
  meth_list[[i]] <- df
}

                                        # Add brnum and feature IDs
brnums <- unique(unlist(lapply(meth_list, function(df) df$brnum)))
meth_matrix <- matrix(NA, nrow = length(meth_files), ncol = length(brnums),
                   dimnames = list(feature_ids, brnums))

for (i in seq_along(meth_list)) {
  df <- meth_list[[i]]
  meth_matrix[i, df$brnum] <- df$meth
}

                                        # Save meth matrix
out_meth   <- file.path(output_path, "normalized_methylation.tsv")
as.data.frame(meth_matrix) |> tibble::rownames_to_column("feature_id") |>
    write.table(file=out_meth, sep="\t", 
                quote=FALSE, row.names=FALSE)

#data.table::fwrite(meth_df, "norm_vmr.tsv", sep='\t', row.names=TRUE)

#### Export annotation
bed_file_path <- here("heritability/hippocampus/_m/vmr.bed")
vmr_bed <- fread(bed_file_path, header = FALSE)
colnames(vmr_bed) <- c("seqnames", "start", "end")

                                        # Format bed file
annot <- tibble::rowid_to_column(vmr_bed, "feature_id") |> 
  mutate(seqnames = paste0("chr", seqnames),
         feature_id = paste0("VMR", row_number())) |>
  select(seqnames, start, end, feature_id) |> 
  mutate(index=feature_id) |> 
  tibble::column_to_rownames("index")

                                        # Save feature bed file
out_vmr <- file.path(output_path, "feature.bed")
data.table::fwrite(annot, out_vmr, sep='\t', row.names=TRUE)

#### Export chr list
out_chr <- file.path(output_path, "vcf_chr_list.txt")
chr_list <- paste0("chr", c(1:22))
write.table(chr_list, file = out_chr, quote = FALSE, 
            row.names = FALSE, col.names = FALSE)


#### R session information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
