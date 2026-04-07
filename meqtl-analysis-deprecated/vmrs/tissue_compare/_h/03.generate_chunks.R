chunk_data <- function(chunk_size, input_dir, output_dir){
    load(file.path(input_dir, "bhat_shat.RData"))
    chunks = split(1:dim(bhat)[1],
                   cut(1:dim(bhat)[1], chunk_size, labels=FALSE))
    for(chunk_num in 1:chunk_size){
        bhat_chunk <- bhat[chunks[[chunk_num]],]
        shat_chunk <- shat[chunks[[chunk_num]],]
        save(bhat_chunk, shat_chunk,
             file=file.path(output_dir, paste0("chunk_", chunk_num, "_bhat_shat.RData"))
             )
    }
}

## Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("
Usage: Rscript 03.generate_chunks.R <input_dir> <output_dir> [chunk_size]

Required:
  input_dir       Directory containing bhat/shat files
  output_dir      Where the chunks will go

Optional:
  chunk_size      Default = 250
")
}

input_dir   <- normalizePath(args[1], mustWork = TRUE)
output_dir  <- normalizePath(args[2], mustWork = TRUE)
chunk_size  <- ifelse(length(args) >= 3, as.integer(args[3]), 250)

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

## Run mashr for specific feature
chunk_data(chunk_size, input_dir, output_dir)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
