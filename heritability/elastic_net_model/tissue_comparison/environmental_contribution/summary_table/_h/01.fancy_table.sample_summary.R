## Summarize results and plot correlation with q-value

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(gtsummary)
})

options(gt.html_tag_check = FALSE)

save_table <- function(pp, fn){
    for(ext in c(".pdf", ".tex", ".rtf")){
        gt::gtsave(as_gt(pp), filename=paste0(fn,ext))
    }
}

load_data <- function(){
    fn <- here("inputs/counts/_m",
               "RNAseq_Collection_postQC_n2844_6publicDataset_2021-01_geneRSE.Rdata")
    load(fn)
    sample_data <- SummarizedExperiment::colData(rse_gene) |> as.data.frame() |>
        filter(Dx %in% c("MDD", "Control"), Region == "DLPFC", Age > 17,
               Race %in% c("AA", "CAUC"), Dataset == "BrainSeq_Phase1") |>
        select(!c("bamFile", "rna_preSwap_BrNum")) |>
        mutate(Race=recode(Race, "AA"="BA", "CAUC"="WA")) |>
        mutate_if(is.character, as.factor)
    return(sample_data)
}

#### MAIN
                                        # Generate phenotype data
pheno_df <- load_data()
data.table::fwrite(pheno_df, file="phenotype_data.tsv", sep="\t",
                   row.names=FALSE, col.names=TRUE)

                                         # Generate pretty tables
fn <- "sample_breakdown.table"
pp <- pheno_df |> select(Dx, Race, Sex, Age, RIN) |>
    mutate(Sex=factor(Sex, labels=c("Female", "Male")),
           Race=factor(forcats::fct_drop(Race), labels=c("BA", "WA")),
           Dx=factor(forcats::fct_drop(Dx), labels=c("CTL", "MDD"))) |>
    tbl_summary(by="Dx", missing="no",
                statistic=all_continuous() ~ c("{mean} ({sd})")) |>
    modify_header(all_stat_cols()~"**{level}**<br>N = {n}") |>
    modify_spanning_header(all_stat_cols()~"**Diagnosis**") |>
    bold_labels() |> italicize_levels()
save_table(pp, fn)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
