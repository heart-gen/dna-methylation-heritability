#### Thread limits for a SLURM allocation ####
##
## data.table sizes its thread pool from the number of cores it can SEE, which
## on a Quest compute node is the whole machine -- not the cgroup SLURM handed
## us. The first full-scale DLPFC run therefore ran with:
##
##     parallel::detectCores()  256
##     data.table::getDTthreads() 128
##
## inside allocations that requested one CPU. That is 128 threads contending for
## one core's worth of scheduling, on every task of every array.
##
## Setting R_DATATABLE_NUM_THREADS in 00_shared/slurm.sh is not enough on its
## own: data.table reads it when the package is FIRST loaded, so any script that
## attaches data.table before the variable is set keeps the node-sized pool. The
## call here is unconditional and runs after the package is available, which is
## the only reliable point. Both are kept -- the environment variable also covers
## OpenMP and BLAS, which never consult this file.
##
## The ceiling is 2x the allocation rather than exactly the allocation, per PI
## direction: hyperthreads make modest oversubscription useful, and the defect
## being fixed is 128-on-1, not 2-on-1.

#' Clamp threaded libraries to this process's CPU allocation.
#'
#' @param cpus CPUs allocated. Defaults to SLURM_CPUS_PER_TASK, then V2_CPUS,
#'   then 1 -- so an interactive run on a login node stays single-threaded.
#' @param factor Ceiling as a multiple of `cpus`.
#' @return The applied thread count, invisibly.
limit_threads <- function(cpus = NULL, factor = 2L) {
    if (is.null(cpus)) {
        cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("V2_CPUS", "1"))
    }
    cpus <- suppressWarnings(as.integer(cpus))
    if (is.na(cpus) || cpus < 1L) cpus <- 1L
    n <- cpus * as.integer(factor)

    if (requireNamespace("data.table", quietly = TRUE)) {
        ## Guard against a pool that is already larger than the allocation:
        ## setDTthreads(throttle=) does not shrink an existing pool.
        data.table::setDTthreads(n, restore_after_fork = TRUE)
        got <- data.table::getDTthreads()
        if (got > n) {
            warning("data.table still reports ", got, " threads after ",
                    "setDTthreads(", n, ")")
        }
    }
    ## RhpcBLASctl is not a hard dependency; the env vars in 00_shared/slurm.sh
    ## cover BLAS when it is absent.
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(cpus)
        RhpcBLASctl::omp_set_num_threads(cpus)
    }
    Sys.setenv(OMP_NUM_THREADS = cpus, R_DATATABLE_NUM_THREADS = cpus)
    invisible(n)
}
