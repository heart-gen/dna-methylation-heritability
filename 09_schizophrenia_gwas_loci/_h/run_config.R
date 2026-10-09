#### Read configuration from a run's immutable snapshot ####
##
## AGENTS.md 5.2 makes a run immutable and 9 makes it reproducible, and the
## submit drivers already snapshot _h/ into runs/{RUN_ID}/code/_h so that
## editing _h/ while a run is queued cannot change what that run executes.
## Configuration was NOT covered by that: every stage re-read config/ from the
## live working tree, so a config edit -- or, as happened on 2026-09-02, a git
## branch switch that removed the file entirely -- changed or broke a run that
## was already in flight, while its manifest attested to the original checksum.
##
## The submit driver now copies config/ into runs/{RUN_ID}/code/config, and this
## loader prefers that snapshot. It falls back to the live config only when a
## run has no snapshot (a hand-run stage), so existing runs keep working.
load_run_config <- function(name, run_dir) {
    f <- file.path(run_dir, "code", "config", paste0(name, ".yml"))
    if (!file.exists(f)) {
        return(load_config(name))
    }
    cfg <- yaml::read_yaml(f)
    attr(cfg, "config_file") <- f
    attr(cfg, "config_sha256") <- file_sha256(f)
    cfg
}
