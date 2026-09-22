## Resource-aware helpers for parallelizing genome scans across loci.
## Kept separate from any one scan function since the locus loop shape
## (independent per-locus fits given a fixed null model) is shared across
## scan.h2lmm(), imputed.snp.scan.h2lmm(), and others.

#' Determine a safe number of cores for parallelizing a genome scan
#'
#' Sizes worker count off the CPU allocation actually available to this
#' process/job -- not just the host machine's total core count -- and,
#' when a per-locus memory estimate is supplied, off available memory too.
#' On a shared SLURM node, \code{parallel::detectCores()} alone reports the
#' node's total cores, not what SLURM allocated to this job; blindly using
#' it can oversubscribe a shared node. SLURM's own environment variables
#' are checked first and take priority whenever present.
#'
#' @param requested.cores DEFAULT: NULL. An explicit user request. Never
#' exceeds the detected CPU/memory budget -- this is a ceiling, not a
#' guaranteed allocation.
#' @param mem.per.locus.mb DEFAULT: NULL. An estimate (e.g. empirically
#' profiled from fitting one locus) of the memory a single locus fit needs,
#' in megabytes. When supplied, available memory is used to additionally
#' cap the number of workers. When NULL, no memory-based cap is applied.
#' @param mem.safety.frac DEFAULT: 0.8. Fraction of detected available
#' memory that may be used for the memory-based cap; leaves headroom for
#' the objects mclapply's forked workers share via copy-on-write (which is
#' not airtight in R -- GC/refcounting can force partial copies over time)
#' plus whatever else is running on the node.
#' @param n.loci DEFAULT: NULL. Total number of loci to be scanned. Worker
#' count is never set higher than is useful for the amount of work present.
#' @param min.loci.per.core DEFAULT: 1. Minimum loci per worker, used with
#' n.loci to avoid over-splitting a small scan.
#' @param verbose DEFAULT: TRUE. Print the resolved core count and how it
#' was derived.
#' @export
#' @examples determine.scan.cores()
determine.scan.cores <- function(requested.cores=NULL,
                                 mem.per.locus.mb=NULL,
                                 mem.safety.frac=0.8,
                                 n.loci=NULL,
                                 min.loci.per.core=1,
                                 verbose=TRUE){
  cpu.budget <- get.allocated.cpus()

  mem.budget.cores <- Inf
  if(!is.null(mem.per.locus.mb)){
    available.mem.mb <- get.available.memory.mb()
    if(!is.na(available.mem.mb)){
      mem.budget.cores <- max(1, floor((available.mem.mb * mem.safety.frac) / mem.per.locus.mb))
    }
    else if(verbose){
      cat("determine.scan.cores(): could not determine available memory on this system; skipping memory-based core limit.\n")
    }
  }

  num.cores <- min(cpu.budget, mem.budget.cores)

  if(!is.null(requested.cores)){
    num.cores <- min(num.cores, requested.cores)
  }
  if(!is.null(n.loci)){
    num.cores <- min(num.cores, max(1, floor(n.loci / min.loci.per.core)))
  }
  num.cores <- as.integer(max(1, floor(num.cores)))

  if(.Platform$OS.type == "windows" & num.cores > 1){
    if(verbose){
      cat("determine.scan.cores(): parallel::mclapply() (fork-based) is not supported on Windows; falling back to num.cores=1.\n")
    }
    num.cores <- 1
  }

  if(verbose){
    mem.msg <- if(is.finite(mem.budget.cores)) paste0(", memory budget: ", mem.budget.cores) else ""
    cat(sep="", "determine.scan.cores(): using ", num.cores, " core(s) (CPU budget: ", cpu.budget, mem.msg, ")\n")
  }
  return(num.cores)
}

## Number of CPUs actually allocated to this process. Checks SLURM's own
## reporting first (in order of specificity/reliability), since
## parallel::detectCores() reports the host's total logical cores, which on
## a shared SLURM node is not the same thing as this job's allocation.
get.allocated.cpus <- function(){
  slurm.cpus.per.task <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if(nzchar(slurm.cpus.per.task)){
    parsed <- suppressWarnings(as.integer(slurm.cpus.per.task))
    if(!is.na(parsed) && parsed >= 1){ return(parsed) }
  }
  slurm.job.cpus <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")
  if(nzchar(slurm.job.cpus)){
    parsed <- parse.slurm.cpu.string(slurm.job.cpus)
    if(!is.na(parsed) && parsed >= 1){ return(parsed) }
  }
  slurm.cpus.on.node <- Sys.getenv("SLURM_CPUS_ON_NODE")
  if(nzchar(slurm.cpus.on.node)){
    parsed <- suppressWarnings(as.integer(slurm.cpus.on.node))
    if(!is.na(parsed) && parsed >= 1){ return(parsed) }
  }
  detected <- tryCatch(parallel::detectCores(logical=FALSE), error=function(e) NA)
  if(is.na(detected) || detected < 1){ detected <- 1L }
  return(as.integer(detected))
}

## SLURM_JOB_CPUS_PER_NODE can be a compressed, comma-separated list like
## "2(x4),3" (2 cpus on each of 4 nodes, then 3 on one more) for multi-node
## jobs. For scan.h2lmm()'s single-process use this job's own allocation is
## the first entry; take its cpu count.
parse.slurm.cpu.string <- function(x){
  first.entry <- strsplit(x, ",", fixed=TRUE)[[1]][1]
  m <- regmatches(first.entry, regexec("^([0-9]+)(\\(x[0-9]+\\))?$", first.entry))[[1]]
  if(length(m) < 2){ return(NA_integer_) }
  as.integer(m[2])
}

## Available memory in megabytes, checked in order of reliability:
## 1) SLURM's own reporting of this job's memory allocation
## 2) cgroup memory limits (catches jobs where --mem was requested but the
##    SLURM_MEM_* env vars weren't exported)
## 3) OS-level available memory (Linux /proc/meminfo, macOS sysctl)
## Returns NA if none of the above resolve -- callers should treat that as
## "unknown" and skip memory-based sizing rather than guessing.
get.available.memory.mb <- function(){
  mem.per.node <- Sys.getenv("SLURM_MEM_PER_NODE")
  if(nzchar(mem.per.node)){
    parsed <- suppressWarnings(as.numeric(mem.per.node))
    if(!is.na(parsed)){ return(parsed) }
  }
  mem.per.cpu <- Sys.getenv("SLURM_MEM_PER_CPU")
  if(nzchar(mem.per.cpu)){
    parsed <- suppressWarnings(as.numeric(mem.per.cpu))
    if(!is.na(parsed)){ return(parsed * get.allocated.cpus()) }
  }

  cgroup.v2 <- "/sys/fs/cgroup/memory.max"
  cgroup.v1 <- "/sys/fs/cgroup/memory/memory.limit_in_bytes"
  if(file.exists(cgroup.v2)){
    val <- suppressWarnings(readLines(cgroup.v2, n=1))
    if(length(val) == 1 && !is.na(val) && val != "max"){
      parsed <- suppressWarnings(as.numeric(val))
      if(!is.na(parsed)){ return(parsed / 1024^2) }
    }
  }
  if(file.exists(cgroup.v1)){
    parsed <- suppressWarnings(as.numeric(readLines(cgroup.v1, n=1)))
    ## cgroup v1 reports a very large sentinel value (~2^63ish) for "no limit"
    if(!is.na(parsed) && parsed < 1e15){
      return(parsed / 1024^2)
    }
  }

  if(file.exists("/proc/meminfo")){
    meminfo <- readLines("/proc/meminfo")
    avail.line <- grep("^MemAvailable:", meminfo, value=TRUE)
    if(length(avail.line) == 1){
      kb <- suppressWarnings(as.numeric(regmatches(avail.line, regexpr("[0-9]+", avail.line))))
      if(!is.na(kb)){ return(kb / 1024) }
    }
  }
  if(identical(Sys.info()[["sysname"]], "Darwin")){
    mem.bytes <- suppressWarnings(tryCatch(as.numeric(system("sysctl -n hw.memsize", intern=TRUE)), error=function(e) NA))
    if(!is.na(mem.bytes)){ return(mem.bytes / 1024^2) }
  }
  return(NA_real_)
}
