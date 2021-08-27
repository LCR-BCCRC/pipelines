#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

threads  <- args[1]
nrml.dir <- args[2]
tbl.out  <- args[3]
pon.path <- args[4]
par.path <- args[5]
build    <- args[6]
choose.samples <- args[7]

library(data.table)
library(dryclean)
library(parallel)

if (!choose.samples %in% c('rand', 'all', 'cluster')){
   stop('Argument 7 must be one of "rand", "all", or "cluster"')
}

nrml.dir.abs <- normalizePath(nrml.dir)
nrml.rds <- Sys.glob(file.path(nrml.dir.abs, "*.rds"))
nrml.smp <- gsub('.cov.rds', '', basename(nrml.rds))

nrml.tbl <- data.table(
    sample = nrml.smp,
    normal_cov = nrml.rds
)

saveRDS(nrml.tbl, tbl.out)

if (choose.samples == 'all') {
    detergent = prepare_detergent(
    	normal.table.path = tbl.out,
    	path.to.save = dirname(pon.path),
    	num.cores = threads,
    	build = build,
    	PAR.file = par.path,
    	save.pon = TRUE,
    	use.all = TRUE
    )
} else if (choose.samples == 'rand') {
    detergent = prepare_detergent(
        normal.table.path = tbl.out,
        path.to.save = dirname(pon.path),
        num.cores = threads,
        build = build,
        PAR.file = par.path,
        save.pon = TRUE,
        use.all = FALSE,
        choose.randomly = TRUE
    )
} else if (choose.samples == 'cluster') {
    detergent = prepare_detergent(
        normal.table.path = tbl.out,
        path.to.save = dirname(pon.path),
        num.cores = threads,
        build = build,
        PAR.file = par.path,
        save.pon = TRUE,
        number.of.samples = min(100,length(nrml.smp)),
        use.all = FALSE,
        choose.by.clustering = TRUE
    )
}


