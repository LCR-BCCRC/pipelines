#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

nrml.tbl.path <- args[1]
drcln.dir     <- args[2]
outdir        <- args[3]
signal.thresh <- args[4]
pct.thresh    <- args[5]
cores         <- args[6]

library(dryclean)
library(parallel)

drcln.dir.abs <- normalizePath(drcln.dir)
drcln.rds     <- Sys.glob(file.path(drcln.dir.abs, "*.rds"))
drcln.smp     <- gsub('.drycleaned.cov.rds', '', basename(drcln.rds))

drcln.tbl <- data.table(
    sample = drcln.smp,
    decomposed_cov = drcln.rds
)

nrml.tbl <- readRDS(nrml.tbl.path)
nrml.tbl <- nrml.tbl[drcln.tbl,,on = 'sample']

saveRDS(nrml.tbl, file.path(dirname(outdir), 'normal_table.rds'))

grm = identify_germline(
    normal.table.path = file.path(dirname(outdir), 'normal_table.rds'),
    path.to.save = dirname(outdir),
    signal.thresh = signal.thresh,
    pct.thresh = pct.thresh,
    save.grm = TRUE,
    num.cores = cores
)

