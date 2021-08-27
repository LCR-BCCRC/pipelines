#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

cov <- normalizePath(args[1])
junc <- normalizePath(args[2])
outdir <- normalizePath(args[3])
cplex_dir <- args[4]
threads <- args[5]

#Sys.setenv(CPLEX_DIR = cplex_dir)

library(JaBbA)
#setwd(outdir)
JaBbA(junc, cov, 
      field = "foreground", mc.cores = as.integer(threads),
      tfield = 'tier',
      outdir = outdir)
