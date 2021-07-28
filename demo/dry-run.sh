#!/bin/bash

# Launches a snakefile of your choice in dry run mode (for debugging)
# Usage: ./dry_run.sh <snakefile.smk> <target_rule> "<snakemake_flags>"
# Example: ./dry_run.sh example.smk example_all 
# snakefile.smk The snakefile you want to run
# target_rule: The name of one of the target rules specified in one of the included Snakefiles
# snakemake_flags: One or more flags for the snakemake to run, specified inside quotation marks


# Default to all targets
snakefile=$1
TARGETS=${2:-all}
snakemake_flags=$3

snakemake --dryrun --cores 24 $snakemake_flags -s $snakefile --printshellcmds --reason --use-conda $TARGETS


