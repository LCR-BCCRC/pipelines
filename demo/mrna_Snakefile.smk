#!/usr/bin/env snakemake

'''
This Snakefile is made to run all the modules compatible with Mrna workflow.
Compatibility of a workflow can be checked by referring to the pairing_config parameter present in a default.yaml file of that module.
'''
##### SETUP #####

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")
MRNA = op.filter_samples(SAMPLES, seq_type = "mrna")


##### REFERENCE_FILES WORKFLOW #####


subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


##### CONFIGURATION FILES #####


# Load module-specific configuration

configfile: "../modules/picard_qc/1.0/config/default.yaml"
configfile: "../modules/salmon/1.1/config/default.yaml"
configfile: "../modules/star/1.4/config/default.yaml"
configfile: "../modules/bam2fastq/1.2/config/default.yaml"
configfile: "../modules/manta/2.3/config/default.yaml"
configfile: "../modules/mixcr/1.1/config/default.yaml"
configfile: "../modules/pathseq/1.0/config/default.yaml"
configfile: "../modules/utils/2.1/config/default.yaml"


# Load project-specific config, which includes the shared 
# configuration and some module-specific config updates
configfile: "mrna_config.yaml"


##### CONFIGURATION UPDATES #####


# Use all samples as a default sample list for each module
config["lcr-modules"]["_shared"]["samples"] = MRNA

##### MODULE SNAKEFILES #####


# Load module-specific snakefiles

include: "../modules/picard_qc/1.0/picard_qc.smk"
include: "../modules/salmon/1.1/salmon.smk"
include: "../modules/star/1.4/star.smk"
include: "../modules/bam2fastq/1.2/bam2fastq.smk"
include: "../modules/manta/2.3/manta.smk"
include: "../modules/mixcr/1.1/mixcr.smk"
include: "../modules/pathseq/1.0/pathseq.smk"
include: "../modules/utils/2.1/utils.smk"

##### TARGETS ######

rule all:
    input:
        rules._picard_qc_all.input,
        rules._salmon_all.input,
        rules._star_all.input,
        rules._bam2fastq_all.input,
        rules._manta_all.input,
        rules._mixcr_all.input,
        rules._pathseq_all.input
