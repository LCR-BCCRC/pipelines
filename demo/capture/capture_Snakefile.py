#!/usr/bin/env snakemake

'''
It will only run through the workflow(eg. capture) sub directory.
'''
##### SETUP #####

import oncopipe as op

SAMPLES = op.load_samples("~/lcr-modules/demo/data/samples.tsv")
CAPTURE = op.filter_samples(SAMPLES, seq_type = "capture")


##### REFERENCE_FILES WORKFLOW #####


subworkflow reference_files:
    workdir:
        "../reference/"
    snakefile:
        "../../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../../workflows/reference_files/2.4/config/default.yaml"


##### CONFIGURATION FILES #####


# Load module-specific configuration
configfile: "../../modules/utils/2.1/config/default.yaml"
configfile: "../../modules/picard_qc/1.0/config/default.yaml"
configfile: "../../modules/bam2fastq/1.2/config/default.yaml"
configfile: "../../modules/manta/2.3/config/default.yaml"
configfile: "../../modules/gridss/1.1/config/default.yaml"
#configfile: "../../modules/vcf2maf/1.2/config/default.yaml"
configfile: "../../modules/sequenza/1.4/config/default.yaml"
configfile: "../../modules/strelka/1.1/config/default.yaml"
configfile: "../../modules/bwa_mem/1.1/config/default.yaml"
configfile: "../../modules/lofreq/1.0/config/default.yaml"
configfile: "../../modules/starfish/2.0/config/default.yaml"


# Load project-specific config, which includes the shared 
# configuration and some module-specific config updates
configfile: "capture_config.yaml"


##### CONFIGURATION UPDATES #####


# Use all samples as a default sample list for each module
config["lcr-modules"]["_shared"]["samples"] = CAPTURE

##### MODULE SNAKEFILES #####


# Load module-specific snakefiles

include: "../../modules/utils/2.1/utils.smk"
include: "../../modules/picard_qc/1.0/picard_qc.smk"
include: "../../modules/manta/2.3/manta.smk"
#include: "../../modules/vcf2maf/1.2/vcf2maf.smk"
include: "../../modules/sequenza/1.4/sequenza.smk"
include: "../../modules/strelka/1.1/strelka.smk"
include: "../../modules/bwa_mem/1.1/bwa_mem.smk"
include: "../../modules/gridss/1.1/gridss.smk"
include: "../../modules/bam2fastq/1.2/bam2fastq.smk"
include: "../../modules/lofreq/1.0/lofreq.smk"
include: "../../modules/starfish/2.0/starfish.smk"



##### TARGETS ######

rule all:
    input:
        rules._picard_qc_all.input,
        rules._bam2fastq_all.input,
        rules._manta_all.input,
        rules._sequenza_all.input,
        rules._lofreq_all.input,
        rules._strelka_all.input,
        rules._bwa_mem_all.input,
        rules._gridss_all.input,
        rules._starfish_all.input,
        #rules._vcf2maf_all.input
