#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####

import os

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["bwa_mem"]`
CFG = op.setup_module(
    name = "bwa_mem",
    version = "1.0",
    subdirectories = ["inputs", "bwa_mem",  "sort_bam", "mark_dups", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _bwa_mem_input_fastq,
    _bwa_mem_symlink_bam,
    _bwa_mem_symlink_sorted_bam,
    _bwa_mem_output_bam,
    _bwa_mem_all,

include: "../../utils/2.0/utils.smk"

sample_ids = list(CFG['samples']['sample_id'])

wildcard_constraints: 
    sample_id = "|".join(sample_ids)

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _bwa_mem_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R2.fastq.gz",
    run:
        op.relative_symlink(input.fastq_1, output.fastq_1)
        op.relative_symlink(input.fastq_2, output.fastq_2)


rule _bwa_mem_run:
    input:
        fastq_1 = str(rules._bwa_mem_input_fastq.output.fastq_1),
        fastq_2 = str(rules._bwa_mem_input_fastq.output.fastq_2),
        fastq1_real = CFG["inputs"]["sample_fastq_1"],
        fastq2_real = CFG["inputs"]["sample_fastq_2"],
        fasta = reference_files("genomes/{genome_build}/bwa_index/bwa-0.7.17/genome.fa")
    output:
        sam = pipe(CFG["dirs"]["bwa_mem"] + "{seq_type}--{genome_build}/{sample_id}_out.sam")
    log:
        stderr = CFG["logs"]["bwa_mem"] + "{seq_type}--{genome_build}/{sample_id}/bwa.stderr.log"
    params:
        opts = CFG["options"]["bwa_mem"]
    conda:
        CFG["conda_envs"]["bwa"]
    threads:
        CFG["threads"]["bwa_mem"]
    resources:
        mem_mb = CFG["mem_mb"]["bwa_mem"]
    shell:
        op.as_one_line("""
        bwa mem -t {threads} 
        {params.opts}
        {input.fasta}
        {input.fastq_1} {input.fastq_2} 
        > {output.sam}
        2> {log.stderr}
        """)


rule _bwa_mem_samtools:
    input:
        sam = str(rules._bwa_mem_run.output.sam)
    output:
        bam = CFG["dirs"]["bwa_mem"] + "{seq_type}--{genome_build}/{sample_id}_out.bam"
    log:
        stderr = CFG["logs"]["bwa_mem"] + "{seq_type}--{genome_build}/{sample_id}/samtools.stderr.log"
    params:
        opts = CFG["options"]["samtools"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["samtools"]
    resources:
        mem_mb = CFG["mem_mb"]["samtools"]
    shell:
        op.as_one_line("""
        samtools view {params.opts}
        {input.sam} > {output.bam} 
        2> {log.stderr}
        """)


rule _bwa_mem_symlink_bam:
    input:
        bam = str(rules._bwa_mem_samtools.output.bam)
    output:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


rule _bwa_mem_symlink_sorted_bam:
    input:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam",
        bwa_mem_bam = str(rules._bwa_mem_samtools.output.bam)
    output:
        bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam"
    run:
        op.relative_symlink(input.bam, output.bam)
        os.remove(input.bwa_mem_bam)
        shell("touch {input.bwa_mem_bam}.deleted")


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _bwa_mem_output_bam:
    input:
        bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}" + CFG["options"]["suffix"] + ".bam",
        bai = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}" + CFG["options"]["suffix"] + ".bam.bai",
        sorted_bam = str(rules._bwa_mem_symlink_sorted_bam.input.bam)
    output:
        bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bam + ".bai")
        os.remove(input.sorted_bam)
        shell("touch {input.sorted_bam}.deleted")


# Generates the target sentinels for each run, which generate the symlinks
rule _bwa_mem_all:
    input:
        expand(str(rules._bwa_mem_output_bam.output.bam),
            zip,
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
