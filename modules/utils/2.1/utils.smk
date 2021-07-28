#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Snakefile Author:    Helena Winata, Bruno Grande
# Module Author:                Bruno Grande
# Additional Contributors:      N/A


##### SETUP #####

import sys, os
from os.path import join

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

_UTILS = config["lcr-modules"]["utils"]
LOG = "/logs/" + op._session.launched_fmt

wildcard_constraints:
    prefix = "[0-9]{2}-.*"


##### RULES #####
# _utils_bam_sort: Sort a BAM file using coordinates
rule _utils_bam_sort:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.sort.bam", 
        prefix = temp(directory("{out_dir}/{prefix}/{suffix}_tmp"))
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_sort.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_sort.stderr.log"
    wildcard_constraints:
        prefix = ".*(sort).*"
    params:
        opts = _UTILS["options"]["bam_sort"],
        memory= lambda wildcards, resources, threads: int(resources.mem_mb/threads/2)
    conda:
        _UTILS["conda_envs"]["samtools"]
    threads:
        _UTILS["threads"]["bam_sort"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["bam_sort"]
    shell:
        op.as_one_line("""
        mkdir -p {output.prefix} &&
        samtools sort 
        -@ {threads} -m $(({params.memory}))M
        {params.opts}
        -T {output.prefix} -o {output.bam} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


# _utils_bam_markdups: Mark duplicates in a BAM file using Picard criteria
rule _utils_bam_markdups:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.mdups.bam"
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_mark_dups.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_mark_dups.stderr.log"
    wildcard_constraints:
        prefix = ".*(mark_dups).*"
    params:
        opts = _UTILS["options"]["bam_markdups"]
    conda:
        _UTILS["conda_envs"]["sambamba"]
    threads:
        _UTILS["threads"]["bam_markdups"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["bam_markdups"]
    shell:
        op.as_one_line("""
        sambamba markdup 
        {params.opts} 
        --nthreads {threads} 
        {input.bam} {output.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


# _utils_bam_index: Index a BAM file
rule _utils_bam_index:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.bam.bai"
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_index.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_index.stderr.log"
    params:
        opts = _UTILS["options"]["bam_index"]
    conda:
        _UTILS["conda_envs"]["samtools"]
    threads:
        _UTILS["threads"]["bam_index"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["bam_index"]
    shell:
        op.as_one_line("""
        samtools index 
        {params.opts} 
        -@ {threads} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


rule _utils_create_intervals: # create_interval_list from bed; default exomes
    input:
        bed = lambda w: _UTILS["inputs"]["bed"][w.genome_build],
        seq_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output: 
        intervals = "reference/exomes/{genome_build}/interval/{id}_intervals.txt"
    log: 
        "reference/exomes/{genome_build}/interval/{id}_intervals.log"
    conda: 
        _UTILS["conda_envs"]["picard"]
    threads:
        _UTILS["threads"]["interval"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["interval"]
    shell:
        op.as_one_line("""
        picard BedToIntervalList
        I={input.bed}
        O={output.intervals}
        SD={input.seq_dict}
        2> {log}
        &&
        chmod a-w {output.intervals}
        """)



