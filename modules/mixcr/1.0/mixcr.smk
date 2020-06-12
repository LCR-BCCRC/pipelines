#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Anita dos Santos
# Module Author:    Anita dos Santos
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "mixcr",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "mixcr", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _mixcr_input_fastq,
    _mixcr_step_2,
    _mixcr_output_txt,
    _mixcr_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _mixcr_input_fastq:
    input:
        fastq = CFG["inputs"]["sample_fastq"]
    output:
        fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq"
    run:
        op.relative_symlink(input.fastq, output.fastq)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _mixcr_step_1:
    input:
        tumour_fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{tumour_id}.fastq",
        normal_fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{normal_id}.fastq",
        fasta = reference_files(CFG["reference"]["genome_fasta"])
    output:
        txt = CFG["dirs"]["mixcr"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.txt"
    log:
        stdout = CFG["logs"]["mixcr"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["mixcr"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --tumour {input.tumour_fastq} --normal {input.normal_fastq}
        --ref-fasta {params.fasta} --output {output.txt} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _mixcr_step_2:
    input:
        txt = rules._mixcr_step_1.output.txt
    output:
        txt = CFG["dirs"]["mixcr"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.txt"
    log:
        stderr = CFG["logs"]["mixcr"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.txt} > {output.txt} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _mixcr_output_txt:
    input:
        txt = rules._mixcr_step_2.output.txt
    output:
        txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.txt"
    run:
        op.relative_symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
rule _mixcr_all:
    input:
        expand(
            [
                rules._mixcr_output_txt.output.txt,
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
