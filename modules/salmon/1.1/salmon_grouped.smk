#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Helena Winata
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd

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

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["salmon"]`
CFG = op.setup_module(
    name = "salmon",
    version = "1.1",
    subdirectories = ["inputs", "salmon", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _salmon_input_fastq,
    _salmon_output,
    _salmon_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _salmon_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"]
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.read1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.read2.fastq.gz"
    group: 
        CFG["group"]["quant"]
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)


rule _salmon_quant:
    input:
        fastq_1 = rules._salmon_input_fastq.output.fastq_1,
        fastq_2 = rules._salmon_input_fastq.output.fastq_2,
        fastq1_real = CFG["inputs"]["sample_fastq_1"], # Placeholders to prevent premature deletion of temp fastqs
        fastq2_real = CFG["inputs"]["sample_fastq_2"],
        index = reference_files(expand("genomes/{genome_build}/salmon_index/salmon-1.3.0/index", genome_build = CFG["transcriptome"]["quant_to"]))
    output:
        quant = expand(CFG["dirs"]["salmon"] + "quant_to_{quant_to}/{{seq_type}}/{{sample_id}}/quant.sf", quant_to=CFG["transcriptome"]["quant_to"])
    log:
        stdout = expand(CFG["logs"]["salmon"] + "quant_to_{quant_to}/{{seq_type}}/{{sample_id}}/quant.stdout.log", quant_to=CFG["transcriptome"]["quant_to"]),
        stderr = expand(CFG["logs"]["salmon"] + "quant_to_{quant_to}/{{seq_type}}/{{sample_id}}/quant.stderr.log", quant_to=CFG["transcriptome"]["quant_to"])
    params:
        opts = CFG["options"]["quant"],
        quant_to = CFG["transcriptome"]["quant_to"]
    conda:
        CFG["conda_envs"]["salmon"]
    threads:
        CFG["threads"]["quant"]
    resources:
        **CFG["resources"]["quant"]
    group: 
        CFG["group"]["quant"]
    shell:
        op.as_one_line("""
        salmon quant -p {threads} 
        {params.opts} 
        -i {input.index}
        -o $(dirname {output.quant})
        -1 {input.fastq_1} -2 {input.fastq_2} 
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _salmon_output:
    input:
        quant = rules._salmon_quant.output.quant
    output:
        quant = CFG["dirs"]["outputs"] + "quant_to_" + CFG["transcriptome"]["quant_to"] + "/{seq_type}/{sample_id}.quant.sf"
    run:
        op.relative_symlink(input.quant, output.quant, in_module=True)


rule export_sample_table:
    input:
        quant = expand(rules._salmon_output.output.quant,
            zip,
            seq_type=CFG["samples"]["seq_type"],
            sample_id=CFG["samples"]["sample_id"])
    output:
        sample_table = CFG["dirs"]["outputs"] + "quant_to_" + CFG["transcriptome"]["quant_to"] + "_matrix/{seq_type}/sample_table.tsv"
    run:
        samples=pd.DataFrame(config["lcr-modules"]["salmon"]["samples"])
        samples.to_csv(output.sample_table, sep='\t', header=True, index=False)

def get_gtf(wildcards):
    if config["lcr-modules"]["salmon"]["transcriptome"]["quant_to"] == "hg38":
        return reference_files("downloads/gencode-33/gencode.annotation.grch38.gtf")
    else:
        return reference_files("downloads/gencode-33/gencode.annotation.grch37.gtf")

rule build_counts_matrix:
    input:
        quant = expand(rules._salmon_output.output.quant,
            zip,
            seq_type=CFG["samples"]["seq_type"],
            sample_id=CFG["samples"]["sample_id"]),
        sample_table = str(rules.export_sample_table.output.sample_table),
        salmon2counts = CFG["inputs"]["salmon2counts"],
        gtf = get_gtf
    output:
        counts_matrix = CFG["dirs"]["outputs"] + "quant_to_" + CFG["transcriptome"]["quant_to"] + "_matrix/{seq_type}/salmon.genes.counts.tsv"
    params:
        path = CFG["dirs"]["outputs"]+"quant_to_" + CFG["transcriptome"]["quant_to"]+"/{seq_type}",
        out_dir = directory(CFG["dirs"]["outputs"] + "quant_to_" + CFG["transcriptome"]["quant_to"] + "_matrix/{seq_type}/")
    conda:
        CFG["conda_envs"]["salmon2counts"]
    resources:
        **CFG["resources"]["matrix"]
    shell:
        op.as_one_line("""
        Rscript {input.salmon2counts}
        {params.path}
        {input.gtf}
        {params.out_dir}
        {input.sample_table}
        sample_id
        """)


# Generates the target sentinels for each run, which generate the symlinks
rule _salmon_all:
    input:
        expand(
            rules.build_counts_matrix.output.counts_matrix,
            seq_type=CFG["samples"]["seq_type"].unique().tolist()
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
