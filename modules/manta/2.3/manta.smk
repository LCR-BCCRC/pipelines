#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Snakefile Author:    Bruno Grande
# Module Author:                Bruno Grande
# Additional Contributors:      N/A


##### SETUP #####

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

# Setup module and store module-specific configuration in `CFG`.
CFG = op.setup_module(
    name = "manta", 
    version = "2.3",
    subdirectories = ["inputs", "chrom_bed", "manta", "augment_vcf", "bedpe", "outputs"]
)

# Define rules to be run locally when using a compute cluster.
localrules: 
    _manta_input_bam,
    _manta_index_bed,
    _manta_configure_paired,
    _manta_configure_unpaired,
    _manta_output_bedpe,
    _manta_output_vcf, 
    _manta_dispatch,
    _manta_all


##### RULES #####


# Symlinks the input BAM files into the module output directory (under '00-inputs/').
rule _manta_input_bam:
    input:
        sample_bam = CFG["inputs"]["sample_bam"],
        sample_bai = CFG["inputs"]["sample_bai"]
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        sample_crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.sample_bam, output.sample_bam)
        op.absolute_symlink(input.sample_bai, output.sample_bai)
        op.absolute_symlink(input.sample_bai, output.sample_crai)


# bgzip-compress and tabix-index the BED file to meet Manta requirement
rule _manta_index_bed:
    input:
        bed = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    output:
        bedz = CFG["dirs"]["chrom_bed"] + "{genome_build}.main_chroms.bed.gz"
    conda:
        CFG["conda_envs"]["tabix"]
    shell:
        op.as_one_line("""
        bgzip -c {input.bed} > {output.bedz}
            &&
        tabix {output.bedz}
        """)


# Configures the manta workflow with the input BAM files and reference FASTA file.
rule _manta_configure_paired:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        config = op.switch_on_wildcard("seq_type", CFG["switches"]["manta_config"]),
        bedz = str(rules._manta_index_bed.output.bedz)
    output:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"]),
        tumour_bam_arg_name = op.switch_on_wildcard("seq_type", CFG["switches"]["tumour_bam_arg_name"])
    wildcard_constraints:
        pair_status = "matched|unmatched"
    conda:
        CFG["conda_envs"]["manta"]
    shell:
        op.as_one_line("""
        configManta.py {params.opts} --referenceFasta {input.fasta} --callRegions {input.bedz}
        --runDir "$(dirname {output.runwf})" {params.tumour_bam_arg_name} {input.tumour_bam}
        --normalBam {input.normal_bam} --config {input.config} > {log.stdout} 2> {log.stderr}
        """)


# Configures the manta workflow with the input BAM files and reference FASTA file.
rule _manta_configure_unpaired:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        config = op.switch_on_wildcard("seq_type", CFG["switches"]["manta_config"]),
        bedz = str(rules._manta_index_bed.output.bedz)
    output:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"]),
        tumour_bam_arg_name = op.switch_on_wildcard("seq_type", CFG["switches"]["tumour_bam_arg_name"])
    wildcard_constraints:
        pair_status = "no_normal"
    conda:
        CFG["conda_envs"]["manta"]
    shell:
        op.as_one_line("""
        configManta.py {params.opts} --referenceFasta {input.fasta} --callRegions {input.bedz}
        --runDir "$(dirname {output.runwf})" {params.tumour_bam_arg_name} {input.tumour_bam}
        --config {input.config} > {log.stdout} 2> {log.stderr}
        """)


# Launches manta workflow in parallel mode and deletes unnecessary files upon success.
rule _manta_run:
    input:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    output:
        variants_dir = directory(CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/"),
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_run.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_run.stderr.log"
    params:
        opts = CFG["options"]["manta"]
    conda:
        CFG["conda_envs"]["manta"]
    threads:
        CFG["threads"]["manta"]
    resources: 
        mem_mb = CFG["mem_mb"]["manta"],
        bam = 1
    shell:
        op.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log.stdout} 2> {log.stderr}
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


# Calculates the tumour and/or normal variant allele fractions (VAF) from the allele counts
# and fixes the sample IDs in the VCF header to match sample IDs used in Snakemake
rule _manta_augment_vcf:
    input:
        variants_dir = str(rules._manta_run.output.variants_dir)        
    output:
        vcf = CFG["dirs"]["augment_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.augmented.vcf"
    log:
        stdout = CFG["logs"]["augment_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_augment_vcf.{vcf_name}.stdout.log",
        stderr = CFG["logs"]["augment_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_augment_vcf.{vcf_name}.stderr.log"
    params:
        opts = CFG["options"]["augment_vcf"], 
        aug_vcf = CFG["scripts"]["augment_manta_vcf"]
    conda:
        CFG["conda_envs"]["augment_manta_vcf"]
    threads:
        CFG["threads"]["augment_vcf"]
    resources: 
        mem_mb = CFG["mem_mb"]["augment_vcf"]
    shell:
        op.as_one_line("""
        {params.aug_vcf} {params.opts} --tumour_id {wildcards.tumour_id} --normal_id {wildcards.normal_id} 
        --vcf_type {wildcards.vcf_name} {input.variants_dir}/{wildcards.vcf_name}.vcf.gz {output.vcf}
        > {log.stdout} 2> {log.stderr}
        """)


# Converts the VCF file into a more tabular BEDPE file, which is easier to handle in R
# and automatically pairs up breakpoints for interchromosomal events.
rule _manta_vcf_to_bedpe:
    input:
        vcf = str(rules._manta_augment_vcf.output.vcf)
    output:
        bedpe = CFG["dirs"]["bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.bedpe"
    log:
        stderr = CFG["logs"]["bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_vcf_to_bedpe.{vcf_name}.stderr.log"
    conda:
        CFG["conda_envs"]["svtools"]
    threads:
        CFG["threads"]["vcf_to_bedpe"]
    resources: 
        mem_mb = CFG["mem_mb"]["vcf_to_bedpe"]
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log.stderr}"


# Symlinks the augmented VCF files
rule _manta_output_vcf:
    input:
        vcf = str(rules._manta_augment_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module=True)


# Symlinks the final BEDPE files
rule _manta_output_bedpe:
    input:
        bedpe = str(rules._manta_vcf_to_bedpe.output.bedpe)
    output:
        bedpe = CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe, in_module=True)


def _manta_predict_output(wildcards):
    """Request symlinks for all Manta VCF/BEDPE files.
    
    This function is required in conjunction with a Snakemake
    checkpoint because Manta produces different files based
    on whether it's run in paired mode or not and based on
    some parameters (like `--rna`). This function dynamically
    generates target symlinks for the raw VCF files and the
    processed BEDPE files based on what was actually produced.
    """

    # Get module config (since `CFG` is temporary)
    CFG = config["lcr-modules"]["manta"]

    # Standard VCF outputs
    vcf_names = ["candidateSV", "candidateSmallIndels"]

    # Create switch to emulate what is done in `_manta_configure`
    configure_switch = op.switch_on_wildcard(
        "seq_type", CFG["options"]["configure"]
    )
    # And use the switch right away to get the command-line parameters used
    configure_params = configure_switch(wildcards)
    # Check if run with '--rna'
    if "--rna" in configure_params:
        vcf_names.append("rnaSV")
    # Otherwise, check it run in `no_normal` mode
    elif wildcards.pair_status == "no_normal":
        vcf_names.append("tumorSV")
    else:
        vcf_names.append("diploidSV")
        vcf_names.append("somaticSV")


    # Some VCF files can't be converted into BEDPE files
    no_bedpe = ["candidateSV", "candidateSmallIndels"]
    vcf_names_with_bedpe = set(vcf_names) - set(no_bedpe)
    vcf_names_without_bedpe = set(vcf_names) & set(no_bedpe)

    # Request the output files based on whether VAF info is available
    outputs_with_bedpe = expand(
        [
            str(rules._manta_output_vcf.output.vcf),
            str(rules._manta_output_bedpe.output.bedpe),
        ],
        vcf_name=vcf_names_with_bedpe,
        **wildcards
    )

    outputs_without_bedpe = expand(
        str(rules._manta_output_vcf.output.vcf),
        vcf_name=vcf_names_without_bedpe,
        **wildcards
    )

    return outputs_with_bedpe + outputs_without_bedpe


# Generates the target symlinks for each run depending on the Manta output VCF files
rule _manta_dispatch:
    input:
        _manta_predict_output
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _manta_all:
    input:
        expand(
            [
                str(rules._manta_dispatch.output.dispatched),
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
