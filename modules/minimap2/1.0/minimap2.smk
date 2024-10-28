#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Haya Shaalan
# Module Author:    Haya Shaalan
# Contributors:     N/A


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
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["minimap2"]`
CFG = op.setup_module(
    name = "minimap2",
    version = "1.0",
    subdirectories = ["inputs", "minimap2", "sort_bam", "mark_dups","outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _minimap2_input_fastq,
    _minimap2_symlink_bam,
    _minimap2_output_bam,
    _minimap2_all,

sample_ids_minimap2 = list(CFG['samples']['sample_id'])
short_read_samples = op.filter_samples(CFG['samples'], seq_type = "genome").sample_id.tolist()
long_read_samples = op.filter_samples(CFG['samples'], seq_type = "promethION").sample_id.tolist()

##### RULES #####


def _input_fastq(wildcards):
    CFG = config["lcr-modules"]["minimap2"]
    if wildcards.seq_type == "promethION":
        fastqs = CFG["inputs"]["sample_fastq"]["promethION"]
    else:
        fastqs = CFG["inputs"]["sample_fastq"]["genome"]
    return(fastqs)
    

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _minimap2_input_fastq:
    input:
        fastq = _input_fastq
    output:
        fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.fastq_{number}.gz"
    run:
        op.absolute_symlink(input.fastq, output.fastq)


def _get_fastq(wildcards):
    CFG = config["lcr-modules"]["minimap2"]
    if wildcards.seq_type == "promethION":
        fastq = expand(str(rules._minimap2_input_fastq.output.fastq), zip,
        seq_type = wildcards.seq_type,
        sample_id = wildcards.sample_id,
        number = "unpaired")
    else:
        fastq = expand(str(rules._minimap2_input_fastq.output.fastq), zip,
        seq_type = wildcards.seq_type,
        sample_id = wildcards.sample_id,
        number = ["1", "2"])
    return(fastq)


rule _minimap2_run:
    input:
        fastq = _get_fastq,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        sam = pipe(CFG["dirs"]["minimap2"] + "{seq_type}--{genome_build}/{sample_id}_out.sam")
    log:
        stdout = CFG["logs"]["minimap2"] + "{seq_type}--{genome_build}/{sample_id}/minimap2.stdout.log",
        stderr = CFG["logs"]["minimap2"] + "{seq_type}--{genome_build}/{sample_id}/minimap2.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["minimap2"])
    conda:
        CFG["conda_envs"]["minimap2"]
    threads:
        CFG["threads"]["minimap2"]
    resources:
        **CFG["resources"]["minimap2"]
    shell:
        op.as_one_line("""
        minimap2 {params.opts}
        -t {threads} 
        {input.fasta} 
        {input.fastq} 
        -o {output.sam} 
        > {log.stdout} 
        2> {log.stderr}
        """)


rule _minimap2_samtools:
    input:
        sam = str(rules._minimap2_run.output.sam)
    output:
        bam = CFG["dirs"]["minimap2"] + "{seq_type}--{genome_build}/{sample_id}_out.bam",
        complete = touch(CFG["dirs"]["minimap2"] + "{seq_type}--{genome_build}/{sample_id}_out.bam.complete")
    log:
        stderr = CFG["logs"]["minimap2"] + "{seq_type}--{genome_build}/{sample_id}/samtools.stderr.log"
    params:
        opts = CFG["options"]["samtools"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["samtools"]
    resources:
        **CFG["resources"]["samtools"]
    shell:
        op.as_one_line("""
        samtools view {params.opts}
        {input.sam} > {output.bam} 
        2> {log.stderr}
        """)


# Create symlink in subdirectory where BAM files will be sorted by the `utils` module
rule _minimap2_symlink_bam:
    input:
        bam = str(rules._minimap2_samtools.output.bam),
        complete = str(rules._minimap2_samtools.output.complete)
    output:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)


rule _minimap2_symlink_sorted_bam:
    input:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam",
        minimap2_bam = str(rules._minimap2_samtools.output.bam)
    output:
        bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam"
    wildcard_constraints: 
        seq_type = "genome"
    run:
        op.absolute_symlink(input.bam, output.bam)
        os.remove(input.minimap2_bam)
        shell("touch {input.minimap2_bam}.deleted")



# This rule will trigger the `utils` rule for sorting the BAM file
# By this point, the sorted BAM file exists, so this rule deletes the original BAM file
rule _minimap2_LR_bam:
    input:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam",
        bai = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam.bai",
        sorted_bam = str(rules._minimap2_symlink_bam.input.bam)
    output:
        bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    wildcard_constraints: 
        seq_type = "promethION"
    run:
        op.relative_symlink(input.bam, output.bam, in_module=True)
        op.relative_symlink(input.bai, output.bam + ".bai", in_module=True)
        os.remove(input.sorted_bam)
        shell("touch {input.sorted_bam}.deleted")


rule _minimap2_SR_bam:
    input:
        bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam",
        bai = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam.bai",
        sorted_bam = str(rules._minimap2_symlink_bam.input.bam)
    output:
        bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    wildcard_constraints: 
        seq_type = "genome"
    run:
        op.relative_symlink(input.bam, output.bam, in_module=True)
        op.relative_symlink(input.bai, output.bam + ".bai", in_module=True)
        os.remove(input.sorted_bam)
        shell("touch {input.sorted_bam}.deleted")

# Generates the target sentinels for each run, which generate the symlinks
rule _minimap2_all:
    input:
        expand(
            [
                str(rules._minimap2_LR_bam.output.bam),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type= "promethION",
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"]),

        expand(
            [
                str(rules._minimap2_SR_bam.output.bam),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type = "genome",
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
