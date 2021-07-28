#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     Ryan Morin


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

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["strelka"]`
CFG = op.setup_module(
    name = "strelka",
    version = "1.1",
    subdirectories = ["inputs", "chrom_bed", "strelka", "filtered", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _strelka_input_bam,
    _strelka_input_vcf,
    _strelka_dummy_vcf,
    _strelka_index_bed,
    _strelka_configure_paired,
    _strelka_configure_unpaired,
    _strelka_filter_combine,
    _strelka_output_filtered_vcf,
    _strelka_all,

wildcard_constraints: 
    var_type = "somatic.snvs|somatic.indels|variants"

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _strelka_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


rule _strelka_dummy_vcf:
    # creates a dummy vcf if users do not specify candidateSmallIndels file
    output:
        touch(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.dummy.tbi")


rule _strelka_input_vcf:
    input:
        manta_vcf = CFG["inputs"]["candidate_small_indels"]
    output:
        vcf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.candidateSmallIndels.vcf.gz",
        tbi = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.candidateSmallIndels.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["tabix"]
    shell:
        op.as_one_line("""
        bgzip -c {input.manta_vcf} > {output.vcf}
            &&
        tabix {output.vcf}
        """)


# bgzip-compress and tabix-index the BED file to meet strelka requirement
rule _strelka_index_bed:
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


def _strelka_get_indel_cli_arg(vcf_in = config["lcr-modules"]["strelka"]["inputs"]["candidate_small_indels"]):
    def _strelka_get_indel_cli_custom(wildcards, input):
        if vcf_in:
            param = f"--indelCandidates={input.indels}"
        else: 
            param = ""
        return param
    return _strelka_get_indel_cli_custom


rule _strelka_configure_paired: # Somatic
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bedz = str(rules._strelka_index_bed.output.bedz),
        indels = str(rules._strelka_input_vcf.output.vcf) if CFG["inputs"]["candidate_small_indels"] else str(rules._strelka_dummy_vcf.output)
    output:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}-{pair_status}/strelka_configure.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stderr.log"
    params:
        indel_arg = _strelka_get_indel_cli_arg(),
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"]),
    wildcard_constraints:
        pair_status = "matched|unmatched"
    conda:
        CFG["conda_envs"]["strelka"]
    shell:
        op.as_one_line("""
        configureStrelkaSomaticWorkflow.py 
        --normalBam={input.normal_bam}
        --tumorBam={input.tumour_bam}
        --referenceFasta={input.fasta}
       --callRegions={input.bedz}
        --runDir=$(dirname {output.runwf})
        {params.indel_arg}
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _strelka_configure_unpaired: # germline
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bedz = str(rules._strelka_index_bed.output.bedz)
    output:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"])
    message:
        "WARNING: {wildcards.seq_type} sample ({wildcards.tumour_id}) is being processed using Strelka Germline workflow. Ensure pairing config for capture is set to run_unpaired_tumours_with: 'unmatched_normal' to run Strelka Somatic workflow"
    wildcard_constraints:
        pair_status = "no_normal"
    conda:
        CFG["conda_envs"]["strelka"]
    shell:
        op.as_one_line("""
        configureStrelkaGermlineWorkflow.py 
        --bam={input.tumour_bam}
        --referenceFasta={input.fasta}
        --callRegions={input.bedz}
        --runDir=$(dirname {output.runwf})
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _strelka_run_unpaired:
    input:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    output:
        vcf_variants = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/variants.vcf.gz"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_run.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_run.stderr.log"
    params:
        opts = CFG["options"]["strelka"]
    conda:
        CFG["conda_envs"]["strelka"]
    threads:
        CFG["threads"]["strelka"]
    resources: 
        mem_mb = op.retry(CFG["mem_mb"]["strelka"], 2),
        bam = 1
    shell:
        #TODO: separate job from cleanup because of errors that can occur during deletion of the tmp directories. This should not cause Strelka to rerun but it does because Snakemake deletes the output files
        op.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log.stdout} 2> {log.stderr}
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)

rule _strelka_run_paired:
    input:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    output:
        vcf_snvs = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somatic.snvs.vcf.gz",
        vcf_indels = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somatic.indels.vcf.gz"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_run.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_run.stderr.log"
    params:
        opts = CFG["options"]["strelka"]
    conda:
        CFG["conda_envs"]["strelka"]
    threads:
        CFG["threads"]["strelka"]
    resources: 
        mem_mb = op.retry(CFG["mem_mb"]["strelka"], 2),
        bam = 1
    shell:
        #TODO: separate job from cleanup because of errors that can occur during deletion of the tmp directories. This should not cause Strelka to rerun but it does because Snakemake deletes the output files
        op.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log.stdout} 2> {log.stderr}
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)

# Combine and filter for PASS variants
rule _strelka_filter_combine:
    input:
        vcf = expand(CFG["dirs"]["strelka"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/results/variants/{var_type}.vcf.gz", 
                    var_type = ["somatic.indels", "somatic.snvs"])
    output:
        vcf = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/combined.passed.vcf.gz",
        vcf_tbi = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/combined.passed.vcf.gz.tbi"
    resources: 
        mem_mb = CFG["mem_mb"]["bcftools_sort"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["bcftools"]
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_filter_combine.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_filter_combine.stderr.log"
    shell:
        op.as_one_line("""
        bcftools concat -a {input.vcf} | 
        bcftools view -f ".,PASS" -Ov | 
        bcftools sort --max-mem {params.mem_mb}M -Oz -o {output.vcf}
        > {log.stdout} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr} 
        """)


# infers name of output files depending on how Strelka was run
def _strelka_get_output(wildcards):
    CFG = config["lcr-modules"]["strelka"]

    if wildcards.pair_status == "no_normal":
        vcf = expand(CFG["dirs"]["filtered"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{var_type}.passed.vcf.gz", 
                    var_type = "variants"
        )
    else:
        vcf = str(rules._strelka_filter_combine.output.vcf)
    return vcf

# Symlinks the final output files into the module results directory (under '99-outputs/'). Links will always use "combined" in the name (dropping odd naming convention used by Strelka in unpaired mode)
rule _strelka_output_filtered_vcf:
    input:
        vcf = _strelka_get_output
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz",
        vcf_tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module=True)
        op.relative_symlink(str(input.vcf) + ".tbi", output.vcf_tbi, in_module=True)



rule _strelka_dispatch:
    input: 
        vcf = str(rules._strelka_output_filtered_vcf.output.vcf)
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _strelka_all:
    input:
        expand(str(rules._strelka_dispatch.output.dispatched), zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"])

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
