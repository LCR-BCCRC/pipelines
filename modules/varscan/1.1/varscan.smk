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
# `CFG` is a shortcut to `config["lcr-modules"]["varscan"]`
CFG = op.setup_module(
    name = "varscan",
    version = "1.1",
    subdirectories = ["inputs", "mpileup", "varscan", "maf", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _varscan_input_bam,
    _varscan_input_chroms,
    _varscan_reheader_vcf,
    _varscan_combine_chroms_vcf,
    _varscan_combine_vcf,
    _varscan_symlink_maf,
    _varscan_output_vcf,
    _varscan_output_maf,
    _varscan_all,

ruleorder: _varscan_reheader_vcf > _varscan_combine_vcf

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _varscan_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bam + ".bai")


# Pulls in list of chromosomes for the genome builds
checkpoint _varscan_input_chroms:
    input:
        txt = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        txt = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)

# generate mpileups for tumour and normal bams separately. 
# If we parallelize this by chromosome we will need 2 * 2 threads per chromosome but this should be a lot more efficient
rule _varscan_bam2mpu:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        mpu = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.mpileup")
    log:
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam2mpu.{chrom}.stderr.log"
    params:
        opts = CFG["options"]["mpileup"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        1    #hardcoded because samtools mpileup does not support more than one thread.
    resources:
        mem_mb = CFG["mem_mb"]["mpileup"]
    shell:
        op.as_one_line("""
        samtools mpileup {params.opts}
        -r {wildcards.chrom}
        -f {input.fasta} {input.bam}
        > {output.mpu}
        2> {log.stderr}
        """)


rule _varscan_somatic:
    input:
        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.{chrom}.mpileup",
        normal_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.{chrom}.mpileup"
    output:
        snp = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_snp.vcf"),
        indel = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_indel.vcf")
    wildcard_constraints:
        pair_status = "matched|unmatched"
    log:
        stdout = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_somatic.stdout.log",
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_somatic.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["somatic"])
    conda:
        CFG["conda_envs"]["varscan"]
    threads:
        CFG["threads"]["somatic"]   #this seems to rarely exceed 300% due to samtools and/or I/O restrictions
    resources:
        mem_mb = CFG["mem_mb"]["somatic"]
    shell:
        op.as_one_line("""
        varscan somatic 
        {input.normal_mpu} {input.tumour_mpu} 
        --output-snp {output.snp} --output-indel {output.indel}
        {params.opts}
        > {log.stdout} 2> {log.stderr} 
        || true
        """)

rule _varscan_unpaired:
    input:
        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.{chrom}.mpileup"
    output:
        vcf = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf")
    wildcard_constraints:
        pair_status = "no_normal",
        vcf_name = "(snp|indel)"
    log:
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_{vcf_name}.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["unpaired"])
    conda:
        CFG["conda_envs"]["varscan"]
    threads:
        CFG["threads"]["unpaired"]
    resources:
        mem_mb = CFG["mem_mb"]["unpaired"]
    shell:
        op.as_one_line("""
        varscan mpileup2{wildcards.vcf_name} 
        {input.tumour_mpu} 
        {params.opts}
        > {output.vcf} 2> {log.stderr}
        """)


rule _varscan_reheader_vcf:
    input:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf",
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf.gz",
        #tbi= CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        contig=$( awk '{{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}}' {input.fai})
            &&
        awk -v var="$contig" '/^#CHROM/ {{ print(var); }} {{print;}}' {input.vcf} |
        bgzip -c > {output.vcf}
            &&
        tabix -p vcf {output.vcf}
        """)

        #bcftools reheader -h {input.header} {input.vcf} | bcftools view -l 1 -o {output.vcf}


def _varscan_request_chrom_vcf(wildcards):
    CFG = config["lcr-modules"]["varscan"]
    with open(checkpoints._varscan_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    vcf_files = expand(str(rules._varscan_reheader_vcf.output.vcf),
        chrom = mains_chroms, **wildcards
    )
    return vcf_files


rule _varscan_combine_chroms_vcf:
    input:
        vcf = _varscan_request_chrom_vcf
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz",
        #tbi= CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}_combine_chroms_vcf.stdout.log",
        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}_combine_chroms_vcf.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    wildcard_constraints:
        vcf_name = "(snp|indel)"
    resources:
        mem_mb = CFG["mem_mb"]["bcftools_sort"]
    shell:
        op.as_one_line(""" 
        bcftools concat -a {input.vcf} | 
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf} 
        > {log.stdout} 2> {log.stderr} 
            &&
        tabix -p vcf {output.vcf} 
        >> {log.stdout} 2>> {log.stderr}
        """)

''' 
bcftools concat -o {output.vcf} {input.vcf} 
    && 
bgzip -c {output.vcf} >> {output.vcf_gz} 
    &&
bcftools index {output.vcf_gz}
'''

def _varscan_request_vcf(wildcards):
    return expand([
        str(rules._varscan_combine_chroms_vcf.output.vcf)
        ],
        vcf_name = ["snp", "indel"], **wildcards)


rule _varscan_combine_vcf:
    input:
        vcf = _varscan_request_vcf
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.vcf.gz",
        tbi= CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged_combine_vcf.stdout.log",
        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged_combine_vcf.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        mem_mb = CFG["mem_mb"]["bcftools_sort"]
    shell:
        op.as_one_line("""
        bcftools concat -a {input.vcf} | 
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf} 
        > {log.stdout} 2> {log.stderr} 
            &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# symlink vcf file to maf directory to run vcf2maf
rule _varscan_symlink_maf:
    input:
        vcf = str(rules._varscan_combine_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.vcf.gz"
    run:
        op.absolute_symlink(input.vcf, output.vcf)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _varscan_output_vcf:
    input:
        vcf = str(rules._varscan_combine_vcf.output.vcf),
        tbi = str(rules._varscan_combine_vcf.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-merged.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-merged.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module=True)
        op.relative_symlink(input.tbi, output.tbi, in_module=True)


rule _varscan_output_maf:
    input:
        #vcf = str(rules._varscan_combine_vcf.output.vcf), # ensure vcf is not deleted before maf is created
        maf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.maf"
    output:
        maf = CFG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-merged.maf"
    run:
        op.relative_symlink(input.maf, output.maf, in_module=True)

"""
def _varscan_get_output(wildcards):
    return expand([
        str(rules._varscan_output_vcf.output.vcf),
        str(rules._varscan_output_maf.output.maf)
        ],
        vcf_name = ["snp", "indel"], **wildcards)


rule _varscan_dispatch:
    input:
        _varscan_get_output
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")
"""

# Generates the target sentinels for each run, which generate the symlinks
rule _varscan_all:
    input:
        expand(
            [
                #str(rules._varscan_dispatch.output.dispatched),
            str(rules._varscan_output_vcf.output.vcf),
            #str(rules._varscan_output_maf.output.maf)
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
