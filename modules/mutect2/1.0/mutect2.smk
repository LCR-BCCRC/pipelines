#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mutect2"]`
CFG = op.setup_module(
    name = "mutect2",
    version = "1.0",
    subdirectories = ["inputs", "mutect2", "filter", "passed", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mutect2_input_bam,
    _mutect2_input_chrs,
    _mutect2_get_sm,
    _mutect2_output_vcf,
    _mutect2_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutect2_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


# Symlink chromosomes used for parallelization
checkpoint _mutect2_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.relative_symlink(input.chrs, output.chrs)


# Retrieves from SM tag from BAM and writes to file
rule _mutect2_get_sm:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    output:
        sm = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}_sm.txt"),
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}_mutect2_get_sm.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    shell:
        "samtools view -H {input.bam} | grep '^@RG' | "
        r"sed 's/.*SM:\([^\t]*\).*/\1/g'"" | uniq > {output.sm} 2> {log.stderr}"


# Launches Mutect2 in matched and unmatched mode
rule _mutect2_run_matched_unmatched:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        tumour_sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sm.txt",
        normal_sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{normal_id}_sm.txt"
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.tbi"),
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.stats")
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stderr.log"
    params:
        opts = CFG["options"]["mutect2_run"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    resources:
        **CFG["resources"]["mutect2_run"]
    wildcard_constraints: 
        pair_status = "matched|unmatched"
    shell:
        op.as_one_line("""
        gatk Mutect2 {params.opts} -I {input.tumour_bam} -I {input.normal_bam}
        -R {input.fasta} -normal $(cat {input.normal_sm}) -O {output.vcf}
        --germline-resource {input.gnomad} -L {wildcards.chrom} 
        > {log.stdout} 2> {log.stderr}
        """)

# Launces Mutect2 in no normal mode
rule _mutect2_run_no_normal:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        tumour_sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sm.txt"
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.tbi"),
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.output.vcf.gz.stats")
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.mutect2_run.stderr.log"
    params:
        opts = CFG["options"]["mutect2_run"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    resources:
        **CFG["resources"]["mutect2_run"]
    wildcard_constraints: 
        pair_status = "no_normal"
    shell:
        op.as_one_line("""
        gatk Mutect2 {params.opts} -I {input.tumour_bam} -R {input.fasta} 
        -O {output.vcf} --germline-resource {input.gnomad} -L {wildcards.chrom}
        > {log.stdout} 2> {log.stderr}
        """)


def _mutect2_get_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.output.vcf.gz",
        chrom = chrs
    )
    return(vcfs)


def _mutect2_get_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    tbis = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.output.vcf.gz.tbi",
        chrom = chrs
    )
    return(tbis)


# Merge chromosome mutect2 VCFs from the same sample
rule _mutect2_merge_vcfs:
    input:
        vcf = _mutect2_get_chr_vcfs,
        tbi = _mutect2_get_chr_tbis
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["mutect2_merge_vcfs"]
    resources:
        **CFG["resources"]["mutect2_merge_vcfs"]
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcf} 2> {log.stderr}
            |
        bcftools sort -O z -o {output.vcf} 2>> {log.stderr} 
            &&
        bcftools index -t --threads {threads} {output.vcf} 2>> {log.stderr}
        """)


def _mutect2_get_chr_stats(wildcards):
    CFG = config["lcr-modules"]["mutect2"]
    chrs = checkpoints._mutect2_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    stats = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.output.vcf.gz.stats",
        chrom = chrs
    )
    return(stats)


# Merge chromosome mutect2 stats for FilterMutectCalls rule
rule _mutect2_merge_stats:
    input:
        stat = _mutect2_get_chr_stats
    output:
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz.stats")
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_merge_stats.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_merge_stats.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    shell:
        op.as_one_line("""
        gatk MergeMutectStats $(for i in {input.stat}; do echo -n "-stats $i "; done)
        -O {output.stat} > {log.stdout} 2> {log.stderr}
        """)


# Marks variants filtered or PASS annotations
rule _mutect2_filter:
    input:
        vcf = str(rules._mutect2_merge_vcfs.output.vcf),
        tbi = str(rules._mutect2_merge_vcfs.output.tbi),
        stat = str(rules._mutect2_merge_stats.output.stat),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.unfilt.vcf.gz"
    log:
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter.stdout.log",
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter.stderr.log"
    params:
        opts = CFG["options"]["mutect2_filter"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_filter"]
    resources:
        **CFG["resources"]["mutect2_filter"]
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls {params.opts} -V {input.vcf} -R {input.fasta}
        -O {output.vcf} > {log.stdout} 2> {log.stderr}
        """)


# Filters for PASS variants
rule _mutect2_filter_passed:
    input:
        vcf = str(rules._mutect2_filter.output.vcf)
    output:
        vcf = CFG["dirs"]["passed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.passed.vcf.gz",
        tbi = CFG["dirs"]["passed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.passed.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["passed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter_passed.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["mutect2_passed"]
    resources:
        **CFG["resources"]["mutect2_passed"]
    shell:
        op.as_one_line(""" 
        bcftools view -f '.,PASS' -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
            &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutect2_output_vcf:
    input:
        vcf = str(rules._mutect2_filter_passed.output.vcf),
        tbi = str(rules._mutect2_filter_passed.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.passed.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.passed.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)


# Generates the target sentinels for each run, which generate the symlinks
rule _mutect2_all:
    input:
        expand(
            [
                str(rules._mutect2_output_vcf.output.vcf),
                str(rules._mutect2_output_vcf.output.tbi)
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
