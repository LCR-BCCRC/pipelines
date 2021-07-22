#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Nicole Thomas
# Module Author:    Jasper Wong
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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["gatk_rnaseq"]`
CFG = op.setup_module(
    name = "gatk_rnaseq",
    version = "1.0",
    subdirectories = ["inputs", "gatk_splitntrim", "base_recal_report", "gatk_applybqsr", "gatk_variant_calling", "merge_vcfs", "gatk_variant_filtration", "passed", "gnomad_filter", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gatk_rnaseq_input_bam,
    _gatk_rnaseq_output_vcf,
    _gatk_rnaseq_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gatk_rnaseq_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


rule _gatk_splitntrim:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["gatk_splitntrim"] +  "bam/{seq_type}--{genome_build}/{sample_id}.split_reassign_mq.bam")
    log:
        stdout = CFG["logs"]["gatk_splitntrim"] + "{seq_type}--{genome_build}/{sample_id}.gatk_splitntrim.stdout.log",
        stderr = CFG["logs"]["gatk_splitntrim"] + "{seq_type}--{genome_build}/{sample_id}.gatk_splitntrim.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["java_opts"],
        gatk_opts = CFG["options"]["gatk_splitntrim"]
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads:
        CFG["threads"]["gatk_splitntrim"]
    resources:
        **CFG["resources"]["gatk_splitntrim"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m {params.opts}" SplitNCigarReads {params.gatk_opts} 
        -R {input.fasta} -I {input.bam} -O {output.bam}
        > {log.stdout} 2> {log.stderr}
        """)

rule _gatk_addRG:
    input:
        bam = str(rules._gatk_splitntrim.output)
    output:
        bam = temp(CFG["dirs"]["gatk_splitntrim"] + "bam_withRG/{seq_type}--{genome_build}/{sample_id}.withRG.bam")
    params:
        sampleName = "{sample_id}",
        platform = CFG["options"]["gatk_addRG"]["platform"],
        unit = CFG["options"]["gatk_addRG"]["unit"],
        stringency = CFG["options"]["gatk_addRG"]["stringency"]
    conda:
        CFG["conda_envs"]["picard"]
    log:
        stdout = CFG["logs"]["gatk_splitntrim"] + "bam_withRG/{seq_type}--{genome_build}/{sample_id}.addRG.stdout.log"
    threads:
        CFG["threads"]["gatk_addRG"]
    resources:
        **CFG["resources"]["gatk_addRG"]
    shell:
        op.as_one_line("""
        picard AddOrReplaceReadGroups VALIDATION_STRINGENCY={params.stringency} I={input.bam} O={output.bam} RGLB={params.sampleName} RGPL={params.platform} RGPU={params.unit} RGSM={params.sampleName} > {log.stdout} 2>&1 && 
        samtools index {output.bam} 
        """)

rule _gatk_base_recalibration:
    input:
        bam = str(rules._gatk_addRG.output),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        table = CFG["dirs"]["base_recal_report"] + "{seq_type}--{genome_build}/{sample_id}.recalibration_report.grp"
    log:
        stdout = CFG["logs"]["base_recal_report"] + "{seq_type}--{genome_build}/{sample_id}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["base_recal_report"] + "{seq_type}--{genome_build}/{sample_id}.gatk_base_recal.stderr.log"
    params:
        dbsnp = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["java_opts"],
        gatk_opts = CFG["options"]["gatk_baserecalibrator"]
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_base_recalibration"]
    resources:
        **CFG["resources"]["gatk_base_recalibration"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m {params.opts}" BaseRecalibrator {params.gatk_opts} 
        -R {input.fasta} -I {input.bam} -known-sites {params.dbsnp} -known-sites {params.gnomad} -O {output.table} > {log.stdout} 2> {log.stderr}
        """)
        
rule _gatk_applybqsr:
    input:
        bam = str(rules._gatk_addRG.output.bam),
        table = str(rules._gatk_base_recalibration.output.table),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}/{sample_id}.recalibrated.bam"),
        bai = temp(CFG["dirs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}/{sample_id}.recalibrated.bai")
    log:
        stdout = CFG["logs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}/{sample_id}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}/{sample_id}.gatk_base_recal.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["java_opts"],
        gatk_opts = CFG["options"]["gatk_applybqsr"]
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_applybqsr"]
    resources:
        **CFG["resources"]["gatk_applybqsr"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m {params.opts}" ApplyBQSR {params.gatk_opts} 
        -R {input.fasta} -I {input.bam} -bqsr-recal-file {input.table} -O {output.bam} > {log.stdout} 2> {log.stderr}
        """)
        
        
rule _gatk_variant_calling:
    input:
        bam = str(rules._gatk_applybqsr.output.bam),
        bai = str(rules._gatk_applybqsr.output.bai),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG["dirs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.gatk_base_recal.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["java_opts"],
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        min_conf_thres = CFG["options"]["gatk_variant_calling"]["min_conf_thres"],
        gatk_opts = CFG["options"]["gatk_variant_calling"]["gatk_opts"]
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_variant_calling"]
    resources:
        **CFG["resources"]["gatk_variant_calling"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m {params.opts}" HaplotypeCaller {params.gatk_opts} 
        -R {input.fasta} -I {input.bam} -dont-use-soft-clipped-bases -stand-call-conf {params.min_conf_thres} -L {wildcards.chrom}
        -dbsnp {params.vcf} -O {output.vcf} > {log.stdout} 2> {log.stderr}
        """)
        
        
def _gatk_rnaseq_get_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["gatk_rnaseq"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["gatk_variant_calling"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}.{chrom}.vcf.gz",
        chrom = chrs
    )
    return(vcfs)


def _gatk_rnaseq_get_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["gatk_rnaseq"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    tbis = expand(
        CFG["dirs"]["gatk_variant_calling"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}.{chrom}.vcf.gz.tbi",
        chrom = chrs
    )
    return(tbis)


# Merge chromosome gatk_rnaseq VCFs from the same sample
rule _gatk_rnaseq_merge_vcfs:
    input:
        vcf = _gatk_rnaseq_get_chr_vcfs,
        tbi = _gatk_rnaseq_get_chr_tbis
    output:
        vcf = CFG["dirs"]["merge_vcfs"] + "{seq_type}--{genome_build}/{sample_id}.output.vcf.gz",
        tbi = CFG["dirs"]["merge_vcfs"] + "{seq_type}--{genome_build}/{sample_id}.output.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["merge_vcfs"] + "{seq_type}--{genome_build}/{sample_id}.merge_vcfs_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["merge_vcfs"]
    resources:
        **CFG["resources"]["merge_vcfs"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcf} 2> {log.stderr}
            |
        bcftools sort -m {params.mem_mb}M -O z -o {output.vcf} 2>> {log.stderr} 
            &&
        bcftools index -t --threads {threads} {output.vcf} 2>> {log.stderr}
        """)
        
        
rule _gatk_variant_filtration:
    input:
        vcf = str(rules._gatk_rnaseq_merge_vcfs.output.vcf),
        tbi = str(rules._gatk_rnaseq_merge_vcfs.output.tbi),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}/{sample_id}.filtered.vcf.gz",
        tbi = CFG["dirs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}/{sample_id}.filtered.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}/{sample_id}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}/{sample_id}.gatk_base_recal.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["java_opts"],
        window = CFG["options"]["gatk_variant_filtration"]["window"],
        cluster_size = CFG["options"]["gatk_variant_filtration"]["cluster_size"],
        filter_expression = CFG["options"]["gatk_variant_filtration"]["filter_expression"]
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_variant_filtration"]
    resources:
        **CFG["resources"]["gatk_variant_filtration"]
    shell:
        # flag to remove: FS - phred score with strand bias > 30; QD - Variant conf/qual by depth < 2; DP - read depth < 5
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m {params.opts}" VariantFiltration -R {input.fasta} -V {input.vcf} 
        -window {params.window} -cluster-size {params.cluster_size}  
        {params.filter_expression}  
        -O {output.vcf} > {log.stdout} 2> {log.stderr}
        """)


# Filters for PASS variants
rule _gatk_rnaseq_filter_passed:
    input:
        vcf = str(rules._gatk_variant_filtration.output.vcf)
    output:
        vcf = CFG["dirs"]["passed"] + "{seq_type}--{genome_build}/temp/{sample_id}.output.passed.vcf.gz",
        tbi = CFG["dirs"]["passed"] + "{seq_type}--{genome_build}/temp/{sample_id}.output.passed.vcf.gz.tbi"
    params:
        opts = CFG["options"]["gatk_rnaseq_filter_passed"]["params"]
    log:
        stderr = CFG["logs"]["passed"] + "{seq_type}--{genome_build}/{sample_id}.mutect2_filter_passed.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["gatk_rnaseq_passed"]
    resources:
        **CFG["resources"]["gatk_rnaseq_passed"]
    shell:
        op.as_one_line(""" 
        bcftools view {params.opts} -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
            &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Remove AF tag in vcfs:
rule _gatk_rnaseq_removeAF:
    input:
        vcf = str(rules._gatk_rnaseq_filter_passed.output.vcf),
        tbi = str(rules._gatk_rnaseq_filter_passed.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = temp(CFG["dirs"]["gnomad_filter"] + "{seq_type}--{genome_build}/{sample_id}.combined.removeAF.vcf.gz"), 
        tbi = temp(CFG["dirs"]["gnomad_filter"] + "{seq_type}--{genome_build}/{sample_id}.combined.removeAF.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["gnomad_filter"] + "{seq_type}--{genome_build}/{sample_id}.removeAF.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["gnomad_filter"]
    resources:
        **CFG["resources"]["gnomad_filter"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} -x INFO/AF {input.vcf} -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Annotate and filter out common gnomad variants
rule _gatk_rnaseq_annotate_gnomad:
    input:
        vcf = str(rules._gatk_rnaseq_removeAF.output.vcf),
        tbi = str(rules._gatk_rnaseq_removeAF.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["gnomad_filter"] + "{seq_type}--{genome_build}/{sample_id}.combined.gnomad.vcf.gz", 
        tbi = CFG["dirs"]["gnomad_filter"] + "{seq_type}--{genome_build}/{sample_id}.combined.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["gnomad_filter"] + "{seq_type}--{genome_build}/{sample_id}.gnomad_filter.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["gnomad_filter"]
    resources:
        **CFG["resources"]["gnomad_filter"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} -a {input.gnomad} -c INFO/AF {input.vcf} | 
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' | 
        bcftools view -i 'INFO/AF < 0.0001' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)
        
        
# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gatk_rnaseq_output_vcf:
    input:
        vcf = str(rules._gatk_rnaseq_annotate_gnomad.output.vcf),
        tbi = str(rules._gatk_rnaseq_annotate_gnomad.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.output.filt.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.output.filt.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _gatk_rnaseq_all:
    input:
        expand(
            [
                str(rules._gatk_rnaseq_output_vcf.output.vcf),
                str(rules._gatk_rnaseq_output_vcf.output.tbi)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            sample_id=CFG["runs"]["tumour_sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
