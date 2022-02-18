#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Jasper Wong
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cnvkit"]`
CFG = op.setup_module(
    name = "cnvkit",
    version = "1.0",
    subdirectories = ["inputs", "coverage", "fix", "cns", "SNPs", "BAF", "plots", "breaks", "geneMetrics", "seg", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cnvkit_input_bam,
    _cnvkit_output_no_capture_space,
    _cnvkit_output,
    _cnvkit_all,


##### RULES #####

# cnvkit reference files
# gene annotation files
rule _get_refFlat:
    output:
        refFlat = CFG["dirs"]["inputs"] + "{genome_build}/refFlat.final.txt"
    params:
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/",
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19",
        txt = CFG["dirs"]["inputs"] + "{genome_build}/refFlat.txt",
        prefix = lambda w: "chr" if "hg" in str({w.genome_build}) else "",
    shell:
        """
            wget {params.url}{params.build}/database/refFlat.txt.gz -O - | gzip -d > {params.txt} &&
            sed 's/chr/{params.prefix}/g' {params.txt} > {output.refFlat}
        """
    

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cnvkit_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)


rule _cnvkit_accessible_regions:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        access = CFG["dirs"]["inputs"] + "reference/access.{genome_build}.bed"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    log:
        stdout = CFG["logs"]["inputs"] + "{genome_build}_access.log"
    shell:
        """
            cnvkit.py access {input.fasta} -o {output.access} &> {log.stdout}
        """


# filters out chrG, chrJ, chrM from bed
rule _cnvkit_filter_main_chrs:
    input:
        access = CFG["dirs"]["inputs"] + "reference/access.{genome_build}.bed"
    output:
        access_main = CFG["dirs"]["inputs"] + "reference/access_main.{genome_build}.bed"
    shell:
        """
            grep -v GL {input.access} | grep -v J | grep -v M > {output.access_main} 
        """


def _cnvkit_get_normals(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    capture_space = CFG["runs"][CFG["runs"]["normal_capture_space"].isin([wildcards.capture_space])]
    capture_space = capture_space[capture_space["normal_seq_type"].isin([wildcards.seq_type])]
    normals = expand(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}.bam", 
                    zip,
                    capture_space = capture_space["normal_capture_space"],
                    seq_type = capture_space["normal_seq_type"],
                    genome_build = capture_space["normal_genome_build"],
                    normal_id = capture_space["normal_sample_id"])
    normals = list(dict.fromkeys(normals))
    return normals


def _cnvkit_get_capspace(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    # default is just taking values from genome_build and capture_space in metadata
    default_bed = reference_files("genomes/" + wildcards.genome_build + "/capture_space/" + wildcards.capture_space + ".bed") 
    # if there is a custom bed in the config that pertains to the capture_space variable, then use that
    if str(wildcards.capture_space) in CFG["options"]["target_bed"].keys():
        bed = CFG["options"]["target_bed"][wildcards.capture_space]
    # if there is "none" instead, then we go towards the default, which is the full exon-space
    elif "none" in str({wildcards.capture_space}):
        if "38" in str({wildcards.genome_build}):
            bed = reference_files("genomes/" + wildcards.genome_build + "/capture_space/exome-utr-grch38.bed")
        else:
            bed = reference_files("genomes/" + wildcards.genome_build + "/capture_space/exome-utr-grch37.bed")
    # if it's anything else, then we can assume it will just use the default, which is what can be taken from the metadata and found via reference_files()
    else:
        bed = default_bed
    return bed
    

rule _cnvkit_build_access_bed:
    input:
        access = CFG["dirs"]["inputs"] + "reference/access_main.{genome_build}.bed",
        bam = _cnvkit_get_normals,
        targets = _cnvkit_get_capspace,
        refFlat = str(rules._get_refFlat.output.refFlat)
    output:
        target = CFG["dirs"]["inputs"] + "reference/{seq_type}--{genome_build}/{capture_space}/target_sites.target.bed",
        antitarget = CFG["dirs"]["inputs"] + "reference/{seq_type}--{genome_build}/{capture_space}/target_sites.antitarget.bed"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    log:
        stdout = CFG["logs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}_autobin.log"
    shell:
        """
            cnvkit.py autobin {input.bam} -t {input.targets} -g {input.access} --annotate {input.refFlat} --short-names --target-output-bed {output.target} --antitarget-output-bed {output.antitarget} &> {log.stdout}
        """


# Coverage for each sample
rule _coverage_target:
    input:
        bam = str(rules._cnvkit_input_bam.output.bam),
        bed = str(rules._cnvkit_build_access_bed.output.target),
    output:
        cov = CFG["dirs"]["coverage"] + "target/{seq_type}--{genome_build}/{capture_space}/{sample_id}.targetcoverage.cnn"
    conda: CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    log:
        stdout = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_target.log"
    shell:
        """
            cnvkit.py coverage {input.bam} {input.bed} -o {output.cov} -p {threads} &> {log.stdout}
        """

rule _coverage_antitarget:
    input:
        bam = str(rules._cnvkit_input_bam.output.bam),
        bed = str(rules._cnvkit_build_access_bed.output.antitarget),
    output:
        cov = CFG["dirs"]["coverage"] + "antitarget/{seq_type}--{genome_build}/{capture_space}/{sample_id}.antitargetcoverage.cnn"
    conda: CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    log:
        stdout = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_antitarget.log"
    shell:
        """
            cnvkit.py coverage {input.bam} {input.bed} -o {output.cov} -p {threads} &> {log.stdout}
        """


# For NORMALS
def get_normals_target(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    runs = CFG["runs"]
    platform = runs[runs['tumour_capture_space'].isin([wildcards.capture_space])]
    normals = expand(CFG["dirs"]["coverage"] +  "target/{seq_type}--{genome_build}/{capture_space}/{normal_id}.targetcoverage.cnn", 
                        zip,
                        seq_type = platform['tumour_seq_type'],
                        genome_build = platform['tumour_genome_build'],
                        normal_id = platform["normal_sample_id"],
                        capture_space=platform["tumour_capture_space"])
    normals = list(dict.fromkeys(normals))
    return normals

def get_normals_anti(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    runs = CFG["runs"]
    platform = runs[runs['tumour_capture_space'].isin([wildcards.capture_space])]
    normals = expand(CFG["dirs"]["coverage"] +  "antitarget/{seq_type}--{genome_build}/{capture_space}/{normal_id}.antitargetcoverage.cnn", 
                        zip,
                        seq_type = platform['tumour_seq_type'],
                        genome_build = platform['tumour_genome_build'],
                        normal_id = platform["normal_sample_id"],
                        capture_space=platform["tumour_capture_space"])
    normals = list(dict.fromkeys(normals))
    return normals


if CFG["options"]["new_normals"] == True:
    rule _create_pon_reference:
        input:
            control_target = get_normals_target,
            control_antitarget = get_normals_anti,
        output:
            pon = CFG["dirs"]["coverage"] +  "normal/{seq_type}--{genome_build}/{capture_space}/normal_reference.cnn"
        params:
            fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
            male_reference = CFG["options"]["male_ref"]
        conda: CFG["conda_envs"]["cnvkit"]
        threads:
            CFG["threads"]["reference"]
        resources:
            **CFG["resources"]["reference"]
        log:
            stdout = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/normal_antitarget.log"
        shell:
            """
                cnvkit.py reference {input.control_target} {input.control_antitarget} --fasta {params.fasta} -o {output.pon} {params.male_reference} &> {log.stdout}
            """


rule _cnvkit_fix:
    input:
        targetcov = CFG["dirs"]["coverage"]  + "target/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.targetcoverage.cnn",
        antitargetcov = CFG["dirs"]["coverage"]  + "antitarget/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.antitargetcoverage.cnn",
        pon_reference = CFG["dirs"]["coverage"]  + "normal/{seq_type}--{genome_build}/{capture_space}/normal_reference.cnn"
    output:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr"
    conda: CFG["conda_envs"]["cnvkit"]
    resources:
        **CFG["resources"]["fix"]
    log:
        stdout = CFG["logs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        """
            cnvkit.py fix {input.targetcov} {input.antitargetcov} {input.pon_reference} -o {output.cnr} &> {log.stdout}
        """

rule _cnvkit_segment:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr"
    output:
        cns = CFG["dirs"]["cns"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cns"
    params:
        method = CFG["options"]["cns"]["method"]
    conda: CFG["conda_envs"]["cnvkit"]
    resources:
        **CFG["resources"]["cns"]
    threads: CFG["threads"]["cns"]
    log: 
        stdout = CFG["logs"]["cns"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.segment.log"
    shell:
        """
            cnvkit.py segment {input.cnr} -o {output.cns} -p {threads} --drop-low-coverage -m {params.method} &> {log.stdout}
        """


# need SNPs not SNVs (i.e. get germline calls using a dbSNP vcf)
rule _cnvkit_dbsnp_to_bed:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        bed = CFG["dirs"]["SNPs"] + "{genome_build}/dbsnp.common_all-151.bed"
    resources: **CFG["resources"]["SNPs"]
    shell:
        op.as_one_line(""" gunzip -c {input.vcf} | awk {{'printf ("%s\\t%s\\t%s\\n", $1,$2-1,$2)'}} | zgrep -v -h "^#" > {output.bed} """)

#### set-up mpileups for BAF calling ####
def _cnvkit_get_chr_mpileups(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["SNPs"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}.{chrom}.vcf.gz", 
        chrom = chrs
    )
    return(mpileups)
    
def _cnvkit_get_chr_mpileups_tbi(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["SNPs"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}.{chrom}.vcf.gz.tbi", 
        chrom = chrs
    )
    return(mpileups)

# vcf needs DP, GT, AD - bcftools -mv calls multiallelic variants (will annotate GT)
# without it, GT will not be annotated
rule _cnvkit_mpileup_per_chrom:
    input:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam",
        fastaFile = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bed = str(rules._cnvkit_dbsnp_to_bed.output.bed)
    output: # creates a temporary file for mpileup
        vcf = temp(CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.{chrom}.vcf.gz.tbi")
    params:
        quality = CFG["options"]["SNPs"]["quality"],
        opts = CFG["options"]["SNPs"]["opts"]
    conda:
        CFG["conda_envs"]["bcftools"]
    resources: 
        **CFG["resources"]["SNPs"]
    group: "cnvkit"
    log:
        stderr = CFG["logs"]["SNPs"] + "{capture_space}/{seq_type}--{genome_build}/{sample_id}/{chrom}.vcf.stderr.log",
    shell:
        """
            bcftools mpileup -T {input.bed} -r {wildcards.chrom} -f {input.fastaFile} -Q {params.quality} {params.opts} -Ou {input.bam} | bcftools call -mv -Oz -o {output.vcf} 2> {log.stderr} && 
            tabix -p vcf {output.vcf}
        """
        

rule _cnvkit_concatenate_vcf:
    input: 
        vcf = _cnvkit_get_chr_mpileups,
        tbi = _cnvkit_get_chr_mpileups_tbi,
    output: 
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.vcf.gz.tbi"
    resources: 
        **CFG["resources"]["SNPs"]
    group: "cnvkit"
    conda:
        CFG["conda_envs"]["bcftools"]
    shell: 
        """
            bcftools concat {input.vcf} -Oz -o {output.vcf} && 
            tabix -p vcf {output.vcf}
        """


# ----------------------------------------------------------------------------------------------- #
# Integrating cnvkit with BAF
# ----------------------------------------------------------------------------------------------- #
# Adds extra columns - One-sample t-test of bin log2 ratios versus 0.0 and ci high and ci low to be able to use filtering by ci in the next step
# Note that the t-test is not used in filtration step, but the ci is
rule _cnvkit_segmetrics_ttest:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns = CFG["dirs"]["cns"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cns"
    output:
        cns = temp(CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.raw.cns")
    params:
        add_col = CFG["options"]["segmetrics"]["add_col"]
    log:
        stdout = CFG["logs"]["BAF"] + "segmetrics/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_segmetrics.log"
    conda: CFG["conda_envs"]["cnvkit"]
    resources: 
        **CFG["resources"]["call"]
    group: "call_cns"
    shell:
        """
            cnvkit.py segmetrics {input.cnr} -s {input.cns} -o {output.cns} {params.add_col} &> {log.stdout}
        """
    

rule _run_cnvkit_call_vcf:
    input:        
        cns = CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.raw.cns",
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.vcf.gz.tbi"
    output: 
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns"
    params:
        rescale = CFG["options"]["BAF"]["rescale"],
        min_depth = CFG["options"]["BAF"]["min_depth"],
        filter_by = CFG["options"]["BAF"]["filter_by"],
        male_ref = CFG["options"]["male_ref"],
        opts = CFG["options"]["BAF"]["opts"]
    log: 
        CFG["logs"]["BAF"] + "call/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_call.log"
    group: "cnvkit"
    conda: 
        CFG["conda_envs"]["cnvkit"]
    resources: 
        **CFG["resources"]["call"]
    shell:
        """
            cnvkit.py call {input.cns} --output {output.cns} -v {input.vcf} --min-variant-depth {params.min_depth} -m {params.rescale} --filter {params.filter_by} {params.male_ref} {params.opts} &> {log}
        """


# plot a scatter plot of amps and dels, also BAF
rule _run_cnvkit_scatter:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns", 
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.vcf.gz.tbi"
    output: 
        png = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_scatter.png"
    params:
        min_depth = CFG["options"]["scatter"]["min_depth"],
        ymax = CFG["options"]["scatter"]["ymax"],
        ymin = CFG["options"]["scatter"]["ymin"]
    conda: 
        CFG["conda_envs"]["cnvkit"]
    group: "cnvkit"
    resources:
        **CFG["resources"]["plots"]
    log:
        stdout = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_scatter.log"
    shell:
        "cnvkit.py scatter {input.cnr} -s {input.cns} --output {output.png} -v {input.vcf} --min-variant-depth {params.min_depth} --y-max {params.ymax} --y-min {params.ymin} &> {log.stdout} "


# plot chromosome diagrams highlighting these amps/dels and also key genes that are located in these CNVs
rule _run_cnvkit_diagram:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns", 
    output:  # only pdf works
        pdf = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_diagram.pdf"
    params:
        threshold = CFG["options"]["diagram"]["threshold"], # to only label genes in high level amps and dels
        male_ref = CFG["options"]["male_ref"]
    conda: 
        CFG["conda_envs"]["cnvkit"]
    group: "cnvkit"
    resources:
        **CFG["resources"]["plots"]
    log:
        stdout = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_diagram.log"
    shell:
        "cnvkit.py diagram {input.cnr} -s {input.cns} --output {output.pdf} -t {params.threshold} {params.male_ref} &> {log.stdout} "


# find potential breakpoints across the CNVs or regions with large CN signal inconsistencies
rule _cnvkit_breaks:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns", 
    output:
        breaks = CFG["dirs"]["breaks"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.genebreaks.txt"
    conda: 
        CFG["conda_envs"]["cnvkit"]
    resources:
        **CFG["resources"]["breaks"]
    log:
        stderr = CFG["logs"]["breaks"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_breaks.log"
    shell:
        """
            cnvkit.py breaks {input.cnr} {input.cns} > {output.breaks} 2> {log.stderr}
        """


# with segments (cns) as input, the minimum probes option defines the segment's bin count
# without cns as input, the gene's weighted bin counts are used instead
rule _cnvkit_genemetrics_seg:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns", 
    output:
        genemetrics = CFG["dirs"]["geneMetrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/segment.gene_cn.txt"
    params:
        threshold = CFG["options"]["geneMetrics"]["threshold"],
        min_segments = CFG["options"]["geneMetrics"]["min_segments"], # to remove false positives that cover a small number of bins
        male_ref = CFG["options"]["male_ref"]
    conda: CFG["conda_envs"]["cnvkit"]
    resources:
        **CFG["resources"]["geneMetrics"]
    shell:
        """
            cnvkit.py genemetrics {input.cnr} -s {input.cns} --threshold {params.threshold} --min-probes {params.min_segments} {params.male_ref} > {output.genemetrics}
        """

rule _cnvkit_genemetrics:
    input:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
    output:
        genemetrics = CFG["dirs"]["geneMetrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/bin.gene_cn.txt"
    params:
        threshold = CFG["options"]["geneMetrics"]["threshold"],
        min_segments = CFG["options"]["geneMetrics"]["min_segments"], # to remove false positives that cover a small number of bins
        male_ref = CFG["options"]["male_ref"]
    conda: CFG["conda_envs"]["cnvkit"]
    resources:
        **CFG["resources"]["geneMetrics"]
    shell:
        """
            cnvkit.py genemetrics {input.cnr} --threshold {params.threshold} --min-probes {params.min_segments} {params.male_ref} > {output.genemetrics}
        """


# can take the intersection of the two methods to filter for a list of genes that confidently have CN change
rule _cnvkit_trusted_genes_cna:
    input:
        segGene = CFG["dirs"]["geneMetrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/segment.gene_cn.txt",
        binGene = CFG["dirs"]["geneMetrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/bin.gene_cn.txt"
    output:
        trusted_genes = CFG["dirs"]["geneMetrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/trusted_genes.txt"
    shell:
        """
            comm -12 <(tail -n+2 {input.segGene} | cut -f1 | sort ) <(tail -n+2 {input.binGene} | cut -f1 | sort ) > {output.trusted_genes}
        """


rule _cnvkit_infer_sex:
    input:
        targetcov = CFG["dirs"]["coverage"] + "target/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.targetcoverage.cnn",
        antitargetcov = CFG["dirs"]["coverage"] + "antitarget/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.antitargetcoverage.cnn",
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns = CFG["dirs"]["cns"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cns",
        call =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns"
    output:
        sex = CFG["dirs"]["geneMetrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/inferred_sex.txt"
    params:
        male_ref = CFG["options"]["male_ref"]
    conda: CFG["conda_envs"]["cnvkit"]
    resources:
        **CFG["resources"]["geneMetrics"]
    shell:
        """
            cnvkit.py sex {input.targetcov} {input.antitargetcov} {input.cnr} {input.cns} {input.call} {params.male_ref} > {output.sex}
        """


rule _cnvkit_to_seg:
    input:
        cns = str(rules._run_cnvkit_call_vcf.output.cns),
    output:
        seg = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["seg"]
    resources:
        **CFG["resources"]["seg"]
    shell:
        """cnvkit.py export seg {input.cns} -o {output.seg} """


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cnvkit_output:
    input:
        call_cns = str(rules._run_cnvkit_call_vcf.output.cns),
        scatter = str(rules._run_cnvkit_scatter.output.png),
        diagram = str(rules._run_cnvkit_diagram.output.pdf),
        breaks = str(rules._cnvkit_breaks.output.breaks),
        gene_seg = str(rules._cnvkit_genemetrics_seg.output.genemetrics),
        geneList = str(rules._cnvkit_trusted_genes_cna.output.trusted_genes),
        sex = str(rules._cnvkit_infer_sex.output.sex),
        seg = str(rules._cnvkit_to_seg.output.seg)
    output:
        call_cns = CFG["dirs"]["outputs"] + "wildcard_capture/BAF_cns/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns",
        scatter = CFG["dirs"]["outputs"] + "wildcard_capture/scatter/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_scatter.png",
        diagram = CFG["dirs"]["outputs"] + "wildcard_capture/diagram/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_diagram.pdf",
        breaks = CFG["dirs"]["outputs"] + "wildcard_capture/breaks/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genebreaks.txt",
        gene_seg = CFG["dirs"]["outputs"] + "wildcard_capture/geneMetrics/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_geneSeg.txt",
        geneList = CFG["dirs"]["outputs"] + "wildcard_capture/geneList/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genebreaks.txt",
        sex = CFG["dirs"]["outputs"] + "wildcard_capture/sex/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genebreaks.txt",
        seg = CFG["dirs"]["outputs"] + "wildcard_capture/seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
    run:
        op.relative_symlink(input.call_cns, output.call_cns, in_module = True)
        op.relative_symlink(input.scatter, output.scatter, in_module = True)
        op.relative_symlink(input.diagram, output.diagram, in_module = True)
        op.relative_symlink(input.breaks, output.breaks, in_module = True)
        op.relative_symlink(input.gene_seg, output.gene_seg, in_module = True)
        op.relative_symlink(input.geneList, output.geneList, in_module = True)
        op.relative_symlink(input.sex, output.sex, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)


# this is to collapse the wildcard {capture_space} to streamline downstream analyses
rule _cnvkit_output_no_capture_space:
    input:
        call_cns = str(rules._run_cnvkit_call_vcf.output.cns),
        scatter = str(rules._run_cnvkit_scatter.output.png),
        diagram = str(rules._run_cnvkit_diagram.output.pdf),
        breaks = str(rules._cnvkit_breaks.output.breaks),
        gene_seg = str(rules._cnvkit_genemetrics_seg.output.genemetrics),
        geneList = str(rules._cnvkit_trusted_genes_cna.output.trusted_genes),
        sex = str(rules._cnvkit_infer_sex.output.sex),
        seg = str(rules._cnvkit_to_seg.output.seg)
    output:
        call_cns = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/BAF_cns/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns"),
        scatter = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/scatter/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_scatter.png"),
        diagram = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/diagram/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_diagram.pdf"),
        breaks = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/breaks/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genebreaks.txt"),
        gene_seg = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/geneMetrics/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_geneSeg.txt"),
        geneList = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/geneList/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genebreaks.txt"),
        sex = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/sex/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genebreaks.txt"),
        seg = touch(CFG["dirs"]["outputs"] + "no_wildcard_capture/seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg")
    params:
        call_cns = CFG["dirs"]["outputs"] + "no_wildcard_capture/BAF_cns/{seq_type}--{genome_build}/{tumour_id}.call.cns",
        scatter = CFG["dirs"]["outputs"] + "no_wildcard_capture/scatter/{seq_type}--{genome_build}/{tumour_id}_scatter.png",
        diagram = CFG["dirs"]["outputs"] + "no_wildcard_capture/diagram/{seq_type}--{genome_build}/{tumour_id}_diagram.pdf",
        breaks = CFG["dirs"]["outputs"] + "no_wildcard_capture/breaks/{seq_type}--{genome_build}/{tumour_id}_genebreaks.txt",
        gene_seg = CFG["dirs"]["outputs"] + "no_wildcard_capture/geneMetrics/{seq_type}--{genome_build}/{tumour_id}_geneSeg.txt",
        geneList = CFG["dirs"]["outputs"] + "no_wildcard_capture/geneList/{seq_type}--{genome_build}/{tumour_id}_genebreaks.txt",
        sex = CFG["dirs"]["outputs"] + "no_wildcard_capture/sex/{seq_type}--{genome_build}/{tumour_id}_genebreaks.txt",
        seg = CFG["dirs"]["outputs"] + "no_wildcard_capture/seg/{seq_type}--{genome_build}/{tumour_id}.seg"
    run:
        op.relative_symlink(input.call_cns, params.call_cns, in_module = True)
        op.relative_symlink(input.scatter, params.scatter, in_module = True)
        op.relative_symlink(input.diagram, params.diagram, in_module = True)
        op.relative_symlink(input.breaks, params.breaks, in_module = True)
        op.relative_symlink(input.gene_seg, params.gene_seg, in_module = True)
        op.relative_symlink(input.geneList, params.geneList, in_module = True)
        op.relative_symlink(input.sex, params.sex, in_module = True)
        op.relative_symlink(input.seg, params.seg, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _cnvkit_all:
    input:
        expand(
            [
                str(rules._cnvkit_output_no_capture_space.output.call_cns),
                str(rules._cnvkit_output_no_capture_space.output.scatter),
                str(rules._cnvkit_output_no_capture_space.output.diagram),
                str(rules._cnvkit_output_no_capture_space.output.breaks),
                str(rules._cnvkit_output_no_capture_space.output.gene_seg),
                str(rules._cnvkit_output_no_capture_space.output.geneList),                
                str(rules._cnvkit_output_no_capture_space.output.sex),
                str(rules._cnvkit_output_no_capture_space.output.seg),
                str(rules._cnvkit_output.output.call_cns),
                str(rules._cnvkit_output.output.scatter),
                str(rules._cnvkit_output.output.diagram),
                str(rules._cnvkit_output.output.breaks),
                str(rules._cnvkit_output.output.gene_seg),
                str(rules._cnvkit_output.output.geneList),                
                str(rules._cnvkit_output.output.sex),
                str(rules._cnvkit_output.output.seg),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            capture_space=CFG["runs"]["tumour_capture_space"],
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
