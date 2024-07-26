#!/usr/bin/env snakemake


# ---------------------------------------------------------------------------- #
##### ATTRIBUTION #####
# ---------------------------------------------------------------------------- #

# Original snakemake author: Jasper Wong
# Module author: Jasper Wong
# Additional contributors: N/A


# ---------------------------------------------------------------------------- #
##### SETUP #####
# ---------------------------------------------------------------------------- #

### Modules ###

import pandas as pd
import numpy as np
import oncopipe as op
import glob
import os

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

### Directories ###
# Setup module and store module-specific configuration in `CFG`.
CFG = op.setup_module(
    name = "ichorcna", 
    version = "1.0",
    subdirectories = ["inputs", "readDepth", "seg", "outputs"]
)

localrules:
    _ichorcna_input_bam,
    _ichorcna_output,
    _ichorcna_all

# ---------------------------------------------------------------------------- #
##### RULES #####
# ---------------------------------------------------------------------------- #

### Set-up dependencies and packages ###
# Download github and all external files for ichorCNA: (needed since their extdata is not complete for all genome builds)
rule _install_ichorcna:
    output:
        complete = CFG["dirs"]["inputs"] + "ichorcna_dependencies_installed.success"
    params:
        outdir = CFG["dirs"]["inputs"] + "ichorCNA/"
    conda:
        CFG["conda_envs"]["ichorcna"]
    shell:
        op.as_one_line("""
        git clone git://github.com/broadinstitute/ichorCNA.git {params.outdir} &&            
        touch {output.complete}""")

# This defines the script/extdata directory used by ichorCNA in the subsequent rules:
ichorDir = CFG["dirs"]["inputs"] + "ichorCNA/inst/extdata/" 

# Symlinks the extdata appropriately
rule _setup_ichorcna_extdata:
    input:
        complete = CFG["dirs"]["inputs"] + "ichorcna_dependencies_installed.success"
    params:
        hg19_1Mb_rds = ichorDir + "HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds",
        hg19_500kb_rds = ichorDir + "HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds",
        hg38_1Mb_rds = ichorDir + "HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds",
        hg38_500kb_rds = ichorDir + "HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds",
        hg19_1000kb_gc = ichorDir + "gc_hg19_1000kb.wig",
        hg19_500kb_gc = ichorDir + "gc_hg19_500kb.wig",
        hg19_50kb_gc = ichorDir + "gc_hg19_50kb.wig",
        hg19_10kb_gc = ichorDir + "gc_hg19_10kb.wig",
        hg38_1000kb_gc = ichorDir + "gc_hg38_1000kb.wig",
        hg38_500kb_gc = ichorDir + "gc_hg38_500kb.wig",
        hg38_50kb_gc = ichorDir + "gc_hg38_50kb.wig",
        hg38_10kb_gc = ichorDir + "gc_hg38_10kb.wig",
        hg19_1000kb_map = ichorDir + "map_hg19_1000kb.wig",
        hg19_500kb_map = ichorDir + "map_hg19_500kb.wig",
        hg19_50kb_map = ichorDir + "map_hg19_50kb.wig",
        hg19_10kb_map = ichorDir + "map_hg19_10kb.wig",
        hg38_1000kb_map = ichorDir + "map_hg38_1000kb.wig",
        hg38_500kb_map = ichorDir + "map_hg38_500kb.wig",
        hg38_50kb_map = ichorDir + "map_hg38_50kb.wig",
        hg38_10kb_map = ichorDir + "map_hg38_10kb.wig",
    output:
        hg19_1Mb_rds = ichorDir + "HD_ULP_PoN_hg19_1Mb_median_normAutosome_median.rds",
        hg19_500kb_rds = ichorDir + "HD_ULP_PoN_hg19_500kb_median_normAutosome_median.rds",
        grch37_1Mb_rds = ichorDir + "HD_ULP_PoN_grch37_1Mb_median_normAutosome_median.rds",
        grch37_500kb_rds = ichorDir + "HD_ULP_PoN_grch37_500kb_normAutosome_median.rds",
        hs37d5_1Mb_rds = ichorDir + "HD_ULP_PoN_hs37d5_1Mb_median_normAutosome_median.rds",
        hs37d5_500kb_rds = ichorDir + "HD_ULP_PoN_hs37d5_500kb_normAutosome_median.rds",
        grch38_1Mb_rds = ichorDir + "HD_ULP_PoN_grch38_1Mb_median_normAutosome_median.rds",
        grch38_500kb_rds = ichorDir + "HD_ULP_PoN_grch38_500kb_median_normAutosome_median.rds",
        grch37_1000kb_gc = ichorDir + "gc_grch37_1000kb.wig",
        grch37_500kb_gc = ichorDir + "gc_grch37_500kb.wig",
        grch37_50kb_gc = ichorDir + "gc_grch37_50kb.wig",
        grch37_10kb_gc = ichorDir + "gc_grch37_10kb.wig",
        hs37d5_1000kb_gc = ichorDir + "gc_hs37d5_1000kb.wig",
        hs37d5_500kb_gc = ichorDir + "gc_hs37d5_500kb.wig",
        hs37d5_50kb_gc = ichorDir + "gc_hs37d5_50kb.wig",
        hs37d5_10kb_gc = ichorDir + "gc_hs37d5_10kb.wig",
        grch38_1000kb_gc = ichorDir + "gc_grch38_1000kb.wig",
        grch38_500kb_gc = ichorDir + "gc_grch38_500kb.wig",
        grch38_50kb_gc = ichorDir + "gc_grch38_50kb.wig",
        grch38_10kb_gc = ichorDir + "gc_grch38_10kb.wig",
        grch37_1000kb_map = ichorDir + "map_grch37_1000kb.wig",
        grch37_500kb_map = ichorDir + "map_grch37_500kb.wig",
        grch37_50kb_map = ichorDir + "map_grch37_50kb.wig",
        grch37_10kb_map = ichorDir + "map_grch37_10kb.wig",
        hs37d5_1000kb_map = ichorDir + "map_hs37d5_1000kb.wig",
        hs37d5_500kb_map = ichorDir + "map_hs37d5_500kb.wig",
        hs37d5_50kb_map = ichorDir + "map_hs37d5_50kb.wig",
        hs37d5_10kb_map = ichorDir + "map_hs37d5_10kb.wig",
        grch38_1000kb_map = ichorDir + "map_grch38_1000kb.wig",
        grch38_500kb_map = ichorDir + "map_grch38_500kb.wig",
        grch38_50kb_map = ichorDir + "map_grch38_50kb.wig",
        grch38_10kb_map = ichorDir + "map_grch38_10kb.wig",
        complete = touch(ichorDir + "symlink.done")
    run:
        op.relative_symlink(params.hg19_1Mb_rds, output.hg19_1Mb_rds)
        op.relative_symlink(params.hg19_500kb_rds, output.hg19_500kb_rds)
        op.relative_symlink(params.hg19_1Mb_rds, output.grch37_1Mb_rds)
        op.relative_symlink(params.hg19_1Mb_rds, output.hs37d5_1Mb_rds)
        op.relative_symlink(params.hg19_500kb_rds, output.grch37_500kb_rds)
        op.relative_symlink(params.hg19_500kb_rds, output.hs37d5_500kb_rds)
        op.relative_symlink(params.hg38_1Mb_rds, output.grch38_1Mb_rds)
        op.relative_symlink(params.hg38_500kb_rds, output.grch38_500kb_rds)
        op.relative_symlink(params.hg19_1000kb_gc, output.grch37_1000kb_gc)
        op.relative_symlink(params.hg19_500kb_gc, output.grch37_500kb_gc)
        op.relative_symlink(params.hg19_50kb_gc, output.grch37_50kb_gc)
        op.relative_symlink(params.hg19_10kb_gc, output.grch37_10kb_gc)
        op.relative_symlink(params.hg19_1000kb_gc, output.hs37d5_1000kb_gc)
        op.relative_symlink(params.hg19_500kb_gc, output.hs37d5_500kb_gc)
        op.relative_symlink(params.hg19_50kb_gc, output.hs37d5_50kb_gc)
        op.relative_symlink(params.hg19_10kb_gc, output.hs37d5_10kb_gc)
        op.relative_symlink(params.hg38_1000kb_gc, output.grch38_1000kb_gc)
        op.relative_symlink(params.hg38_500kb_gc, output.grch38_500kb_gc)
        op.relative_symlink(params.hg38_50kb_gc, output.grch38_50kb_gc)
        op.relative_symlink(params.hg38_10kb_gc, output.grch38_10kb_gc)
        op.relative_symlink(params.hg19_1000kb_map, output.grch37_1000kb_map)
        op.relative_symlink(params.hg19_500kb_map, output.grch37_500kb_map)
        op.relative_symlink(params.hg19_50kb_map, output.grch37_50kb_map)
        op.relative_symlink(params.hg19_10kb_map, output.grch37_10kb_map)
        op.relative_symlink(params.hg19_1000kb_map, output.hs37d5_1000kb_map)
        op.relative_symlink(params.hg19_500kb_map, output.hs37d5_500kb_map)
        op.relative_symlink(params.hg19_50kb_map, output.hs37d5_50kb_map)
        op.relative_symlink(params.hg19_10kb_map, output.hs37d5_10kb_map)
        op.relative_symlink(params.hg38_1000kb_map, output.grch38_1000kb_map)
        op.relative_symlink(params.hg38_500kb_map, output.grch38_500kb_map)
        op.relative_symlink(params.hg38_50kb_map, output.grch38_50kb_map)
        op.relative_symlink(params.hg38_10kb_map, output.grch38_10kb_map)

### Run ichorCNA ###
# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ichorcna_input_bam:
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
     

# set-up for CRAM files (readCounter does not work with CRAM)
# deeptools to get .bw from .bam and .cram
rule _ichorcna_bamCoverage:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai",
        ichorcna_package = CFG["dirs"]["inputs"] + "ichorcna_dependencies_installed.success",
        symlink_complete = ichorDir + "symlink.done"
    output:
        bw = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/bw/{sample_id}.bin{binSize}.bw"
    params:
        binSize = CFG["options"]["deeptools"]["binSize"],
        qual = CFG["options"]["deeptools"]["qual"],
        excludeFlag = CFG["options"]["deeptools"]["flagExclude"],
        opt = CFG["options"]["deeptools"]["opt"],
        dirOut = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/bw/"
    conda: CFG["conda_envs"]["deeptools"]
    threads: CFG["threads"]["deeptools"]
    resources:
        **CFG["resources"]["deeptools"]
    log:
        CFG["logs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/bw/{sample_id}.bin{binSize}.log"
    shell:
        """
            mkdir -p {params.dirOut}; 
            bamCoverage -b {input.bam} --binSize {params.binSize} --minMappingQuality {params.qual} --samFlagExclude {params.excludeFlag} {params.opt} -o {output.bw} -p {threads}  
        """


# Converts bigWig to Wig
rule _ichorcna_bigwigToWig:
    input:
        bw = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/bw/{sample_id}.bin{binSize}.bw"
    output:
        wig_int = temp(CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{sample_id}.bin{binSize}{chrom}.wig"),
    conda: CFG["conda_envs"]["ucsc-bigwigtowig"]
    threads: CFG["threads"]["ucsc"]
    resources:
        **CFG["resources"]["ucsc"]
    wildcard_constraints:
        chrom = ".+(?<!--fixed)"
    log:
        CFG["logs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{sample_id}.bin{binSize}.{chrom}.log"
    shell:
        """
            bigWigToWig {input.bw} {output.wig_int} -chrom={wildcards.chrom} 
        """


# This function will reformat the wig file to one that can be used for ichorCNA
rule _ichorcna_spread_centromeres:
    input:
        wig_int = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{sample_id}.bin{binSize}{chrom}.wig",
    output:
        wig = temp(CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{sample_id}.bin{binSize}.{chrom}--fixed.wig")
    conda: CFG["conda_envs"]["bedops_tools"]
    threads: CFG["threads"]["ucsc"]
    resources:
        **CFG["resources"]["ucsc"]
    wildcard_constraints:
        chrom = ".+(?<!--fixed)"
    log:
        CFG["logs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{sample_id}.bin{binSize}.{chrom}.bedops.log"
    shell:
        """
            echo -e "fixedStep chrom={wildcards.chrom} start=1 step={wildcards.binSize} span={wildcards.binSize} " > {output.wig} &&
            intersectBed -a <( bedops --chop {wildcards.binSize} --header  {input.wig_int} ) -b {input.wig_int} -wa -wb | awk '{{print $7}}' >> {output.wig} 2>> {log}
        """


# This function is used to get the wigs of the main chromosomes, which will be stitched together
def get_chrom_wigs(wildcards):
    CFG = config["lcr-modules"]["ichorcna"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    wig = expand(
        CFG["dirs"]["readDepth"] + "{{seq_type}}--{{genome_build}}/{{binSize}}/wig/{{sample_id}}.bin{{binSize}}.{chrom}--fixed.wig", 
        chrom = chrs
    )
    return(wig)


rule _ichorcna_wigCompile:
    input:
        get_chrom_wigs
    output:
        wig = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{sample_id}.bin{binSize}.wig"
    shell:
        """
            cat {input} > {output.wig}
        """



# This function will return a comma-separated list of chromosomes to include in runIchorCNA
def get_chromosomes_R(wildcards):
    chromosomesR=[]
    stringStart="c('"
    for i in range(1,23):
        chromosomesR.append(str(i))
    chromosomesR.append("X")
    if "38" in str(wildcards.genome_build):
        chromosomesR = ["chr" + x for x in chromosomesR]
    chromosomesR= "','".join(chromosomesR)
    stringEnd="')"
    return stringStart + chromosomesR + stringEnd

# New functions just for specific alt-genome builds    
def _which_gcwig(wildcards):
    CFG = config['lcr-modules']['ichorcna']
    this_genome_build = str(wildcards.genome_build)
    try:
        wigs = CFG["options"]["run"]["ichorCNA_gcWig_custom"][this_genome_build]
    except NameError:
        wigs = None
    try:
        wigs = CFG["options"]["run"]["ichorCNA_gcWig_custom"][wildcards.genome_build]
    except NameError:
        wigs = None
    if wigs is not None and wigs != "":
        return wigs
    elif "38" in str({wildcards.genome_build}):
        if "1000000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg38_1000kb.wig"
        elif "500000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg38_500kb.wig"
        elif "50000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg38_50kb.wig"
        elif "10000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg38_10kb.wig"
        else:
            wig = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_gcWig"])
            return wig
    else:
        if "1000000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg19_1000kb.wig"
        elif "500000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg19_500kb.wig"
        elif "50000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg19_50kb.wig"
        elif "10000" in str({wildcards.binSize}):
            return "inst/extdata/gc_hg19_10kb.wig"
        else:
            wig = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_gcWig"])
            return wig

def _which_mapwig(wildcards):
    CFG = config['lcr-modules']['ichorcna']
    this_genome_build = str(wildcards.genome_build)
    try:
        wigs = CFG["options"]["run"]["ichorCNA_mapWig_custom"][this_genome_build]
    except NameError:
        wigs = None
    try:
        wigs = CFG["options"]["run"]["ichorCNA_mapWig_custom"][wildcards.genome_build]
    except NameError:
        wigs = None
    if wigs is not None and wigs != "":
        return wigs
    elif "38" in str({wildcards.genome_build}):
        if "1000000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg38_1000kb.wig"
        elif "500000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg38_500kb.wig"
        elif "50000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg38_50kb.wig"
        elif "10000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg38_10kb.wig"
        else:
            wig = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_mapWig"])
            return wig
    else:
        if "1000000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg19_1000kb.wig"
        elif "500000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg19_500kb.wig"
        elif "50000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg19_50kb.wig"
        elif "10000" in str({wildcards.binSize}):
            return "inst/extdata/map_hg19_10kb.wig"
        else:
            wig = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_mapWig"])
            return wig

def _which_normPanel(wildcards):
    CFG = config['lcr-modules']['ichorcna']
    this_genome_build = str(wildcards.genome_build)
    try:
        rds = CFG["options"]["run"]["ichorCNA_normalPanel_custom"][this_genome_build]
    except NameError:
        rds = None
    try:
        rds = CFG["options"]["run"]["ichorCNA_normalPanel_custom"][wildcards.genome_build]
    except NameError:
        rds = None
    if rds is not None and rds != "":
        return rds
    elif "38" in str({wildcards.genome_build}):
        if "1000000" in str({wildcards.binSize}):
            return "inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds"
        elif "500000" in str({wildcards.binSize}):
            return "inst/extdata/HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds"
        else:
            rds = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_normalPanel"])
            return rds
    else:
        if "1000000" in str({wildcards.binSize}):
            return "inst/extdata/HD_ULP_PoN_hg19_1Mb_median_normAutosome_median.rds"
        elif "500000" in str({wildcards.binSize}):
            return "inst/extdata/HD_ULP_PoN_hg19_1Mb_median_normAutosome_median.rds"
        else:
            rds = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_normalPanel"])
            return rds

rule _run_ichorcna:
    input:
        tum = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/wig/{tumour_id}.bin{binSize}.wig",
    output:
        corrDepth = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.correctedDepth.txt",
        param = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.params.txt",
        cna = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cna.seg",
        segTxt = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.seg.txt",
        seg = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.seg",
        plot = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}/{tumour_id}_genomeWide.pdf",
    params:
        ichorDir = CFG["dirs"]["inputs"] + "ichorCNA/",
        outDir = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/",
        rscript = CFG["options"]["run"]["ichorCNA_rscript"],
        name = "{tumour_id}",
        ploidy = CFG["options"]["run"]["ichorCNA_ploidy"],
        normal = CFG["options"]["run"]["ichorCNA_normal"],
        gcwig = _which_gcwig,
        mapwig = _which_mapwig,
        normalpanel = _which_normPanel,
        estimateNormal = CFG["options"]["run"]["ichorCNA_estimateNormal"],
        estimatePloidy = CFG["options"]["run"]["ichorCNA_estimatePloidy"],
        estimateClonality = CFG["options"]["run"]["ichorCNA_estimateClonality"],
        scStates = CFG["options"]["run"]["ichorCNA_scStates"],
        maxCN = CFG["options"]["run"]["ichorCNA_maxCN"],
        includeHOMD = CFG["options"]["run"]["ichorCNA_includeHOMD"],
        chrs = get_chromosomes_R,
        chrTrain = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_chrTrain"]),
        genomeBuild = "{genome_build}",
        genomeStyle = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_genomeStyle"]),
        centromere = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_centromere"]),
        fracReadsChrYMale = CFG["options"]["run"]["ichorCNA_fracReadsInChrYForMale"],
        minMapScore = CFG["options"]["run"]["ichorCNA_minMapScore"],
        maxFracGenomeSubclone = CFG["options"]["run"]["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone = CFG["options"]["run"]["ichorCNA_maxFracCNASubclone"],
        exons = CFG["options"]["run"]["ichorCNA_exons"],
        txnE = CFG["options"]["run"]["ichorCNA_txnE"],
        txnStrength = CFG["options"]["run"]["ichorCNA_txnStrength"],
        plotFileType = CFG["options"]["run"]["ichorCNA_plotFileType"],
        plotYlim = CFG["options"]["run"]["ichorCNA_plotYlim"],
        libdir = CFG["dirs"]["inputs"] + "ichorCNA/" + CFG["options"]["run"]["ichorCNA_libdir"]
    conda: CFG["conda_envs"]["ichorcna"]
    threads: CFG["threads"]["run"]
    resources:
        **CFG["resources"]["run"]
    log:
        stdout = CFG["logs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}.stdout.log",
        stderr = CFG["logs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}.stderr.log"
    shell:
         op.as_one_line("""
            Rscript --vanilla {params.rscript} 
            --id {params.name} 
            --libdir {params.libdir} 
            --WIG {input.tum} 
            --gcWig {params.ichorDir}{params.gcwig} 
            --mapWig {params.ichorDir}{params.mapwig} 
            --normalPanel {params.ichorDir}{params.normalpanel} 
            --ploidy \"{params.ploidy}\" 
            --normal \"{params.normal}\" 
            --maxCN {params.maxCN} 
            --includeHOMD {params.includeHOMD} 
            --chrs \"{params.chrs}\" 
            --chrTrain \"{params.chrTrain}\" 
            --genomeStyle {params.genomeStyle} 
            --genomeBuild {params.genomeBuild} 
            --estimateNormal {params.estimateNormal} 
            --estimatePloidy {params.estimatePloidy} 
            --estimateScPrevalence {params.estimateClonality} 
            --scStates \"{params.scStates}\" 
            --centromere {params.ichorDir}{params.centromere} 
            --exons.bed {params.exons} 
            --txnE {params.txnE} 
            --txnStrength {params.txnStrength} 
            --minMapScore {params.minMapScore} 
            --fracReadsInChrYForMale {params.fracReadsChrYMale} 
            --maxFracGenomeSubclone {params.maxFracGenomeSubclone} 
            --maxFracCNASubclone {params.maxFracCNASubclone} 
            --plotFileType {params.plotFileType} 
            --plotYLim \"{params.plotYlim}\" 
            --outDir {params.outDir} > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ichorcna_output:
    input:
        corrDepth = str(rules._run_ichorcna.output.corrDepth),
        param = str(rules._run_ichorcna.output.param),
        cna = str(rules._run_ichorcna.output.cna),
        segTxt = str(rules._run_ichorcna.output.segTxt),
        seg = str(rules._run_ichorcna.output.seg),
        plot = str(rules._run_ichorcna.output.plot)
    output:
        corrDepth = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/corrDepth/{binSize}/{tumour_id}--{normal_id}--{pair_status}.corrDepth.txt",
        param = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/param/{binSize}/{tumour_id}--{normal_id}--{pair_status}.param.txt",
        cna = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/binCNA/{binSize}/{tumour_id}--{normal_id}--{pair_status}.cna.seg",
        segTxt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg_txt/{binSize}/{tumour_id}--{normal_id}--{pair_status}.seg.txt",
        seg = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg/{binSize}/{tumour_id}--{normal_id}--{pair_status}.seg",
        plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plot/{binSize}/{tumour_id}--{normal_id}--{pair_status}_genomeWide.pdf"
    run:
        op.relative_symlink(input.corrDepth, output.corrDepth, in_module=True)
        op.relative_symlink(input.param, output.param, in_module=True)
        op.relative_symlink(input.cna, output.cna, in_module=True)
        op.relative_symlink(input.segTxt, output.segTxt, in_module=True)
        op.relative_symlink(input.seg, output.seg, in_module=True)
        op.relative_symlink(input.plot, output.plot, in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _ichorcna_all:
    input:
        expand(
            [
                str(rules._ichorcna_output.output.corrDepth),
                str(rules._ichorcna_output.output.param),
                str(rules._ichorcna_output.output.cna),
                str(rules._ichorcna_output.output.segTxt),
                str(rules._ichorcna_output.output.seg),
                str(rules._ichorcna_output.output.plot)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            pair_status=CFG["runs"]["pair_status"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            binSize=[CFG["options"]["deeptools"]["binSize"]] * len(CFG["runs"]["tumour_sample_id"]))



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
