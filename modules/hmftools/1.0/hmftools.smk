#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
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
# `CFG` is a shortcut to `config["lcr-modules"]["hmftools"]`
CFG = op.setup_module(
    name = "hmftools",
    version = "1.0",
    subdirectories = ["inputs", "prepare_strelka", "amber", "cobalt", "purple", "linx", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _hmftools_input_bam,
    _hmftools_input_strelka,
    _hmftools_strelka_sample_names, 
    _hmftools_input_gridss, 
    _hmftools_input_references,
    _hmftools_get_cobalt_gc,
    _hmftools_get_amber_snps, 
    _hmftools_get_purple_drivers, 
    _hmftools_get_linx_db,
    _hmftools_linx_prepare_ensembl,
    _hmftools_purple_output,
    _hmftools_purple_plots,  
    _hmftools_all


VERSION_MAP = {
    "grch37": "hg19",
    "hs37d5": "hg19",
    "hg38": "hg38"
}

possible_genome_builds = VERSION_MAP.keys()
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        "Samples table includes genome builds not yet compatible with this module. "
        "This module is currently only compatible with {possible_genome_builds}. "
    )

wildcard_constraints: 
    genome_build = "|".join(VERSION_MAP.keys()), 
    pair_status = "matched|unmatched"


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _hmftools_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"], 
        bai = CFG["inputs"]["sample_bai"], 
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai", 
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

rule _hmftools_input_strelka: 
    input: 
        strelka_vcf = CFG["inputs"]["strelka_vcf"], 
    output: 
        strelka_vcf = CFG["dirs"]["inputs"] + "strelka_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic.combined.vcf.gz" 
    run: 
        op.absolute_symlink(input.strelka_vcf, output.strelka_vcf)

rule _hmftools_input_gridss: 
    input: 
        gridss_somatic_vcf = CFG["inputs"]["gridss_somatic"], 
        gridss_somatic_tbi = CFG["inputs"]["gridss_somatic_tbi"], 
        gridss_filtered_vcf = CFG["inputs"]["gridss_somatic_filtered"], 
        gridss_filtered_tbi = CFG["inputs"]["gridss_somatic_filtered_tbi"]
    output: 
        gridss_somatic_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz",
        gridss_somatic_tbi = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz.tbi",
        gridss_filtered_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz", 
        gridss_filtered_tbi = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz.tbi"
    run: 
        op.absolute_symlink(input.gridss_somatic_vcf, output.gridss_somatic_vcf)
        op.absolute_symlink(input.gridss_somatic_tbi, output.gridss_somatic_tbi)
        op.absolute_symlink(input.gridss_filtered_vcf, output.gridss_filtered_vcf)
        op.absolute_symlink(input.gridss_filtered_tbi, output.gridss_filtered_tbi)

# Rules to download and setup reference files

rule _hmftools_input_references: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        genome_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa", 
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa.fai", 
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.dict"
    shell: 
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fai} {output.genome_fai} &&
        ln -s {input.genome_dict} {output.genome_dict}
        """)

rule _hmftools_get_cobalt_gc: 
    output: 
        gc = CFG["dirs"]["inputs"] + "references/{genome_build}/cobalt/GC_profile.1000bp.cnp"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/cobalt",
        alt_build = lambda w: VERSION_MAP[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.gc} {params.url}/GC_profile.{params.alt_build}.1000bp.cnp'

rule _hmftools_get_amber_snps: 
    output: 
        vcf = CFG["dirs"]["inputs"] + "references/{genome_build}/amber/GermlineHetPon.vcf.gz", 
        snpcheck = CFG["dirs"]["inputs"] + "references/{genome_build}/amber/GermlineHetPon.snpcheck.vcf.gz"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/amber",
        alt_build = lambda w: VERSION_MAP[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.vcf} {params.url}/GermlineHetPon.{params.alt_build}.vcf.gz; '
        'wget -O {output.snpcheck} {params.url}/GermlineHetPon.{params.alt_build}.snpcheck.vcf.gz'

rule _hmftools_get_purple_drivers: 
    output: 
        hotspots = CFG["dirs"]["inputs"] + "references/{genome_build}/purple/KnownHotspots.vcf.gz", 
        gene_panel = CFG["dirs"]["inputs"] + "references/{genome_build}/purple/DriverGenePanel.tsv"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/purple",
        alt_build = lambda w: VERSION_MAP[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.hotspots} {params.url}/KnownHotspots.{params.alt_build}.vcf.gz && '
        'wget -O {output.gene_panel} {params.url}/DriverGenePanel.{params.alt_build}.tsv'

rule _hmftools_get_linx_db: 
    output: 
        directory(CFG["dirs"]["inputs"] + "references/linx_db/")
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/linx/"
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -r -np -nd -P {output} -A .bed,.csv {params.url} '

# Prepare Strelka VCF files for use with PURPLE
# SnpEff annotation enables driver discovery logic

rule _hmftools_strelka_sample_names: 
    input: 
        vcf = rules._hmftools_input_strelka.output.strelka_vcf
    output: 
        vcf = temp(CFG["dirs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/tmp.strelka.vcf")
    log: CFG["dirs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_sample_names.log"
    conda: 
        CFG["conda_envs"]["bcftools"]
    threads: CFG["threads"]["strelka_sample_names"]
    resources: 
        **CFG["resources"]["strelka_sample_names"]
    shell: 
        op.as_one_line("""
        bcftools view -Ov {input.vcf} | 
        sed 's/TUMOR/{wildcards.tumour_id}/g' | 
        sed 's/NORMAL/{wildcards.normal_id}/g'  
        > {output.vcf}
        """)

rule _hmftools_snpeff_strelka: 
    input: 
        vcf = str(rules._hmftools_strelka_sample_names.output.vcf)
    output: 
        sample_key = temp(CFG["dirs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sample_key.txt"),
        vcf = temp(CFG["dirs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.snpeff.vcf.gz")
    resources: 
        **CFG["resources"]["snpeff"]
    params: 
        snpeff_build = lambda w: {
            "grch37": "GRCh37.75", 
            "hs37d5": "GRCh37.75", 
            "hg38": "hg38"
        }[w.genome_build], 
        config = "$(readlink -e $(which snpEff)).config",
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    log: 
        CFG["logs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/snpeff_strelka.log"
    conda: 
        CFG["conda_envs"]["snpeff"]
    threads: 
        CFG["threads"]["snpeff"]
    shell: 
        op.as_one_line("""
        printf "{wildcards.normal_id}\t{wildcards.tumour_id}\n" > {output.sample_key} && 
        snpEff -Xmx{params.mem_mb}m   
        -c {params.config} -noStats
        -cancer -cancerSamples {output.sample_key} 
        {params.snpeff_build} {input.vcf} | 
        bcftools view -Oz -o {output.vcf} - && 
        bcftools index -t {output.vcf}
        """)


rule _hmftools_annotate_strelka: 
    input: 
        vcf = str(rules._hmftools_snpeff_strelka.output.vcf)
    output: 
        vcf = temp(CFG["dirs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.snpeff.annotated.vcf.gz")
    params: 
        purple_jar = "$(dirname $(readlink -e $(which PURPLE)))/purple.jar"
    log: CFG["logs"]["prepare_strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/annotate_strelka.log"
    conda: CFG["conda_envs"]["purple"]
    threads: 
        CFG["threads"]["annotate_strelka"]
    resources: 
        **CFG["resources"]["annotate_strelka"]
    shell: 
        op.as_one_line("""
        java -Xmx{resources.mem_mb}m 
            -cp {params.purple_jar} com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth
            -in {input.vcf} -out {output.vcf}
            2>&1 | tee -a {log}
        """)

# Run AMBER to calculate BAFs
rule _hmftools_amber_matched: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        snps = str(rules._hmftools_get_amber_snps.output.vcf)
    output: 
        vcf = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status, matched}/{tumour_id}.amber.baf.vcf.gz"
    resources: 
        **CFG["resources"]["amber"]
    params:
        options = CFG["options"]["amber"], 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8) 
    log: CFG["logs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/amber.log"
    conda: 
        CFG["conda_envs"]["amber"]
    threads: 
        CFG["threads"]["amber"]
    shell: 
        op.as_one_line("""
        AMBER -Xmx{params.jvmheap}m
        -reference {wildcards.normal_id} -reference_bam {input.normal_bam}
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam}
        -output_dir `dirname {output.vcf}`
        -threads {threads}
        -loci {input.snps}
        {params.options}
        2>&1 | tee -a {log} 
        """)

rule _hmftools_amber_unmatched: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        snps = str(rules._hmftools_get_amber_snps.output.vcf)
    output: 
        vcf = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status, unmatched}/{tumour_id}.amber.baf.vcf.gz"
    resources: 
        **CFG["resources"]["amber"]
    params:
        options = CFG["options"]["amber"], 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    log: CFG["logs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/amber.log"
    conda: 
        CFG["conda_envs"]["amber"]
    threads: 
        CFG["threads"]["amber"]
    shell: 
        op.as_one_line("""
        AMBER -Xmx{params.jvmheap}m
        -tumor_only 
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam} 
        -output_dir `dirname {output.vcf}`
        -threads {threads}
        -loci {input.snps}
        {params.options}
        2>&1 | tee -a {log} 
        """)

# Run COBALT to estimate depth across the genome
rule _hmftools_cobalt: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        gc_profile = str(rules._hmftools_get_cobalt_gc.output.gc)
    output: 
        ratio = CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cobalt.ratio.tsv"
    log: ratio = CFG["logs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cobalt.log"
    resources: 
        **CFG["resources"]["cobalt"]
    params:
        options = CFG["options"]["cobalt"], 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: 
        CFG["conda_envs"]["cobalt"]
    threads: 
        CFG["threads"]["cobalt"]
    shell: 
        op.as_one_line("""
        COBALT -Xmx{params.jvmheap}m
        -reference {wildcards.normal_id} -reference_bam {input.normal_bam} 
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam} 
        -output_dir `dirname {output.ratio}` 
        -threads {threads} 
        -gc_profile {input.gc_profile} 
        {params.options} 
        2>&1 | tee -a {log} 
        """)

# Run PURPLE for final CNV calling 

# Define variables for output file names
purple_out = [
    "purity.tsv", 
    "purity.range.tsv", 
    "cnv.gene.tsv", 
    "sv.vcf.gz", 
]
purple_plots = [
        "circos",
        "input",
        "map",
        "purity.range",
        "segment"
    ]

rule _hmftools_purple_matched:
    input: 
        amber = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.amber.baf.vcf.gz",
        cobalt = CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cobalt.ratio.tsv",
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        strelka_vcf = rules._hmftools_annotate_strelka.output.vcf, 
        gridss_somatic_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz",
        gridss_filtered_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz",
        reference_fa = str(rules._hmftools_input_references.output.genome_fa), 
        gene_panel = str(rules._hmftools_get_purple_drivers.output.gene_panel), 
        hotspots = str(rules._hmftools_get_purple_drivers.output.hotspots), 
        gc_profile = str(rules._hmftools_get_cobalt_gc.output.gc)
    output: 
        files = expand(CFG["dirs"]["purple"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{{tumour_id}}.purple.{out_file}", 
            out_file = purple_out), 
        plots = expand(CFG["dirs"]["purple"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/plot/{{tumour_id}}.{plot_name}.png", 
            plot_name = purple_plots)
    log: CFG["logs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/purple.log"
    resources: 
        **CFG["resources"]["purple"]
    params: 
        outdir = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}", 
        options = CFG["options"]["purple"],
        circos = "`which circos`", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    wildcard_constraints: 
        pair_status = "matched"
    conda: 
        CFG["conda_envs"]["purple"]
    threads: 
        CFG["threads"]["purple"]
    shell: 
        op.as_one_line("""
        PURPLE -Xmx{params.jvmheap}m -driver_catalog 
            -reference {wildcards.normal_id} 
            -tumor {wildcards.tumour_id} 
            -output_dir {params.outdir} 
            -amber `dirname {input.amber}` 
            -cobalt `dirname {input.cobalt}` 
            -gc_profile {input.gc_profile} 
            -ref_genome {input.reference_fa}
            -hotspots {input.hotspots} 
            -driver_gene_panel {input.gene_panel}  
            -somatic_vcf {input.strelka_vcf} 
            -structural_vcf {input.gridss_filtered_vcf} 
            -sv_recovery_vcf {input.gridss_somatic_vcf} 
            -circos {params.circos} 
            {params.options}
            2>&1 | tee -a {log}
        """)

rule _hmftools_purple_unmatched:
    input: 
        amber = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.amber.baf.vcf.gz",
        cobalt = CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cobalt.ratio.tsv",
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        gridss_somatic_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz",
        gridss_filtered_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz",
        reference_fa = str(rules._hmftools_input_references.output.genome_fa), 
        gc_profile = str(rules._hmftools_get_cobalt_gc.output.gc)
    output: 
        files = expand(CFG["dirs"]["purple"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{{tumour_id}}.purple.{out_file}", 
            out_file = purple_out), 
        plots = expand(CFG["dirs"]["purple"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/plot/{{tumour_id}}.{plot_name}.png", 
            plot_name = purple_plots)
    log: CFG["logs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/purple.log"
    resources: 
        **CFG["resources"]["purple"]
    params: 
        outdir = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}", 
        options = CFG["options"]["purple"],
        circos = "`which circos`", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    wildcard_constraints: 
        pair_status = "unmatched"
    conda: 
        CFG["conda_envs"]["purple"]
    threads: 
        CFG["threads"]["purple"]
    shell: 
        op.as_one_line("""
        PURPLE -Xmx{params.jvmheap}m 
            -reference {wildcards.normal_id} 
            -tumor {wildcards.tumour_id} 
            -output_dir {params.outdir} 
            -amber `dirname {input.amber}` 
            -cobalt `dirname {input.cobalt}` 
            -gc_profile {input.gc_profile} 
            -ref_genome {input.reference_fa} 
            -structural_vcf {input.gridss_filtered_vcf} 
            -sv_recovery_vcf {input.gridss_somatic_vcf} 
            -circos {params.circos} 
            -tumor_only  
            {params.options}
            2>&1 | tee -a {log}
        """)

# Prepare LINX Ensembl data cache

rule _hmftools_linx_prepare_ensembl: 
    output: 
        gene_data = CFG["dirs"]["inputs"] + "references/ensembl/HG{ensembl_build}/ensembl_gene_data.csv", 
        transcript_data = CFG["dirs"]["inputs"] + "references/ensembl/HG{ensembl_build}/ensembl_trans_exon_data.csv", 
        protein_data = CFG["dirs"]["inputs"] + "references/ensembl/HG{ensembl_build}/ensembl_protein_features.csv", 
        splice_data = CFG["dirs"]["inputs"] + "references/ensembl/HG{ensembl_build}/ensembl_trans_splice_data.csv"
    params: 
        url = op.switch_on_wildcard("ensembl_build", CFG["switches"]["ensembl_url"]),
        jar = "$(dirname $(readlink -e $(which linx)))/sv-linx.jar"
    conda: 
        CFG["conda_envs"]["linx"]
    shell: 
        op.as_one_line("""
        java -cp {params.jar} com.hartwig.hmftools.linx.gene.GenerateEnsemblDataCache
            -ensembl_db {params.url} -ensembl_user "anonymous" -ensembl_pass ""
            -output_dir `dirname {output.splice_data}` -ref_genome_version HG{wildcards.ensembl_build}
        """)
    
def _hmftools_get_ensembl_build(wildcards): 

    # CFG = config["lcr-modules"]["gridss"]

    if wildcards.genome_build == "hg38":
        ensembl_build = "38"
    else:
        ensembl_build = "37"

    ensembl_cache = expand(str(rules._hmftools_linx_prepare_ensembl.output.splice_data), ensembl_build = ensembl_build)

    return ensembl_cache



# Run LINX to cluster and visualize CNV and SV data
rule _hmftools_linx: 
    input: 
        purple_vcf = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.purple.sv.vcf.gz", 
        ensembl_cache = _hmftools_get_ensembl_build, 
        linx_db = str(rules._hmftools_get_linx_db.output)
    output: 
        clusters = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.linx.clusters.tsv", 
        svs = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.linx.svs.tsv"
    log: CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/linx.log"
    resources: 
        **CFG["resources"]["linx"]
    params: 
      alt_build = lambda w: VERSION_MAP[w.genome_build], 
      ensembl_build = lambda w: {
          "grch37": "HG37",
          "hs37d5": "HG37", 
          "hg38": "HG38"
      }[w.genome_build],
      jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
      options = CFG["options"]["linx"]
    conda: 
        CFG["conda_envs"]["linx"]
    threads: 
        CFG["threads"]["linx"]
    shell: 
        op.as_one_line("""
        linx -Xmx{params.jvmheap}m 
            -sample {wildcards.tumour_id} 
            -ref_genome_version {params.ensembl_build} 
            -sv_vcf {input.purple_vcf} 
            -purple_dir `dirname {input.purple_vcf}` 
            -output_dir `dirname {output.clusters}`  
            -gene_transcripts_dir `dirname {input.ensembl_cache}` 
            -fragile_site_file {input.linx_db}/fragile_sites_hmf.{params.alt_build}.csv 
            -line_element_file {input.linx_db}/line_elements.{params.alt_build}.csv 
            -replication_origins_file {input.linx_db}/heli_rep_origins.{params.alt_build}.bed 
            -viral_hosts_file {input.linx_db}/viral_host_ref.csv 
            -known_fusion_file {input.linx_db}/known_fusion_data.{params.alt_build}.csv 
            -check_fusions 
            -check_drivers 
            -write_vis_data 
            {params.options} 
            2>&1 | tee -a {log}
        """)

rule _hmftools_linx_viz: 
    input:
        clusters = rules._hmftools_linx.output.clusters,
        ensembl_cache = _hmftools_get_ensembl_build
    output:
        plots = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot"),
        data = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/data")
    log: CFG["logs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/linx_viz.log"
    resources:
        **CFG["resources"]["linx_viz"]
    params: 
        linx_jar = "$(ls $(dirname $(readlink -e $(which linx)))/*.jar)", 
        circos = "$(which circos)", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        options = CFG["options"]["linx_viz"]
    conda:
        CFG["conda_envs"]["linx"]
    threads: 
        CFG["threads"]["linx_viz"]
    
    shell:
        op.as_one_line("""
        java -Xmx{params.jvmheap}m -cp {params.linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser 
            -sample {wildcards.tumour_id} 
            -gene_transcripts_dir `dirname {input.ensembl_cache}` 
            -plot_out {output.plots} 
            -data_out {output.data} 
            -segment `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_segments.tsv 
            -link `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_sv_data.tsv 
            -exon `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_gene_exon.tsv 
            -cna `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_copy_number.tsv 
            -protein_domain `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_protein_domain.tsv 
            -fusion `dirname {input.clusters}`/{wildcards.tumour_id}.linx.fusions.tsv 
            -circos {params.circos} 
            -threads {threads} 
            {params.options}
            2>&1 | tee -a {log}
        """) 

rule _hmftools_linx_viz_annotate: 
    input:
        plots = rules._hmftools_linx_viz.output.plots, 
        clusters = rules._hmftools_linx.output.clusters,
        svs = rules._hmftools_linx.output.svs, 
        ensembl_cache = _hmftools_get_ensembl_build
    output:
        annotated = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plots_annotated.done",
    log: CFG["logs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/linx_viz_annotate.log"
    resources:
        **CFG["resources"]["linx_viz"]
    params: 
        linx_jar = "$(ls $(dirname $(readlink -e $(which linx)))/*.jar)", 
        circos = "$(which circos)", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        options = CFG["options"]["linx_viz_annotate"]
    conda:
        CFG["conda_envs"]["linx"]
    threads: 
        CFG["threads"]["linx_viz"]
    
    shell:
        op.as_one_line("""
        touch {log};
        tail -n +2 {input.clusters} | grep -v "ARTIFACT" | grep -v "INCOMPLETE" | cut -f 1 | while read cluster; do 
            genes=$(cat {input.svs} | 
                        awk -v c=$cluster 'BEGIN {{FS="\\t"}} {{OFS=";"}} $3 == c {{print $12,$13}}' | 
                        tr ';' '\\n' | sort | uniq | tr '\\n' ','); 
            if [[ $genes =~ [A-Z] ]]; then
                echo $cluster $genes >> {log}; 

                java -Xmx{params.jvmheap}m -cp {params.linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser 
                    -sample {wildcards.tumour_id} 
                    -gene_transcripts_dir `dirname {input.ensembl_cache}` 
                    -plot_out `dirname {input.clusters}`/plot 
                    -data_out `dirname {input.clusters}`/data 
                    -segment `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_segments.tsv 
                    -link `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_sv_data.tsv 
                    -exon `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_gene_exon.tsv 
                    -cna `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_copy_number.tsv 
                    -protein_domain `dirname {input.clusters}`/{wildcards.tumour_id}.linx.vis_protein_domain.tsv 
                    -fusion `dirname {input.clusters}`/{wildcards.tumour_id}.linx.fusions.tsv 
                    -circos {params.circos} 
                    -threads {threads} 
                    -gene $genes
                    -clusterId $cluster 
                    {params.options}
                    >> {log} 2>&1;
            fi;           
        done 
        && touch {output.annotated}
        """)

        

# Symlinks the final output files into the module results directory (under '99-outputs/')

rule _hmftools_purple_output:
    input:
        files = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.purple.{out_file}" 
    output:
        files = CFG["dirs"]["outputs"] + "purple_output/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purple.{out_file}" 
    run:
        op.relative_symlink(input.files, output.files, in_module=True)

rule _hmftools_purple_plots: 
    input:
        plots = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot/{tumour_id}.{plot_name}.png"
    output: 
        plots = CFG["dirs"]["outputs"] + "purple_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{plot_name}.png"
    run: 
        op.relative_symlink(input.plots, output.plots, in_module=True)


rule _hmftools_linx_plots: 
    input:
        plots = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot", 
        annotated = rules._hmftools_linx_viz_annotate.output.annotated
    output: 
        plots = CFG["dirs"]["outputs"] + "linx_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.symlinked"
    shell: 
        op.as_one_line("""
        workdir=$PWD &&
        cd `dirname {output.plots}` && 
        find $workdir/{input.plots} -type f -name "*.png" -exec cp -s {{}} . \; && 
        touch $workdir/{output.plots} && 
        cd $workdir
        """)

rule _hmftools_dispatch: 
    input: 
        files = expand(CFG["dirs"]["outputs"] + "purple_output/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.purple.{out_file}", 
            out_file = purple_out), 
        plots = expand(CFG["dirs"]["outputs"] + "purple_plots/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{plot_name}.png", 
            plot_name = purple_plots),
        linx = rules._hmftools_linx_plots.output.plots
    output: 
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _hmftools_all:
    input:
        expand(
            [
                str(rules._hmftools_dispatch.output.dispatched),
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
