#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     Bruno Grande


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
# `CFG` is a shortcut to `config["lcr-modules"]["sequenza"]`
CFG = op.setup_module(
    name = "sequenza",
    version = "1.4",
    subdirectories = ["inputs", "seqz", "sequenza", "igv_seg", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _sequenza_input_bam,
    _sequenza_input_chroms,
    _sequenza_input_dbsnp_pos,
    _sequenza_cnv2igv,
    _sequenza_output_seg,
    _sequenza_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _sequenza_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Pulls in list of chromosomes for the genome builds
checkpoint _sequenza_input_chroms:
    input:
        txt = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        txt = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)


# Pulls in list of chromosomes for the genome builds
rule _sequenza_input_dbsnp_pos:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        pos = CFG["dirs"]["inputs"] + "dbsnp/{genome_build}/dbsnp.common_all-151.pos"
    log:
        stderr = CFG["logs"]["inputs"] + "{genome_build}/sequenza_input_dbsnp_pos.stderr.log"
    resources:
        **CFG["resources"]["vcf_sort"]
    shell:
        op.as_one_line("""
        gzip -dc {input.vcf}
            |
        awk 'BEGIN {{FS="\t"}} $0 !~ /^#/ {{print $1 ":" $2}}' 2>> {log.stderr}
            |
        LC_ALL=C sort -S {resources.mem_mb}M > {output.pos} 2>> {log.stderr}
        """)


rule _sequenza_bam2seqz:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        gc_wiggle = reference_files("genomes/{genome_build}/annotations/gc_wiggle.window_50.wig.gz"),
        genome = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        seqz = temp(CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.binned.seqz.gz")
    log:
        stderr = CFG["logs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_bam2seqz.{chrom}.stderr.log"
    params:
        bam2seqz_opts = CFG["options"]["bam2seqz"],
        seqz_binning_opts = CFG["options"]["seqz_binning"]
    conda:
        CFG["conda_envs"]["sequenza-utils"]
    threads:
        CFG["threads"]["bam2seqz"]
    resources:
        **CFG["resources"]["bam2seqz"]    
    shell:
        op.as_one_line("""
        sequenza-utils bam2seqz {params.bam2seqz_opts} -gc {input.gc_wiggle} --fasta {input.genome} 
        --normal {input.normal_bam} --tumor {input.tumour_bam} --chromosome {wildcards.chrom} 2>> {log.stderr}
            | 
        sequenza-utils seqz_binning {params.seqz_binning_opts} --seqz - 2>> {log.stderr}
            |
        gzip > {output} 2>> {log.stderr}
        """)


def _sequenza_request_chrom_seqz_files(wildcards):
    CFG = config["lcr-modules"]["sequenza"]
    with open(checkpoints._sequenza_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    seqz_files = expand(
        CFG["dirs"]["seqz"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.binned.seqz.gz", 
        chrom=mains_chroms
    )
    return seqz_files


rule _sequenza_merge_seqz:
    input:
        seqz = _sequenza_request_chrom_seqz_files,
        merge_seqz = CFG["inputs"]["merge_seqz"],
        gc = reference_files("genomes/{genome_build}/annotations/gc_wiggle.window_50.wig.gz")
    output:
        seqz = CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.binned.unfiltered.seqz.gz"
    log:
        stderr = CFG["logs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_merge_seqz.stderr.log"
    threads:
        CFG["threads"]["merge_seqz"]
    resources: 
        **CFG["resources"]["merge_seqz"]  
    shell:
        op.as_one_line("""
        bash {input.merge_seqz} {input.seqz} 2>> {log.stderr}
            | 
        gzip > {output.seqz} 2>> {log.stderr}
        """)


rule _sequenza_filter_seqz:
    input:
        seqz = str(rules._sequenza_merge_seqz.output.seqz),
        filter_seqz = CFG["inputs"]["filter_seqz"],
        dbsnp_pos = str(rules._sequenza_input_dbsnp_pos.output.pos),
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed")
    output:
        seqz = CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.binned.filtered.seqz.gz"
    log:
        stderr = CFG["logs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_filter_seqz.stderr.log"
    conda:
        CFG["conda_envs"]["bedtools"]
    threads: 
        CFG["threads"]["filter_seqz"]
    resources: 
        **CFG["resources"]["filter_seqz"] 
    shell:
        op.as_one_line("""
        SEQZ_BLACKLIST_BED_FILES='{input.blacklist}'
        {input.filter_seqz} {input.seqz} {input.dbsnp_pos} 2>> {log.stderr}
            |
        gzip > {output.seqz} 2>> {log.stderr}
        """)


rule _sequenza_run:
    input:
        seqz = CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.binned.{filter_status}.seqz.gz",
        run_sequenza = CFG["inputs"]["run_sequenza"],
        assembly = reference_files("genomes/{genome_build}/version.txt"),
        chroms = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt"),
        x_chrom = reference_files("genomes/{genome_build}/genome_fasta/chromosome_x.txt")
    output:
        segments = CFG["dirs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{filter_status}/sequenza_segments.txt"
    log:
        stdout = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_run.{filter_status}.stdout.log",
        stderr = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_run.{filter_status}.stderr.log"
    conda:
        CFG["conda_envs"]["r-sequenza"]
    threads:
        CFG["threads"]["sequenza"]
    resources: 
        **CFG["resources"]["sequenza"] 
    shell:
        op.as_one_line("""
        Rscript {input.run_sequenza} {input.seqz} {input.assembly} {input.chroms} {input.x_chrom} 
        $(dirname {output.segments}) {threads} > {log.stdout} 2> {log.stderr}
        """)


rule _sequenza_cnv2igv:
    input:
        segments = str(rules._sequenza_run.output.segments),
        cnv2igv =  CFG["inputs"]["cnv2igv"]
    output:
        igv = CFG["dirs"]["igv_seg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{filter_status}/sequenza_segments.igv.seg"
    log:
        stderr = CFG["logs"]["igv_seg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cnv2igv.{filter_status}.stderr.log"
    params:
        opts = CFG["options"]["cnv2igv"]
    conda:
        CFG["conda_envs"]["cnv2igv"]
    shell:
        op.as_one_line("""
        python {input.cnv2igv} {params.opts} --sample {wildcards.tumour_id} 
        {input.segments} > {output} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sequenza_output_seg:
    input:
        seg = str(rules._sequenza_cnv2igv.output.igv)
    output:
        seg = CFG["dirs"]["outputs"] + "{filter_status}_seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.igv.seg"
    run:
        op.relative_symlink(input.seg, output.seg, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _sequenza_all:
    input:
        expand(
            [
                str(rules._sequenza_output_seg.output.seg).replace("{filter_status}", "filtered"),
                str(rules._sequenza_output_seg.output.seg).replace("{filter_status}", "unfiltered"),
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
