#!/usr/bin/env snakemake


##### SETUP #####


# Import package with useful functions for developing analysis modules.
import modutils as md

# Make sure the `CFG` variable doesn't exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

# Setup module and store module-specific configuration in `CFG`.
CFG = md.setup_module(
    config = config, 
    name = "manta", 
    version = "1.0",
    subdirs = ["inputs", "manta", "bedpe", "outputs"]
)

# Define rules to be run locally when using a compute cluster.
localrules: _manta_input_bam, _manta_configure, _manta_output_bedpe, _manta_all


##### RULES #####


# Symlinks the input BAM files into the module output directory (under '00-inputs/').
rule _manta_input_bam:
    input:
        CFG["inputs"].get("sample_bam") or unpack(md.locate_bam(CFG.get("bam_directory")))
    output:
        sample_bam = CFG["dirs"]["inputs"] + "{seq_type}/{sample_id}.bam"
    run:
        md.symlink(input.sample_bam, output.sample_bam)
        md.symlink(input.sample_bam + ".bai", output.sample_bam + ".bai")


# Configures the manta workflow with the input BAM files and reference FASTA file.
rule _manta_configure:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "{seq_type}/{normal_id}.bam",
        config = CFG["inputs"]["manta_config"]
    output:
        runwf = CFG["dirs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        CFG["logs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.log"
    params:
        opts   = md.make_seqtype_specific(CFG["options"]["configure"]),
        fasta  = config["reference"]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("manta") or "envs/manta.yaml"
    shell:
        md.as_one_line("""
        configManta.py {params.opts} --referenceFasta {params.fasta} 
        --runDir "$(dirname {output.runwf})" --tumourBam {input.tumour_bam}
        --normalBam {input.normal_bam} > {log} 2>&1
        """)


# Launches manta workflow in parallel mode and deletes unnecessary files upon success.
rule _manta_run:
    input:
        runwf = rules._manta_configure.output.runwf
    output:
        vcf = CFG["dirs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somaticSV.vcf.gz"
    log:
        CFG["logs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/manta_run.log"
    params:
        opts   = CFG["options"]["manta"]
    conda:
        CFG["conda_envs"].get("manta") or "envs/manta.yaml"
    threads:
        CFG["threads"].get("manta") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("manta") or 1000
    shell:
        md.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log} 2>&1
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


# Fixes sample IDs in VCF header for compatibility with svtools vcftobedpe Otherwise, 
# manta uses the sample name from the BAM read groups, which may not be useful.
rule _manta_fix_vcf_ids:
    input:
        vcf  = rules._manta_run.output.vcf
    output:
        vcf = pipe(CFG["dirs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somaticSV.with_ids.vcf")
    log:
        CFG["logs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/manta_fix_vcf_ids.log"
    shell:
        md.as_one_line("""
        gzip -dc {input.vcf}
            |
        awk 'BEGIN {{FS=OFS="\\t"}}
        $1 == "#CHROM" {{$10="{wildcards.normal_id}"; $11="{wildcards.tumour_id}"}}
        {{print $0}}' > {output.vcf} 2> {log}
        """)


# Calculates the tumour and normal variant allele fraction (VAF) from the allele counts
# and creates new fields in the INFO column for convenience.
rule _manta_calc_vaf:
    input:
        vcf  = rules._manta_fix_vcf_ids.output.vcf,
        cvaf = CFG["inputs"]["calc_manta_vaf"]
    output:
        vcf = pipe(CFG["dirs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somaticSV.with_ids.with_vaf.vcf")
    log:
        CFG["logs"]["manta"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/manta_calc_vaf.log"
    conda:
        CFG["conda_envs"].get("manta") or "envs/manta.yaml"
    shell:
        "{input.cvaf} {input.vcf} > {output.vcf} 2> {log}"


# Converts the VCF file into a more tabular BEDPE file, which is easier to handle in R
# and automatically pairs up breakpoints for interchromosomal events.
rule _manta_vcf_to_bedpe:
    input:
        vcf  = rules._manta_calc_vaf.output.vcf
    output:
        bedpe = CFG["dirs"]["bedpe"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bedpe"
    log:
        CFG["logs"]["bedpe"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}/manta_vcf_to_bedpe.log"
    conda:
        CFG["conda_envs"].get("manta") or "envs/manta.yaml"
    threads:
        CFG["threads"].get("vcf_to_bedpe") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("vcf_to_bedpe") or 1000
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log}"


# Symlinks the final BEDPE file
rule _manta_output_bedpe:
    input:
        bedpe = rules._manta_vcf_to_bedpe.output.bedpe
    output:
        bedpe = CFG["dirs"]["outputs"] + "{seq_type}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    run:
        md.symlink(input.bedpe, output.bedpe)


rule _manta_all:
    input:
        vcfs = expand(rules._manta_vcf_to_bedpe.output.bedpe, zip,
                      seq_type=CFG["paired_runs"]["tumour_seq_type"],
                      tumour_id=CFG["paired_runs"]["tumour_sample_id"],
                      normal_id=CFG["paired_runs"]["normal_sample_id"],
                      pair_status=CFG["paired_runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete the CFG variable to avoid interfering with other code
del CFG
