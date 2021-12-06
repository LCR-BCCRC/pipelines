#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     Ryan Morin


##### SETUP #####

import sys, os
from os.path import join

import oncopipe as op

# Setup module and store module-specific configuration in `CONFIG`
CFG = op.setup_module(
    name = "vcf2maf",
    version = "1.3",
    subdirectories = ["inputs","decompressed","vcf2maf","crossmap","outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _vcf2maf_input_vcf,
    _vcf2maf_gnomad_filter_maf,
    _vcf2maf_output_maf,
    _vcf2maf_crossmap,
    _vcf2maf_all

VCF2MAF_GENOME_VERSION_MAP = {
    "grch37": "GRCh37",
    "hg38": "GRCh38",
    "hs37d5": "GRCh37"
}

#set variable for prepending to PATH based on config
VCF2MAF_SCRIPT_PATH = CFG['inputs']['src_dir']

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _vcf2maf_input_vcf:
    input:
        vcf_gz = CFG["inputs"]["sample_vcf_gz"]
    output:
        vcf_gz = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf.gz",
        index = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf_gz, output.vcf_gz)
        op.relative_symlink(input.vcf_gz + ".tbi", output.index)

rule _vcf2maf_annotate_gnomad:
    input:
        vcf = str(rules._vcf2maf_input_vcf.output.vcf_gz),
        normalized_gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.normalized.{genome_build}.vcf.gz")
    output:
        vcf = temp(CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.annotated.vcf")
    conda:
        CFG["conda_envs"]["bcftools"]
    resources: 
        **CFG["resources"]["annotate"]
    threads: 
        CFG["threads"]["annotate"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} -a {input.normalized_gnomad} {input.vcf} -c "INFO/gnomADg_AF:=INFO/AF" -o {output.vcf}
        """)

rule _vcf2maf_run:
    input:
        vcf = str(rules._vcf2maf_annotate_gnomad.output.vcf),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vep_cache = CFG["inputs"]["vep_cache"]
    output:
        maf = temp(CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.maf"),
        vep = temp(CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.annotated.vep.vcf")
    log:
        stdout = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stdout.log",
        stderr = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stderr.log",
    params:
        opts = CFG["options"]["vcf2maf"],
        build = lambda w: VCF2MAF_GENOME_VERSION_MAP[w.genome_build],
        custom_enst = op.switch_on_wildcard("genome_build", CFG["switches"]["custom_enst"])
    conda:
        CFG["conda_envs"]["vcf2maf"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        **CFG["resources"]["vcf2maf"]
    shell:
        op.as_one_line("""
        VCF2MAF_SCRIPT_PATH={VCF2MAF_SCRIPT_PATH};
        PATH=$VCF2MAF_SCRIPT_PATH:$PATH;
        VCF2MAF_SCRIPT="$VCF2MAF_SCRIPT_PATH/vcf2maf.pl";
        if [[ -e {output.maf} ]]; then rm -f {output.maf}; fi;
        if [[ -e {output.vep} ]]; then rm -f {output.vep}; fi;
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*;
        if [[ $(which vcf2maf.pl) =~ $VCF2MAF_SCRIPT ]]; then
            echo "using bundled patched script $VCF2MAF_SCRIPT";
            echo "Using $VCF2MAF_SCRIPT to run {rule} for {wildcards.tumour_id} on $(hostname) at $(date)" > {log.stderr};
            vcf2maf.pl
            --input-vcf {input.vcf}
            --output-maf {output.maf}
            --tumor-id {wildcards.tumour_id}
            --normal-id {wildcards.normal_id}
            --ref-fasta {input.fasta}
            --ncbi-build {params.build}
            --vep-data {input.vep_cache}
            --vep-path $vepPATH
            {params.opts}
            --custom-enst {params.custom_enst}
            --retain-info gnomADg_AF
            >> {log.stdout} 2>> {log.stderr};
        else echo "WARNING: PATH is not set properly, using $(which vcf2maf.pl) will result in error during execution. Please ensure $VCF2MAF_SCRIPT exists." > {log.stderr};fi
        """)

rule _vcf2maf_gnomad_filter_maf:
    input:
        maf = str(rules._vcf2maf_run.output.maf)
    output:
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.maf",
        dropped_maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.dropped.maf.gz"
    params:
        opts = CFG["options"]["gnomAD_cutoff"],
        temp_file = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.dropped.maf"
    shell:
        op.as_one_line("""
        cat {input.maf} | perl -lane 'next if /^(!?#)/; my @cols = split /\t/; @AF_all =split/,/, $cols[114]; $skip=0; for(@AF_all){{$skip++ if $_ > {params.opts}}} if ($skip) {{print STDERR;}} else {{print;}};' > {output.maf} 2>{params.temp_file}
            &&
        gzip {params.temp_file}
            &&
        touch {output.dropped_maf}
        """)

def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

rule _vcf2maf_crossmap:
    input:
        maf = rules._vcf2maf_gnomad_filter_maf.output.maf,
        convert_coord = CFG["inputs"]["convert_coord"],
        chains = get_chain
    output:
        dispatched =  CFG["dirs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.converted"
    log:
        stdout = CFG["logs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.crossmap.stdout.log",
        stderr = CFG["logs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.crossmap.stderr.log"
    conda:
        CFG["conda_envs"]["crossmap"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        **CFG["resources"]["crossmap"]
    params:
        out_name = CFG["dirs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.converted_",
        chain = lambda w: "hg38ToHg19" if "38" in str({w.genome_build}) else "hg19ToHg38",
        file = ".maf"
    shell:
        op.as_one_line("""
        {input.convert_coord}
        {input.maf}
        {input.chains}
        {params.out_name}{params.chain}{params.file}
        crossmap
        > {log.stdout} 2> {log.stderr}
        && touch {output.dispatched}
        """)


rule _vcf2maf_output_maf:
    input:
        maf = str(rules._vcf2maf_gnomad_filter_maf.output.maf),
        maf_converted = str(rules._vcf2maf_crossmap.output.dispatched)
    output:
        maf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_{base_name}.maf"
    params:
        chain = lambda w: "hg38ToHg19" if "38" in str({w.genome_build}) else "hg19ToHg38"
    run:
        op.relative_symlink(input.maf, output.maf)
        op.relative_symlink((input.maf_converted+str("_")+str(params.chain)+str(".maf")), (output.maf[:-4]+str(".converted_")+str(params.chain)+str(".maf")))

# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(str(rules._vcf2maf_output_maf.output.maf), zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"],
            base_name = [CFG["vcf_base_name"]] * len(CFG["runs"]["tumour_sample_id"]))

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)