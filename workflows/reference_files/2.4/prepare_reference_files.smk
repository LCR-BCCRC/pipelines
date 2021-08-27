#!/usr/bin/env snakemake


##### SETUP #####


import os


##### REFERENCE_FILES MODULE #####


# Must load configfile before including `reference_files.smk`
# Also used for wildcard values and tool versions
configfile: "config/default.yaml"

# Check that the reference_directory is provided
assert "reference_directory" in config, (
    "Re-run snakemake with `--config reference_directory=/path/to/reference_files` in your command."
)

# Confirm that current working directory matched config value
original_dir = os.getcwd()
reference_dir = config["reference_directory"]
os.chdir(reference_dir)


# Include the `reference_files` module
include: os.path.join(original_dir, "reference_files.smk")


##### REFERENCE_FILES TARGETS #####


# Ask for all copies of all output files in `reference_files.smk`
rule all:
    input:
        expand(
            [
                rules.get_genome_fasta_download.output.fasta,
                rules.index_genome_fasta.output.fai,
                rules.get_main_chromosomes_download.output.txt,
                rules.get_main_chromosomes_download.output.bed,
                rules.get_main_chromosomes_download.output.chrx,
                rules.get_main_chromosomes_withY_download.output.txt,
                rules.get_main_chromosomes_withY_download.output.bed,
                rules.store_genome_build_info.output.version,
                rules.store_genome_build_info.output.provider,
                rules.create_bwa_index.output.prefix,
                rules.create_gatk_dict.output.dict, 
                rules.get_gencode_download.output.gtf,
                rules.get_dbsnp_download.output.vcf,
                rules.create_star_index.output.index,
                rules.calc_gc_content.output.wig,
                rules.get_blacklist_download.output.bed, 
                rules.get_repeatmasker_download.output.bed, 
                rules.get_mutect2_pon.output.vcf, 
                rules.get_af_only_gnomad_vcf.output.vcf,
                rules.get_jabba_gc_rds.output.rds,
                rules.get_jabba_map_rds.output.rds,
                rules.get_par_rds.output.rds,
                rules.jabba_pon_make_germline_filter.output.germline,
                rules.jabba_pon_make_pon.output.pon
            ],
            genome_build=config["genome_builds"].keys(),
            bwa_version=config["tools"]["bwa"]["version"],
            gencode_release=config["wildcard_values"]["gencode_release"],
            dbsnp_build=config["wildcard_values"]["dbsnp_build"],
            star_version=config["tools"]["star"]["version"],
            star_overhang=config["wildcard_values"]["star_overhang"],
            gc_window_size=config["wildcard_values"]["gc_window_size"]
        )
