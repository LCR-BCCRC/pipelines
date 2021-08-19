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

# Added by Chris. Pair the capture space with the appropriate genome build
# If we don't specify this seperately, snakemake tries to pair all genome builds
# with all capture spaces
def generate_capture_targets():
    genome_builds = config["genome_builds"].keys() # For sanity checking later
    cap_targets = []
    for capspace in config["capture_space"].keys():
        genome_version = config["capture_space"][capspace]["genome"]
        for genome_build in GENOME_VERSION_GROUPS[genome_version]:

            # Sanity check to make sure this genome build actually exists
            if genome_build not in genome_builds:
                raise AttributeError("When processing the capspace \'%s\', reference \'%s\' was specified but does not exist as a genome build" % (capspace, genome_build))
            # Add all rules for the capture space processing here
            tabix_target = expand(rules.compress_index_capspace_bed.output.tabix, capture_space=capspace, genome_build=genome_build)
            interval_list = expand(rules.create_interval_list.output.interval_list, capture_space=capspace, genome_build=genome_build)
            contig_log = expand(rules.check_capspace_contigs.output.contig_log, capture_space=capspace, genome_build=genome_build)

            cap_targets.append(tabix_target)
            cap_targets.append(interval_list)
            cap_targets.append(contig_log)

    return cap_targets

##### REFERENCE_FILES TARGETS #####

rule all:
    input:
        expand(
            [
                # DO NOT PUT CAPTURE RELATED RULES HERE! Put them in generate_capture_targets()
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
                rules.prepare_mutect2_pon.output.vcf, 
                rules.prepare_af_only_gnomad_vcf.output.vcf
            ],
            genome_build=config["genome_builds"].keys(),
            bwa_version=config["tools"]["bwa"]["version"],
            gencode_release=config["wildcard_values"]["gencode_release"],
            dbsnp_build=config["wildcard_values"]["dbsnp_build"],
            star_version=config["tools"]["star"]["version"],
            star_overhang=config["wildcard_values"]["star_overhang"],
            gc_window_size=config["wildcard_values"]["gc_window_size"],
        ),
        generate_capture_targets()
