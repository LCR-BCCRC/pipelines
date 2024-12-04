#!/usr/bin/env python3

import os
import warnings
import numpy as np
import pandas as pd
import oncopipe as op
import copy
import sys
import logging
import traceback

def log_exceptions(exctype, value, tb):
    logging.critical(''.join(traceback.format_tb(tb)))
    logging.critical('{0}: {1}'.format(exctype, value))

sys.excepthook = log_exceptions

def main():
    with open(snakemake.log[0], "w") as stdout:
        # Set up logging
        sys.stdout = stdout

        try:
            input_bam = snakemake.input["bam_file"]
            input_bai = snakemake.input["bai_file"]
            
            maf = snakemake.input["filter_maf"]

            batch_options = snakemake.params["batch_options"]

            # Print run info for logging
            print(f"Setting up batch scripts using the following inputs:\n\
            Bam files:\t{input_bam}\n\
            Bai files:\t{input_bai}\n\
            Filtered maf:\t{maf}\n\
            Batch options:\t{batch_options}")

            if not isinstance(maf, list):
                maf = [maf]

            empty_mafs = []

            for m in maf:
                # Skip if no variants in outfile
                input_maf = open(m, "r")

                line_count = 0
                for line in input_maf:
                    line_count += 1
                    if line_count > 1:
                        # Return to top of MAF 
                        input_maf.seek(0)
                        break
                if line_count < 2:
                    input_maf.close()
                    empty_mafs.append(m)

            if len(empty_mafs) != 0:
                if all(m in empty_mafs for m in maf):
                    touch_outputs(
                        output_dir = snakemake.params["batch_dir"],
                        seq_type = snakemake.wildcards["seq_type"],
                        genome_build = snakemake.wildcards["genome_build"],
                        presets = snakemake.params["igv_presets"],
                        sample_id = snakemake.wildcards["sample_id"],
                        suffix = snakemake.params["suffix"],
                        finished_file = snakemake.output["finished"]
                    )
                    exit()
                for e in empty_mafs:
                    maf.remove(e)

            # Read MAF file and create dataframe
            regions = get_regions_df(
                maf,
                padding=batch_options["padding"]
            )

            # Create the batch scripts
            generate_igv_batches(
                regions = regions,
                bam = input_bam,
                bai = input_bai,
                output_dir = snakemake.params["batch_dir"],
                snapshot_dir = snakemake.params["snapshot_dir"],
                genome_build = snakemake.params["genome_build"],
                seq_type = snakemake.params["seq_type"],
                suffix = snakemake.params["suffix"],
                igv_presets = snakemake.params["igv_presets"],
                igv_options = batch_options["igv_options"],
                max_height = batch_options["max_height"],
                tissue_status = snakemake.params["tissue_status"],
                sleep_timer = batch_options["sleep_timer"]
            )

            touch_output = open(snakemake.output[0], "w")
            touch_output.close()
        
        except Exception as e:
            logging.error(e, exc_info=1)
            raise

def touch_outputs(output_dir, seq_type, genome_build, presets, sample_id, suffix, finished_file):
    sample_suffix = sample_id + suffix + ".batch"
    for preset in presets:
        os.makedirs(os.path.join(output_dir, "merged_batch_scripts", "--".join([seq_type, genome_build]), preset), exist_ok = True)
        merged_batch = os.path.join(output_dir, "merged_batch_scripts", "--".join([seq_type, genome_build]), preset, sample_suffix)
        merged_file = open(merged_batch, "w")
        merged_file.close()
    touch_finished = open(finished_file, "w")
    touch_finished.close()

def get_regions_df(input_maf, padding):
    # Read MAF as dataframe
    if len(input_maf) > 1:
        maf = pd.concat([pd.read_table(file, comment="#", sep="\t") for file in input_maf])
    else:
        maf = pd.read_table(input_maf[0], comment="#", sep="\t")

    chrom = (maf["Chromosome"].astype(str)).apply(lambda x: x.replace("chr",""))

    # Specify regions that will be captured by IGV based on variant positions
    region_start = (maf["Start_Position"]).astype(str)
    snapshot_start = (maf["Start_Position"] - padding).astype(str)
    snapshot_end = (maf["End_Position"] + padding).astype(str)
    snapshot_coordinates = "chr" + chrom + ":" + snapshot_start + "-" + snapshot_end
    regions = "chr" + chrom + ":" + region_start

    regions_df = pd.DataFrame(
        {"chromosome": "chr" + chrom,
        "region": regions,
        "region_name": maf.Hugo_Symbol,
        "tumour_id": maf.Tumor_Sample_Barcode,
        "normal_id":  maf.Matched_Norm_Sample_Barcode,
        "ref_allele": maf.Reference_Allele,
        "alt_allele": maf.Tumor_Seq_Allele2,
        "snapshot_coordinates": snapshot_coordinates,
        "padding": padding,
        "pair_status": maf.pair_status
        }
    )

    samp_id = snakemake.wildcards["sample_id"]

    assert len(regions_df["normal_id"].drop_duplicates()) == 1, f"More than one normal ID found within the MAF files' `Matched_Norm_Sample_Barcode` column for this sample: {samp_id}. Please double check MAF files: {input_maf}"

    return regions_df

def output_lines(lines, batch_output):
    output = open(batch_output, "w")
    lines.append("")
    text = "\n".join(lines)
    output.write(text)
    output.close()

def generate_igv_batch_per_row(sleep_interval, preset, options, coordinates, directory, child_dir, seq_build, chrom_directory, snapshot_filename):
    lines = []

    lines.append(f"goto {coordinates}")

    snapshot_regions_dir = os.path.join(directory, seq_build, child_dir, preset, chrom_directory, "")

    lines.append(f"snapshotDirectory {snapshot_regions_dir}")
    lines.append("collapse")
    for igv_option in options[preset]:
        lines.append(igv_option)
    lines.append(f"setSleepInterval {sleep_interval}")
    lines.append(f"snapshot {snapshot_filename}")
    
    return lines

def generate_igv_batch_header(bam, index, max_height, genome_build):
    lines = []

    genome_build = genome_build.replace("grch37","hg19")

    lines.append(f"load {bam} index={index}")
    lines.append(f"maxPanelHeight {max_height}")
    lines.append(f"genome {genome_build}")

    return lines

def generate_igv_batches(regions, bam, bai, output_dir, snapshot_dir, genome_build, seq_type, suffix, igv_presets, igv_options, max_height, tissue_status, sleep_timer=2000):
    for preset in igv_presets:

        merged_batch_suffix = snakemake.wildcards["sample_id"] + suffix + ".batch"
        
        grouped_regions = regions.groupby("region")

        # Group by genomic coordinate
        for coordinate, variants in grouped_regions:
            all_lines = []

            header = generate_igv_batch_header(bam=bam, index=bai, max_height=max_height, genome_build=genome_build)
            all_lines.extend(header)

            seq_type_build = f"{seq_type}--{genome_build}"
            chrom_dir = variants["chromosome"].unique()[0]

            filename = []
            filename.append(variants["region"].unique()[0])
            filename.append(variants["region_name"].unique()[0])

            batch_filename = filename.copy()
            batch_filename.append(snakemake.wildcards["sample_id"])
            batch_filename = "--".join(batch_filename) + suffix + ".batch"

            if tissue_status == "tumour":
                # Iterate over variants at same position so that instructions to create multiple snapshots 
                # with different alleles in their filenames are added to the batch script.

                for _, row in variants.iterrows():

                    snap_filename = filename.copy()
                    snap_filename.append(f"{row.ref_allele}_{row.alt_allele}")
                    snap_filename.append(snakemake.wildcards["sample_id"])
                    snap_filename = "--".join(snap_filename) + suffix + ".png"

                    lines = generate_igv_batch_per_row(
                        sleep_interval = sleep_timer,
                        preset = preset,
                        options = igv_options,
                        coordinates = row.snapshot_coordinates,
                        directory = snapshot_dir,
                        child_dir = tissue_status,
                        seq_build = seq_type_build,
                        chrom_directory = chrom_dir,
                        snapshot_filename = snap_filename
                    )

                    all_lines.extend(lines)

            elif tissue_status == "normal":

                # Only need to create one snapshot since only one allele will be expected from the ref

                snap_filename = filename.copy()
                ref_allele = variants["ref_allele"].unique()[0]
                snap_filename.append(f"{ref_allele}_{ref_allele}")
                snap_filename.append(snakemake.wildcards["sample_id"])
                snap_filename = "--".join(snap_filename) + suffix + ".png"

                lines = generate_igv_batch_per_row(
                    sleep_interval = sleep_timer,
                    preset = preset,
                    options = igv_options,
                    coordinates = variants["snapshot_coordinates"].unique()[0],
                    directory = snapshot_dir,
                    child_dir = tissue_status,
                    seq_build = seq_type_build,
                    chrom_directory = chrom_dir,
                    snapshot_filename = snap_filename
                )

                all_lines.extend(lines)

            # Make subdirectories if necessary because snakemake won't make them since rule is a checkpoint
            os.makedirs(os.path.join(output_dir, "single_batch_scripts", seq_type_build, preset), exist_ok=True)

            batch_file_path = os.path.join(output_dir, "single_batch_scripts", seq_type_build, preset, batch_filename)

            output_lines(all_lines, batch_file_path)
        
        os.makedirs(os.path.join(output_dir, "merged_batch_scripts", seq_type_build, preset), exist_ok=True)

        merged_preset_path = os.path.join(output_dir, "merged_batch_scripts", seq_type_build, preset, merged_batch_suffix)

        merged_preset_touch = open(merged_preset_path, "w")
        merged_preset_touch.close()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        filename=snakemake.log[1],
        filemode='w'
    )
    main()
