#!/usr/bin/env python

"""
This script reformats MAF, HotMAPS, OncodriveCLUSTL or GENOMIC_POSITION files into BED files
"""

import os
import sys
import pandas as pd
import oncopipe as op
import shutil
import logging
import traceback

def log_exceptions(exctype, value, tb):
    logging.critical(''.join(traceback.format_tb(tb)))
    logging.critical('{0}: {1}'.format(exctype, value))

sys.excepthook = log_exceptions

def main():

    with open(snakemake.log[0], "w") as stdout:
        sys.stdout = stdout

        try:
            regions_file = snakemake.input["regions"]
            regions_format = snakemake.params["regions_format"]

            output_file = snakemake.output["regions"]

            line_count = 0
            with open(regions_file, "r") as handle:
                for line in handle:
                    line_count += 1
                    if line_count > 1:
                        break
            if line_count < 2:
                touch_output = open(output_file, "w")
                touch_output.close()
                exit()

            if regions_format == "mutation_id":
                global REGIONS_BUILD
                REGIONS_BUILD = snakemake.params["regions_build"]
                REGIONS_BUILD = REGIONS_BUILD.lower()

            if regions_format == "bed":
                # Do not need to reformat for liftover
                shutil.copy(regions_file, output_file)
                exit()

            # Reformat for liftover based on regions format
            regions_formatted = format_regions(regions_file, regions_format)
            regions_formatted = regions_formatted.drop_duplicates()

            # Add empty column to end for liftover
            regions_formatted.insert(3, '', '')

            # Output regions file
            regions_formatted.to_csv(output_file, sep="\t", index=False)
        
        except Exception as e:
            logging.error(e, exc_info=1)
            raise

def format_mutation_id(mutation_id):
    # Read regions into dataframe
    mutation_id = pd.read_table(mutation_id, comment="#", sep="\t")

    # Create columns required for liftover in BED format
    genomic_pos_col = f"mutation_id_{REGIONS_BUILD}"

    for col, idx in {"chr_std": 0, "start": 1, "end": 1}.items():
        mutation_id[col] = mutation_id.apply(lambda x: str(x[genomic_pos_col]).split(":")[idx].replace("chr",""), axis=1)
    end = mutation_id["end"].apply(lambda x: int(x) + 1)

    mutation_id_reformatted = pd.DataFrame(
        {
            "chrom": "chr" + mutation_id["chr_std"],
            "start": mutation_id["start"],
            "end": end
        }
    )

    # Remove duplicate rows
    mutation_id_reformatted = mutation_id_reformatted.drop_duplicates(keep='first')

    return mutation_id_reformatted

def format_hotmaps(hotmaps_regions):
    # Read regions into dataframe
    hotmaps_regions = pd.read_table(hotmaps_regions, comment="#", sep="\t")

    # Convert HotMAPS coordinates to BED format

    hotmaps_regions["chr_std"] = hotmaps_regions.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
    chr_std = "chr" + hotmaps_regions["chr_std"].map(str)
    end = hotmaps_regions["Start_Position"].apply(lambda x: x + 1)

    hotmaps_reformatted = pd.DataFrame(
        {
            "chrom": chr_std,
            "start": hotmaps_regions["Start_Position"],
            "end": end
        }
    )
    return hotmaps_reformatted

def format_clustl(clustl_regions):
    # Read regions into dataframe
    clustl_regions = pd.read_table(clustl_regions, comment="#", sep="\t")

    # Create columnns required for BED format
    chr_std = "chr" + clustl_regions["Chromosome"].map(str)
    end = clustl_regions["Start_Position"].apply(lambda x: x + 1)
    clustl_reformatted = pd.DataFrame(
        {
            "chrom": chr_std,
            "start": clustl_regions["Start_Position"],
            "end": end
        }
    )
    return clustl_reformatted

def format_maf(maf):
    # Read regions into dataframe
    maf_regions = pd.read_table(maf, comment="#", sep="\t")

    # Create dataframe in BED format
    chr_std = "chr" + maf_regions["Chromosome"].map(str).replace("chr","")
    end = maf_regions["Start_Position"].apply(lambda x: x + 1)

    maf_reformatted = pd.DataFrame(
        {
            "chrom": chr_std,
            "start": maf_regions["Start_Position"],
            "end": end
        }
    )

    return maf_reformatted

def format_regions(regions, regions_format):
    format_functions = {
        "oncodriveclustl": format_clustl,
        "hotmaps": format_hotmaps,
        "mutation_id": format_mutation_id,
        "maf": format_maf
    }

    return format_functions[regions_format](regions)

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        filename=snakemake.log[1],
        filemode='w'
    )
    main()
