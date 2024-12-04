#!/usr/bin/env python

import os
import sys
import logging
import traceback
import pandas as pd
import oncopipe as op


def log_exceptions(exctype, value, tb):
    logging.critical(''.join(traceback.format_tb(tb)))
    logging.critical('{0}: {1}'.format(exctype, value))

sys.excepthook = log_exceptions

def main():

    with open(snakemake.log[0], "w") as stdout:
        # Set up logging
        sys.stdout = stdout
        
        try:

            maf_file = snakemake.input[0]

            regions_file = snakemake.input[1]
            regions_format = snakemake.params[0]

            metadata = snakemake.params[1]

            output_file = snakemake.output[0]

            # Return empty dataframe if no lines in MAF
            line_count = count_lines(maf_file)
            if line_count == 1:
                empty_maf = pd.read_table(maf_file, comment="#", sep="\t")
                # Add columns required by workflow
                required_columns = ["seq_type","genome_build","chr_std"]
                empty_maf = empty_maf.assign(**{col:None for col in required_columns if col not in empty_maf.columns})
                write_output(empty_maf, output_file)
                exit()

            maf = maf_add_columns(maf=maf_file, metadata=metadata, wildcards=snakemake.wildcards)

            # Perform filtering
            filtered_maf = maf_filter(
                maf=maf, 
                regions=regions_file,
                regions_format=regions_format
                )

            write_output(filtered_maf, output_file)
        
        except Exception as e:
            logging.error(e, exc_info=1)
            raise

def count_lines(maf):
    with open(maf, "r") as handle:
        total_lines = len(handle.readlines())
    return total_lines

def filter_by_bed(maf, regions):

    # Remove row containing column names
    regions = regions[regions[0].str.contains("chrom")==False]

    # Create common columns between BED and MAF
    regions["chr_std"] = regions.apply(lambda x: "chr" + str(x[0]).replace("chr",""), axis=1)
    regions["genomic_pos_std"] = regions["chr_std"] + ":" + regions[1].map(str)

    maf["chr_std"] = maf.apply(lambda x: "chr" + str(x["Chromosome"]).replace("chr",""), axis=1)
    maf["genomic_pos_std"] = maf["chr_std"] + ":" + maf["Start_Position"].map(str)

    filtered_maf = maf[maf["genomic_pos_std"].isin(regions["genomic_pos_std"])]
    return filtered_maf

def filter_by_maf(maf, regions):

    # Create common column by which to subset MAF
    for df in [maf, regions]:
        df["chr_std"] = df.apply(lambda x: "chr" + str(x["Chromosome"]).replace("chr",""), axis=1)
        df["genomic_pos_std"] = df["chr_std"] + ":" + df["Start_Position"].map(str)

    # Subset the MAF
    filtered_maf = maf[maf["genomic_pos_std"].isin(regions["genomic_pos_std"])]
    return filtered_maf

def maf_filter(maf, regions, regions_format):
    
    if regions_format != "bed":
        regions_df = pd.read_table(regions, comment="#", sep="\t")
    else:
        regions_df = pd.read_table(regions, comment="#", sep="\t", header=None)

    # Return empty dataframe without filtering if df is empty
    if len(maf)==0:
        return maf

    filter_functions = {
        "maf": filter_by_maf,
        "bed": filter_by_bed
        }
    
    return filter_functions[regions_format](maf, regions_df)

def maf_add_columns(maf, metadata, wildcards):
    # Read input MAF as df
    maf = pd.read_table(maf, comment="#", sep="\t")

    sample_id = snakemake.wildcards["tumour_id"]
    seq_type = snakemake.wildcards["seq_type"]
    genome_build = snakemake.wildcards["genome_build"]
    normal_sample_id = snakemake.wildcards["normal_sample_id"]
    pair_status = snakemake.wildcards["pair_status"]

    maf["seq_type"] = seq_type
    maf["genome_build"] = genome_build
    maf["normal_sample_id"] = normal_sample_id
    maf["pair_status"] = pair_status
    
    return maf

def write_output(maf, outfile):
    maf.to_csv(outfile, sep="\t", na_rep="NA", index=False)

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        filename=snakemake.log[1],
        filemode='w'
    )
    
    main()
