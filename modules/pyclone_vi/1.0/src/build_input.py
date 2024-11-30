from operator import index
import numpy as np
import pandas as pd


def main(args):
    df = build_input_df(
        args.cnv_file,
        args.snv_file,
        sample_id=args.sample_id,
        cnv_caller=args.cnv_caller,
        sex=args.sex,
        tumour_content=args.tumour_content,
        pyclone_tool=args.pyclone_tool
    )

    assert np.all(df["major_cn"] >= df["minor_cn"])

    df = df[df["major_cn"] > 0]

    df.to_csv(args.out_file, index=False, sep="\t")


def build_input_df(cnv_file, snv_file, sample_id="tumour", cnv_caller="battenberg", sex=None, tumour_content=1.0, pyclone_tool="pyclone-vi"):
    if cnv_caller == "battenberg":
        cnv_df = load_battenberg_cnv_df(cnv_file, sample_id=sample_id, sex=sex)
    else:
        cnv_df = load_sequenza_cnv_df(cnv_file, sex=sex)
    snv_df = load_snv_df(snv_file, sample_id=sample_id)
    df = merge_files(cnv_df, snv_df)
    df["tumour_content"] = tumour_content
    if pyclone_tool == "pyclone":
        df = df[["mutation_id", "ref_counts", "alt_counts", "normal_cn", "major_cn", "minor_cn", "tumour_content", "sample_id"]].rename(columns = {"alt_counts": "var_counts"})
    return df


def load_battenberg_cnv_df(file_name, sample_id="tumour", sex=None, solution="A"):
    def get_dominant_clone(x, solution="A"):
        frac_1_key = "frac1_{}".format(solution)
        frac_2_key = "frac2_{}".format(solution)
        if np.isnan(x[frac_2_key]):
            clone = 1
        else:
            if x[frac_1_key] >= x[frac_2_key]:
                clone = 1
            else:
                clone = 2
        return clone

    def get_major_cn(x, solution="A"):
        clone = get_dominant_clone(x, solution=solution)
        return x["nMaj{clone}_{solution}".format(clone=clone, solution=solution)]

    def get_minor_cn(x, solution="A"):
        clone = get_dominant_clone(x, solution=solution)
        return x["nMin{clone}_{solution}".format(clone=clone, solution=solution)]

    df = pd.read_csv(file_name, sep="\t")
    df = df.rename(columns={"chr": "chrom", "startpos": "start", "endpos": "end"})
    df["major_cn"] = df.apply(lambda row: get_major_cn(row, solution=solution), axis=1)
    df["minor_cn"] = df.apply(lambda row: get_minor_cn(row, solution=solution), axis=1)
    df["normal_cn"] = df["chrom"].apply(lambda row: get_normal_cn(row, sex=sex))
    df["sample_id"] = sample_id
    df = df[["sample_id", "chrom", "start", "end", "normal_cn", "major_cn", "minor_cn"]]
    if sex is None:
        df = df[~df["chrom"].isin(["X", "Y"])]
    df.to_csv("loaded_battenberg.tsv", sep="\t", index=False)
    return df


def load_sequenza_cnv_df(file_name, sex=None):
    df = pd.read_csv(file_name, sep="\t")
    df = df.rename(columns={"chromosome": "chrom", "start.pos": "start", "end.pos": "end"})
    df["major_cn"] = df.apply(lambda row: max(row["A"], row["B"]), axis=1)
    df["minor_cn"] = df.apply(lambda row: min(row["A"], row["B"]), axis=1)
    df["normal_cn"] = df["chrom"].apply(lambda row: get_normal_cn(row, sex=sex))
    df = df[["chrom", "start", "end", "normal_cn", "major_cn", "minor_cn"]]
    if sex is None:
        df = df[~df["chrom"].isin(["X", "Y"])]
    df["sample_id"] = "tumour"
    df = df.dropna()
    return df


def get_normal_cn(chrom, sex):
    if sex == "M":
        if chrom in ["X", "Y"]:
            cn = 1
        else:
            cn = 2
    else:
        cn = 2
    return cn


def load_snv_df(file_name, sample_id="tumour"):
    df = pd.read_csv(file_name, sep="\t")
    # PyClone only works on SNPs, not InDels
    # df = df[df["Variant_Type"].isin(["SNP"])]
    # Ignore intergenic mutations (IGR)
    df = df[~df["Variant_Classification"].isin(["IGR"])]
    # Acutally I can't do this sub-sampling here because the mutations in all files for all tumours
    # per patient need to be the same! Have to do this at the merging step.
    # Except that by the time this file is generated, info about coding status of
    # mutations is lost. Need to think about this.
    # Separate the df into coding and non-coding mutations first
    # df_coding = df[~df["Variant_Classification"].isin(["3'Flank", "5'Flank", "Intron"])]
    # df_noncoding = df[df["Variant_Classification"].isin(["3'Flank", "5'Flank", "Intron"])]
    # # To get up to 5000 total mutations, first check how many coding mutations there are:
    # if len(df_coding.index) >= 5000:
    #     df = df_coding.sample(n = 5000)
    # else:
    #     to_add = min([(5000 - len(df_coding.index)), len(df_noncoding.index)])
    #     df_noncoding = df_noncoding.sample(n = to_add)
    #     df = pd.concat([df_coding, df_noncoding])
    df = df.rename(columns={
        "Chromosome": "chrom",
        "Start_Position": "coord",
        "Reference_Allele": "ref",
        "Tumor_Seq_Allele2": "alt",
        "t_ref_count": "ref_counts",
        "t_alt_count": "alt_counts"
    })
    df = df[["chrom", "coord", "ref", "alt", "ref_counts", "alt_counts"]]
    df["mutation_id"] = df.apply(lambda row: "{chrom}:{coord}:{ref}:{alt}".format(**row.to_dict()), axis=1)
    df["sample_id"] = sample_id
    return df


def position_segment_merge(positions, segments):
    """
    Merge positions with segments that contain them

    Args:
        positions (pandas.DataFrame): ['chrom', 'coord'] columns required
        segments (pandas.DataFrame): ['chrom', 'start', 'end'] columns required

    Returns:
        pandas.DataFrame: merged table with ['chrom', 'coord', 'start', 'end'] columns


    Assuming a set of non-overlapping segments, merge a set of positions so that
    each entry in the new table provides coord and containing segment start/end
    """

    positions = positions[['chrom', 'coord']].copy().sort_values(by=['chrom','coord'])
    positions["chrom"] = positions["chrom"].astype(str)
    positions["chrom"] = positions["chrom"].str.strip()
    segments = segments[['chrom', 'start', 'end']].copy().sort_values(by=['chrom','start'])
    segments["chrom"] = segments["chrom"].astype(str)
    segments["chrom"] = segments["chrom"].str.strip()

    merged = positions.merge(segments, left_on='chrom', right_on='chrom', how="left")\
               .sort_values(by=['chrom', 'coord'])


    merged['start'] = merged['start'].fillna(method='ffill')
    merged['end'] = merged['end'].fillna(method='ffill')


    merged = merged[(merged['coord'] >= merged['start']) &
                    (merged['coord'] <= merged['end'])]

    return merged


def merge_files(cnv_df, snv_df):
    df = []
    for sample_id in cnv_df['sample_id'].unique():
        df_1 = snv_df[snv_df['sample_id'] == sample_id]
        df_2 = cnv_df[cnv_df['sample_id'] == sample_id]
        df_2["chrom"] = df_2["chrom"].astype(str).str.strip()
        merged = position_segment_merge(df_1, df_2)
        merged["chrom"] = merged["chrom"].astype(str).str.strip()
        merged = pd.merge(merged, df_2, on=['chrom', 'start', 'end'])
        merged = merged[['chrom', 'coord', 'major_cn', 'minor_cn', 'normal_cn']]
        # df_1.to_csv("snv_df.tsv", sep="\t", index=False)
        df_1["chrom"] = df_1["chrom"].astype(str).str.strip()
        merged["chrom"] = merged["chrom"].astype(str).str.strip()
        merged = pd.merge(merged, df_1, on=['chrom', 'coord'])
        # merged.to_csv("merged.tsv", sep="\t", index=False)
        merged = merged[['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn']]
        df.append(merged)
    df = pd.concat(df)
    df.drop_duplicates(subset=['mutation_id', 'sample_id'], inplace=True)
    df['major_cn'], df['minor_cn'] = df['major_cn'].astype('Int64'), df['minor_cn'].astype('Int64')
    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ic", "--cnv-file", required=True)

    parser.add_argument("-is", "--snv-file", required=True)

    parser.add_argument("-id", "--sample-id", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-c", "--cnv-caller", default="battenberg", choices=["battenberg", "sequenza"])

    parser.add_argument("-s", "--sex", default=None, choices=["M", "F"])

    parser.add_argument("-t", "--tumour-content", default=1.0, type=float)

    parser.add_argument("-p", "--pyclone_tool", default = "pyclone-vi", choices = ["pyclone-vi", "pyclone"])

    cli_args = parser.parse_args()

    main(cli_args)
