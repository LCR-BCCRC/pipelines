import pandas as pd


def main(args):
    df = pd.read_csv(args.in_file, sep="\t")

    out_df = pd.merge(
        df.groupby("cluster_id")["mutation_id"].nunique().reset_index(),
        df[["cluster_id", "cellular_prevalence"]].drop_duplicates(),
        on="cluster_id"
    )

    out_df.insert(0, "patient_id", args.patient_id)

    out_df = out_df.rename(columns={"mutation_id": "num_snvs"})

    out_df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-p", "--patient-id", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
