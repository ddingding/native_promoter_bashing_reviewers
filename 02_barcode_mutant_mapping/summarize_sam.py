import sys

sys.path.append("../src/")
from pacbio import summarize_sam, get_bc_from_muts, clean_muts
import argparse
import pandas as pd

s_to_bc_loc = {
    "2A": (140, 180),
    "2B": (140, 180),
    "2C": (140, 180),
    "3A": (1420, 1460),
    "3B": (645, 685),
    "3C": (1330, 1370),
}

parser = argparse.ArgumentParser()

parser.add_argument("--samfile", dest="samfile", type=str)
parser.add_argument("--file_out_name", dest="file_out_name", type=str)
parser.add_argument("--df_out_name", dest="df_out_name", type=str)
parser.add_argument("--sample", dest="sample", type=str)

args = parser.parse_args()


if __name__ == "__main__":
    sample = args.sample
    ##### 1. summarize sam file
    summarize_sam(
        args.samfile, args.file_out_name
    )  # writes a tab separated file per sam file

    ##### 2. cleaning code
    df = pd.read_csv(args.file_out_name, delimiter="\t", header=None)
    df = df.rename(columns={0: "bc", 1: "muts"})

    # drop nan values
    print("before dropping nan:", len(df))
    df = df.dropna()
    print("after dropping nan:", len(df))

    # flag all the erros
    df["valid_mut"] = (df.muts != "typeError") & (df.muts != "noMDtag")
    df_muts = df.loc[df.valid_mut]
    print("valid muts", len(df_muts))

    # get bc
    df_muts["bc_from_mut"] = df_muts.muts.apply(
        lambda x: get_bc_from_muts(x, s_to_bc_loc[sample][0], s_to_bc_loc[sample][1])[0]
    )
    df_muts["muts_wo_bc"] = df_muts.muts.apply(
        lambda x: get_bc_from_muts(x, s_to_bc_loc[sample][0], s_to_bc_loc[sample][1])[1]
    )

    # df_muts["bc"] = df_muts.muts.apply(lambda x: get_bc_from_muts(x)[0])
    # df_muts["muts_wo_bc"] = df_muts.muts.apply(lambda x: get_bc_from_muts(x)[1])

    df_muts = df_muts.dropna()
    # collecting consecutive deletions into a single position, flag wild-type,
    # get rid of position 0 isnertions

    df_muts["muts_clean"] = df_muts.apply(lambda r: clean_muts(r.muts_wo_bc), axis=1)

    df_muts.to_csv(args.df_out_name, index=False)
    print("**************************************")