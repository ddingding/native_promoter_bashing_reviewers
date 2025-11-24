import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import statsmodels.api as sm
import sys
import os

sys.path.append("/data/davidding/dms/dms_plants/src/")
import nucleotides as nt
import plotting as pt
import linreg as lr
from constants import sample_to_positions

import pickle
from os import listdir

plt.style.use("/data/davidding/dms/dms_plants/src/paper_style1.mplstyle")

from matplotlib import font_manager

font_path = "/data/davidding/dms/dms_plants/src/Arial.ttf"  # Your font path goes here
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)


# input data for most samples
# din_df = "/data/davidding/dms/illumina_data/data/clean/illumina/lrr_rcs/"
# din_lr = "/data/davidding/dms/illumina_data/linreg_all/" # output data

# to process illumina run 3, which has a specific focus on the 3A samples
din_df = "/data/davidding/dms/illumina_data/data/clean/illumina_run3/lrr_rcs/"
din_lr = "/data/davidding/dms/illumina_data/data/clean/illumina_run3/linreg/"

percentile = 99


def ohe_all_muts(df, possible_muts, mut_col="muts_oprc"):
    # one hot encode all mutations from dataframe
    ohe = np.zeros((df.shape[0], len(possible_muts)))
    for i in range(df.shape[0]):
        for mut in df[mut_col].values[i].split(":"):
            ohe[i, possible_muts.index(mut)] = 1
    return ohe


def is_of_interest(mutstr, len_5utr=200):
    if mutstr == "wt":
        return True
    # only returns mutations if all mutations fall in the -2kb and 5'utr region, or contain the very prevalent -2179 mutation that is widely distributed
    ind_muts = nt.get_ind_muts(mutstr)
    p_to_mutstr = dict(zip([m[1] for m in ind_muts], ["".join(map(str, m)) for m in ind_muts]))
    m_to_keep = [p for p in p_to_mutstr if (p < len_5utr and p > -2000) or p == -2179]

    if len(m_to_keep) == len(p_to_mutstr):
        return True
    else:
        return False


def get_Xb(df, chpt_fout, mut_col="muts_oprc"):

    counter_muts = nt.count_ind_muts(df, mut_col=mut_col)
    possible_muts = sorted(counter_muts.keys())

    X = ohe_all_muts(df, possible_muts, mut_col=mut_col)
    X_b = np.hstack((X, np.ones((X.shape[0], 1))))
    np.save(chpt_fout, X_b)

    return X_b, possible_muts


def clean_design_matrix(X_b, threshold=1e-6):
    """Clean design matrix by removing problematic columns"""
    # Remove constant columns
    non_constant = ~np.all(X_b == X_b[0, :], axis=0)

    # Check correlation matrix
    corr = np.corrcoef(X_b[:, non_constant].T)

    # Find eigenvalues to identify linear dependencies
    eigenvals = np.linalg.eigvals(corr)
    good_cols = np.abs(eigenvals) > threshold

    # Get column indices to keep
    cols_to_keep = np.where(non_constant)[0][good_cols]
    # reattach the bias term
    cols_to_keep = list(cols_to_keep) + [X_b.shape[1] - 1]
    print(f"Removed {X_b.shape[1] - len(cols_to_keep)} columns due to linear dependencies")

    return X_b[:, cols_to_keep], cols_to_keep


def fit_linreg(y, X_b_clean, fout="./test.pickle"):
    print(y.shape)
    # using statsmodels to fit the model
    print("fitting model ...")
    try:
        res = sm.OLS(y, X_b_clean).fit()
    except np.linalg.LinAlgError:
        print("falling back to QR decomposition")
        # If still failing, try with more robust solver
        res = sm.OLS(y, X_b_clean).fit(method="qr")

    res.save(fout)
    return res


def fit_collect_betas(df, f_name, fit_col, folder_name=None, gene_n=None):
    if folder_name is None:
        folder_name = f_name + f"_{fit_col}"
    else:
        folder_name = folder_name

    print("fitting model for", fit_col)
    # lin_reg for mean_lrr
    os.makedirs(din_lr + folder_name, exist_ok=True)
    # get the X_b
    X_b, poss_muts = get_Xb(
        df, chpt_fout=din_lr + folder_name + f_name + f"_{fit_col}.npy", mut_col="mut_tss0"
    )

    X_b_clean = X_b

    pickle.dump(
        poss_muts,
        open(din_lr + folder_name + f"_poss_muts_{fit_col}.pickle", "wb"),
    )
    print(X_b_clean.shape)

    y = df[fit_col].values

    res = fit_linreg(y, X_b_clean, fout=din_lr + folder_name + f"_sm_ols_{fit_col}.pickle")

    sample_n = f_name
    # collect betas.
    # fetch the data

    X_b = X_b_clean

    # get the possible mutations
    counter_muts = nt.count_ind_muts(df, mut_col="mut_tss0")
    poss_muts = sorted(counter_muts.keys())

    # plot the distribution
    df_wt, df_mut, df_del = nt.get_df_wt_mut_del(
        df,
        "mut_tss0",
        mut_str="TCGAGCAGCTGGCTT",
        mut_strict=False,  # ie. any string that contains this sequence
        n_del=300,
    )

    pt.plot_reproduce(
        df,
        df_wt=df_wt,
        df_mut=df_mut,
        mut_label="enhancer",
        df_del=df_del,
        del_label=">300 del",
        fout_plot=din_lr + folder_name + "_enhancer_del",
    )

    # extract the coefficients
    df_coef, bonf = lr.extract_betas(res, poss_muts, counter_muts, sample_n, gene_n=gene_n)

    # add whether the mean_lrr is in the top 0.5% percentile fo the distribution
    df_coef = lr.add_is_top_percentile(df_coef, df, percentile=percentile, mut_col="mut_tss0")

    df_coef.to_csv(din_lr + folder_name + "_df_coef_0.csv", index=False)

    # plot the volcano
    lr.plot_volcano(df_coef, bonf, din_lr + folder_name + "_volcano")
    lr.check_linreg(res, y, fout=din_lr + folder_name + "_check_linreg")


i = 0
for f_name in listdir(din_df):
    f_name = f_name.split(".")[0]
    print("fitting mode for", f_name, "iteration", i)

    # hack to pass the gene name on
    if "run3" in din_df:
        gene_n = f_name.split("_")[0][-2:]
    else:
        gene_n = None
    print("gene_n:", gene_n)

    df = pd.read_csv(din_df + f_name + ".csv")
    print(f"shape of df: {df.shape}")

    # running into a statsmodels init gesdd failed problem.
    # solution is to filter out mutations outside of the promoter and 5'UTR region, then it works.
    if "3A" in f_name:

        len_5utr = sample_to_positions["3A"][-1] - sample_to_positions["3A"][-2]

        df["mut_oi"] = df["mut_tss0"].apply(
            lambda mutstr: is_of_interest(mutstr, len_5utr=len_5utr)
        )
        df = df[df["mut_oi"] == True]

    print(f"shape of df after filtering: {df.shape}")

    for fit_col in ["mean_lrr", "lrr_rep1", "lrr_rep2"]:
        fit_collect_betas(df, f_name, fit_col, folder_name=f_name + f"_{fit_col}/", gene_n=gene_n)
