import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

import plotting as pt
import nucleotides as nt
import constants as co

import statsmodels.api as sm
from scipy.stats import skew
from statsmodels.stats.stattools import durbin_watson, jarque_bera


def get_Xb_lr(df, chpt_fout=None, mut_col="mut_tss0"):
    # fetch the X matrix with a 1 bias term
    counter_muts = nt.count_ind_muts(df, mut_col=mut_col)
    possible_muts = sorted(counter_muts.keys())

    X = ohe_all_muts(df, possible_muts, mut_col=mut_col)
    X_b = np.hstack((X, np.ones((X.shape[0], 1))))
    if chpt_fout != None:
        np.save(chpt_fout, X_b)

    return X_b, possible_muts


def load_data(model_path, df_path, x_path=None, poss_mut_path=None):
    # fetch
    # read the raw mutation file in.
    df = pd.read_csv(df_path)
    y = df["mean_lrr"].values

    # load the X matrix
    if x_path != None:
        X_b = np.load(x_path)
    else:
        X_b, poss_muts = get_Xb_lr(df)

    # load the linear regression model
    res = sm.load(model_path)

    # load_the possible mutations used to fit model
    if poss_mut_path is not None:
        poss_muts = pickle.load(open(poss_mut_path, "rb"))
    else:
        poss_muts = None

    return X_b, y, df, res, poss_muts


def check_linreg(res, y, fout=None):

    # try to verify some assumptions of linear regression:
    # 1. residuals should be normally distributed (looking at jarque bera statistic)
    # 2. residuals should be homoscedastic ()
    # 3. residuals should not have autocorrlation with lag 1 (looking at durbin watson)

    print(skew(res.resid))

    print("durbin watson: [~2]", durbin_watson(res.resid))
    print(
        "jarque bera: [~0]", jarque_bera(res.resid)
    )  # skew of < 0: left skewed, kurtosis >3: heavier tails than a normal

    plt.figure(figsize=(4, 4))
    plt.scatter(res.fittedvalues, res.resid, s=0.1, c="black")
    plt.xlabel("fitted values")
    plt.ylabel("residuals")
    if fout != None:
        plt.savefig(fout + "_residuals.png", format="png", dpi=300)
    plt.show()

    # verify that fit is good
    pt.plot_corr_marginal(
        y, res.fittedvalues, xlabel="True LRR", ylabel="Predicted LRR", figsize=(4, 4), fout=fout
    )


def extract_betas(res, possible_muts, counter_ind_muts, f_name, sig_level=0.05, gene_n=None):
    # get the coefficients and p-values associated with each mutation
    df_coef = pd.DataFrame(
        {"mut": possible_muts, "beta": res.params[:-1], "pval": res.pvalues[:-1]}
    )

    bonferoni = sig_level / len(df_coef)
    df_coef["sig"] = df_coef.pval < bonferoni
    print("found ", df_coef.sig.sum(), " significant mutations")

    print(df_coef.shape)

    # get the count of how often each mutation appears in dataset
    df_count = pd.DataFrame(
        {"mut": list(counter_ind_muts.keys()), "n": list(counter_ind_muts.values())}
    )

    df_coef = df_coef.merge(df_count, on="mut")

    # get rid of mutants that are N's, ie. sequencing artefacts.
    df_coef = df_coef[~df_coef["mut"].str.contains("N")]

    # DEprecated: should already be done in 08_clean_processing. shift mutations
    # shift the positions of the mutations s.t. position 0 is the beginning of the TSS
    # get rid of mutants that are outside of the -2kb promoter and the 5'UTR

    # add a column whether mutation is in promoter or not
    if gene_n is None:
        gene_n = f_name.split("_")[0]
    tss_pos = co.sample_to_positions[gene_n][1] - 1
    len_5utr = co.sample_to_positions[gene_n][-1] - tss_pos + 1
    df_coef["mut_in_prom_utr"] = df_coef.mut.apply(
        lambda m: nt.is_in_promoter_utr(m, -2000, len_5utr)
    )

    # df_coef = df_coef.sort_values(by='beta', ascending=False)
    return df_coef, bonferoni


def extract_betas_reg(res, possible_muts, counter_ind_muts, f_name):
    # get the coefficients and p-values associated with each mutation
    df_coef = pd.DataFrame({"mut": possible_muts, "beta": res.params[:-1]})

    print(df_coef.shape)

    # get the count of how often each mutation appears in dataset
    df_count = pd.DataFrame(
        {"mut": list(counter_ind_muts.keys()), "n": list(counter_ind_muts.values())}
    )

    df_coef = df_coef.merge(df_count, on="mut")

    # get rid of mutants that are N's, ie. sequencing artefacts.
    df_coef = df_coef[~df_coef["mut"].str.contains("N")]

    # DEprecated: should already be done in 08_clean_processing. shift mutations
    # shift the positions of the mutations s.t. position 0 is the beginning of the TSS
    # get rid of mutants that are outside of the -2kb promoter and the 5'UTR

    # add a column whether mutation is in promoter or not
    gene_n = f_name.split("_")[0]
    tss_pos = co.sample_to_positions[gene_n][1] - 1
    len_5utr = co.sample_to_positions[gene_n][-1] - tss_pos + 1
    df_coef["mut_in_prom_utr"] = df_coef.mut.apply(
        lambda m: nt.is_in_promoter_utr(m, -2000, len_5utr)
    )

    # df_coef = df_coef.sort_values(by='beta', ascending=False)
    return df_coef


def find_cutoff(data, percentile=95):
    """
    Finds the cutoff above which a certain percentage of the data lie.

    Parameters:
    - data: list of numbers
    - percentile: the percentile above which the data lies (default is 95 for the top 5%)

    Returns:
    - cutoff: the value above which the top percentile of data lie
    """
    data = list(data)
    if not data:
        raise ValueError("Data list cannot be empty")

    # Calculate the percentile cutoff
    cutoff = np.percentile(data, percentile)

    return cutoff


def add_is_top_percentile(df_coef, df, percentile=99.5, mut_col="mut_tss0"):
    # adds whether the mean effect of a mutation across its observed reads is in the top percentile of the LRR distribution
    cutoff_lrr = find_cutoff(list(df.mean_lrr.values), percentile)
    print("cutoff:", cutoff_lrr)
    plt.figure(figsize=(4, 4))
    plt.hist(df.mean_lrr, bins=100, color="black")
    plt.axvline(cutoff_lrr, color="red")
    plt.xlabel("mean LRR")
    plt.ylabel("count")
    plt.title(f"where the percentile {percentile} cutoff is:")
    plt.show()

    df_coef["mean_lrr_observed_reads"] = df_coef.mut.apply(
        lambda m: np.mean(nt.get_df_with_muts(df, m, mut_col=mut_col, strict=True).mean_lrr)
    )

    df_coef["is_top_p"] = df_coef.mean_lrr_observed_reads > cutoff_lrr

    return df_coef.sort_values(by="beta", ascending=False)


def plot_volcano(df_coef, bonferoni, fout=None):
    plt.figure(figsize=(4, 4))
    plt.scatter(df_coef["beta"], -np.log(df_coef["pval"]), s=0.5, c="black")
    plt.axhline(-np.log(bonferoni), color="red")
    plt.xlabel("beta")
    plt.ylabel("-log(pval)")
    if fout != None:
        plt.savefig(fout + ".png", format="png", dpi=300)
    plt.show()


def inspect_top_muts(df_coef_top, df, f_name, mut_col="mut_tss0", n_del=300):
    # inspect the top positive significant mutations

    for i, r in df_coef_top.iterrows():
        m = r.mut
        # fetch the wt, enhancer and deletion dataframes.
        df_wt, df_mut, df_del = nt.get_df_wt_mut_del(df, mut_col, mut_str=m, n_del=n_del)

        pt.plot_reproduce(
            df,
            df_wt=df_wt,
            df_mut=df_mut,
            mut_label=m,
            df_del=df_del,
            del_label=f">{n_del} del",
            legend_fs=6,
        )

        #### find any co-occuring mutations
        # find all barcodes with this mutation
        df_hit = nt.get_df_with_muts(df, m, mut_col=mut_col, strict=True)

        # print the mutations
        for i, m in enumerate(df_hit[mut_col]):
            print(f"read{i+1}:", m)

        #
        pt.plot_cooccuring_lr(df_hit, df, f_name, mut_col=mut_col)


def ohe_all_muts(df, possible_muts, mut_col="mut_tss0"):
    # one hot encode all mutations from dataframe
    ohe = np.zeros((df.shape[0], len(possible_muts)))
    for i in range(df.shape[0]):
        for mut in df[mut_col][i].split(":"):
            ohe[i, possible_muts.index(mut)] = 1
    return ohe


def get_Xb(df, chpt_fout, mut_col="mut_tss0"):
    # fetch the X matrix with a 1 bias term
    counter_muts = nt.count_ind_muts(df, mut_col=mut_col)
    possible_muts = sorted(counter_muts.keys())

    X = ohe_all_muts(df, possible_muts, mut_col=mut_col)
    X_b = np.hstack((X, np.ones((X.shape[0], 1))))
    np.save(chpt_fout, X_b)

    return X_b, possible_muts
