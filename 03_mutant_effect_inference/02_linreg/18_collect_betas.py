import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib

import statsmodels.api as sm
import sys

sys.path.append("/data/davidding/dms/dms_plants/src/")
import nucleotides as nt
import linreg as lr
import plotting as pt

plt.style.use("/data/davidding/dms/dms_plants/src/paper_style1.mplstyle")


from matplotlib import font_manager
font_path = "/data/davidding/dms/dms_plants/src/Arial.ttf"  # Your font path goes here
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)


din_lr = "/data/davidding/dms/illumina_data/linreg_all/"

din_df = "/data/davidding/dms/illumina_data/data/clean/illumina/lrr_rcs/"

percentile = 99

list_samples = ["2A_rice",
                "2B_sorghum",
                "2C_sorghum",
                "2C_rice",
                "3A_sorghum",
                "3B_sorghum",
                "3C_sorghum"]

for sample_n in list_samples:
    print("running", sample_n)
    for fit_col in ["mean_lrr", "lrr_rep1", "lrr_rep2"]:
        # fetch the data
        X_b, y, df, res, _ = lr.load_data(
            din_lr + sample_n + f"_{fit_col}/" + f"_sm_ols_{fit_col}.pickle",  # path to model
            din_df + sample_n + ".csv",  # path to dataframe
            din_lr + sample_n + f"_{fit_col}/" + sample_n + f"_{fit_col}.npy",  # path to X_b
        )

        # get the possible mutations
        counter_muts = nt.count_ind_muts(df, mut_col="mut_tss0")
        poss_muts = sorted(counter_muts.keys())

        # plot the distribution
        importlib.reload(nt)
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
            fout_plot=din_lr + sample_n + f"_{fit_col}/" + "_enhancer_del",
        )

        # extract the coefficients
        df_coef, bonf = lr.extract_betas(res, poss_muts, counter_muts, sample_n)

        # add whether the mean_lrr is in the top 0.5% percentile fo the distribution
        df_coef = lr.add_is_top_percentile(df_coef, df, percentile=percentile, mut_col="mut_tss0")

        df_coef.to_csv(din_lr + sample_n + f"_{fit_col}/" + "_df_coef_0.csv", index=False)

        # plot the volcano
        lr.plot_volcano(df_coef, bonf, din_lr + sample_n + f"_{fit_col}/" + "_volcano")
        lr.check_linreg(res, y, fout=din_lr + sample_n + f"_{fit_col}/" + "_check_linreg")
