import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from typing import Optional
import pandas as pd
import numpy as np
import nucleotides as nt


def scatter_with_marginals(
    x,
    y,
    x1=[],
    y1=[],
    title="scatter",
    save_path=None,
    fig_size=(8, 8),
    xlabel="lrr_rep1 Beta",
    ylabel="lrr_rep2 Beta",
    wt_col="dodgerblue",
    xlim=None,
    ylim=None,
    n_bins=40,
    n_bins_y=None,
    s=0.1,
    align_text_bottom_right=True,
    fontsize=10,
    tick_linewidth=0.5,
):
    from matplotlib import gridspec
    from matplotlib.ticker import MaxNLocator

    # Set up the axes with gridspec
    fig = plt.figure(figsize=fig_size)
    gs = gridspec.GridSpec(4, 4, wspace=0.4, hspace=0.4)  # Decrease spacing between subplots
    ax_main = plt.subplot(gs[1:4, 0:3])
    ax_xhist = plt.subplot(gs[0, 0:3], sharex=ax_main)
    ax_yhist = plt.subplot(gs[1:4, 3], sharey=ax_main)

    # Scatterplot
    ax_main.scatter(x, y, alpha=0.5, color="black", s=s)
    if len(x1) != 0 and len(y1) != 0:
        ax_main.scatter(x1, y1, color=wt_col, s=3)
    ax_main.set_xlabel(xlabel, fontsize=fontsize)
    ax_main.set_ylabel(ylabel, fontsize=fontsize)
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    ax_main.set_xlim([min_val - 0.1, max_val + 0.1])
    ax_main.set_ylim([min_val - 0.1, max_val + 0.1])
    for spine in ax_main.spines.values():
        spine.set_linewidth(tick_linewidth)  # Decrease axis linewidth
    print(title)

    # Set tick label font sizes for main plot
    ax_main.tick_params(axis="both", which="major", labelsize=fontsize, width=tick_linewidth)

    # Pearson correlation
    r, _ = pearsonr(x, y)
    r_text = f"r: {r:.2f}\nn: {len(x)}"
    # Add textbox in the top left corner of the main scatter plot

    if align_text_bottom_right:
        verticalalignment = "bottom"
        horizontalalignment = "right"
        x_pos = 0.97
        y_pos = 0.02
    else:
        x_pos = 0.05
        y_pos = 0.95
        verticalalignment = "top"
        horizontalalignment = "left"

    ax_main.text(
        x_pos,
        y_pos,
        r_text,
        transform=ax_main.transAxes,
        fontsize=fontsize,
        verticalalignment=verticalalignment,
        horizontalalignment=horizontalalignment,
    )

    # Marginal histograms
    n, bins, patches = ax_xhist.hist(
        x, bins=n_bins, color="black", alpha=0.7, density=True, log=True
    )
    if len(x1) != 0:
        ax_xhist.hist(x1, bins=bins, color=wt_col, alpha=0.7, density=True, log=True)

    if n_bins_y is None:
        n_bins_y = n_bins

    n, bins, patches = ax_yhist.hist(
        y, bins=n_bins, orientation="horizontal", color="black", alpha=0.7, density=True, log=True
    )
    if len(y1) != 0:
        if n_bins_y is None:
            n_bins_y = n_bins
        ax_yhist.hist(
            y1,
            bins=n_bins_y,
            color=wt_col,
            orientation="horizontal",
            alpha=0.7,
            density=True,
            log=True,
        )

    # Set exactly one tick for each marginal histogram

    # ax_xhist.yaxis.set_major_locator(MaxNLocator(nbins=1))
    # ax_yhist.xaxis.set_major_locator(MaxNLocator(nbins=1))

    # Set tick label font sizes for marginal plots
    ax_xhist.tick_params(axis="both", which="major", labelsize=fontsize, width=tick_linewidth)
    ax_yhist.tick_params(axis="both", which="major", labelsize=fontsize, width=tick_linewidth)

    # Turn off frequency axis spines for marginal plots
    ax_xhist.spines["left"].set_visible(
        False
    )  # Turn off y-axis (frequency) spine for top histogram
    ax_xhist.spines["right"].set_visible(False)
    ax_xhist.spines["top"].set_visible(False)
    ax_xhist.spines["bottom"].set_linewidth(tick_linewidth)  # Keep shared x-axis

    ax_yhist.spines["bottom"].set_visible(
        False
    )  # Turn off x-axis (frequency) spine for right histogram
    ax_yhist.spines["top"].set_visible(False)
    ax_yhist.spines["right"].set_visible(False)
    ax_yhist.spines["left"].set_linewidth(tick_linewidth)  # Keep shared y-axis

    # Format the tick labels to be more readable
    ax_xhist.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    ax_yhist.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    ax_xhist.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax_yhist.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    ax_xhist.set_yticks([])
    ax_yhist.set_xticks([])

    # Hide tick labels on marginals
    plt.setp(ax_xhist.get_xticklabels(), visible=False)
    plt.setp(ax_yhist.get_yticklabels(), visible=False)

    # Remove axis labels on marginals
    ax_xhist.set_ylabel("")
    ax_yhist.set_xlabel("")
    if xlim != None:
        ax_xhist.set_xlim(xlim)
    if ylim != None:
        ax_yhist.set_ylim(ylim)

    # Adjust layout with minimal padding
    plt.tight_layout(pad=0.1)  # Decrease padding around the entire figure
    if save_path:
        plt.savefig(save_path, dpi=500, bbox_inches="tight")
        plt.show()
        plt.close()
    else:
        plt.show()


def plot_num_muts(df, saveout=None):
    fig = plt.figure(figsize=(2, 1))
    ax = plt.gca()
    ax.hist(
        df.num_muts,
        bins=max(df.num_muts) + 1,
        color="black",
        alpha=0.9,
        log=True,
        density=True,
        histtype="step",
    )
    plt.xlabel("num. of mutations\n per variant")
    plt.ylabel("frequency")
    plt.xlim([0, 20])
    ax.set_xticks([0, 5, 10, 15, 20])
    plt.setp(ax.get_xticklabels(), rotation=45)

    if saveout:
        plt.savefig(saveout, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close()
    else:
        plt.show()


def plot_corr_marginal(
    x,
    y,
    x2=[],
    y2=[],
    plot_margin=0.04,
    figsize=(2, 2),
    s=0.3,
    c="black",
    alpha=1,
    ticksize=7,
    fout=None,
    plot_log_hist=False,
    plot_diag=True,
    diag_alpha=0.5,
    diag_col="orange",
    # xticks=[0, 200],
    # yticks=[0, 200],
    xlim=None,
    ylim=None,
    xlabel=None,
    ylabel=None,
    plt_title=None,
):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 5])

    ax2 = fig.add_subplot(gs[2])
    ax2.scatter(x, y, s=s, alpha=alpha, c="black")
    if len(x2) != 0 and len(y2) != 0:
        ax2.scatter(x2, y2, s=s, alpha=alpha, c="red")
    """
    if plot_diag:
        ax2.plot(
            [min(xticks) - plot_margin, max(xticks) + plot_margin],
            [min(xticks) - plot_margin, max(xticks) + plot_margin],
            color=diag_col,
            alpha=diag_alpha,
        )
    ax2.set_xticks(xticks)
    ax2.set_yticks(yticks)
    ax2.set_xticklabels(xticks, fontsize=ticksize)
    ax2.set_yticklabels(yticks, fontsize=ticksize)
    """
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)

    if ylabel != None:
        ax2.set_ylabel(ylabel)
    if xlabel != None:
        ax2.set_xlabel(xlabel)

    ax4 = fig.add_subplot(gs[0], sharex=ax2)
    ax4.hist(x, color=c, density=True, alpha=0.5, log=plot_log_hist)
    if len(x2) != 0 and len(y2) != 0:
        ax4.hist(x2, color="red", density=True, alpha=0.5, log=plot_log_hist)
    ax4.axis("off")

    # plot the distance marginal on the right
    ax1 = fig.add_subplot(gs[3], sharey=ax2)
    ax1.hist(y, color=c, density=True, alpha=0.5, orientation="horizontal", log=plot_log_hist)
    if len(x2) != 0 and len(y2) != 0:
        ax1.hist(
            y2, color="red", density=True, alpha=0.5, orientation="horizontal", log=plot_log_hist
        )
    ax1.axis("off")

    fig.patch.set_visible(False)
    ax1.patch.set_visible(False)
    ax2.patch.set_visible(False)
    ax4.patch.set_visible(False)

    # try not to cutoff xticklabels
    fig.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0)  # has to come after fig.tight_layout

    if plt_title != None:
        plt.title(plt_title)

    if fout != None:
        plt.savefig(fout + ".svg", format="svg")
        plt.savefig(fout + ".pdf", format="pdf")
        plt.savefig(fout + ".png", format="png", dpi=800)
    # print("var explained {}".format(get_var_explain(x, y)))
    # print("mse {}".format(get_mse(x, y)))
    print("pearsonr {}".format(pearsonr(x, y)[0]))
    print("spearmanr {}".format(spearmanr(x, y)[0]))
    plt.show()


def plot_reproduce(
    df,
    df_wt: Optional[pd.DataFrame] = None,
    df_mut: Optional[pd.DataFrame] = None,
    mut_label=None,
    df_del: Optional[pd.DataFrame] = None,
    del_label=None,
    fout_plot=None,
    bins=100,
    hist_ylim=[0, 2],
    figsize=(6, 3),
    legend_fs=8,
):

    fig, axs = plt.subplots(1, 2, figsize=figsize)
    axs[0].scatter(df.lrr_rep1, df.lrr_rep2, label="all", color="grey", s=0.5)
    if isinstance(df_mut, pd.DataFrame):
        axs[0].scatter(df_mut.lrr_rep1, df_mut.lrr_rep2, label=mut_label, color="pink", s=3)
    if isinstance(df_del, pd.DataFrame):
        axs[0].scatter(df_del.lrr_rep1, df_del.lrr_rep2, label=del_label, color="blue", s=1)
    if isinstance(df_wt, pd.DataFrame):
        axs[0].scatter(df_wt.lrr_rep1, df_wt.lrr_rep2, label="wt", color="red", s=1)
    axs[0].legend(
        loc="upper left",
        fontsize=legend_fs,
        handlelength=2,
        handletextpad=0.5,
        borderpad=0.5,
        markerscale=0.8,
    )
    axs[0].set_xlabel("lrr replicate 1")
    axs[0].set_ylabel("lrr replicate 2")
    print("pearson:", pearsonr(df.lrr_rep1, df.lrr_rep2))
    print("spearmann:", spearmanr(df.lrr_rep1, df.lrr_rep2))

    n, bins, patches = axs[1].hist(
        df.mean_lrr, bins=bins, label="all", color="grey", density=True, alpha=0.5
    )
    if isinstance(df_mut, pd.DataFrame):
        axs[1].hist(
            df_mut.mean_lrr, bins=bins, label=mut_label, color="pink", density=True, alpha=0.5
        )
    if isinstance(df_del, pd.DataFrame):
        axs[1].hist(
            df_del.mean_lrr, bins=bins, label=del_label, color="blue", density=True, alpha=0.5
        )
    if isinstance(df_wt, pd.DataFrame):
        axs[1].hist(df_wt.mean_lrr, bins=bins, label="wt", color="red", density=True, alpha=0.5)
    axs[1].set_ylim(hist_ylim)
    axs[1].legend(
        loc="upper left",
        fontsize=legend_fs,
        handlelength=2,
        handletextpad=0.5,
        borderpad=0.5,
        markerscale=0.8,
    )
    axs[1].set_xlabel("mean lrr")

    if fout_plot != None:
        plt.savefig(fout_plot + ".png", dpi=300)
    plt.show()

    """
    # plot the scatter plot
    plt.figure(figsize=(2,2))
    plt.scatter(df.lrr_rep1, df.lrr_rep2, label='all', color='grey', s=0.5)
    
    plt.scatter(df_mut.lrr_rep1, df_mut.lrr_rep2, label=mut_label, color='pink', s=3)
    plt.scatter(df_del.lrr_rep1, df_del.lrr_rep2, label=f'>{n_del} deletion', color='blue', s=1)
    plt.scatter(df_wt.lrr_rep1, df_wt.lrr_rep2, label='wt', color='red', s=1)
    plt.legend(loc='upper left', fontsize = 3,handlelength=2, handletextpad=0.5, borderpad=0.5, markerscale=0.8)
    plt.savefig(fout_plot + '.png', dpi=300)
    plt.savefig(fout_plot + '.svg')
    plt.show()
    # plot the histogram
    plt.figure(figsize=(2,2))
    bins=bins
    n, bins, patches = plt.hist(df.mean_lrr, bins=bins,label='all', color='grey', density=True, alpha=0.5)
    plt.hist(df_mut.mean_lrr, bins=bins, label=mut_label, color='pink', density=True, alpha=0.5)
    plt.hist(df_del.mean_lrr, bins=bins, label=f'>{n_del} deletion', color='blue', density=True, alpha=0.5)
    plt.hist(df_wt.mean_lrr, bins=bins, label='wt', color='red', density=True, alpha=0.5)
    plt.ylim(hist_ylim)
    plt.legend(loc='upper left', fontsize = 3,handlelength=2, handletextpad=0.5, borderpad=0.5, markerscale=0.8)
    plt.savefig(fout_plot + '_mean_lrr.png', dpi=300, bbox_inches='tight')
    plt.savefig(fout_plot + '_mean_lrr.svg', bbox_inches='tight')
    plt.show()
    """


"""
def plot_reproducibility(df_high,df_high_wt, plot_title, label_2nd_df = 'wt', figsize=(3,3)):
    plt.figure(figsize=figsize)
    plt.scatter(df_high.lrr_rep1, df_high.lrr_rep2, s=0.5, label='all')
    plt.scatter(df_high_wt.lrr_rep1, df_high_wt.lrr_rep2, s=0.5, c='red', label=label_2nd_df)
    plt.title(plot_title)
    plt.legend(loc='upper left', fontsize = 'small',handlelength=2, handletextpad=0.5, borderpad=0.5, markerscale=0.8)
    plt.xlabel('lrr replicate 1')
    plt.ylabel('lrr replicate 2')
    plt.show()
"""


def plot_multiple_muts(
    list_muts,
    df_high_reads,
    sample_n,
    wt_col="muts_clean",
    mut_col="muts_clean_shifted_oi",
    size_mut=5,
    plotout=None,
):

    # same as above but for multiple mutants.
    n_muts = len(list_muts)
    if n_muts > 0:
        nrows = int(np.ceil(np.sqrt(n_muts)))
        ncols = nrows

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3 * ncols, 3 * nrows))
        if n_muts > 1:
            axes = axes.flatten()
        df_wt = df_high_reads.loc[df_high_reads[wt_col] == "wt"]
        print("found wt reads:", len(df_wt))
        df_plot = df_high_reads.dropna(subset=[mut_col])  # so no nan error when selecting
        # Iterate through the list of mutations and plot each in a subplot
        for i, m in enumerate(list_muts):
            print(m)
            df_with_mut = nt.get_df_with_muts(df_high_reads, m, mut_col=mut_col, strict=True)
            if n_muts > 1:
                ax = axes[i]
            else:
                ax = axes
            ax.scatter(df_plot.lrr_rep1, df_plot.lrr_rep2, s=0.2, color="grey")  # , label='all')
            ax.scatter(df_wt.lrr_rep1, df_wt.lrr_rep2, s=0.4, c="red")  # , label='wt')
            ax.scatter(df_with_mut.lrr_rep1, df_with_mut.lrr_rep2, s=size_mut, c="pink", label=m)

            ax.set_title(f"{sample_n}: {m}, {len(df_with_mut)}/{len(df_high_reads)}", fontsize=4)
            ax.legend(
                loc="lower left",
                fontsize=7,
                handlelength=2,
                handletextpad=0.5,
                borderpad=0.5,
                markerscale=0.8,
            )
            # ax.set_xlabel('lrr replicate 1')
            # ax.set_ylabel('lrr replicate 2')

        # Remove any unused subplots
        if n_muts > 1:
            for j in range(i + 1, len(axes)):
                fig.delaxes(axes[j])

        # plt.tight_layout(pad=1, w_pad=1, h_pad=1.0)
        if plotout != None:
            plt.savefig(plotout + ".png", dpi=400)
            plt.savefig(plotout + ".svg")
        plt.show()


def plot_cooccuring_lr(df_top_hit, df, sample_n, mut_col="muts_oprc", plot_out=None):
    # get all the individual mutants that co-occur in any barcode with the mutant of interest
    count_muts = nt.count_ind_muts(df_top_hit, mut_col=mut_col)

    # plot how many other mutations occur with this top hit mutations
    plt.figure(figsize=(2, 2))
    plt.hist(count_muts.values())
    plt.axvline(len(df_top_hit) / 2, color="red")
    plt.xlabel("# mut co-occurs")
    if plot_out != None:
        plt.savefig(plot_out + "_cooccur.png", dpi=300)
        plt.savefig(plot_out + "_cooccur.svg")
    plt.show()

    # select the co-occuring mutations that occur at least half of the time
    count_muts_co_occur = {k: v for k, v in count_muts.items() if v >= len(df_top_hit) / 2}
    print(
        f"top hit co-occurs with {len(count_muts)} other mutations, {len(count_muts_co_occur)} of which are found in at least half of the occurences"
    )

    if plot_out:
        plot_multi_out = plot_out + "_cooccur"
    else:
        plot_multi_out = None
    plot_multiple_muts(
        list(count_muts_co_occur.keys()),
        df,
        sample_n,
        wt_col=mut_col,
        mut_col=mut_col,
        size_mut=5,
        plotout=plot_multi_out,
    )


def plot_position_count(
    df,
    plot_title,
    mut_col="muts_clean_shifted_oi",
    high_lrr_threshold=2,
    ylim=[0, 0.01],
    xlim=[0, 2100],
    bins=100,
):

    df = df.loc[df[mut_col] != "wt"]

    df["positions"] = df[mut_col].apply(lambda x: [t[1] for t in pacbio.get_ind_muts(x)])

    # flatten the list of positions
    all_mut_locs = df.positions.tolist()
    all_mut_locs = [int(x) for xs in all_mut_locs for x in xs if x != np.nan]

    high_lrr_locs = df.loc[df.mean_lrr > high_lrr_threshold].positions.tolist()
    high_lrr_locs = [int(x) for xs in high_lrr_locs for x in xs if x != np.nan]

    plt.figure(figsize=(20, 5))
    n, bins, patches = plt.hist(all_mut_locs, bins=bins, density=True, label="all")
    plt.hist(high_lrr_locs, bins=bins, alpha=0.5, density=True, label=f"lrr >{high_lrr_threshold}")
    plt.legend(loc="upper center")
    plt.title(plot_title)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel("""position""")
    plt.ylabel("density")
    plt.show()


def plot_distribution_num_muts(
    df, plot_title, num_mut_col="num_muts", bins=1000, plotout=None, xlim=[0, 10]
):
    # plotting number of mutations
    plt.figure(figsize=(2, 2))
    plt.hist(df[num_mut_col], bins=1000)
    plt.xlim(xlim)
    plt.xlabel("number of mutations")
    plt.title(plot_title)
    if plotout != None:
        plt.savefig(plotout + ".png", dpi=300, bbox_inches="tight")
        plt.savefig(plotout + ".svg", bbox_inches="tight")
    plt.tight_layout()
    plt.show()


def plot_dist_mean_lrr(df, df_wt, plot_title, col_plot="mean_lrr", bins=100):
    plt.figure(figsize=(3, 3))
    n, bins, patches = plt.hist(df[col_plot], bins=bins, density=True, label="all", alpha=0.5)
    plt.hist(df_wt.mean_lrr, bins=bins, color="red", density=True, alpha=0.5, label="wt")
    plt.title(plot_title)
    plt.show()


def plot_mutation_type_distribution(
    gene_data_dict,
    mutation_type_col="mutation_type",
    gene_names=None,
    colors=None,
    figsize=(8, 6),
    fontsize=12,
    title="Mutation Type Distribution by Gene",
    ylabel="Proportion",
    save_path=None,
):
    """
    Create a normalized stacked barplot showing the distribution of mutation types for genes.

    Parameters:
    -----------
    gene_data_dict : dict
        Dictionary with gene names as keys and DataFrames as values.
        Each DataFrame should contain a column with mutation types.
    mutation_type_col : str
        Column name containing mutation types ('insertion', 'deletion', 'substitution')
    gene_names : list, optional
        List of gene names to use as x-axis labels. If None, uses dictionary keys.
    colors : dict, optional
        Dictionary mapping mutation types to colors. If None, uses default colors.
    figsize : tuple
        Figure size (width, height)
    fontsize : int
        Font size for all text elements
    title : str
        Plot title
    ylabel : str
        Y-axis label
    save_path : str, optional
        Path to save the plot

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects

    Example:
    --------
    gene_data = {
        'Gene1': df1,
        'Gene2': df2,
        'Gene3': df3
    }
    fig, ax = plot_mutation_type_distribution(gene_data, mutation_type_col='mut_type')
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    # Default colors for mutation types
    if colors is None:
        colors = {
            "insertion": "#FF6B6B",  # Red
            "deletion": "#4ECDC4",  # Teal
            "substitution": "#45B7D1",  # Blue
        }

    # Use provided gene names or dictionary keys
    if gene_names is None:
        gene_names = list(gene_data_dict.keys())

    # Ensure we have exactly 3 genes
    if len(gene_names) != 3:
        raise ValueError("This function is designed for exactly 3 genes")

    # Calculate mutation type counts for each gene
    mutation_counts = {}
    mutation_types = ["insertion", "deletion", "substitution"]

    for gene in gene_names:
        if gene not in gene_data_dict:
            raise ValueError(f"Gene '{gene}' not found in gene_data_dict")

        df = gene_data_dict[gene]
        if mutation_type_col not in df.columns:
            raise ValueError(
                f"Column '{mutation_type_col}' not found in DataFrame for gene '{gene}'"
            )

        # Count each mutation type
        counts = df[mutation_type_col].value_counts()
        gene_counts = {}

        for mut_type in mutation_types:
            gene_counts[mut_type] = counts.get(mut_type, 0)

        mutation_counts[gene] = gene_counts

    # Convert to normalized proportions
    normalized_data = {}
    for gene in gene_names:
        total = sum(mutation_counts[gene].values())
        if total > 0:
            normalized_data[gene] = {
                mut_type: count / total for mut_type, count in mutation_counts[gene].items()
            }
        else:
            normalized_data[gene] = {mut_type: 0 for mut_type in mutation_types}

    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for stacked bar plot
    x_pos = np.arange(len(gene_names))
    bottoms = np.zeros(len(gene_names))

    # Plot each mutation type as a stack
    for mut_type in mutation_types:
        values = [normalized_data[gene][mut_type] for gene in gene_names]
        ax.bar(
            x_pos,
            values,
            bottom=bottoms,
            label=mut_type.capitalize(),
            color=colors[mut_type],
            alpha=0.8,
        )
        bottoms += values

    # Customize the plot
    ax.set_xlabel("Gene", fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize + 2)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(gene_names, fontsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)

    # Add legend
    ax.legend(fontsize=fontsize, loc="upper right")

    # Set y-axis to show proportions (0 to 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])

    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis="y")

    # Tight layout
    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax
