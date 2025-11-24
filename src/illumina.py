# tools for illumina analysis
import matplotlib.pyplot as plt
import pacbio
import nucleotides as nt
import plotting as pt
import numpy as np
from collections import Counter
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from constants import dic_sample_to_muts_prevalent


def create_lrr_file(f_dna, f_rna, dout):
    """
    Create log read ratio file by merging DNA and RNA count files and calculating LRR.

    Args:
        n_dna (str): Name of DNA sample
        n_rna (str): Name of RNA sample
        f_to_path (dict): Dictionary mapping sample names to file pathåså
        dout (str): Output directory path
    """

    df_dna = pd.read_csv(f_dna)
    df_rna = pd.read_csv(f_rna)

    # making sure to merge on outer to not miss mutants.
    df_merge = pd.merge(df_dna, df_rna, on="bc", how="outer", suffixes=("_dna", "_rna"))
    print(len(df_dna), len(df_rna), len(df_merge))

    # adding a pseudocount to dna
    df_merge["bc_count_dna"] = df_merge["bc_count_dna"].fillna(0)
    df_merge["bc_count_dna"] = df_merge["bc_count_dna"] + 1

    # adding a pseudocount to rna
    df_merge["bc_count_dna"] = df_merge["bc_count_dna"].fillna(0)
    df_merge["bc_count_dna"] = df_merge["bc_count_dna"] + 1

    # add a log read ratio
    n_reads_dna = df_merge.bc_count_dna.sum()
    n_reads_rna = df_merge.bc_count_rna.sum()
    df_merge["lrr"] = df_merge.apply(
        lambda r: np.log(r.bc_count_rna / n_reads_rna) - np.log(r.bc_count_dna / n_reads_dna),
        axis=1,
    )

    f_name_out = f_rna.split("/")[-1][:-4] + "_lrr.csv"
    df_merge.to_csv(dout + f_name_out, index=False)


def drop_bc_nan_filter_len(df, rep):
    # filtering only for barcodes of the right length of 15 nucleotides.
    len_raw = len(df)
    df = df.dropna(subset=["bc"])
    len_dropna_bc = len(df)
    df["len_bc"] = df.bc.apply(lambda x: len(x))
    df = df.loc[df.len_bc.isin([15, 3])]
    len_bc15 = len(df)
    print(
        f"for rep {rep}: \n droppedn from {len_raw}: \n  after dropna: {len_dropna_bc} \n after bc length filtering: {len_bc15})"
    )
    return df


def merge_replicates(sample, list_reps, dout, use_old_format=True):
    """
    Merge log read ratio data from two replicate samples.

    Args:
        sample: Sample name (e.g. '2B_sorghum')
        list_reps: List of replicate names to merge
        dout: Output directory path

    Returns:
        Merged DataFrame containing data from both replicates
    """
    print(sample, list_reps)

    if use_old_format:
        df1 = pd.read_csv(dout + list_reps[0] + "_lrr.csv")
    else:
        df1 = pd.read_csv(list_reps[0])

    df1 = drop_bc_nan_filter_len(df1, list_reps[0])

    # rename columns to keep track in merged dataframe
    dic_col_rename1 = dict(zip(df1.columns, [c + "_rep1" for c in df1.columns]))
    df1 = df1.rename(columns=dic_col_rename1)

    if use_old_format:
        df2 = pd.read_csv(dout + list_reps[1] + "_lrr.csv")
    else:
        df2 = pd.read_csv(list_reps[1])

    df2 = drop_bc_nan_filter_len(df2, list_reps[1])

    # rename columns to keep track in merged dataframe
    dic_col_rename2 = dict(zip(df2.columns, [c + "_rep2" for c in df2.columns]))
    df2 = df2.rename(columns=dic_col_rename2)

    print(len(df1), len(df2))

    display(df1)
    display(df2)

    df_merged = df1.merge(df2, left_on="bc_rep1", right_on="bc_rep2", how="outer")

    display(df_merged)
    df_merged.to_csv(dout + "reps/" + sample + "_lrr.csv", index=False)
    print(len(df_merged))

    return df_merged


def drop_bc_nan_filter_len(df, rep, list_allowed_bc_lens=[15, 3]):
    # filtering only for barcodes of the right length of 15 nucleotides (or 3 for the 'XXX' barcode), and dropping nan barcodes.
    # Reason:
    # for full length samples, could only define the barcode in illumina by the 5' end, and then taking the remaining 15 nucleotides.
    len_raw = len(df)
    df = df.dropna(subset=["bc"])
    len_dropna_bc = len(df)
    df["len_bc"] = df.bc.apply(lambda x: len(x))
    df = df.loc[df.len_bc.isin(list_allowed_bc_lens)]
    len_bc15 = len(df)
    print(
        f"for rep {rep}: \n droppedn from {len_raw}: \n  after dropna: {len_dropna_bc} \n after bc length filtering: {len_bc15}"
    )
    return df


def get_forbidden_df(
    df,
    sample_n,
    list_forbidden,
    dout="/data/davidding/dms/dms_plants/data/illumina/regression/",
    mut_col="muts_wo_bc",
):
    # filter out reads that have mutations in list_forbidden, and write out

    df = df.dropna(subset=["mean_lrr"])  # drop the ones that don't have a mean_lrr
    df["list_pos_mut"] = df.apply(
        lambda r: [x[1] for x in nt.get_ind_muts(str(r[mut_col]))], axis=1
    )

    # exclude mutations in list_forbidden
    df["mut_forbidden"] = df.apply(
        lambda r: any([x in list_forbidden for x in r.list_pos_mut]), axis=1
    )
    df_no_forbidden = df.loc[df.mut_forbidden == False]
    print(
        f"filtering truseq mutations out filters out this many barcodes {len(df) - len(df_no_forbidden)} out of {len(df)}"
    )

    # write these to file
    df.to_csv(dout + sample_n + ".csv", index=False)
    df_no_forbidden.to_csv(dout + sample_n + "_no_forbidden.csv", index=False)

    return df, df_no_forbidden


def test_correct_prevalent_mutations():
    list_muts_remove = ["T3575A", "C5576T", "A5574G", "T5575A"]
    # test that all mutations are removed and can make wild-types
    assert correct_prevalent_mutations("T3575A:C5576T:A5574G:T5575A", list_muts_remove) == "wt"
    # test that if a mutation from the list of remove is not found, the alternate mutation is introduced.
    assert correct_prevalent_mutations("C5576T:A5574G:T5575A", list_muts_remove) == "A3575T"
    # test that a mutation to a different base at one of the positions to remove is preserved
    assert correct_prevalent_mutations("T3575C:C5576T:A5574G:T5575A", list_muts_remove) == "T3575C"

    assert correct_prevalent_mutations("T2964C:G4018T:G4424A", []) == "T2964C:G4018T:G4424A"
    # test_correct_prevalent_mutations()


def correct_prevalent_mutations(mutstr, list_muts_remove):
    # takes a mut string, like 'A3570G:T3575A:G4128A:G4976_:G4977_'
    # and remove mutations that are found in list_muts_remove, ie. useful when almost all the reads have this particular mutation
    # effectively getting a new reference sequence with these mutations in it.

    # Note: this only works for list_muts_remove that are coded as single base, position and mutated base, ie. 'A3570G', not 'A2570__'

    if str(mutstr).lower() == "nan":
        return np.nan
    if str(mutstr).lower() == "wt":
        return "wt"

    list_muts_in = mutstr.split(":")
    # remove mutations from that list if they are in list_muts_remove
    list_muts_out = [m for m in list_muts_in if m not in list_muts_remove]

    # the new reference should contain these mutations if they are not found in the mutstr
    # search for positions that are not found in the mutstr, but are found in the list_muts_remove, if so, add the reversed mutation to the list (ie. A5574G -> G5574A based on the new ref)
    # the positions that are mutated in the mutstr
    mut_pos_to_mutstr = dict(zip([int(t[1]) for t in nt.get_ind_muts(mutstr)], list_muts_in))

    mut_pos_to_mutstr_remove = dict(
        zip([int(nt.break_mut(t)[1]) for t in list_muts_remove], list_muts_remove)
    )

    #### adding the reverse mutation, ie. from the new reference, for reads that don't have mutations in these positions
    # find the positions to reverse the mutation given a new reference.
    pos_to_reverse = [p for p in mut_pos_to_mutstr_remove.keys() if p not in mut_pos_to_mutstr]

    reversed_muts_to_add = []
    for p in pos_to_reverse:
        orig_mutstr = mut_pos_to_mutstr_remove[p]
        orig_wt, orig_pos, orig_mut = nt.break_mut(orig_mutstr)

        reversed_muts_to_add.append(orig_mut + str(orig_pos) + orig_wt)
        """
        # do not deal with insertions/deletions for now: something like so
        # if the new reference is a deletion at a particular position
        # then the new mutation from that reference should b
        if orig_mut == '_':
            reversed_mut = 'I' + str(orig_pos) + orig_wt

        # if the new reference has an insertion at a particular site, 
        # then flag the reverse mutation as a deletion the length of the insertion
        
        # reversed string
        reversed_mut = ''.join(map(str,[::-1]))
        # if some of these are reversing A4987_ to _4987A, then the underscore should be replaced with an I signifying insertion.
        if reversed_mut[0] == '_':
            reversed_mut = 'I' + reversed_mut[1:]
        
        if reversed_mut[-1] == 'I':
        """

    list_muts_out += reversed_muts_to_add

    # if there are no mutations, return wt
    if len(list_muts_out) == 0:
        return "wt"
    # otherwise, return the list joined by colons
    return ":".join(list_muts_out)


def add_reref_clean_hr_cols(f_name, dout_no_truseq, dout_oprc, gene_n=None, old_format=True):
    # add a column with mutations where a new reference (based on extremely prevalent mutations) is used
    # this column cleaned up, and a column with the number of mutations is added.
    # finally filter for high read barcodes.

    if old_format:
        gene_n = f_name.split("_")[0]

        df = pd.read_csv(dout_no_truseq + f_name + ".csv")
    else:
        df = pd.read_csv(f_name)

    # new column with mutations with the Original Position, new Reference.
    df["muts_opr"] = df.muts_orig_pos.apply(
        lambda x: correct_prevalent_mutations(x, dic_sample_to_muts_prevalent[gene_n])
    )

    # cleaning conseuctive deletions up etc., and remove a prevalent mutation in 3C
    if gene_n != "3C":
        # new column with mutations with the Original Position, new Reference, Cleaned, ie. consecutive deletions are merged together.
        df["muts_oprc"] = df.muts_opr.apply(lambda x: nt.clean_muts(x))

    # for 3C: remove the I1592AA mutations as they are found next to a 10bp poly-A tract, also checked that linear regression doesn't pick this up: it's +/- 0.2 in coefficient.
    elif gene_n == "3C":

        def remove_mut(mutstr, mut_remove="I1592A"):
            removed_mut_list = [m for m in mutstr.split(":") if not m.startswith(mut_remove)]
            if len(removed_mut_list) == 0:
                return "wt"
            else:
                return ":".join(removed_mut_list)

        df["muts_opr_no_i1592a"] = df.muts_opr.apply(lambda x: remove_mut(x, mut_remove="I1592A"))
        df["muts_oprc"] = df.muts_opr_no_i1592a.apply(lambda x: nt.clean_muts(x))

    # add mutcount
    df["num_muts"] = df.muts_oprc.apply(
        lambda x: len(x.split(":")) if str(x).lower() != "wt" else 0
    )

    # fetch only bcs with high reads and seen in both RNA samples at least once.
    df_h = nt.fetch_high_reads(df, min_reads=20, min_reads_rna=1)

    if old_format:
        df_h.to_csv(dout_oprc + f_name + ".csv", index=False)
    else:
        df_h.to_csv(dout_oprc + f_name.split("/")[-1], index=False)

    print(len(df_h))
    print("# of wt:", len(df.loc[df.num_muts == 0]))

    pt.plot_distribution_num_muts(
        df_h, plot_title=f_name, num_mut_col="num_muts", bins=1000, plotout=None, xlim=[0, 10]
    )


def calc_df_enrichment(df, mean_lrr_cutoff, mut_col="muts_wo_bc"):
    # sort individual mutations based on their enrichment in a quadrant,
    # as well as calculate their mean and median lrr values in reads that contain them

    # split into high mean_lrr and low mean_lrr
    df_high_lrr = df.loc[df.mean_lrr > mean_lrr_cutoff]
    df_low_lrr = df.loc[df.mean_lrr <= mean_lrr_cutoff]

    # print(len(df_high_lrr)+ len(df_low_lrr), len(df))

    # counting in the high lrr df
    counter_ind_muts_high_lrr = nt.count_ind_muts(df_high_lrr, mut_col=mut_col)
    df_ind_muts_high_lrr = pd.DataFrame(counter_ind_muts_high_lrr.items(), columns=["mut", "count"])

    # counting mutation occurence in the low lrr df
    counter_ind_muts_low_lrr = nt.count_ind_muts(df_low_lrr, mut_col=mut_col)
    df_ind_muts_low_lrr = pd.DataFrame(counter_ind_muts_low_lrr.items(), columns=["mut", "count"])

    # sanity to count int the total df
    counter_ind_muts_total = nt.count_ind_muts(df, mut_col=mut_col)
    df_ind_muts_total = pd.DataFrame(counter_ind_muts_total.items(), columns=["mut", "count"])

    # merging the dataframes together
    df_high_low = df_ind_muts_high_lrr.merge(
        df_ind_muts_low_lrr, on="mut", how="outer", suffixes=("_high_lrr", "_low_lrr")
    ).fillna(0)
    df_all = df_high_low.merge(df_ind_muts_total, on="mut", how="outer").fillna(0)
    # df_all['sum_low_high'] = df_all.count_high_lrr + df_all.count_low_lrr

    total_c_high_lrr = sum(df_all.count_high_lrr)
    total_c_low_lrr = sum(df_all.count_low_lrr)
    print(total_c_high_lrr, total_c_low_lrr)
    df_all["frac_high_lrr"] = (df_all.count_high_lrr + 1) / total_c_high_lrr
    df_all["frac_low_lrr"] = (df_all.count_low_lrr + 1) / total_c_low_lrr
    df_all["enrichment"] = df_all.frac_high_lrr / df_all.frac_low_lrr

    # adding the mutation position
    df_all["pos"] = df_all.mut.apply(lambda x: nt.get_ind_muts(x)[0][1])

    # adding the mean and median lrr values
    df_all["mean_mean_lrr"] = df_all.mut.apply(
        lambda x: np.mean(nt.get_df_with_muts(df, x, mut_col=mut_col).mean_lrr)
    )
    df_all["median_lrr"] = df_all.mut.apply(
        lambda x: np.median(nt.get_df_with_muts(df, x, mut_col=mut_col).mean_lrr)
    )

    df_all = df_all.sort_values(by="enrichment", ascending=False)
    return df_all


def look_into_cooccurence(df_top_hit, df, df_enrich_min_count, sample_n, plot_out=None):
    print("LegacyError: This function is now called plot_cooccuring_lr in plotting...")
    # get all the individual mutants that co-occur in any barcode with the mutant of interest
    # for dealing with df_enrich
    count_muts = count_ind_muts(df_top_hit)

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

    print("their enrichment values are")
    display(df_enrich_min_count.loc[df_enrich_min_count.mut.isin(list(count_muts_co_occur.keys()))])

    plot_multiple_muts(
        list(count_muts_co_occur.keys()),
        df,
        sample_n,
        mut_col="muts_wo_bc",
        size_mut=1,
        plotout=plot_out + "_cooccur",
    )
