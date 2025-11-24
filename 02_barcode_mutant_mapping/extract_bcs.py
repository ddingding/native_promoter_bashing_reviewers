import argparse
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq

S_TO_FLANK_5PR_FWD = {
    "2A": "GAACTATACAAATAA",
    "2B": "GAACTATACAAATAA",
    "2C": "GAACTATACAAATAA",
    "3A": "TTACTATATATAATCCAATC",
    "3B": "GTTTTGCTTGGTTTTGGAGA",
    "3C": "TGATGAGTAGTACTGTGCGC",
}
S_TO_FLANK_5PR_REV = {
    "2A": "TTATTTGTATAGTTC",
    "2B": "TTATTTGTATAGTTC",
    "2C": "TTATTTGTATAGTTC",
    "3A": "GATTGGATTATATATAGTAA",
    "3B": "TCTCCAAAACCAAGCAAAAC",
    "3C": "GCGCACAGTACTACTCATCA",
}

S_TO_FLANK_3PR_FWD = {
    "2A": "CTCGAGGTTCGAGTA",
    "2B": "CTCGAGGTTCGAGTA",
    "2C": "CTCGAGGTTCGAGTA",
    "3A": "AGATCGGAAGAGCGTCGTGT",
    "3B": "AGATCGGAAGAGCGTCGTGT",
    "3C": "AGATCGGAAGAGCGTCGTGT",
}
S_TO_FLANK_3PR_REV = {
    "2A": "TACTCGAACCTCGAG",
    "2B": "TACTCGAACCTCGAG",
    "2C": "TACTCGAACCTCGAG",
    "3A": "ACACGACGCTCTTCCGATCT",
    "3B": "ACACGACGCTCTTCCGATCT",
    "3C": "ACACGACGCTCTTCCGATCT",
}


def extract_bcs(fqin, sample, fiveprime_offset=5):
    """Extract barcodes using flanking reference sequences"""
    r_to_bc = {}
    ref_seq_5prime_fwd = S_TO_FLANK_5PR_FWD[sample]
    ref_seq_5prime_rev = S_TO_FLANK_5PR_REV[sample]
    ref_seq_3prime_fwd = S_TO_FLANK_3PR_FWD[sample]
    ref_seq_3prime_rev = S_TO_FLANK_3PR_REV[sample]

    unmatched_bc_str = "XXXXXXXXXXXXXXX"  # filler barcode for reads that dont find a bc
    for r in SeqIO.parse(fqin, "fastq"):
        id = r.id
        seq = r.seq
        # Extract 5' flank
        left_pos_fwd = seq.find(ref_seq_5prime_fwd)
        left_pos_rev = seq.find(ref_seq_5prime_rev)
        # Extract 3' flank
        right_pos_fwd = seq.find(ref_seq_3prime_fwd)
        right_pos_rev = seq.find(ref_seq_3prime_rev)
        if (left_pos_fwd != -1 and right_pos_fwd != -1):
            bc_start = left_pos_fwd + len(ref_seq_5prime_fwd) - fiveprime_offset  # index of first letter of bc minus 5bp of 5' sequence
            bc_end = right_pos_fwd  # +1 of index of last letter of bc
            if (bc_end - bc_start) <= 200:
                bc = seq[bc_start:bc_end]
                r_to_bc[id] = bc
            else:
                r_to_bc[id] = unmatched_bc_str
        elif (left_pos_rev != -1 and right_pos_rev != -1): # if read is reversed
            bc_start = right_pos_rev + len(ref_seq_3prime_rev)
            bc_end = left_pos_rev + fiveprime_offset # add 5 bp of 5' sequence
            if (bc_end - bc_start) <= 200:
                bc = seq[bc_start:bc_end]
                bc = Seq(bc).reverse_complement() # reverse the barcode
                r_to_bc[id] = str(bc)
            else:
                r_to_bc[id] = unmatched_bc_str
        else:
            r_to_bc[id] = unmatched_bc_str
    return pd.DataFrame(r_to_bc.items(), columns=["read_n", "bc"])


def main():
    # Snakemake shell: python extract_bcs.py --input {input} --reads_df {output.readsBCDF} --unique_bc_df {output.uniqueBCDF} --sample {wildcards.sample}
    parser = argparse.ArgumentParser(
        description="extracts 1) [read:barcode mapping] and 2) [unique barcode:reads] mapping --- then writes to csv"
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="path toinput fastq file"
    )
    parser.add_argument(
        "-r",
        "--reads_csv",
        type=str,
        required=True,
        help="path to output csv mapping read:barcode",
    )
    parser.add_argument(
        "-u",
        "--unique_bc_csv",
        type=str,
        required=True,
        help="path to output csv mapping unique bc:reads",
    )

    parser.add_argument(
        "-s",
        "--sample",
        type=str,
        required=True,
        help="name of sample",
    )

    args = parser.parse_args()
    fqin = args.input
    rdfout = args.reads_csv
    udfout = args.unique_bc_csv
    sample = args.sample

    readsBCDF = extract_bcs(fqin=fqin, sample=sample)
    print(f"Length of readsBCDF: {readsBCDF.shape[0]}")
    readsBCDF.to_csv(rdfout, index=False)
    # uniqueBCDF = pd.DataFrame(readsBCDF["bc"].value_counts())
    # # Remove any sequence with fewer than 10 bases in the barcode position
    # mask = uniqueBCDF["bc"].str.len() > 10
    # uniqueBCDF = uniqueBCDF.loc[mask]
    # uniqueBCDF.reset_index(inplace=True)

    uniqueBCDF = (
        readsBCDF.groupby("bc")["read_n"].apply(lambda x: ":".join(x)).reset_index()
    )
    mask = uniqueBCDF["bc"].str.len() > 10
    uniqueBCDF = uniqueBCDF.loc[mask]
    uniqueBCDF.reset_index(drop=True, inplace=True)
    uniqueBCDF.columns = ["bc", "read_ids"]
    print(f"Length of uniqueBCDF: {uniqueBCDF.shape[0]}")
    uniqueBCDF.to_csv(udfout, index=False)


if __name__ == "__main__":
    main()
