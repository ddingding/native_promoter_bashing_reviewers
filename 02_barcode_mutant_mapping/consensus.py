import argparse
import os
import tempfile
import warnings
import Bio.SeqIO
from Bio import SeqIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
import os, io, random
import pandas as pd
import numpy as np
import pysam
from pysam import bcftools
import math
import time

TEMP_DIR = "/data/szchen/dms/snakemake/consensus-calling/temp-bc-bams"


def sam_to_bam(samfile, bamfile_path, skip_repeat=False):
    """Converts sam to a bamfile and returns index object for fast read access"""
    # Check if the BAM file already exists
    if skip_repeat:
        if os.path.exists(bamfile_path):
            print(f"{bamfile_path} already exists. Skipping conversion.")
    else:
        # Convert to bam
        pysam.sort("-o", bamfile_path, samfile)

    # Index bamfile
    alignmentBAMFile = bamfile_path
    bamfile = pysam.AlignmentFile(alignmentBAMFile)
    indexedBAMFile = pysam.IndexedReads(bamfile)
    indexedBAMFile.build()
    return indexedBAMFile


def bam_to_consensus(bamfile, bamindex, unique_bc_df, fqout, temp_dir):
    """
    produces consensus sequences from a bamfile and a list of unique barcodes
    fqout: file to write consensus sequence for each barcode
    """
    startTime = time.time()
    sample = os.path.basename(bamfile).split("_")[1]
    temp_dir = os.path.join(temp_dir, sample)
    with open(fqout, "w") as f:
        template = pysam.AlignmentFile(bamfile)

        for i, row in unique_bc_df.iterrows():
            if (i % 100000 == 0) or (i == 200) or (i == 20) or (i == 1000):
                print(
                    str(100 * i / unique_bc_df.shape[0])
                    + "% finished in "
                    + str(time.time() - startTime)
                    + " seconds"
                )
            BC = row["bc"]
            read_ids = row["read_ids"].split(":")
            # temp_bam_unsort = f"{temp_dir}/{BC}_{sample}.bam"
            temp_bam_unsort = os.path.join(temp_dir, f"{BC}_{sample}.bam")
            # temp_bam_sort = f"{temp_dir}tempReadssorted_{sample}.bam"
            temp_bam_sort = os.path.join(temp_dir, f"{BC}_{sample}_sort.bam")
            BCReads = pysam.AlignmentFile(temp_bam_unsort, "wb", template=template)
            for read in read_ids:
                iterator = bamindex.find(read)
                for x in iterator:  # x is an alignment
                    BCReads.write(x)
            BCReads.close()
            pysam.sort("-o", temp_bam_sort, temp_bam_unsort)

            consensus = ""
            f.write(">" + BC + "\n")
            try:
                #         consensus = pysam.consensus('-l 10000',
                # #                                     '--config hifi',
                #                                     'tempReadssorted.bam')
                consensus = os.popen(
                    f"samtools consensus -f fastq -l 100000 --config hifi {temp_bam_sort}"
                ).read()
            except Exception as e:
                print(e)
            f.write(consensus.split()[1] + "\n")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generates consensus sequences for each barcode"
    )
    parser.add_argument(
        "-s", "--samfile", type=str, required=True, help="path to input .sam file"
    )
    parser.add_argument(
        "-u",
        "--unique_bc_csv",
        type=str,
        required=True,
        help="path to CSV mapping unique barcodes to reads",
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="name of output FASTQ file"
    )
    return parser.parse_args()


def main():
    # Snakemake usage: 'python consensus.py --samfile {input.samfile} --unique_bc_csv {input.unique_bcs} --o {output}'
    startTime = time.time()

    # load arguments
    args = parse_args()
    samfile = args.samfile
    fqout = args.output
    # reads_bc_df = pd.read_csv(args.reads_csv)
    unique_bc_df = pd.read_csv(args.unique_bc_csv)
    # filter df of unique bcs to bcs with more than one read
    unique_bc_df_filt = unique_bc_df[unique_bc_df["read_ids"].str.contains(":")]
    # remove XXXXX unmatched barcode
    unique_bc_df_filt = unique_bc_df_filt[
        ~unique_bc_df_filt["bc"].str.contains("XXXXXX")
    ]
    unique_bc_df_filt.reset_index(drop=True, inplace=True)
    bamfile_path = f"{samfile.split('.')[0]}.bam"
    bamindex = sam_to_bam(
        samfile, bamfile_path=bamfile_path
    )  # create indexed bamfile from samfile

    bam_to_consensus(
        bamfile=bamfile_path,
        bamindex=bamindex,
        unique_bc_df=unique_bc_df_filt,
        fqout=fqout,
        temp_dir=TEMP_DIR,
    )

    # report elapsed time
    seconds = time.time() - startTime
    mins, secs = divmod(seconds, 60)
    hrs, mins = divmod(mins, 60)
    print(f"processed {samfile} in {hrs} hours {mins} minutes {secs} seconds")


if __name__ == "__main__":
    main()
