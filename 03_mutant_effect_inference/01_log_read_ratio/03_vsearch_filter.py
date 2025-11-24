# 2 filter by quality

from glob import glob
import os

# find all the combined fragment .fastq files

# illumina_dir = '/data/davidding/dms/illumina_data/PROTOPLASTSUMMER24-425876467/BCLConvert_07_28_2024_02_06_09Z-753277525/'
# for missing 2A samples
illumina_dir = "/data/davidding/dms/illumina_data/missing2a/"


proc_fastq_files = [f for f in glob(illumina_dir + "**/*extendedFrags.fastq", recursive=True)]
print(proc_fastq_files)

for f in proc_fastq_files:
    fout = f[: -len(".fastq")] + "_filtered.fasta"
    vsearchCmd = (
        f"/home/davidding/apps/vsearch-2.28.1/bin/vsearch --fastq_filter {f} --fastq_truncqual 20 "
        f"--fastq_maxns 3 --fastq_maxee 0.5 --fastq_ascii 33 --fastaout {fout}"
    )
    print(vsearchCmd)
    os.system(vsearchCmd)
