# for merging protoplast summer 24 reads

from glob import glob
import pandas as pd
import os

# parent directory
# illumina_dir = '/data/davidding/dms/illumina_data/PROTOPLASTSUMMER24-425876467/BCLConvert_07_28_2024_02_06_09Z-753277525/'

# for missing 2A samples
illumina_dir = "/data/davidding/dms/illumina_data/missing2a/"

# max number of bases that can overlap
max_n_base_overlap = 150


"""
# read in the metadata
df_meta = pd.read_csv("./metadata/eg_naming_pp_summer_24.csv")
df_meta["f_name_prefix"] = df_meta.SAMPLE_SUBMISSION_NAME.apply(lambda x: "EG_" + x.lower()[3:])
df_meta["new_name"] = df_meta.apply(lambda r: r.SAMPLE.lower() + "_" + r.SPECIES.lower(), axis=1)

for all the samples
# find all fastq_files
fastq_files = [
    f
    for f in glob(illumina_dir + "**/*.fastq", recursive=True)
    if f.endswith(".fastq") and "EG_protoplast" in f
]
"""
# for 2A
fastq_files = [f for f in glob(illumina_dir + "**/*.fastq", recursive=True) if f.endswith(".fastq")]
for f in fastq_files:
    print(f)
print(f"found {len(fastq_files)} fastq files")

# find the pairs of reads corresponding to each other.
# creating a dictionary of absolute direcotry with the prefix of the file
dir_pref_to_f_list = {}
for f in sorted(fastq_files):
    pref = f[: -len("R1_001.fastq")]
    if pref not in dir_pref_to_f_list:
        dir_pref_to_f_list[pref] = [f]
    else:
        dir_pref_to_f_list[pref].append(f)
print("dir_pref_to_f_list", dir_pref_to_f_list)


# flash merge fastq files
for dir_pref, f_list in dir_pref_to_f_list.items():

    """
    # only for the original, non2A samples since they are named differently
    # finding the human readable sample description and use as the new name
    f_name_prefix = "_".join(dir_pref.split("/")[-1].split("_")[:3])
    new_name = df_meta.loc[df_meta.f_name_prefix == f_name_prefix].new_name.values[0]
    """
    # for 2A
    new_name = dir_pref.split("/")[-1]

    if len(f_list) != 2:
        print(f"error: more than 2 files found for the same prefix {f_name_prefix}, {f_list}")
        continue
    else:
        # merging into the same directory as the input fastq file
        output_dir = "/".join(dir_pref.split("/")[:-1]) + "/"
        print("output_dir", output_dir)
        flash_cmd = f"/home/davidding/apps/FLASH-1.2.11/flash {f_list[0]} {f_list[1]} -M {max_n_base_overlap} -o {new_name} -d {output_dir}"
        print(flash_cmd)
        os.system(flash_cmd)
