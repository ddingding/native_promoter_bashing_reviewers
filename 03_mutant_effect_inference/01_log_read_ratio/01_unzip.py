# gunzip all fastq.gz files in subdirectories

import pandas as pd
import os
from os import listdir
import importlib

from glob import glob

# most illumina files
# illumina_dir = '/data/davidding/dms/illumina_data/PROTOPLASTSUMMER24-425876467/BCLConvert_07_28_2024_02_06_09Z-753277525/'

# for missing 2A samples of rice and sorghum
illumina_dir = "/data/davidding/dms/illumina_data/missing2a/"

# 0 gunzip files
print(illumina_dir)
gz_files = [f for f in glob(illumina_dir + "**/*.gz", recursive=True) if f.endswith(".fastq.gz")]
for f in gz_files:
    cmd = f"gunzip {f}"
    os.system(cmd)
