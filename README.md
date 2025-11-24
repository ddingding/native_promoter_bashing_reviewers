# Native Promoter Bashing Analysis

Analysis code and data for plant promoter bashing experiments for three sorghum genes: PsbS, Raf1, and SBPase.

## Overview

This repository contains the complete computational pipeline for analyzing promoter mutation libraries, from library design to figure generation. The analysis includes:

- Library design and oligonucleotide ordering
- Barcode-mutant mapping from sequencing data
- Statistical inference of mutational effects
- Figure generation for publication

3. **External tools required:**
   - [minimap2](https://github.com/lh3/minimap2) for sequence alignment
   - [vsearch](https://github.com/torognes/vsearch) for sequence processing
   - [FLASH](https://ccb.jhu.edu/software/FLASH/) for read merging

## Repository Structure

```
├── 01_library_creation/          # Library design and oligonucleotide generation
├── 02_barcode_mutant_mapping/    # Process sequencing data to map barcodes to mutants
├── 03_mutant_effect_inference/   # Statistical analysis of mutational effects, starting with Illumina high-throughput counts.
├── 04_figures/                   # Generate figures
├── data/                         # Input data and reference files
├── src/                          # Shared utility functions and constants
└── README.md
```

## Analysis Pipeline

Note that the filepaths in all notebooks correspond to our compute environment and will have to respecified by you if you rerun this analysis.


### 1. Library Creation (`01_library_creation/`)
- Design mutational libraries for PsbS, Raf1, and SBPase promoters
- Generate oligonucleotides for synthesis and subcloning
- **Main notebook:** `01_order_library.ipynb`

### 2. Barcode-Mutant Mapping (`02_barcode_mutant_mapping/`)
- Process PacBio sequencing data
- Map barcodes to specific mutations
- **Pipeline:** Snakemake workflow in `snakefile`

### 3. Effect Inference (`03_mutant_effect_inference/`)
- Calculate log read ratios from Illumina data for each mutant: `01_log_read_ratio/`
- Fit linear regression models to infer individual mutational effects: `02_linreg/`

### 4. Figure Generation (`04_figures/`)
- Generate all publication figures

## File Naming Conventions

### Sample Codes
- **2A, 2B, 2C**: Synthetic GFP reporter constructs
- **3A, 3B, 3C**: Native gene constructs
- **A**: SBPase gene
- **B**: Raf1 gene  
- **C**: PsbS gene

## Mutation Nomenclature

- `wt`: Wildtype (no mutation)
- `I{pos}{seq}`: Insertion at position `pos` of sequence `seq`
- `{wt_base}{pos}{each deletion denoted by a single _}`: Deletion starting at position `pos` of length `len`
- `{wt_base}{pos}{nt}`: Substitution at position `pos` to from `wt_base` nucleotide `nt`

Position numbering:
- Negative positions: upstream of transcription start site (TSS)
- Position 0: transcription start site
- Positive positions: downstream of TSS

## Key Results

The analysis identifies:
- Critical regulatory regions in plant promoters
- Quantitative effects of individual mutations
- Reproducible patterns across replicates
- Comparison with computational predictions

## Data

Raw sequencing data and processed results are stored in the `data/` directory, organized by:
- `library_order/`: Synthesized library sequences
- `library_results/`: Processed effect measurements
   - the mut column refers to the mutation relative to the transcriptional start site, and the relevant column for the inferred effects is beta (in log2 values).
- `ref/`: Reference sequences for the various wild-type plasmids
- `nanoluc/`: nanoluciferase measurements in rice and sorghum for specific mutations.
- and other various datasets used in the paper


## Citation

If you use this code, please cite our paper.

## Contact

For questions about the analysis or code, please contact davidding [at] berkeley [dot] edu
