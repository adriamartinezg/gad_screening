# Genomic landscape and ecological distribution of the GAD system across global biomes
This repository contains the complete bioinformatic framework developed to identify, calibrate, and map the glutamate decarboxylase system (GAD, genes gadA, gadB, and gadC) across diverse microbial biomes. 
It is associated to the Final Master's Thesis in Bioinformatics wrote by Adrià Martínez García.

By combining profile Hidden Markov Models (pHMMs) with custom Trust Thresholds (TCs), this pipeline ensures high-precision discrimination against non-functional paralogs in large-scale datasets (RefSeq, cFMD, and MGnify).

## Repository Structure

```{plaintext}
├── model/                      # Trained .hmm files for gadA, gadB, and gadC
├── markers/                    # Concatenated marker protein sequences (.faa) for phylogenetic alignment
├── thresholds.csv              # Calculated Trusted Cutoffs (TCs) for validation
├── install_dependencies.sh     # Script to install mandatory tools and environment libraries
│
├── calibration/                # Snakemake pipeline for model validation (LOOCV)
│   ├── Snakemake
│   ├── config.yaml
│   ├── negatives.faa           # Decoy/negative control protein sequences
│   └── scripts/
│       ├── clustering.sh
│       ├── loocv_hmmer.sh
│       ├── parser2.sh
│       ├── tmhmm_deeper_loocv.sh
│       └── aggregate_loocv.py
│
├── screening/                  # Core classifier script and database-specific logic
│   ├── classifier.sh           # Main execution script for targeted mining
│   ├── cfmd/
│   │   ├── process_cfmd.sh
│   │   ├── parse_cfmd.sh
│   │   └── taxonomy_cfmd.sh
│   ├── mgnify/
│   │   ├── mgnify.sh
│   │   └── parse_mgnify.sh
│   └── refseq/
│       ├── counter.sh
│       ├── genomes_status.R
│       └── get_taxid.sh
│
└── plotting/                   # Downstream data analysis and R statistical visualization
    ├── refseq_plots.R
    ├── mgnify_plots.R
    ├── cfmd_plots.R
    └── combined_plots.R        # Relative abundance, Fisher's Exact Test, and enrichment plots
```

## Setup & Prerequisites

### Software Dependencies
To set up the initial environment and compile necessary dependencies, run:

```{bash}
bash install_dependencies.sh
```
NOTE: To perform any phylogenetic analysis with GTDB, the database in force must be installed

## Configuration Warning for Screening Scripts
Before running any screening pipeline (screening/), you must manually update the hardcoded file paths and subfolder links inside the bash scripts (`classifier.sh`, database processors, etc.) to match your local infrastructure and folder layout.

## Model calibration
Move the markers directory inside the calibration directory. Then, adjust directories in `config.yaml`and launch the `Snakefile`.

## Metagenomic Screening (/screening)
Houses the main algorithm (`classifier.sh`) optimized for three distinct repositories:

### cFMD (Curated Food Metagenomic Data): 
Processes, parses, and taxonomically maps food-associated MAGs (`process_cfmd`, `parse_cfmd` and `taxonomy_cfmd`).

External Requirement: Users must download `download_mags.sh` and `cFMD_metadata.tsv` directly from [SegataLab Repository](https://github.com/SegataLab/cFMD).

### MGnify: 
Processes and extracts ecological data from diverse environmental and host-associated datasets (`mgnify.sh`, `parse_mgnify.sh`).

### RefSeq: 
Filters, clean-ups, and screens reference genomes using structural filters and taxonomical matching tools (`get_taxid.sh`, `genomes_status.R` and `counter.sh`).

## Statistical Analysis & Plotting (/plotting)
Contains specialized R scripts to ingest parsed outputs from your screening files. These scripts automatically compute:

Relative Abundance Metrics: To map raw occurrences.

Fisher's Exact Tests: Implemented as a bias-correction step to balance out repository sampling skews.

Enrichment Plots: To visually pinpoint niche-specific trends (e.g., sudden pH drop environments like dairy fermentations vs. stable buffers).

NOTE: Get sure to install R packages indicated at the top of the Rscripts.
