Contextualized_metabolic_models

Python pipeline for contextualizing and sampling E. coli metabolic models using transcriptomic data. Integrates Riptide for expression-guided optimization and COBRApy for flux sampling, enabling parallel analyses and saving condition-specific models and results for comparative metabolic studies.

Overview

This repository provides a complete computational pipeline for contextualizing and sampling metabolic models of E. coli based on transcriptomic data.

The workflow integrates:

- Riptide for gene expression–guided model contextualization
- COBRApy for flux sampling of context-specific models
- Parallel processing for high-throughput execution across multiple experimental conditions

The goal is to generate condition-specific models (raw and normalized) and extract their metabolic flux distributions for downstream comparative analysis.

Workflow Summary
1. Data Preparation

Inputs:

- Experimental metadata (41467_2019_13483_MOESM4_ESM.xlsx)
- Gene expression matrix (log₂ TPM) in the same file
- Gene copy number / multiplicity data (raw_molteplicity.xlsx)

The pipeline filters valid samples based on:

- Positive growth rate
- E. coli grown in M9 minimal medium
- Specific carbon sources (e.g., acetate or others)

2. Contextual Model Optimization (Riptide)

Each filtered condition is used to generate a contextualized model with Riptide.

The script:

- Builds a condition-specific model from iML1515
- Sets up the M9 medium composition and adjusts carbon uptake to match experimental growth rate
- Runs riptide.maxfit() with the corresponding transcriptomic profile
- Automatically kills any process exceeding a time limit (default: 1500 s)

Outputs are stored as SBML models under:

raw_exp_maxfit_riptide_runs_iML1515/
normalized_exp_maxfit_riptide_runs_iML1515/

3. Flux Sampling (COBRApy)

For each contextualized model (raw and normalized), the pipeline performs:

- Flux sampling using OptGPSampler (1000 samples, thinning = 100)
- Parallel execution for raw and normalized models of the same condition

Results are saved as:

cobrapy_samples/<condition>/
├── raw_flux_sampling.csv
└── normalized_flux_sampling.csv

4. Integration and Analysis

Each run can be associated with metadata such as:

- Growth rate
- Condition name
- Additional metrics (e.g., Tao times, distances between raw and normalized flux states)

Repository Structure
├── riptide_pipeline.py
├── flux_sampling_pipeline.py
├── iML1515.xml
├── raw_molteplicity.xlsx
├── 41467_2019_13483_MOESM4_ESM.xlsx
├── cobrapy_samples/
├── raw_exp_maxfit_riptide_runs_iML1515/
└── normalized_exp_maxfit_riptide_runs_iML1515/
├── plots/  # includes example visualizations etc.
