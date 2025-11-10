Contextualized_metabolic_models

Python pipeline for contextualizing and sampling E. coli metabolic models using transcriptomic data. Integrates Riptide for expression-guided optimization and COBRApy for flux sampling, enabling parallel analyses and saving condition-specific models and results for comparative metabolic studies.

Overview

This repository provides a complete computational pipeline for contextualizing and sampling metabolic models of E. coli based on transcriptomic data.

The workflow integrates:

Riptide for gene expression–guided model contextualization

COBRApy for flux sampling of context-specific models

Parallel processing for high-throughput execution across multiple experimental conditions

The goal is to generate condition-specific models (raw and normalized) and extract their metabolic flux distributions for downstream comparative analysis.

Workflow Summary
1. Data Preparation

Inputs:

Experimental metadata (41467_2019_13483_MOESM4_ESM.xlsx)

Gene expression matrix (log₂ TPM) in the same file

Gene copy number / multiplicity data (raw_molteplicity.xlsx)

The pipeline filters valid samples based on:

Positive growth rate

E. coli grown in M9 minimal medium

Specific carbon sources (e.g., acetate or others)

2. Contextual Model Optimization (Riptide)

Each filtered condition is used to generate a contextualized model with Riptide.

The script:

Builds a condition-specific model from iML1515

Sets up the M9 medium composition and adjusts carbon uptake to match experimental growth rate

Runs riptide.maxfit() with the corresponding transcriptomic profile

Automatically kills any process exceeding a time limit (default: 1500 s)

Outputs are stored as SBML models under:

raw_exp_maxfit_riptide_runs_iML1515/
normalized_exp_maxfit_riptide_runs_iML1515/

3. Flux Sampling (COBRApy)

For each contextualized model (raw and normalized), the pipeline performs:

Flux sampling using OptGPSampler (1000 samples, thinning = 100)

Parallel execution for raw and normalized models of the same condition

Results are saved as:

cobrapy_samples/<condition>/
│
├── raw_flux_sampling.csv
└── normalized_flux_sampling.csv

4. Integration and Analysis

Each run can be associated with metadata such as:

Growth rate

Condition name

Additional metrics (e.g., Tao times, distances between raw and normalized flux states)

Repository Structure
├── riptide_pipeline.py              # Step 1-2: Riptide contextualization with multiprocessing
├── flux_sampling_pipeline.py        # Step 3: Flux sampling for contextualized models
├── iML1515.xml                      # E. coli base model (SBML format)
├── raw_molteplicity.xlsx            # Multiplicity and Tao data
├── 41467_2019_13483_MOESM4_ESM.xlsx # Metadata + expression dataset
├── cobrapy_samples/                 # Output directory for flux sampling results
├── raw_exp_maxfit_riptide_runs_iML1515/      # Output: raw-contextualized models
└── normalized_exp_maxfit_riptide_runs_iML1515/ # Output: normalized-contextualized models

Installation & Requirements
1. Create a Python environment
python3 -m venv ecoli_env
source ecoli_env/bin/activate

2. Install dependencies
pip install cobra riptide pandas numpy openpyxl tqdm gurobipy


Note:
Gurobi requires a valid license (academic or commercial) to be used as solver in COBRApy.
You can alternatively switch to optlang.glpk_interface if you don't have one.

3. Prepare input data

Ensure the following files are present in the root directory:

iML1515.xml

41467_2019_13483_MOESM4_ESM.xlsx

raw_molteplicity.xlsx

Running the Pipeline
Step 1 — Generate contextualized models
python3 riptide_pipeline.py


This step:

Loads and filters metadata

Runs Riptide optimizations in parallel (configurable threads)

Saves optimized SBML models for each condition

Step 2 — Perform flux sampling
python3 flux_sampling_pipeline.py


This step:

Loads the Riptide-optimized models (raw + normalized)

Performs flux sampling using COBRApy’s OptGPSampler

Saves sampled flux distributions for each condition

Parallelization Control

Both scripts are designed for parallel execution:

Parameter	Script	Description	Default
threads_number	riptide_pipeline.py	Number of Riptide processes to run simultaneously	2
timeout_value	riptide_pipeline.py	Max seconds before a Riptide process is killed	1500
processes	flux_sampling_pipeline.py	Number of parallel sampling jobs (raw + normalized per condition)	2
Output Overview
File	Description
raw_flux_sampling.csv	Flux distribution matrix from raw contextual model
normalized_flux_sampling.csv	Flux distribution matrix from normalized contextual model
terminated.txt	List of Riptide processes that were terminated due to timeout
report.csv (optional)	Summary of per-condition metrics (growth, tao, flux distances)
Notes & Tips

The filtering functions (filter_dataset) automatically merge biological replicates by averaging log₂ TPM values.

Expression values are converted back from log₂ TPM to TPM before being passed to Riptide.

If the number of conditions is large, consider increasing threads_number and processes depending on your CPU core availability.

To visualize sampling results, use tools such as:

import pandas as pd
df = pd.read_csv('cobrapy_samples/glucose/raw_flux_sampling.csv')
df.describe()

Example: Typical Workflow
# 1. Generate Riptide contextual models
python3 riptide_pipeline.py

# 2. Run flux sampling
python3 flux_sampling_pipeline.py

# 3. Inspect outputs
ls cobrapy_samples/


Output:

acetate/
glucose/
...


Each folder contains the flux sampling CSVs for raw and normalized contextual models.

Citation

If you use this pipeline, please cite:

Jenior, M.L. et al., “Transcriptome-guided parsimonious flux analysis improves context-specific metabolic models”, Nature Communications (2019).

Ebrahim, A. et al., “COBRApy: COnstraints-Based Reconstruction and Analysis for Python”, BMC Systems Biology (2013).
