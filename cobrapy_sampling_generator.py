#!/usr/bin/env python3
"""
Flux Sampling Generator for Contextualized E. coli Models
----------------------------------------------------------
This script performs flux sampling on context-specific models generated
via Riptide (both raw and normalized expression-based models).
Each condition is sampled in parallel for raw and normalized datasets,
and results are saved in dedicated output folders.

Author: (Your Name)
"""

import os
import pandas as pd
import cobra
from multiprocessing import Pool


# -------------------------------------------------------------
# Save flux sampling results to CSV files
# -------------------------------------------------------------
def save_flux_sampling(raw_sampling, normalized_sampling, output_dir, condition):
    """
    Saves raw and normalized flux sampling matrices into a dedicated directory.
    """
    output_path = os.path.join(output_dir, condition)
    os.makedirs(output_path, exist_ok=True)

    print(f'[INFO] Saving flux sampling results for: {condition}')

    raw_sampling.to_csv(os.path.join(output_path, 'raw_flux_sampling.csv'))
    normalized_sampling.to_csv(os.path.join(output_path, 'normalized_flux_sampling.csv'))


# -------------------------------------------------------------
# Generate flux sampling for a given condition and dataset path
# -------------------------------------------------------------
def generate_sampling(condition, dataset_path):
    """
    Loads the SBML model corresponding to the given condition and performs
    flux sampling using the OptGPSampler from COBRApy.
    """
    model_path = None

    # Locate the folder containing the model for the condition
    for folder_name in os.listdir(dataset_path):
        if condition in folder_name:
            model_path = os.path.join(dataset_path, folder_name, 'model.sbml')
            break

    if model_path is None or not os.path.exists(model_path):
        raise FileNotFoundError(f'[ERROR] Could not find model for condition: {condition}')

    # Load the SBML model
    ecoli_model = cobra.io.read_sbml_model(model_path)

    # Perform flux sampling (1000 samples, thinning factor 100)
    sampler = cobra.sampling.OptGPSampler(ecoli_model, thinning=100, processes=1)
    flux_samples = sampler.sample(1000)

    return flux_samples


# -------------------------------------------------------------
# Filter metadata for M9 medium and specific carbon sources
# -------------------------------------------------------------
def filter_dataset():
    """
    Reads metadata and TAO values, filters metadata for valid M9 medium
    conditions with positive growth rates, and keeps only selected carbon sources.
    """
    tao_df = pd.read_excel('raw_molteplicity.xlsx', sheet_name='tao', index_col='conditions')
    metadata_df = pd.read_excel('41467_2019_13483_MOESM4_ESM.xlsx', sheet_name='Metadata', index_col='Sample ID')

    # Define carbon sources to include (can be extended as needed)
    carbon_sources = ['acetate']
    '''
    carbon_sources = [
        'glucose', 'fumarate', 'pyruvate', 'fructose', 'glycerol', 'xylose',
        'sorbitol', 'D-ribose', 'glucarate', 'N-acetylglucosamine',
        'galactose', 'gluconate', 'D-lyxose', 'D-arabinose', 'm-tartrate'
    ]
    '''

    # Filter metadata
    metadata_df = metadata_df[
        (metadata_df['Growth Rate (1/hr)'] > 0) &
        (metadata_df['Base Media'] == 'M9') &
        (metadata_df['Carbon Source (g/L)'].str.extract(r'([A-Za-z0-9\-]+)', expand=False).isin(carbon_sources))
    ]

    # Keep only first replicate per condition
    metadata_df = metadata_df.drop([cond for cond in metadata_df.index if int(cond.split('_')[-1]) > 1])
    metadata_df.index = [idx.rsplit('_', 1)[0].rstrip('_') for idx in metadata_df.index]

    return metadata_df, tao_df


# -------------------------------------------------------------
# Main execution
# -------------------------------------------------------------
if __name__ == '__main__':
    # Define paths and initialize counter
    output_dir = 'cobrapy_samples'
    os.makedirs(output_dir, exist_ok=True)
    completed_count = 0

    # Paths to Riptide-optimized model directories
    dataset_paths = [
        '/home/mason/Desktop/Riptide_optimization/new_dataset/raw_exp_maxfit_dataset/raw_exp_maxfit_riptide_runs_iML1515/',
        '/home/mason/Desktop/Riptide_optimization/new_dataset/normalized_exp_maxfit_dataset/normalized_exp_maxfit_riptide_runs_iML1515/'
    ]

    # Load filtered metadata and TAO sheet
    metadata_df, tao_df = filter_dataset()

    # Initialize report dataframe
    report_df = pd.DataFrame(index=metadata_df.index,
                             columns=['raw--normalized distance', 'Growth Rate (1/hr)', 'tao(min)'])

    # Iterate through each condition in metadata
    for condition in metadata_df.index:
        raw_path, normalized_path = dataset_paths

        # Prepare arguments for parallel processing
        tasks = [(condition, raw_path), (condition, normalized_path)]

        print(f'[INFO] Starting flux sampling for condition: {condition}')

        # Perform sampling for raw and normalized models in parallel
        with Pool(processes=2) as pool:
            raw_sampling, normalized_sampling = pool.starmap(generate_sampling, tasks)

        completed_count += 1
        print(f'[INFO] Completed sampling for {condition} ({completed_count}/{len(metadata_df.index)})')

        # Save results
        save_flux_sampling(raw_sampling, normalized_sampling, output_dir, condition)
