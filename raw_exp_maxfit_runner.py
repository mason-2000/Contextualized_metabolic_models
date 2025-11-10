#!/usr/bin/env python3
"""
Parallel Riptide MaxFit Runner for E. coli iML1515 model
---------------------------------------------------------
This script reads transcriptomic and metadata information from Excel files,
filters and averages replicates, sets up the iML1515 model environment for
each condition, and runs Riptide's `maxfit()` optimization in parallel.
Each run is terminated if it exceeds a specified time limit.
"""

import pandas as pd
import numpy as np
import multiprocessing as mp
import cobra
import riptide
import time
import math


# -------------------------------------------------------------
# Utility function to initialize the "terminated.txt" log file
# -------------------------------------------------------------
def initialize_termination_log():
    with open('terminated.txt', 'w') as f:
        f.write('')  # clear file contents
    return open('terminated.txt', 'a')  # return open file handle for appending


# -------------------------------------------------------------
# Riptide worker: runs a single optimization and saves output
# -------------------------------------------------------------
def riptide_worker(metadata_df, expression_dict, sample_id,
                   frac_min=0.5, frac_max=0.95, threshold=math.exp(-8),
                   samples=1000, gpr=True):
    """
    Runs Riptide maxfit optimization for one sample.
    Saves the optimized model in SBML format.
    """
    print(f'[INFO] Optimizing model for sample: {sample_id}...')

    # Build condition-specific E. coli model
    ecoli_model = build_ecoli_model(metadata_df, sample_id)

    # Run Riptide optimization
    result = riptide.maxfit(
        model=ecoli_model,
        transcriptome=expression_dict,
        frac_min=frac_min,
        frac_max=frac_max,
        gpr=gpr,
        samples=samples,
        frac_step=0.05
    )

    # Compare objective function before and after optimization
    initial_obj = ecoli_model.optimize().objective_value
    final_obj = result.model.optimize().objective_value
    print(f'[INFO] Objective value changed from {initial_obj:.4f} to {final_obj:.4f}')

    # Save the resulting model
    riptide.save_output(
        riptide_obj=result,
        path=f'raw_exp_maxfit_riptide_runs_iML1515/{sample_id}',
        file_type='sbml'
    )


# -------------------------------------------------------------
# Runs multiple Riptide processes in parallel with timeout
# -------------------------------------------------------------
def run_riptide_batch(metadata_df, expression_df, timeout_sec=1500, n_threads=2):
    """
    Launches Riptide runs in parallel using multiprocessing.
    Each batch runs up to n_threads samples at a time.
    Kills any process that exceeds the time limit.
    """
    terminated_log = initialize_termination_log()

    start_idx = 0
    end_idx = n_threads

    while start_idx < len(expression_df.columns):
        active_processes = []

        # Launch processes for the current batch
        for col_name in expression_df.columns[start_idx:end_idx]:
            sample_id = col_name
            expression_dict = expression_df[sample_id].to_dict()

            process = mp.Process(
                target=riptide_worker,
                args=(metadata_df, expression_dict, sample_id),
                name=sample_id
            )
            process.start()
            active_processes.append(process)

        batch_start_time = time.monotonic()

        # Monitor processes and terminate if timeout exceeded
        for proc in active_processes:
            elapsed = time.monotonic() - batch_start_time
            remaining_time = max(0, timeout_sec - elapsed)

            proc.join(timeout=remaining_time)

            if proc.is_alive():
                proc.terminate()
                proc.join()
                print(f'[WARNING] Terminated hung process for sample {proc.name}')
                terminated_log.write(f'KILLED: {proc.name}\n')

        start_idx += n_threads
        end_idx += n_threads


# -------------------------------------------------------------
# Builds E. coli model from metadata for a given condition
# -------------------------------------------------------------
def build_ecoli_model(metadata_df, sample_id):
    """
    Loads the iML1515 model, sets the growth medium and carbon source uptake,
    and adjusts the uptake rate to match the observed growth rate.
    """
    ecoli_model = cobra.io.read_sbml_model('iML1515.xml')

    # Extract carbon source name
    carbon_source = metadata_df.loc[sample_id, 'Carbon Source (g/L)'].split('(')[0]

    # Mapping from carbon source to exchange reaction ID
    carbon_exchanges = {
        'glucose': 'EX_glc__D_e',
        'fumarate': 'EX_fum_e',
        'pyruvate': 'EX_pyr_e',
        'fructose': 'EX_fru_e',
        'glycerol': 'EX_glyc_e',
        'xylose': 'EX_xyl__D_e',
        'sorbitol': 'EX_sbt__D_e',
        'D-ribose': 'EX_rib__D_e',
        'glucarate': 'EX_glcr_e',
        'N-acetylglucosamine': 'EX_acgam_e',
        'galactose': 'EX_gal_e',
        'gluconate': 'EX_glcn_e',
        'D-lyxose': 'EX_lyx__L_e',
        'D-arabinose': 'EX_arab__L_e',
        'm-tartrate': 'EX_tartr__L_e',
        'acetate': 'EX_ac_e'
    }

    # Define M9 medium composition
    m9_medium = {
        carbon_exchanges[carbon_source]: 4.0,
        'EX_cobalt2_e': 1000.0, 'EX_h_e': 1000.0, 'EX_h2o_e': 1000.0,
        'EX_k_e': 1000.0, 'EX_cu2_e': 1000.0, 'EX_mg2_e': 1000.0,
        'EX_mn2_e': 1000.0, 'EX_mobd_e': 1000.0, 'EX_na1_e': 1000.0,
        'EX_nh4_e': 1000.0, 'EX_ca2_e': 1000.0, 'EX_ni2_e': 1000.0,
        'EX_o2_e': 1000.0, 'EX_cl_e': 1000.0, 'EX_pi_e': 1000.0,
        'EX_zn2_e': 1000.0, 'EX_so4_e': 1000.0, 'EX_fe3_e': 1000.0,
        'EX_fe2_e': 1000.0, 'EX_cbl1_e': 1000.0
    }

    # Reset all exchange bounds to 0 and apply M9 medium uptake rates
    for ex in ecoli_model.exchanges:
        if ex.id in m9_medium:
            ex.lower_bound = -float(m9_medium[ex.id])
        else:
            ex.lower_bound = 0.0

    # Set objective to WT biomass
    ecoli_model.reactions.BIOMASS_Ec_iML1515_core_75p37M.objective_coefficient = 0
    ecoli_model.reactions.BIOMASS_Ec_iML1515_WT_75p37M.objective_coefficient = 1
    ecoli_model.solver = 'gurobi'

    # Adjust uptake rate to match observed growth rate
    target_growth = metadata_df.loc[sample_id, 'Growth Rate (1/hr)']
    exchange_rxn = ecoli_model.reactions.get_by_id(carbon_exchanges[carbon_source])

    while ecoli_model.optimize().objective_value < target_growth and exchange_rxn.lower_bound > -100:
        exchange_rxn.lower_bound -= 0.01

    print(f'[INFO] {carbon_source} uptake for {sample_id}: {exchange_rxn.lower_bound:.2f}')
    return ecoli_model


# -------------------------------------------------------------
# Reads dataset from Excel
# -------------------------------------------------------------
def read_dataset():
    """
    Reads metadata and gene expression data from the provided Excel file.
    """
    metadata = pd.read_excel(
        '41467_2019_13483_MOESM4_ESM.xlsx',
        sheet_name='Metadata',
        index_col='Sample ID'
    )
    expression = pd.read_excel(
        '41467_2019_13483_MOESM4_ESM.xlsx',
        sheet_name='Expression Data',
        index_col='log-TPM'
    )
    return metadata, expression


# -------------------------------------------------------------
# Checks if two expression columns belong to the same condition
# -------------------------------------------------------------
def same_condition(exp_columns, i, j):
    base_i = exp_columns[i].rsplit('_', 1)[0].rstrip('_')
    base_j = exp_columns[j].rsplit('_', 1)[0].rstrip('_')
    return base_i == base_j


# -------------------------------------------------------------
# Filters dataset and averages replicates
# -------------------------------------------------------------
def filter_dataset():
    """
    Filters metadata for valid M9 media conditions and positive growth,
    keeps only supported carbon sources, and averages replicate expression data.
    """
    valid_sources = [
        'glucose', 'fumarate', 'pyruvate', 'fructose', 'glycerol', 'xylose', 'sorbitol',
        'D-ribose', 'glucarate', 'N-acetylglucosamine', 'galactose', 'gluconate',
        'D-lyxose', 'D-arabinose', 'm-tartrate', 'acetate'
    ]

    metadata, expression = read_dataset()

    # Filter metadata
    metadata = metadata[
        (metadata['Growth Rate (1/hr)'] > 0) &
        (metadata['Base Media'] == 'M9') &
        (metadata['Carbon Source (g/L)'].str.extract(r'([A-Za-z0-9\-]+)', expand=False).isin(valid_sources))
    ]

    # Keep only columns corresponding to filtered samples
    expression = expression.loc[:, [col for col in expression.columns if col in metadata.index]]

    # Keep only one replicate per condition in metadata
    filtered_metadata = metadata.drop(
        [cond for cond in metadata.index if int(cond.split('_')[-1]) > 1]
    )
    filtered_metadata.index = [idx.rsplit('_', 1)[0].rstrip('_') for idx in filtered_metadata.index]

    # Average replicates in expression data
    averaged_expression = pd.DataFrame(index=expression.index)
    exp_cols = list(expression.columns)

    i = 0
    while i < len(exp_cols):
        replicate_group = [exp_cols[i]]
        j = i + 1
        while j < len(exp_cols) and same_condition(exp_cols, i, j):
            replicate_group.append(exp_cols[j])
            j += 1

        group_name = exp_cols[i].rsplit('_', 1)[0].rstrip('_')
        averaged_expression[group_name] = expression[replicate_group].mean(axis=1)
        i = j

    # Convert from log2(TPM) to TPM
    averaged_expression = np.power(2, averaged_expression)

    return filtered_metadata, averaged_expression


# -------------------------------------------------------------
# Reads raw multiplicity data
# -------------------------------------------------------------
def read_raw_multiplicity():
    return pd.read_excel('raw_molteplicity.xlsx', sheet_name='molteplicity', index_col='locus_tag')


# -------------------------------------------------------------
# Main execution
# -------------------------------------------------------------
if __name__ == '__main__':
    metadata_df, expression_df = filter_dataset()
    multiplicity_df = read_raw_multiplicity()

    # Keep only genes common to both expression and multiplicity data
    common_genes = [gene for gene in multiplicity_df.index if gene in expression_df.index]
    expression_df = expression_df.loc[common_genes]
    multiplicity_df = multiplicity_df.loc[common_genes]

    # Run Riptide batch optimization
    run_riptide_batch(metadata_df, expression_df)

