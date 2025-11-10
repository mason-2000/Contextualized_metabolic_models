#!/usr/bin/env python3
"""
This script integrates transcriptomic data (raw and normalized for gene multiplicity)
with the iML1515 genome-scale metabolic model of E. coli using RIPTiDe.
It contextualizes the model under different growth conditions described in Sastry et al. (2019),
by setting the correct media composition and carbon source prior to model fitting.
"""

import pandas as pd
import cobra
import riptide
import multiprocessing as mp
import time
import math


# ------------------------------------------------------------
# Utility function to initialize the log file for terminated runs
# ------------------------------------------------------------
def initialize_terminated_log():
    with open('terminated.txt', 'w') as f:
        f.write('')
    log_file = open('terminated.txt', 'a')
    return log_file


# ------------------------------------------------------------
# Worker process: runs RIPTiDe optimization (maxfit)
# ------------------------------------------------------------
def worker(metadata, expression_dict, condition_id,
           frac_min=0.5, frac_max=0.95, threshold=math.exp(-8),
           samples=1000, gpr=True):
    """
    Runs RIPTiDe maxfit for a given experimental condition.
    Saves the resulting context-specific model.
    """
    print(f'[INFO] Optimizing model for condition: {condition_id}...')

    ecoli_model = prepare_ecoli_model(metadata, condition_id)
    result = riptide.maxfit(
        model=ecoli_model,
        transcriptome=expression_dict,
        frac_min=frac_min,
        frac_max=frac_max,
        gpr=gpr,
        samples=samples,
        frac_step=0.05
    )

    print(f'Objective value changed from: {ecoli_model.optimize().objective_value} '
          f'to: {result.model.optimize().objective_value}')

    riptide.save_output(
        riptide_obj=result,
        path=f'normalized_exp_maxfit_riptide_runs_iML1515/{condition_id}',
        file_type='sbml'
    )


# ------------------------------------------------------------
# Launches RIPTiDe runs in parallel and terminates hung processes
# ------------------------------------------------------------
def run_riptide(metadata, normalized_expression, timeout=1500, num_threads=1):
    """
    Runs RIPTiDe maxfit in parallel for all experimental conditions.
    Terminates any process exceeding the timeout limit.
    """
    terminated_log = initialize_terminated_log()
    start = 0
    stop = num_threads

    while start < len(normalized_expression.columns):
        processes = []

        # Create a process for each sample
        for sample_id in normalized_expression.columns[start:stop]:
            exp_dict = normalized_expression[sample_id].to_dict()
            p = mp.Process(target=worker, args=(metadata, exp_dict, sample_id), name=sample_id)
            p.start()
            processes.append(p)

        # Wait for completion or kill after timeout
        t1 = time.monotonic()
        for p in processes:
            elapsed = time.monotonic() - t1
            remaining = timeout - elapsed

            p.join(timeout=remaining)

            if p.is_alive():
                p.terminate()
                p.join()
                print(f'[WARNING] Terminated hung process for sample: {p.name}')
                terminated_log.write(f'KILLED: {p.name}\n')

        start += num_threads
        stop += num_threads


# ------------------------------------------------------------
# Prepares the E. coli model for the given condition
# ------------------------------------------------------------
def prepare_ecoli_model(metadata, condition_id):
    """
    Loads the iML1515 model, applies M9 medium,
    and adjusts carbon source uptake to match target growth rate.
    """
    ecoli_model = cobra.io.read_sbml_model('iML1515.xml')

    # Extract carbon source information from metadata
    carbon_source_name = metadata.loc[condition_id, 'Carbon Source (g/L)'].split('(')[0]

    carbon_source_rxns = {
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

    # Define M9 minimal medium exchange reactions
    m9_medium = {
        carbon_source_rxns[carbon_source_name]: 4.0,
        'EX_cobalt2_e': 1000.0,
        'EX_h_e': 1000.0,
        'EX_h2o_e': 1000.0,
        'EX_k_e': 1000.0,
        'EX_cu2_e': 1000.0,
        'EX_mg2_e': 1000.0,
        'EX_mn2_e': 1000.0,
        'EX_mobd_e': 1000.0,
        'EX_na1_e': 1000.0,
        'EX_nh4_e': 1000.0,
        'EX_ca2_e': 1000.0,
        'EX_ni2_e': 1000.0,
        'EX_o2_e': 1000.0,
        'EX_cl_e': 1000.0,
        'EX_pi_e': 1000.0,
        'EX_zn2_e': 1000.0,
        'EX_so4_e': 1000.0,
        'EX_fe3_e': 1000.0,
        'EX_fe2_e': 1000.0,
        'EX_cbl1_e': 1000.0
    }

    # Reset all exchange bounds to 0 and apply M9 medium uptake rates
    for ex in ecoli_model.exchanges:
        if ex.id in m9_medium:
            ex.lower_bound = -float(m9_medium[ex.id])
        else:
            ex.lower_bound = 0.0

    # Set biomass objective to WT
    ecoli_model.reactions.BIOMASS_Ec_iML1515_core_75p37M.objective_coefficient = 0
    ecoli_model.reactions.BIOMASS_Ec_iML1515_WT_75p37M.objective_coefficient = 1
    ecoli_model.solver = 'gurobi'

    # Adjust carbon source uptake to match target growth rate
    target_growth = metadata.loc[condition_id, 'Growth Rate (1/hr)']

    while ecoli_model.optimize().objective_value < target_growth and \
          ecoli_model.reactions.get_by_id(carbon_source_rxns[carbon_source_name]).lower_bound > -100:
        ecoli_model.reactions.get_by_id(carbon_source_rxns[carbon_source_name]).lower_bound -= 0.01

    print(f'[INFO] {carbon_source_name} uptake for model {condition_id}: '
          f'{ecoli_model.reactions.get_by_id(carbon_source_rxns[carbon_source_name]).lower_bound}')

    return ecoli_model


# ------------------------------------------------------------
# Normalize expression data by gene multiplicity
# ------------------------------------------------------------
def normalize_expression(raw_multiplicity, expression):
    """
    Converts log2(TPM) values to linear scale and normalizes them
    by the corresponding raw multiplicity values.
    """
    normalized_expression = pd.DataFrame(index=raw_multiplicity.index, columns=raw_multiplicity.columns)
    conditions = list(expression.columns)

    for gene, _ in expression.iterrows():
        for condition in conditions:
            log_exp = float(expression.loc[gene, condition])
            exp_value = math.pow(2, log_exp)
            multiplicity = float(raw_multiplicity.loc[gene, condition])
            normalized_expression.loc[gene, condition] = exp_value / multiplicity

    return normalized_expression


# ------------------------------------------------------------
# Read multiplicity values from Excel file
# ------------------------------------------------------------
def read_raw_multiplicity():
    return pd.read_excel('raw_molteplicity.xlsx', sheet_name='molteplicity', index_col='locus_tag')


# ------------------------------------------------------------
# Read metadata and expression data from Excel file
# ------------------------------------------------------------
def read_expression_dataset():
    metadata = pd.read_excel('41467_2019_13483_MOESM4_ESM.xlsx', sheet_name='Metadata', index_col='Sample ID')
    expression = pd.read_excel('41467_2019_13483_MOESM4_ESM.xlsx', sheet_name='Expression Data', index_col='log-TPM')
    return metadata, expression


# ------------------------------------------------------------
# Check if two expression column names correspond to the same condition
# ------------------------------------------------------------
def same_condition(exp_col, i, j):
    condition_i = exp_col[i].rsplit('_', 1)[0].rstrip('_')
    condition_j = exp_col[j].rsplit('_', 1)[0].rstrip('_')
    return condition_i == condition_j


# ------------------------------------------------------------
# Filter dataset for valid M9 conditions and merge replicates
# ------------------------------------------------------------
def filter_dataset():
    """
    Filters the Sastry et al. dataset for valid M9 media and known carbon sources.
    Averages biological replicates of the same condition.
    """
    carbon_sources = [
        'glucose', 'fumarate', 'pyruvate', 'fructose', 'glycerol', 'xylose',
        'sorbitol', 'D-ribose', 'glucarate', 'N-acetylglucosamine',
        'galactose', 'gluconate', 'D-lyxose', 'D-arabinose',
        'm-tartrate', 'acetate'
    ]

    metadata, expression = read_expression_dataset()

    # Filter metadata for valid growth and media conditions
    metadata = metadata[
        (metadata['Growth Rate (1/hr)'] > 0) &
        (metadata['Base Media'] == 'M9') &
        (metadata['Carbon Source (g/L)'].str.extract(r'([A-Za-z0-9\-]+)', expand=False).isin(carbon_sources))
    ]

    # Keep only expression columns present in filtered metadata
    expression = expression.drop([col for col in expression.columns if col not in metadata.index], axis='columns')

    # Keep one replicate set per condition (remove "_2", "_3"...)
    filtered_metadata = metadata.drop([cond for cond in metadata.index if int(cond.split('_')[-1]) > 1])
    filtered_metadata.index = [idx.rsplit('_', 1)[0].rstrip('_') for idx in filtered_metadata.index]

    filtered_expression = pd.DataFrame(index=expression.index)
    exp_cols = list(expression.columns)

    # Merge replicates by averaging their expression values
    i = 0
    while i < len(exp_cols):
        replicates = [exp_cols[i]]
        j = i + 1
        while j < len(exp_cols) and same_condition(exp_cols, i, j):
            replicates.append(exp_cols[j])
            j += 1

        condition_root = exp_cols[i].rsplit('_', 1)[0].rstrip('_')
        filtered_expression[condition_root] = expression[replicates].mean(axis=1)
        i = j

    return filtered_metadata, filtered_expression


# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------
if __name__ == '__main__':
    metadata, expression = filter_dataset()
    raw_multiplicity = read_raw_multiplicity()

    # Keep only shared genes between datasets
    shared_genes = [gene for gene in raw_multiplicity.index if gene in expression.index]
    expression = expression.loc[shared_genes]
    raw_multiplicity = raw_multiplicity.loc[shared_genes]

    # Normalize expression values by multiplicity
    normalized_expression = normalize_expression(raw_multiplicity, expression)

    # Run RIPTiDe for all experimental conditions
    run_riptide(metadata, normalized_expression)

