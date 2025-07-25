import os
import pandas as pd
import glob
import re

def aggregate_true_interaction_results(base_dir="true_interaction_results"):
    """
    Parses Spectre output directories and merges them with simulation hit
    data to create a single results file.
    """
    print(f"INFO: Starting aggregation in base directory: '{base_dir}'")

    print("INFO: Loading all simulation hit data...")
    hits_files = glob.glob(os.path.join(base_dir, 'freq_*_dist_*', 'true_interaction_hits.tsv'))

    if not hits_files:
        print(f"ERROR: No 'true_interaction_hits.tsv' files found in any subdirectories of '{base_dir}'.")
        return

    # Read all TSV files into a list of DataFrames
    all_hits_dfs = [pd.read_csv(f, sep='\t') for f in hits_files]
    # Combine them into a single DataFrame
    master_hits_df = pd.concat(all_hits_dfs, ignore_index=True)

    # Create a unique key for merging, ensuring positions are integers
    master_hits_df['lookup_key'] = master_hits_df.apply(
        lambda row: f"{int(row['pos_A'])}_{int(row['pos_B'])}", axis=1
    )
    # Handle any potential duplicates from re-runs
    master_hits_df.drop_duplicates(subset='lookup_key', keep='first', inplace=True)
    master_hits_df.set_index('lookup_key', inplace=True)
    
    print(f"INFO: Loaded {len(master_hits_df)} unique simulation hits.")

    print("\nINFO: Loading all Spectre analysis results...")
    spectre_test_variables = [
        'test1a_b_clade_min_p0.5', 'test1a_b_clade_min_p0.8', 'test1a_r2_clade_min_p0.5',
        'test1a_b_clade_min_p0.9', 'test1a_b_clade_min_p0.95',
        'test1b_b_min', 'test2_overall_min_b'
    ]
    spectre_info_variables = [
        'snp1', 'snp2', 'snp3',
        'freq_snp1', 'freq_snp2', 'freq_S',
        'freq_K1', 'freq_K2', 'freq_O'
    ]
    target_spectre_variables = spectre_test_variables + spectre_info_variables
    
    all_spectre_results = []
    
    # Find all Spectre result files recursively
    search_pattern = os.path.join(base_dir, '**', '*_results.csv')
    result_files = glob.glob(search_pattern, recursive=True)

    if not result_files:
        print("WARN: No Spectre result files were found.")
        # Create an empty DataFrame to avoid crashing the merge step
        spectre_results_df = pd.DataFrame(columns=['lookup_key'] + target_spectre_variables).set_index('lookup_key')
    else:
        print(f"INFO: Found {len(result_files)} Spectre result files to process.")
        for f_path in result_files:
            try:
                # Find the parent directory that matches the position-based name
                run_id_dir = None
                current_path = os.path.dirname(f_path)
                while current_path.startswith(base_dir) and current_path != base_dir:
                    potential_run_id = os.path.basename(current_path)
                    if re.match(r'\d+_\d+_\d+_\d+', potential_run_id):
                        run_id_dir = potential_run_id
                        break
                    current_path = os.path.dirname(current_path)

                if not run_id_dir:
                    continue

                match = re.match(r'\d+_(\d+)_\d+_(\d+)', run_id_dir)
                if not match:
                    continue
                
                pos_a, pos_b = map(int, match.groups())
                lookup_key = f"{pos_a}_{pos_b}"
                
                spectre_df = pd.read_csv(f_path, index_col=0, header=0)
                
                result_dict = {'lookup_key': lookup_key}
                for var in target_spectre_variables:
                    result_dict[var] = spectre_df.loc[var].iloc[0] if var in spectre_df.index else None
                
                all_spectre_results.append(result_dict)
            except Exception as e:
                print(f"WARN: Could not process file '{f_path}'. Reason: {e}")

        spectre_results_df = pd.DataFrame(all_spectre_results)
        if not spectre_results_df.empty:
            spectre_results_df.set_index('lookup_key', inplace=True)

    print("\nINFO: Merging simulation data with Spectre results...")
    final_df = master_hits_df.merge(spectre_results_df, left_index=True, right_index=True, how='left')
    final_df.reset_index(drop=True, inplace=True) # Clean up index

    output_filename = "aggregated_true_interaction_results.csv"
    final_df.to_csv(output_filename, index=False)

    print("-" * 50)
    print(f"INFO: Done")
    print(f"INFO: Results for {len(final_df)} runs saved to '{output_filename}'")
    print("-" * 50)


if __name__ == '__main__':
    aggregate_true_interaction_results()
