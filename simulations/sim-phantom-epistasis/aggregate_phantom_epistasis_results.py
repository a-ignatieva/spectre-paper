import os
import pandas as pd
import glob
import re

def aggregate_spectre_results(base_dir="phantom_epistasis_results"):
    """
    Parses Spectre output directories and merges them with simulation hit
    data to create a single results file.
    """
    print(f"INFO: Starting aggregation in base directory: '{base_dir}'")

    # Define all variables to extract from each Spectre summary CSV.
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

    all_results = []

    # Find all the frequency directories (e.g., 'results/01_pct').
    freq_dirs = glob.glob(os.path.join(base_dir, '*_pct'))

    if not freq_dirs:
        print(f"FATAL ERROR: No frequency directories found in '{base_dir}'.")
        return

    for freq_dir in freq_dirs:
        try:
            frequency_category = os.path.basename(freq_dir)
            print(f"\n--- Processing frequency category: {frequency_category} ---")

            # Load the simulation hits data for this frequency
            hits_file_path = os.path.join(freq_dir, 'significant_hits.tsv')
            if not os.path.exists(hits_file_path):
                print(f"WARN: No 'significant_hits.tsv' found in {freq_dir}. Skipping.")
                continue
            
            hits_df = pd.read_csv(hits_file_path, sep='\t')
            hits_df['lookup_key'] = hits_df.apply(
                lambda row: f"{int(row['pos_A'])}_{int(row['pos_B'])}", axis=1
            )
            
            if hits_df['lookup_key'].duplicated().any():
                print(f"WARN: Duplicate interacting pairs found in {hits_file_path}. Keeping first instance.")
                hits_df = hits_df.drop_duplicates(subset='lookup_key', keep='first')

            hits_df.set_index('lookup_key', inplace=True)

            # Find and process all Spectre results for this frequency
            search_pattern = os.path.join(freq_dir, 'relate_trees', '*', '*_results.csv')
            result_files = glob.glob(search_pattern)

            if not result_files:
                print(f"INFO: No Spectre results found for {frequency_category}.")
                continue

            print(f"INFO: Found {len(result_files)} Spectre result files to process.")

            for f_path in result_files:
                # Extract info from Spectre run and merge with hits data
                run_id_dir = os.path.basename(os.path.dirname(f_path))
                
                match = re.match(r'\d+_(\d+)_\d+_(\d+)', run_id_dir)
                if not match:
                    print(f"WARN: Could not parse run_id from directory '{run_id_dir}'. Skipping file {f_path}")
                    continue
                
                pos_a, pos_b = map(int, match.groups())
                lookup_key = f"{pos_a}_{pos_b}"

                if lookup_key not in hits_df.index:
                    print(f"WARN: No matching hit found for run '{run_id_dir}'. Skipping.")
                    continue
                
                # Get all the original simulation data for this hit.
                hit_data = hits_df.loc[lookup_key].to_dict()
                
                hit_data['frequency_category'] = frequency_category
                hit_data['run_id'] = run_id_dir

                # Read the Spectre results CSV.
                spectre_df = pd.read_csv(f_path, index_col=0, header=0)
                
                # Extract all target Spectre variables.
                for var in target_spectre_variables:
                    value = spectre_df.loc[var].iloc[0] if var in spectre_df.index else None
                    hit_data[var] = value
                
                # Helper function to safely parse position from a 'chr:pos' string.
                def get_pos(snp_str):
                    if pd.isna(snp_str) or ':' not in str(snp_str):
                        return None
                    try:
                        return int(str(snp_str).split(':')[-1])
                    except (ValueError, IndexError):
                        return None
                
                # Add parsed position columns for Spectre's SNPs
                hit_data['spectre_pos_snp1'] = get_pos(hit_data.get('snp1'))
                hit_data['spectre_pos_snp2'] = get_pos(hit_data.get('snp2'))
                hit_data['spectre_pos_snp3'] = get_pos(hit_data.get('snp3'))
                
                all_results.append(hit_data)

        except Exception as e:
            print(f"ERROR: An unexpected error occurred while processing {freq_dir}. Reason: {e}")

    # Create and save the final aggregated DataFrame
    if not all_results:
        print("ERROR: No data was successfully aggregated.")
        return
        
    final_df = pd.DataFrame(all_results)

    # Define the combined list of columns for the final output.
    # Start with original hit info, then add Spectre freqs/positions, then Spectre test results.
    original_hit_cols = list(pd.read_csv(hits_file_path, sep='\t', nrows=0).columns)
    
    spectre_parsed_cols = [
        'spectre_pos_snp1', 'spectre_pos_snp2', 'spectre_pos_snp3',
        'freq_snp1', 'freq_snp2', 'freq_S', 'freq_K1', 'freq_K2', 'freq_O'
    ]

    output_columns = ['frequency_category', 'run_id'] + original_hit_cols + spectre_parsed_cols + spectre_test_variables
    
    # Ensure all columns exist in the dataframe before trying to reorder/select them.
    final_df_cols = [col for col in output_columns if col in final_df.columns]
    final_df = final_df[final_df_cols]

    # Save the aggregated results to a new CSV file.
    output_filename = "aggregated_phantom_epistasis_results_full.csv"
    final_df.to_csv(output_filename, index=False)

    print("-" * 50)
    print(f"INFO: Done")
    print(f"INFO: Results for {len(final_df)} runs saved to '{output_filename}'")
    print("-" * 50)


if __name__ == '__main__':
    aggregate_spectre_results()
