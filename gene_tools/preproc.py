from collections import defaultdict
import os
import pandas as pd

def foldersLoad(root_path, folders):
    """
    Load tab-separated data files from multiple subfolders into a dictionary.

    Args:
        root_path (str): Root directory containing the subfolders.
        folders (list): List of subfolder names to load data from.

    Returns:
        dict: Dictionary mapping trait names to lists of dataframes.
    """
    all_data = defaultdict(list)
    for folder in folders:
        folder_path = os.path.join(root_path, folder)
        for file_name in os.listdir(folder_path):
            trait = file_name.split('_')[0]
            file_path = os.path.join(folder_path, file_name)
            df = pd.read_csv(file_path, sep='\t')
            all_data[trait].append(df)
    return all_data



def cleanDic(raw_dict, method_names):
    """
    Aggregates multiple method-specific dataframes into a single dataframe per trait.
    
    Args:
        raw_dict (dict): Dictionary where each key is a trait and the value is a list of dataframes
                         (one for each method: e.g., GWAS, eQTL, etc.).
        method_names (list of str): List of method names in the same order as the dataframes in the list.

    Returns:
        dict: A dictionary {trait: merged dataframe} where each dataframe contains:
              - EnsemblId
              - Trait
              - One p-value column per method (e.g., "GWAS_pvalue")
              - One beta column per method (e.g., "eQTL_b"), if available
    """
    result = {}

    for trait, df_list in raw_dict.items():
        method_dfs = []
        for df, method in zip(df_list, method_names):
            df_temp = df.copy()

            # Rename columns first
            rename_dict = {}
            if "p_value" in df_temp.columns:
                rename_dict["p_value"] = f"{method}_pvalue"
            if "b_ivw" in df_temp.columns:
                rename_dict["b_ivw"] = f"{method}_b"
            df_temp = df_temp.rename(columns=rename_dict)

            # Keep only relevant columns
            cols_to_keep = ["EnsemblId"] + list(rename_dict.values())
            df_temp = df_temp[cols_to_keep]

            method_dfs.append(df_temp)

        # Merge all method-specific dataframes on EnsemblId
        merged_df = method_dfs[0]
        for next_df in method_dfs[1:]:
            merged_df = pd.merge(merged_df, next_df, on="EnsemblId", how="outer")

        # Add the Trait column
        merged_df.insert(1, "Trait", trait)
        result[trait] = merged_df

    return result