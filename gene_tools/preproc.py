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



def cleanDic(trait_data, trait_name):
    """
    Clean and merge dataframes related to a specific trait.

    Args:
        trait_data (list): List of dataframes for the same trait.
        trait_name (str): Name of the trait.

    Returns:
        pd.DataFrame: Merged and cleaned dataframe with scores and ranks.
    """
    cleaned_dict = {}

    for df in trait_data:
        method_name = df['Method'][0].split('-')[0]
        main_data = {
            'EnsemblId': df['EnsemblId'],
            f'{method_name}_p_value': df['p_value']
        }

        if 'b_ivw' in df.columns:
            main_data[f'{method_name}_b_ivw'] = df['b_ivw']

        id_trait_df = pd.DataFrame({
            'EnsemblId': df['EnsemblId'],
            'Trait': trait_name
        })
        method_df = pd.DataFrame(main_data)

        cleaned_dict[df['Method'][0]] = [method_df, id_trait_df]

    keys = list(cleaned_dict.keys())
    merged_df = pd.merge(cleaned_dict[keys[0]][1], cleaned_dict[keys[0]][0])

    for key in keys[1:]:
        merged_df = pd.merge(merged_df, cleaned_dict[key][0], how='outer', on='EnsemblId')

    return merged_df