from collections import defaultdict
import os
import pandas as pd


def foldersLoad(root_path, folders):
    """
    Load tab-separated data files from multiple subfolders into a dictionary.
    Skips files that cannot be read.
    """
    all_data = defaultdict(list)
    for folder in folders:
        folder_path = os.path.join(root_path, folder)
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            print(f"📄 Attempting to load: {file_path}")
            try:
                df = pd.read_csv(file_path, sep='\t')
                trait = file_name.split('_')[0]
                all_data[trait].append(df)
            except Exception as e:
                print(f"❌ Failed to load {file_path}: {e}")
    return all_data




def cleanDic(raw_dict, method_names = ['eQTL', 'Exome', 'GWAS', 'pQTL'], output = 'percentile'):
    """
    Aggregates multiple method-specific dataframes into a single dataframe per trait.
    
    Args:
        raw_dict (dict): Dictionary where each key is a trait and the value is a list of dataframes
                         (one for each method: e.g., GWAS, eQTL, etc.).
        method_names (list of str): List of method names in the same order as the dataframes in the list. 
                        Default is ['eQTL', 'Exome', 'GWAS', 'pQTL']
        output (str): weather the output should be in percentiles per method ('percentile', default) 
                        or in p_values and betas ('stats')

    Returns:
        dict: A dictionary {trait: merged dataframe} where each dataframe contains:
              - EnsemblId
              - Trait
              - One p-value column per method (e.g., "GWAS_pvalue")
              - One beta column per method (e.g., "eQTL_b"), if available
    """
    if output not in ['percentile', 'stats']:
        raise KeyError(f'"{output}" not an argument. Try "percentile" (default) or "stats"')
    

    result = {}

    for trait, df_list in raw_dict.items():
        method_dfs = []
        for df, method in zip(df_list, method_names):
            df_temp = df.copy()
            
            #If output is percentile
            if output == 'percentile' :
                #compute ranks

                beta = False
                threshold = 0.05 / len(df_temp)

                if "b_ivw" in df_temp.columns:
                    beta = True
                    my_betas = df_temp["p_value"] <= threshold
                    my_pvals = df_temp["p_value"] > threshold

                    ranking_betas = df_temp.loc[my_betas,"b_ivw"].abs().rank(method="min",ascending=False)
                    ranking_pvalues = df_temp.loc[my_pvals,"p_value"].rank(method="min",ascending=True)
                    ranking_pvalues +=len(ranking_betas)
                    bothranks =pd.concat([ranking_betas,ranking_pvalues])
                    df_temp["Score"] =bothranks.sort_index() / len(df_temp)
                else :
                    df_temp["Score"] = df_temp['p_value']

                df_temp['Rank'] = df_temp['Score'].rank()
                df_temp['Percentile'] = (df_temp['Rank'] / len(df_temp)) * 100
                # Rename columns
                rename_dict = {}
                rename_dict['Percentile'] = f"{method}_percentile"

                df_temp = df_temp.rename(columns=rename_dict)
            
            elif output == 'stats':

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

