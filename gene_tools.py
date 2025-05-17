import os
import pandas as pd
from collections import defaultdict

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
            f'{method_name}_percentile': df['Percentile']
        }

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


def compute_scores(df):
    """
    Compute summary score, rank and percentile based on available p-values.

    Args:
        df (pd.DataFrame): Merged dataframe for a trait.

    Returns:
        pd.DataFrame: Dataframe with added Score, Rank, and Percentile columns.
    """
    if "Method" in df.columns and "p_value" in df.columns and "b_ivw" in df.columns():
        method_name= df["Method"].iloc[0]
        if method_name in ["eQTL_GWAS_blood","pQTL-GWAS"]:
            threshold = 0.05/len(df.rows())
            my_betas = df["p_value"] <= threshold
            my_pvals = df["p_value"] > threshold
            ranking_betas = df.loc[my_betas,"b_ivw"].abs().rank(method="min",ascending=False)
            ranking_pvalues = df.loc[my_pvals,"p_value"].rank(method="min",ascending=True)
            ranking_pvalues +=len(ranking_betas)
            
            bothranks =pd.concat([ranking_betas,ranking_pvalues])
            df["Score"] =bothranks.sort_index() / len(df)
        else:
            perc_cols = [col for col in df.columns if col.endswith('_percentile')]
            df['Score'] = df[perc_cols].mean(axis=1)

    df['Rank'] = df['Score'].rank()
    df['Percentile'] = (df['Rank'] / len(df)) * 100
    return df


def geneScores(all_data):
    """
    Process all traits: clean, merge and compute scores for each trait.

    Args:
        all_data (dict): Dictionary mapping traits to lists of dataframes.

    Returns:
        dict: Dictionary mapping traits to cleaned and scored dataframes.
    """
    result = {}
    for trait, dfs in all_data.items():
        merged = cleanDic(dfs, trait)
        scored = compute_scores(merged)
        result[trait] = scored
    return result

def NaCount(dataframe, show=False):
    """
    Inspect the presence of NaN values in a gene scoring dataframe.

    Args:
        dataframe (pd.DataFrame): A dataframe containing gene prioritization scores.
        show (bool): If True, prints the results to the console. Default is False.

    Returns:
        dict: Summary statistics with counts of NaNs and relevant columns.
    """
    traitsNa = dataframe['Trait'].isna().sum()
    ensemblNa = dataframe['EnsemblId'].isna().sum()
    lendf = len(dataframe)

    # Find score columns
    score_cols = [col for col in dataframe.columns if col.endswith('_percentile')]
    
    # NaN count per score column
    colsNan = dataframe[score_cols].isna().sum()

    # Count rows where all score columns are NaN
    all_score_nan = dataframe[score_cols].isna().all(axis=1).sum()

    if show:
        print("Trait NaNs:", traitsNa)
        print("EnsemblId NaNs:", ensemblNa)
        print("Total rows:", lendf)
        print("NaNs per score column:")
        print(colsNan)
        print("Rows with all score columns as NaN:", all_score_nan)

    return {
        "Trait_NaNs": traitsNa,
        "EnsemblId_NaNs": ensemblNa,
        "Total_rows": lendf,
        "Score_columns": score_cols,
        "NaNs_per_score_column": colsNan.to_dict(),
        "Rows_with_all_score_NaNs": all_score_nan
    }


#Main execution block

if __name__ == "__main__":
    root_path = '00_BPStart/00_gene_prioritization'
    folders = ['eQTL_GWAS_blood', 'Exome', 'GWAS', 'pQTL_GWAS']

    all_dat = foldersLoad(root_path, folders)
    cleaned_data = geneScores(all_dat)