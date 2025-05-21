import pandas as pd
from .preproc import cleanDic

def compute_scores(df):
    """
    Compute summary score, rank and percentile based on available p-values.

    Args:
        df (pd.DataFrame): Merged dataframe for a trait.

    Returns:
        pd.DataFrame: Dataframe with added Score, Rank, and Percentile columns.
    """

    perc_cols = [col for col in df.columns if col.endswith('_percentile')]

    df["NA_Count"] = df[perc_cols].isna().sum(axis=1)

    if "Method" in df.columns and "p_value" in df.columns and "b_ivw" in df.columns:
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
            
    
            df['Score'] = df[perc_cols].mean(axis=1)

    df['Rank'] = df['Score'].rank()
    df['Percentile'] = (df['Rank'] / len(df)) * 100
    return df



def NA_filtering(df, max_na=0):
    """
    Filter out rows with more than `max_na` missing scores.

    Argsument:
        df (pd.DataFrame): Scored dataframe containing 'NA_Count'.
        max_na (int): Maximum allowed number of NAs per row.

    Returns:
        pd.DataFrame: Filtered dataframe.
    """

    return df[df.isna().sum(axis=1) <= 0]



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
