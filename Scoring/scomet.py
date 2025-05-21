from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np


def Mean(df):

    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]

    df["Prioscore_mean"] = df[percentile_cols].mean(axis=1)

    return df



def Max(df):

    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]

    df["Prioscore_max"] = df[percentile_cols].max(axis=1)

    return df



def Min(df):

    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]

    df["Prioscore_min"] = df[percentile_cols].min(axis=1)

    return df


def Median(df):

    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]

    df["Prioscore_median"] = df[percentile_cols].median(axis=1)

    return df


def pca(df):
    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]
    filtered = df[percentile_cols].dropna()

    scaler = StandardScaler()
    scaled_filtered = scaler.fit_transform(filtered)

    pca_model = PCA(n_components=2)
    df_pca = pca_model.fit_transform(scaled_filtered)

    # Create a Series with PCA scores indexed like `filtered`
    pca_scores = pd.Series(df_pca[:, 0], index=filtered.index, name="PCA")

    # Copy original df and add PCA column where possible
    df_with_pca = df.copy()
    df_with_pca.loc[pca_scores.index, "PCA"] = pca_scores

    return df_with_pca

def Product(df) :
    """
    Compute a 'product'  score from log-transformed percentiles.

    Args:
        df (pd.DataFrame): Input dataframe with the apropriate columns.
        
    Returns:
        DataFrame: Modified dataframe with an added product-based score column.
    """
    # Identify percentile columns
    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]

    # Compute log of percentiles (must exclude NA to avoid warnings)
    log_vals = np.log(df[percentile_cols])

    # Compute the mean of log values across available (non-NA) columns
    df["Prioscore_product"] = log_vals.sum(axis=1) / log_vals.notna().sum(axis=1)

    return df
