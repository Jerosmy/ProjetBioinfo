from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


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


def pca(stat_df, 
        wrk_df = None, 
        filling = "fill", 
        explicit = False,
        plot_cor = False):
    """
    Computes the Principal component analysis of the different traits p_values.

    Args:
        stat_df (pd.DataFrame): A trait dataframe containing p_values of different techniques
        wrk_df (pd.DataFrame): A trait dataframe to which append the pca percentile ranking ('stat_df' used by default)
        filling (str): Whether to fill in Nan ('fill', default) or exclude them ('drop')
        explicit (boolean): Whether to output explicit pca results (False by default)
        plot_cor (boolean): Show the scatter plot and correlation value between pca rankings 
    Returns:
        A dictionary containning the initial df with a new pca column ('full-df'),
        zscores for each column ('zscores') and
        the pca results ('pca'), if explicit == True

    """
    if wrk_df is None :
        wrk_df = stat_df
    
    p_valcols = [col for col in stat_df.columns if col.endswith("_pvalue")]

    if filling == "fill":
        stat_df[p_valcols] = stat_df[p_valcols].T.fillna(stat_df[p_valcols].median(axis=1)).T
    elif filling == 'drop' :
        stat_df = stat_df.dropna()

    #Filtered and borned data (.clip to avoid 0 or 1 that would give ±inf)
    filtered = stat_df[p_valcols].clip(lower=1e-10, upper=1 - 1e-10)

    #Normalize
    zscores = pd.DataFrame(norm.ppf(1 - filtered), columns=filtered.columns, index=filtered.index)

    #Run PCA
    pca_model = PCA(n_components=2)
    df_pca = pca_model.fit_transform(zscores)

    # Create a Series with PCA scores indexed like `filtered`
    pca_scores = pd.Series(df_pca[:, 0], index=filtered.index, name="PCA")


    #PC axes have arbitrary orientation, see how it correlates with the pvals and invert if needed
    corr, _ = spearmanr(pca_scores, stat_df[p_valcols].mean(axis=1))
    if corr > 0: # type: ignore (turns off unesceary warnings)
        pca_scores *= -1  # Inverts PCA values that are inverted (i.e. if big 1stPC val is also with big pval)
    print(corr)

    if plot_cor:
        mean_pval = stat_df[p_valcols].median(axis=1)
        plt.scatter(mean_pval, pca_scores)
        plt.xlabel("Medianne des p-values")
        plt.ylabel("Score PCA (PC1)")
        plt.title("PC1 vs median des p-values")
        plt.grid(True)
        plt.show()
        print(corr)

    #Sort it and transform in percentile
    sorted_pca = pd.Series(100 * (1 - (pca_scores.rank(method='min') / pca_scores.shape[0])), index=pca_scores.index, name="Prioscore_PCA")

    # Copy original df and add PCA column where possible
    df_with_pca = wrk_df.copy()
    df_with_pca.loc[pca_scores.index, "Prioscore_PCA"] = sorted_pca

    if explicit :
        results = {}
        #Isolate principal components
        pca_df = pd.DataFrame(df_pca, columns=[f"PC{i+1}" for i in range(df_pca.shape[1])], index=zscores.index)
        results['pca'] = pca_df
        results['zscores'] = zscores
        results['full_df'] = df_with_pca
    else :
        results = df_with_pca

    return results

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
