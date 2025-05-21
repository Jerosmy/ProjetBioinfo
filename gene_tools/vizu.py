import matplotlib.pyplot as plt

def exploPlot(all_dat, trait="LDL", method="GWAS", save_path="expl_analysis.png"):
    """
    Generates exploratory graphs (scatter + histo) for a given trait and method

    Args:
        all_dat (dict): Dictionary containing the dataframes by trait
        trait (str): Name of the trait (as it appears in all_dat.keys())
        method (str): Chosen method to display
        save_path (str): Path where to save file ('None' to discard)
    """
    # Find the right DataFrame
    dataplot = all_dat.get(trait, [])
    df = None
    for d in dataplot:
        if d["Method"].iloc[0] == method:
            df = d
            break

    if df is None:
        raise ValueError(f"No method '{method}' found for trait '{trait}'.")

    # Create Figure
    plt.figure(figsize=(18, 5))

    # 1. Scatter P-value vs Ranking
    plt.subplot(1, 3, 1)
    plt.scatter(df["p_value"][::15], df["Ranking"][::15], s=8, alpha=0.8,
                color="purple", edgecolors="none")
    plt.xlabel("P-Values")
    plt.ylabel("Ranking")
    plt.title(f"{method}: P-value vs Ranking")
    plt.grid(True, linestyle="-", alpha=0.3)

    # 2. p-vals histograms
    plt.subplot(1, 3, 2)
    plt.hist(df["p_value"].dropna(), bins=20, color="skyblue", edgecolor="black")
    plt.xlabel("P-values")
    plt.ylabel("Frequency")
    plt.title(f"{method}: P-values distribution")
    plt.grid(True, linestyle="-", alpha=0.15)

    # 3. Rankings histograms
    plt.subplot(1, 3, 3)
    plt.hist(df["Ranking"].dropna(), bins=9, color="#ADEBB3", edgecolor='black')
    plt.xlabel("Rankings")
    plt.ylabel("Frequency")
    plt.title(f"{method}: Rankings distribution")
    plt.grid(True, linestyle='-', alpha=0.15)

    plt.tight_layout()

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

    plt.close()

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import pandas as pd

def compare_percentiles_boxplot(
    df,
    trait,
    method,             # e.g. "mean", "median", "max"
    reference_method,   # e.g. "GWAS", "eQTL"
    cutoff=0.01
):
    """
    Compare percentile values from:
    - top X% genes based on Prioscore_{method}
    - top X% genes based on {reference_method}_percentile

    Performs one-sided t-test (H1: Prioscore > reference).

    Parameters:
    - df: DataFrame
    - trait: str, for plot title
    - method: str, e.g. 'mean', 'median', 'product'
    - reference_method: str, e.g. 'GWAS', 'eQTL'
    - cutoff: float, top X% (e.g. 0.01 for top 1%)

    Returns:
    - t-stat, p-value
    - boxplot
    """

    method_col = f"Prioscore_{method}"
    ref_col = f"{reference_method}_percentile"

    if method_col not in df.columns or ref_col not in df.columns:
        raise ValueError(f"Missing columns: '{method_col}' or '{ref_col}'")

    df_clean = df[[method_col, ref_col]].dropna()
    n = int(len(df_clean) * cutoff)

    top_method_vals = df_clean.nsmallest(n, method_col)[method_col]
    top_reference_vals = df_clean.nsmallest(n, ref_col)[ref_col]

    t_stat, p_val = ttest_ind(top_method_vals, top_reference_vals, alternative="greater")

    plot_df = pd.DataFrame({
        f"Prioscore_{method}": top_method_vals,
        f"{reference_method}_percentile": top_reference_vals
    }).melt(var_name="Method", value_name="Percentile")

    plt.figure(figsize=(8, 6))
    sns.boxplot(x="Method", y="Percentile", data=plot_df)
    plt.title(f"{trait} — Top {int(cutoff*100)}% comparison\nT = {t_stat:.2f}, p = {p_val:.2e}")
    plt.ylabel("Percentile")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.tight_layout()
    plt.show()

    return t_stat, p_val
