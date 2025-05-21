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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

def compare_drug_target_percentiles(
    df,
    merged,
    trait,
    reference_method="GWAS",
    method="min",
    targetthreshold=3
):
    """
    For a given trait, compare the percentile distributions of your drug targets
    between the GWAS reference and Prioscore_min, via:
      - one-sided t-test (H1: Prioscore_min < GWAS)
      - side-by-side boxplots

    Parameters
    ----------
    df : pd.DataFrame
        Trait-specific table. Must contain columns:
          - "Trait"
          - "EnsemblId"
          - "{reference_method}_percentile"  (e.g. "GWAS_percentile")
          - "Prioscore_{method}"             (e.g. "Prioscore_min")
    merged : pd.DataFrame
        Full drug-reference table. Must contain columns:
          - "Trait" or "trait"
          - "EnsemblId"
          - "Sum"
    trait : str
        Which trait to pull (matches df["Trait"] and merged["Trait"]).
    reference_method : str, default "GWAS"
        Prefix of the percentile column in `df` (e.g. "GWAS" → "GWAS_percentile").
    method : str, default "min"
        Suffix for the Prioscore column (e.g. "min" → "Prioscore_min").
    targetthreshold : int, default 3
        Minimum `Sum` in `merged` to call a row a drug target.

    Returns
    -------
    t_stat : float
        t‐statistic for the one‐sided test.
    p_val : float
        p‐value for the one‐sided test.
    """
    # 1) Filter to just that trait in both tables
    df_t = df[df["Trait"] == trait]
    key_trait = "Trait" if "Trait" in merged.columns else "trait"
    dr_t = merged[(merged[key_trait] == trait) & (merged["Sum"] >= targetthreshold)]

    # 2) Restrict to drug‐target genes
    target_ids = set(dr_t["EnsemblId"])
    df_t  = df_t[df_t["EnsemblId"].isin(target_ids)]
    if df_t.empty:
        raise ValueError(f"No drug targets found for trait {trait} with Sum>={targetthreshold}")

    # 3) Grab the two columns
    ref_col    = f"{reference_method}_percentile"
    prio_col   = f"Prioscore_{method}"
    for col in (ref_col, prio_col):
        if col not in df_t.columns:
            raise KeyError(f"Column {col!r} not in df for trait {trait}")

    ref_vals  = df_t[ref_col].dropna()
    prio_vals = df_t[prio_col].dropna()
    

    # 4) One‐sided t-test: H1 = Prioscore_min < GWAS
    t_stat, p_val = ttest_ind(prio_vals, ref_vals, alternative="less",nan_policy="omit")

    # 5) Build the boxplot
    plot_df = pd.DataFrame({
        reference_method:   ref_vals,
        f"Prioscore_{method}": prio_vals
    }).melt(var_name="Method", value_name="Percentile")

    plt.figure(figsize=(6,6))
    sns.boxplot(x="Method", y="Percentile", data=plot_df)
    plt.title(f"{trait} drug‐targets: {reference_method} vs Prioscore_{method}\n"
              f"T = {t_stat:.2f}, p = {p_val:.2e}")
    plt.ylabel("Percentile")
    plt.grid(linestyle="--", alpha=0.3)
    plt.tight_layout()
    plt.show()

    return t_stat, p_val
