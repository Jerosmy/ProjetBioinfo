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

# this keeps only rows where *both* values are present
    df_clean = df_t[[ref_col, prio_col]].dropna()


    ref_vals  = df_clean[ref_col]
    prio_vals = df_clean[prio_col]
    

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

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from gene_tools.analysis import evaluate_OR, evaluate_trait_scores, run_full_trait_pipeline
from Scoring.scomet import Mean, Max, Min, Median, Product

def plot_baseline_vs_prioscore(
    dic_clean,
    merged,
    paper_cutoff=0.05,
    our_cutoff=0.01,
    targetthreshold=3,
    figsize=(10,14),
    palette="tab10"
):
    """
    For each trait in dic_clean:
      • Load the paper’s best method & OR at paper_cutoff (hard-coded table).
      • Run your Prioscore pipeline at our_cutoff to get ORs for Mean/Max/Min/Median/Product.
      • Pick your top 3 Prioscore methods.
      • Plot a grouped horizontal bar chart: Baseline vs your top-3 for every trait.

    Parameters
    ----------
    dic_clean : dict
        trait → trait-specific DataFrame (must contain “Trait” & EnsemblId).
    merged : pd.DataFrame
        full drug reference (must contain “Trait”/“trait”, “Sum”, “EnsemblId”).
    paper_cutoff : float, default 0.05
        cutoff used in the paper for baseline ORs.
    our_cutoff : float, default 0.01
        cutoff you use for your Prioscore methods.
    targetthreshold : int, default 3
        minimum Sum to label a row as a drug target.
    figsize : tuple, default (10,14)
        figure size for the plot.
    palette : str or list, default "tab10"
        seaborn color palette.
    """
    # 1) paper baseline dict (at paper_cutoff)
    baseline_or = {
        "LDL":                  ["eQTL-GWAS",   7.5],
        "total_cholesterol":    ["Exome-GWAS",  8.0],
        "T2D":                  ["GWAS",        4.0],
        "T1D":                  ["Exome-GWAS",  3.0],
        "CKD":                  ["GWAS",        5.0],
        "Endometriosis":        ["Exome-GWAS",  1.0],
        "VTE":                  ["GWAS",        4.0],
        "CAD":                  ["GWAS",        4.0],
        "Systolic_BP":          ["eQTL-GWAS",   5.0],
        "Diastolic_BP":         ["GWAS",        3.5],
        "Atrial_fibrillation":  ["GWAS",        3.0],
        "Stroke":               ["Exome-GWAS",  3.0],
        "Rheumatoid_arthritis": ["eQTL-GWAS",   4.0],
        "Osteoporosis":         ["eQTL-GWAS",   3.5],
        "Osteoarthritis":       ["eQTL-GWAS",   1.5],
        "Atopic_dermatitis":    ["GWAS",        3.75],
        "Psoriasis":            ["eQTL-GWAS",   2.5],
        "IBD":                  ["GWAS",        4.0],
        "IBS":                  ["pQTL-GWAS",   8.0],
        "Asthma":               ["GWAS",        4.0],
        "COPD":                 ["eQTL-GWAS",   2.0],
        "Pneumonia":            ["eQTL-GWAS",   6.0],
        "Multiple_sclerosis":   ["eQTL-GWAS",   3.75],
        "Epilepsy":             ["eQTL-GWAS",   3.0],
        "Parkinsons_disease":   ["eQTL-GWAS",   3.0],
        "Alzheimer":            ["eQTL-GWAS",   2.0],
        "Schizophrenia":        ["GWAS",        3.0],
        "MDD":                  ["Exome-GWAS",  3.0],
        "Bipolar":              ["GWAS",        1.0],
        "Glaucoma":             ["Exome-GWAS",  2.0],
    }

    # 2) run your pipeline at our_cutoff
    finaldic = run_full_trait_pipeline(
        dic_clean,
        drug_reference=merged,
        scoring_functions=[Mean, Max, Min, Median, Product],
        evaluate_or_fn=evaluate_OR,
        evaluate_score_fn=evaluate_trait_scores,
        score_kwargs={"cutoff": our_cutoff, "overlap_method": "topvalues"},
        or_kwargs={
            "score_cols": [
                "Prioscore_mean", "Prioscore_max",
                "Prioscore_min", "Prioscore_median",
                "Prioscore_product"
            ],
            "method": "topvalues",
            "cutoff": our_cutoff
        },
        verbose=False
    )

    # 3) assemble plotting DataFrame
    rows = []
    for trait, info in finaldic.items():
        # your top-3 Prioscore methods
        our_items = info["OR_eval"].items()
        top3 = sorted(our_items, key=lambda kv: kv[1][0], reverse=True)[:3]
        for method_col, (orval, pval, n) in top3:
            rows.append({
                "Trait": trait,
                "Method": method_col.replace("Prioscore_",""),
                "OR": orval,
                "Source": f"Ours ({int(our_cutoff*100)}%)"
            })
        # add paper baseline
        bm, bo = baseline_or[trait]
        rows.append({
            "Trait": trait,
            "Method": bm,
            "OR": bo,
            "Source": f"Baseline ({int(paper_cutoff*100)}%)"
        })

    df_plot = pd.DataFrame(rows)

    # 4) plot
    plt.figure(figsize=figsize)
    sns.set_style("whitegrid")
    ax = sns.barplot(
        data=df_plot,
        y="Trait",
        x="OR",
        hue="Method",
        dodge=True,
        order=list(baseline_or.keys()),
        palette=palette
    )
    ax.axvline(1, ls="--", color="gray")
    ax.set_title(f"Baseline {paper_cutoff*100:.0f}% vs. Top‐3 Prioscore ({our_cutoff*100:.0f}%)")
    ax.set_xlabel("Odds‐Ratio")
    ax.set_ylabel("")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.show()
