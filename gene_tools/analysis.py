import numpy as np
from scipy.stats import fisher_exact
import pandas as pd
from scipy.stats import fisher_exact

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
    score_cols = [col for col in dataframe.columns if (col.endswith('_percentile') or col.endswith('_pvalue') or col.endswith('_b'))]
    
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



from scipy.stats import fisher_exact
from scipy.stats import fisher_exact

def evaluate_OR(
    reference,
    mydf,
    scoring_functions=None,
    score_cols=None,
    sum_threshold=3,
    cutoff=0.01,
    method="topvalues",  # or 'lessthan'
    printer=True
):
    """
    Compute odds ratios and Fisher's exact test p-values for enrichment of drug targets
    among top-ranked genes based on prioritization scores.

    Now:
      - auto-filters `reference` to the trait in `mydf`
      - applies `scoring_functions` to `mydf` to generate Prioscore_* columns
      - if `score_cols` is None, auto-detects any column starting with "Prioscore_"
    """
    # 1) Infer trait name
    trait_vals = mydf["Trait"].dropna().unique()
    if len(trait_vals) != 1:
        raise ValueError(f"Expected exactly one unique 'Trait' in df, found: {trait_vals}")
    trait_name = trait_vals[0]

    # 2) Filter reference to this trait (works whether column is 'Trait' or 'trait')
    ref = reference
    if "Trait" in ref.columns:
        ref = ref[ref["Trait"] == trait_name]
    elif "trait" in ref.columns:
        ref = ref[ref["trait"] == trait_name]

    # 3) Copy & apply scoring functions to build Prioscore_* columns
    df_proc = mydf.copy()
    if scoring_functions:
        for fn in scoring_functions:
            df_proc = fn(df_proc)

    # 4) Decide which columns to evaluate
    if score_cols is None:
        score_cols = [c for c in df_proc.columns if c.startswith("Prioscore_")]
    else:
        missing = [c for c in score_cols if c not in df_proc.columns]
        if missing:
            raise KeyError(f"Score columns not found: {missing}")

    # 5) Build sets for Fisher’s test
    drug_targets = set(ref.loc[ref["Sum"] >= sum_threshold, "EnsemblId"])
    all_ids      = set(df_proc["EnsemblId"])

    results = {}
    trait_result = {}

    # 6) Compute OR + p-value per score column
    for col in score_cols:
        if method == "topvalues":
            top_n = int(cutoff * len(df_proc))
            subset = df_proc.sort_values(by=col, ascending=True).head(top_n)
            desc   = f"top {cutoff*100:.1f}%"
        elif method == "lessthan":
            subset = df_proc[df_proc[col] < cutoff * 100]
            desc   = f"score < {cutoff*100}"
        else:
            raise ValueError("Invalid method: choose 'topvalues' or 'lessthan'")

        top_ids  = set(subset["EnsemblId"])
        rest_ids = all_ids - top_ids

        A = len(top_ids & drug_targets)
        B = len(top_ids - drug_targets)
        C = len(rest_ids & drug_targets)
        D = len(rest_ids - drug_targets)

        oddsratio, p_value = fisher_exact([[A, B], [C, D]])

        if printer:
            print(f"[{trait_name}] {col}: OR = {oddsratio:.2f}, p = {p_value:.4e}, drug targets in {desc} = {A}")

        trait_result[col] = [round(oddsratio, 3), round(p_value, 3), A]

    results[trait_name] = trait_result
    return results





def evaluate_trait_scores(
    df,
    drug_reference,
    scoring_functions=None,
    score_columns=None,
    overlap_method="topvalues",
    cutoff=0.01,
    printer=True,
    targetthreshold=3
):
    """
    Evaluate score overlaps with drug targets for a single trait.

    Parameters:
        df (pd.DataFrame): Trait-specific DataFrame (must contain exactly one 'Trait').
        drug_reference (pd.DataFrame): DataFrame with 'trait', 'Sum', and 'EnsemblId'.
        scoring_functions (list of callables): each fn(df) -> df with new Prioscore_* columns.
        score_columns (list of str, optional): Which columns to evaluate. 
            If None, we'll auto-detect any column starting with "Prioscore_".
        overlap_method (str): "topvalues" or "lessthan".
        cutoff (float): fraction (for 'topvalues') or percentile (for 'lessthan').
        printer (bool): whether to print results.
        targetthreshold (int): minimum Sum to call a row a drug target.

    Returns:
        dict: { trait_name: { score_col: percent_overlap } }
    """
    # 1) Infer trait name
    trait_vals = df["Trait"].dropna().unique()
    if len(trait_vals) != 1:
        raise ValueError(f"Expected exactly one unique 'Trait' in df, found: {trait_vals}")
    trait_name = trait_vals[0]

    # 2) Make a working copy & apply scoring functions
    df_proc = df.copy()
    if scoring_functions:
        for fn in scoring_functions:
            df_proc = fn(df_proc)

    # 3) Decide which columns to evaluate
    if score_columns is None:
        # pick up anything starting with "Prioscore_"
        score_columns = [c for c in df_proc.columns if c.startswith("Prioscore_")]
    else:
        # sanity‐check that they exist
        missing = [c for c in score_columns if c not in df_proc.columns]
        if missing:
            raise KeyError(f"Score columns not found: {missing}")

    # 4) Build drug target set for this trait
    drug_targets = set(
        drug_reference.loc[
            (drug_reference["trait"] == trait_name) &
            (drug_reference["Sum"] >= targetthreshold),
            "EnsemblId"
        ]
    )

    n = len(df_proc)
    trait_result = {}

    # 5) Compute overlap for each score column
    for col in score_columns:
        if overlap_method == "topvalues":
            top_n = int(cutoff * n)
            subset = df_proc.sort_values(by=col, ascending=True).head(top_n)
        elif overlap_method == "lessthan":
            subset = df_proc[df_proc[col] < cutoff * 100]
        else:
            raise ValueError("Invalid overlap_method: choose 'topvalues' or 'lessthan'")

        top_ids = set(subset["EnsemblId"])
        overlap = top_ids & drug_targets
        percent = 100 * len(overlap) / len(top_ids) if top_ids else 0.0

        trait_result[col] = round(percent, 3)
        if printer:
            print(f"[{trait_name}] {col}: {percent:.2f}% overlap")

    return {trait_name: trait_result}


def run_full_trait_pipeline(
    trait_dict,
    drug_reference,
    scoring_functions,
    evaluate_score_fn = evaluate_trait_scores,
    evaluate_or_fn  = evaluate_OR,
    score_kwargs=None,
    or_kwargs=None,
    verbose=True
):
    """
    Run scoring + evaluation across all traits in a dictionary of DataFrames.


    Args:
    - trait_dict (dict): {trait_name: df}
    - drug_reference (pd.DataFrame): merged drug reference with 'trait', 'Sum', 'EnsemblId'
    - scoring_functions (list): list of functions like [Mean, Max, Min, Median, Product]
    - evaluate_score_fn (func): function like evaluate_trait_scores
    - evaluate_or_fn (func): function like evaluate_OR
    - score_kwargs (dict): optional kwargs for evaluate_score_fn
    - or_kwargs (dict): optional kwargs for evaluate_or_fn

    Returns:
    - dict: {trait: {'data': df_with_scores, 'score_eval': %, 'OR_eval': [OR, p, n]}}
    """
    if score_kwargs is None:
        score_kwargs = {}
    if or_kwargs is None:
        or_kwargs = {}

    results = {}

    for trait, df in trait_dict.items():
        df = df.copy()
        curdrugref = drug_reference[drug_reference['trait'] == trait]
        # Apply all scoring functions
        for func in scoring_functions:
            df = func(df)

        # Evaluate % of drug targets in top scores
        score_result = evaluate_score_fn(df=df, drug_reference=curdrugref, printer=False, **score_kwargs)

        # Evaluate OR, p-value, n
        or_result = evaluate_or_fn(reference=curdrugref, mydf=df, printer=False, **or_kwargs)

        if verbose:
            print(f"Successfully processed {trait}, target: {curdrugref.shape}")

        # Package everything
        results[trait] = {
            "data": df,
            "score_eval": list(score_result.values())[0], 
            "OR_eval": list(or_result.values())[0]        
        }

    return results
