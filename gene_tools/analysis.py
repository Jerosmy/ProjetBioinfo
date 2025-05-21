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




def evaluate_OR(
    reference,
    mydf,
    score_cols=["Prioscore_mean", "Prioscore_max", "Prioscore_min", "Prioscore_median", "Prioscore_product"],
    sum_threshold=3,
    cutoff=0.01,
    method="topvalues",  # or 'lessthan'
    printer=True
):
    """
    Compute odds ratios and Fisher's exact test p-values for enrichment of drug targets
    among top-ranked genes based on prioritization scores.

    Parameters:
        reference (pd.DataFrame): Drug target info with 'Sum' and 'EnsemblId'.
        mydf (pd.DataFrame): Trait-specific DataFrame (must contain 'Trait').
        score_cols (list): Score columns to evaluate.
        sum_threshold (int): Minimum 'Sum' to define drug targets.
        cutoff (float): Top X% or percentile threshold.
        method (str): "topvalues" or "lessthan".
        printer (bool): Whether to print results.

    Returns:
        dict: {trait: {score_col: [OR, p-value, drug target count in top]}}
    """
    # Extract trait name
    trait_col_values = mydf["Trait"].dropna().unique()
    if len(trait_col_values) != 1:
        raise ValueError(f"Expected exactly one unique 'Trait' in df, found: {trait_col_values}")
    trait_name = trait_col_values[0]

    # Filter to shared genes only
    valid_ids = set(mydf["EnsemblId"]) & set(reference["EnsemblId"])
    mydf = mydf[mydf["EnsemblId"].isin(valid_ids)]
    reference = reference[reference["EnsemblId"].isin(valid_ids)]

    results = {}
    drug_targets = set(reference.loc[reference["Sum"] >= sum_threshold, "EnsemblId"])
    all_ids = set(mydf["EnsemblId"])

    trait_result = {}

    for col in score_cols:
        if method == "topvalues":
            subset = mydf.sort_values(by=col, ascending=True).head(int(cutoff * len(mydf)))
            desc = f"top {int(cutoff * 100)}%"
        elif method == "lessthan":
            subset = mydf[mydf[col] < cutoff * 100]
            desc = f"score < {cutoff * 100}"
        else:
            raise ValueError("Invalid method: choose 'topvalues' or 'lessthan'")

        top_ids = set(subset["EnsemblId"])
        rest_ids = all_ids - top_ids

        A = len(top_ids & drug_targets)
        B = len(top_ids - drug_targets)
        C = len(rest_ids & drug_targets)
        D = len(rest_ids - drug_targets)

        oddsratio, p_value = fisher_exact([[A, B], [C, D]])

        if printer:
            print(f"[{trait_name}] {col}: OR = {oddsratio:.2f}, p = {p_value:.4e}, drug targets in {desc} = {A}")

        trait_result[col] = [round(oddsratio,3), round(p_value,3), round(A,3)] # type: ignore


    results[trait_name] = trait_result
    return results




        
def evaluate_trait_scores(
    df,
    drug_reference,
    score_columns=["Prioscore_mean", "Prioscore_max", "Prioscore_min", "Prioscore_median", "Prioscore_product"],
    overlap_method="topvalues",
    cutoff=0.01,
    printer=True,
    targetthreshold=3
):
    """
    Evaluate score overlaps with drug targets for a single trait.

    Returns:
    - dict: {trait: {score_column: float (percent overlap)}}
    """
    trait_col_values = df["Trait"].dropna().unique()
    if len(trait_col_values) != 1:
        raise ValueError(f"Expected exactly one unique 'Trait' in df, found: {trait_col_values}")
    trait_name = trait_col_values[0]

    results = {}

    drug_targets = set(
        drug_reference[
            (drug_reference["trait"] == trait_name) & (drug_reference["Sum"] >= targetthreshold)
        ]["EnsemblId"]
    )

    n = len(df)
    trait_result = {}

    for col in score_columns:
        if overlap_method == "topvalues":
            subset = df.sort_values(by=col, ascending=True).head(int(cutoff * n))
        elif overlap_method == "lessthan":
            subset = df[df[col] < cutoff * 100]
        else:
            raise ValueError("Invalid overlap_method: choose 'topvalues' or 'lessthan'")

        top_ids = set(subset["EnsemblId"])
        overlap = top_ids & drug_targets
        percent = 100 * len(overlap) / len(top_ids) if top_ids else 0.0

        trait_result[col] = round(percent,3)

        if printer:
            print(f"[{trait_name}] {col}: {percent:.2f}% overlap")

    results[trait_name] = trait_result
    return results


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
