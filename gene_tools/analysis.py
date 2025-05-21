
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




from scipy.stats import fisher_exact

def evaluate_prioritization(reference, mydf, score_cols = ["Prioscore_mean", "Prioscore_max", "Prioscore_min", "Prioscore_median"]
, sum_threshold=3, top_percent=0.01):
    """
    Compute odds ratios and Fisher's exact test p-values for enrichment of drug targets
    among top-ranked genes based on given prioritization scores.

    Parameters:
        reference (pd.DataFrame): DataFrame containing 'Sum' and 'EnsemblId' columns used to define drug targets.
        mydf (pd.DataFrame): DataFrame containing the same 'EnsemblId' and prioritization score columns.
        score_cols (list, optional): Columns to evaluate (default: mean, max, min, median prioritizations).
        sum_threshold (int, optional): Minimum 'Sum' to define drug targets.
        top_percent (float, optional): Fraction of top genes to consider (e.g., 0.01 for top 1%).

    Returns:
        Nothing (prints OR, p-value, and count of overlapping drug targets).
    """
    if score_cols is None:
        score_cols = ["Prioscore_mean", "Prioscore_max", "Prioscore_min", "Prioscore_median","PCA"]

    drug_targets = set(reference.loc[reference["Sum"] >= sum_threshold, "EnsemblId"])
    all_ids = set(mydf["EnsemblId"])

    for col in score_cols:
        ranked = mydf.sort_values(by=col, ascending=True)
        top = ranked.head(int(top_percent * len(mydf)))
        top_ids = set(top["EnsemblId"])
        rest_ids = all_ids - top_ids

        A = len(top_ids & drug_targets)
        B = len(top_ids - drug_targets)
        C = len(rest_ids & drug_targets)
        D = len(rest_ids - drug_targets)

        oddsratio, p_value = fisher_exact([[A, B], [C, D]])

        print(f"{col}: OR = {oddsratio:.2f}, p = {p_value:.4e}, drug targets in top {float(top_percent * 100)}% = {A}")
        


