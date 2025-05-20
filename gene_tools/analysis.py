
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


from Scoring.scomet import Mean, Max, Min, Median, pca

Max(ldl_clean["LDL"])
Mean(ldl_clean["LDL"])
Median(ldl_clean["LDL"])
Min(ldl_clean["LDL"])


    

drug_targets = set(df.loc[df["Sum"] >= 3, "EnsemblId"])
score_cols = ["Prioscore_mean", "Prioscore_max", "Prioscore_min", "Prioscore_median"]

for col in score_cols:
    top = ldl_clean["LDL"].sort_values(by=col, ascending=True).head(int(0.01 * len(ldl_clean["LDL"])))
    top_ids = set(top["EnsemblId"])
    overlap = top_ids & drug_targets
    percent = 100 * len(overlap) / len(top_ids)
    print(f"{col}: {percent:.2f}% of top 1% are drug targets")




