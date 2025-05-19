from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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


def PCA(df):

    percentile_cols = [col for col in df.columns if "_" in col and col.split("_")[1] == "percentile"]

    filtered = df[percentile_cols].dropna()

    scaletool = StandardScaler()

    scaled_filtered = scaletool.fit_transform(filtered)

    pca = PCA(n_components = 2)

    df_pca = pca.fit_transform(scaled_filtered)

    filtered["Prioscore_PCA"] = df_pca[:,0]

    return filtered