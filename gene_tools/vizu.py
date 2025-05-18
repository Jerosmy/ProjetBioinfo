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