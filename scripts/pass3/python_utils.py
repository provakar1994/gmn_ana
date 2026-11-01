import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_cut_scan_trend_light(
    df,
    ylim,
    title="",
    bg_models=("R1", "R2", "R3", "R5")  # default active models
):
    """
    Generalized cut scan plot.

    Parameters
    ----------
    df : pandas.DataFrame
    ylim : tuple
        (ymin, ymax) for first subplot.
    title : str
    bg_models : tuple/list
        Background models to include.
        Options: "R1","R2","R3","R4","R5","R6"
    """

    sns.set_context('talk')

    # Map model names to plotting properties
    bg_config = {
        "R1": {"color": "b",     "alpha": 0.8, "label": "Anti-dy bg"},
        "R2": {"color": "g",     "alpha": 0.7, "label": "Inel-MC bg"},
        "R3": {"color": "r",     "alpha": 1.0, "label": "Poly2 bg"},
        "R4": {"color": "m",     "alpha": 0.6, "label": "Gaussian bg"},
        "R5": {"color": "gold",  "alpha": 0.5, "label": "Poly3 bg"},
        "R6": {"color": "cyan",  "alpha": 0.4, "label": "Anti-dt bg"},
        "R7": {"color": "#E76F51",  "alpha": 0.3, "label": "Randoms+Inel-MC bg"},
    }

    fig, ax = plt.subplots(
        1, 1,
        figsize=(14, 7),
    )

    for model in bg_models:
        if model in bg_config:
            ax.errorbar(
                df['cut'],
                df[model],
                yerr=df[f"{model}err"],
                marker='o',
                linestyle='',
                color=bg_config[model]["color"],
                alpha=bg_config[model]["alpha"],
                label=bg_config[model]["label"]
            )

    ax.set_ylim(*ylim)
    ax.set_ylabel(r'$R_{n/p}^{sf}$')
    ax.grid(axis="both")
    ax.legend(loc='lower center')
    ax.set_title(title)

    plt.xticks(rotation=90)
    fig.tight_layout()
    plt.show()

def plot_cut_scan_trend(
    df,
    ylim1,
    ylim2=None,
    ylim3=None,
    title="",
    bg_models=("R1", "R2", "R3", "R5")  # default active models
):
    """
    Generalized W^2 cut scan plot.

    Parameters
    ----------
    df : pandas.DataFrame
    ylim1 : tuple
        (ymin, ymax) for first subplot.
    ylim2 : tuple, optional
        (ymin, ymax) for second subplot.
    ylim3 : tuple, optional
        (ymin, ymax) for third subplot.
    title : str
    bg_models : tuple/list
        Background models to include.
        Options: "R1","R2","R3","R4","R5","R6"
    """

    # Map model names to plotting properties
    bg_config = {
        "R1": {"color": "b",     "alpha": 0.8, "label": "Anti-dy bg"},
        "R2": {"color": "g",     "alpha": 0.7, "label": "Inel-MC bg"},
        "R3": {"color": "r",     "alpha": 1.0, "label": "Poly2 bg"},
        "R4": {"color": "m",     "alpha": 0.6, "label": "Gaussian bg"},
        "R5": {"color": "gold",  "alpha": 0.5, "label": "Poly3 bg"},
        "R6": {"color": "cyan",  "alpha": 0.4, "label": "Anti-dt bg"},
        "R7": {"color": "#E76F51",  "alpha": 0.3, "label": "Randoms+Inel-MC bg"},
    }

    fig, ax = plt.subplots(
        4, 1,
        figsize=(5.5, 7),
        dpi=100,
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1, 1, 1]}
    )

    # ---- First subplot ----
    for model in bg_models:
        if model in bg_config:
            ax[0].errorbar(
                df['cut'],
                df[model],
                yerr=df[f"{model}err"],
                marker='o',
                linestyle='',
                color=bg_config[model]["color"],
                alpha=bg_config[model]["alpha"],
                label=bg_config[model]["label"]
            )

    ax[0].set_ylim(*ylim1)
    ax[0].set_ylabel(r'$R_{n/p}^{sf}$')
    ax[0].grid(axis="both")
    ax[0].legend(loc='lower center')
    ax[0].set_title(title)

    # ---- Second subplot ----
    ax[1].errorbar(
        df['cut'], df['RMCnf'], yerr=df['RMCnferr'],
        marker='o', linestyle='', color="#2A9D8F", alpha=0.8,
        label='Pure MC Signal'
    )

    if ylim2 is not None:
        ax[1].set_ylim(*ylim2)

    ax[1].set_ylabel(r'$R_{n/p}^{MC}$')
    ax[1].grid(axis="both")

   # ---- Third subplot ----
    ax[2].scatter(
        df['cut'], df['chi27']/df['NDF7'],
        marker='o', linestyle='', color="#2A9D8F", alpha=0.8,
        label='Poly 2 Bg Fit'
    )

    if ylim3 is not None:
        ax[2].set_ylim(*ylim3)

    ax[2].set_ylabel(r'$\chi^2$/NDF')
    ax[2].grid(axis="both")

    # ---- Fourth subplot ----
    total = df['Yp3'] + df['Yn3'] + df['Ybg3']

    ax[3].plot(df['cut'], (df['Yp3'] + df['Yn3']) / total,
               marker='o', linestyle='', color="#264653", alpha=0.8, label='Sig')

    ax[3].plot(df['cut'], df['Ybg3'] / total,
               marker='o', linestyle='', color="#E76F51", alpha=0.8, label='Bg')

    ax[3].set_ylabel('Norm. Counts')
    ax[3].grid(axis="both")
    ax[3].legend()
 

    plt.xticks(rotation=90)
    fig.tight_layout()
    plt.show()

def plot_cut_scan_histo(
    df,
    nbins=10,
    title="",
    bg_model="R3"  # default active models
):
    """
    Generalized cut scan plot.

    Parameters
    ----------
    df : pandas.DataFrame
    nbins : int
        Number of bins for the histogram.
    title : str
    bg_model : str
        Background model to include.
        Options: "R1","R2","R3","R4","R5","R6"
    """

    sns.set_context('talk')

    df[bg_model].hist(
        bins=nbins,
        color='green',
        edgecolor='black',
        alpha=0.7
    )
    plt.title(title)
    plt.xlabel(r'$R_{n/p}^{sf}$')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.show()

    print(f"Mean: {df[bg_model].mean():.4f}")
    print(f"Std Dev: {df[bg_model].std():.4f}")