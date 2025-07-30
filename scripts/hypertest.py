# Import libraries
import argparse
import sys
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import to_rgb

# Define hypergeometric test function
def hypergeom_test(state, class_label, matrix_df):
    M = matrix_df[class_label].sum() + matrix_df.drop(columns=[class_label]).sum().sum()
    n = matrix_df[class_label].sum()
    N = matrix_df.loc[state].sum()
    x = matrix_df.loc[state, class_label]
    if N == 0 or n == 0:
        return np.nan
    return hypergeom.sf(x - 1, M, N, n)

# Parse arguments
def main():
    parser = argparse.ArgumentParser(description="Run hypergeometric tests on a chromatin state matrix.")
    parser.add_argument("path", help="Path to matrix .tsv file")
    parser.add_argument("--output", help="Output file for raw p-values (.tsv)", default=None)
    parser.add_argument("--plot", action="store_true", help="Display heatmap of -log10(p-values)")
    parser.add_argument("--plot-output", help="Optional file to save the heatmap (e.g., 'plot.png')")
    parser.add_argument("--rangecap", help="Upper limit for -log(p-value) in the heatmap color scale; cells containing higher values will appear with the same color", default=15)
    parser.add_argument("--fdr", action="store_true", help="Apply FDR correction (Benjamini-Hochberg)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Load matrix
    matrix_df = pd.read_csv(args.path, sep="\t", index_col=0)

    # Perform hypergeometric test
    results = []
    for state in matrix_df.index:
        for class_label in matrix_df.columns:
            pval = hypergeom_test(state, class_label, matrix_df)
            results.append((state, class_label, pval))

    pval_df = pd.DataFrame(results, columns=["State", "Class", "PValue"])

    # Apply FDR correction if requested
    if args.fdr:
        pvals = pval_df["PValue"].values
        _, fdr_corrected, _, _ = multipletests(pvals, method='fdr_bh')
        pval_df["AdjPValue"] = fdr_corrected
        value_column = "AdjPValue"
    else:
        value_column = "PValue"

    # Save to file
    if args.output:
        pval_df.to_csv(args.output, sep="\t", index=False)
        print(f"Results saved to {args.output}")
    else:
        print(pval_df)

    # Plot if requested
    if args.plot:
        pivot = pval_df.pivot(index="State", columns="Class", values=value_column)

        min_pval = 1e-300  # avoid log10(0)
        log_p = -np.log10(pivot.clip(lower=min_pval))  # Contains the values that will be shown in the cells

        # Dynamic range to show color scale in significant p-value cells.
        min_display = 0
        max_display = args.rangecap
        norm = plt.Normalize(vmin=min_display, vmax=max_display)


        # Create figure
        fig = plt.figure(figsize=(15, 8))

        # Add background gradient
        ax_bg = fig.add_axes([0, 0, 1, 1], zorder=0)
        ax_bg.axis("off")
        top_rgb = np.array(to_rgb("palegoldenrod"))
        bottom_rgb = np.array(to_rgb("midnightblue"))
        height = 1000
        gradient = np.linspace(0, 1, height)[:, None]
        gradient_rgb = bottom_rgb + (top_rgb - bottom_rgb) * (1 - gradient)
        gradient_rgb = np.tile(gradient_rgb[:, None, :], (1, 1000, 1))
        ax_bg.imshow(
            gradient_rgb,
            aspect="auto",
            extent=[0, 1, 0, 1],
            transform=fig.transFigure,
            zorder=0,
            origin="lower"
        )

        # Main heatmap
        ax = fig.add_subplot(111, zorder=1)
        cmap = plt.get_cmap("viridis")
        norm = plt.Normalize(log_p.min().min(), log_p.max().max())

        heatmap = sns.heatmap(
            log_p,
            cmap=cmap,
            annot=False,
            fmt=".1f",
            linewidths=0.5,
            cbar_kws={"label": f"-log10({value_column})"},
            ax=ax,
            norm=norm  # Aplica el rango de color entre 1.3 y 10
        )

        # Manual annotation with contrast-aware text
        for y in range(log_p.shape[0]):
            for x in range(log_p.shape[1]):
                val = log_p.iloc[y, x]
                if pd.isna(val):
                    text = "NaN"
                    text_color = "black"
                else:
                    bg_color = cmap(norm(val))[:3]
                    brightness = 0.299 * bg_color[0] + 0.587 * bg_color[1] + 0.114 * bg_color[2]
                    text_color = 'black' if brightness > 0.6 else 'white'
                    text = f"{val:.1f}"
                ax.text(x + 0.5, y + 0.5, text, ha='center', va='center', color=text_color)

        # Y-axis label colors based on chromatin state
        state_colors_rgb = {
            "TssA": (255, 0, 0),
            "TssFlnk": (255, 69, 0),
            "TssFlnkU": (255, 69, 0),
            "TssFlnkD": (255, 69, 0),
            "Tx": (0, 128, 0),
            "TxWk": (0, 100, 0),
            "EnhG1": (194, 225, 5),
            "EnhG2": (194, 225, 5),
            "EnhA1": (255, 195, 77),
            "EnhA2": (255, 195, 77),
            "EnhWk": (255, 255, 0),
            "ZNF/Rpts": (102, 205, 170),
            "Het": (138, 145, 208),
            "TssBiv": (205, 92, 92),
            "EnhBiv": (189, 183, 107),
            "ReprPC": (128, 128, 128),
            "ReprPCWk": (192, 192, 192),
            "Quies": (255, 255, 255),
        }

        state_colors_hex = {
            k: '#%02x%02x%02x' % rgb for k, rgb in state_colors_rgb.items()
        }

        for label in ax.get_yticklabels():
            text = label.get_text()
            if text in state_colors_hex:
                label.set_color(state_colors_hex[text])

        plt.title(f"State Enrichment: -log10({value_column})", color='white')
        plt.tight_layout()

        if args.plot_output:
            plt.savefig(args.plot_output, dpi=300, bbox_inches='tight')
            print(f"Heatmap saved to {args.plot_output}")
        else:
            plt.show()

if __name__ == "__main__":
    main()