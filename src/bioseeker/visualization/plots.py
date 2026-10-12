"""Visualization of conservation analysis."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_conservation_regression(
    bicodon_df: pd.DataFrame,
    slope: float,
    r_squared: float,
    sigma: float,
    output_path: str,
) -> None:
    """Plot Observed vs Expected conservation rates with regression and outlier bands."""
    fig, ax = plt.subplots(figsize=(10, 8))

    x = bicodon_df["ExpectedRate"].to_numpy(dtype=float)
    y = bicodon_df["ConservationRate"].to_numpy(dtype=float)
    outlier = bicodon_df["Outlier"].to_numpy()

    normal = outlier == ""
    high = outlier == "high"
    low = outlier == "low"

    ax.scatter(x[normal], y[normal], color="gray", alpha=0.5, label="Normal", zorder=2)
    if high.any():
        ax.scatter(x[high], y[high], color="magenta", alpha=0.8, label="High outliers", zorder=3)
    if low.any():
        ax.scatter(x[low], y[low], color="blue", alpha=0.8, label="Low outliers", zorder=3)

    max_val = max(np.nanmax(x), np.nanmax(y)) if len(x) > 0 else 1.0
    line_x = np.array([0, max_val])
    line_y = slope * line_x

    ax.plot(
        line_x, line_y, color="black", label=f"Fit (m={slope:.2f}, R²={r_squared:.3f})", zorder=1
    )

    if np.isfinite(sigma) and sigma > 0:
        ax.fill_between(
            line_x,
            line_y - 3 * sigma,
            line_y + 3 * sigma,
            color="gray",
            alpha=0.15,
            label="±3σ band",
            zorder=0,
        )

    ax.plot(line_x, line_x, "k--", alpha=0.3, label="Expected=Observed", zorder=1)

    ax.set_xlabel("Expected Conservation Rate")
    ax.set_ylabel("Observed Conservation Rate")
    ax.set_title("Bicodon Conservation: Observed vs Expected")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
