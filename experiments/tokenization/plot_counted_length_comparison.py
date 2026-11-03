"""Plot the table's BRICS SAFE-text counts versus GraphBPE fixed-motif counts.

This is a comparison of different counting units, not a full-encoding
compression measurement. Motif counts exclude attachment/connectivity syntax.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = json.loads(args.input_report.read_text(encoding="utf-8"))
    specs = [
        ("BRICS+SAFE: text tokens", report["text"]["BRICS+SAFE"], "#2e86a6"),
        ("GraphBPE+SAFE: motif IDs", report["fragments"]["GraphBPE"], "#a8427d"),
    ]
    series = []
    for name, summary, color in specs:
        frequencies = np.asarray(summary["frequencies"], dtype=int)
        values = np.repeat(frequencies[:, 0], frequencies[:, 1])
        if len(values) != report["matched_n"] or int(values.sum()) != summary["total"]:
            raise ValueError("Frequency totals do not match the reported cohort")
        if not np.isclose(values.mean(), summary["mean"]):
            raise ValueError("Mean does not match source summary")
        series.append((name, values, color))
    bins = np.geomspace(0.9, max(v.max() for _, v, _ in series) + 0.5, 55)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, ax = plt.subplots(figsize=(10, 6.2))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.84, bottom=0.21)
    fig.suptitle("BRICS+SAFE text tokens vs GraphBPE+SAFE motif IDs", fontsize=15, y=0.97)
    ax.set_title(f"{report['matched_n']:,} matched molecules", fontsize=12, pad=12)
    for name, values, color in series:
        density, edges = np.histogram(values, bins=bins, density=True)
        assert np.isclose(np.sum(density * np.diff(edges)), 1)
        ax.hist(
            values,
            bins=bins,
            density=True,
            color=color,
            alpha=0.58,
            edgecolor=color,
            linewidth=0.5,
            label=f"{name} (mean {values.mean():.2f})",
        )
        ax.axvline(values.mean(), color=color, ls="--", lw=1.8)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Number of counted tokens / IDs per molecule",
        ylabel="Probability density",
    )
    ax.grid(alpha=0.18)
    ax.set_axisbelow(True)
    ax.legend(fontsize=9, framealpha=0.96)
    fig.text(
        0.5,
        0.105,
        "Dashed lines: mean counted length. Full length range retained.",
        ha="center",
        fontsize=10,
    )
    fig.text(
        0.5,
        0.035,
        "Different units: text counts include connectivity; motif counts exclude it.\n"
        "This is not a complete-encoding compression ratio.",
        ha="center",
        fontsize=10,
        color="#6b4331",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)
    print(
        json.dumps(
            {
                name: {"n": len(v), "total": int(v.sum()), "mean": float(v.mean())}
                for name, v, _ in series
            }
        )
    )


if __name__ == "__main__":
    main()
