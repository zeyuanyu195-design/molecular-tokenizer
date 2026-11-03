"""BRICS versus GraphBPE distributions on the existing matched molecular cohort.

Layout follows DemoDiff Appendix B.2 Figure 9, arXiv:2510.08744v1.
Figures show this experiment's measurements, not digitized paper data.
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from safe.tokenizer import SAFETokenizer


def rows(path):
    with path.open(encoding="utf-8-sig") as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return dict(filename=path.name, bytes=path.stat().st_size, sha256=digest.hexdigest())


def collect(graph_path, brics_path, demo_path):
    graph, brics, demo = {}, {}, {}
    for path, target, method in (
        (graph_path, graph, "graph"),
        (brics_path, brics, "brics"),
        (demo_path, demo, "demo"),
    ):
        for row in rows(path):
            key = row["record_index"]
            if key in target:
                raise ValueError(f"Duplicate record {key} in {path}")
            item = dict(
                chembl_id=row["chembl_id"], smiles=row.get("canonical_smiles", row.get("smiles"))
            )
            if method == "demo":
                if "smiles" not in row or "safe" in row:
                    raise ValueError("Expected original-SMILES DemoDiff output")
                item.update(ok=row["status"] == "ok", count=row.get("num_graph_tokens"))
            else:
                prefix = "bpe" if method == "graph" else "brics"
                count = row[f"num_{prefix}_fragments"]
                if count != len(row[f"{prefix}_groups"]) or count < 1:
                    raise ValueError(f"Invalid fragment count at {key}")
                exact = row[f"{prefix}_boundary_exact"]
                if exact and sorted(map(sorted, row[f"{prefix}_groups"])) != sorted(
                    map(sorted, row["safe_groups"])
                ):
                    raise ValueError(f"Partition mismatch at {key}")
                if len(row["safe"].split(".")) != row["num_safe_fragments"]:
                    raise ValueError(f"SAFE fragment mismatch at {key}")
                if row["roundtrip_ok"] is not True:
                    raise ValueError(f"Unverified molecular round trip at {key}")
                item.update(
                    count=count, exact=exact, safe=row["safe"], safe_count=row["num_safe_fragments"]
                )
                if method == "graph" and count != len(row["bpe_node_ids"]):
                    raise ValueError(f"GraphBPE node count mismatch at {key}")
            target[key] = item
    keys = sorted(
        key
        for key in graph.keys() & brics.keys() & demo.keys()
        if graph[key]["exact"] and demo[key]["ok"]
    )
    if not keys:
        raise ValueError("Empty matched cohort")
    for key in keys:
        identities = {
            (source[key]["chembl_id"], source[key]["smiles"]) for source in (graph, brics, demo)
        }
        if len(identities) != 1 or graph[key]["count"] != demo[key]["count"]:
            raise ValueError(f"Matched-cohort identity/node-count mismatch at {key}")
    return keys, graph, brics


def describe(values):
    frequencies = Counter(map(int, values))
    return dict(
        n=len(values),
        total=int(sum(values)),
        mean=float(np.mean(values)),
        median=float(np.median(values)),
        p95=float(np.percentile(values, 95)),
        p99=float(np.percentile(values, 99)),
        maximum=int(max(values)),
        frequencies=[[x, frequencies[x]] for x in sorted(frequencies)],
    )


def plot_pair(series, title, xlabel, output):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "svg.fonttype": "none",
            "svg.hashsalt": "brics-graphbpe-20261103",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    colors = ["#2e86a6", "#a8427d"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.7), layout="constrained")
    fig.suptitle(title, fontsize=16)
    max_count = max(max(values) for _, values in series)
    # Identical log-spaced edges for both series, covering every observation.
    log_bins = np.geomspace(0.9, max_count + 0.5, 55)
    linear_bins = np.arange(0.5, max_count + 1.5, 1)
    percentiles = []
    for (name, values), color in zip(series, colors):
        values = np.asarray(values)
        n = len(values)
        x, counts = np.unique(values, return_counts=True)
        a = axs[0, 0]
        density, edges = np.histogram(values, bins=log_bins, density=True)
        assert np.isclose(np.sum(density * np.diff(edges)), 1)
        a.hist(
            values,
            bins=log_bins,
            density=True,
            alpha=0.60,
            color=color,
            edgecolor=color,
            linewidth=0.5,
            label=name,
        )
        a.axvline(
            np.median(values),
            color=color,
            ls="--",
            lw=1.6,
            label=f"{name} median: {np.median(values):g}",
        )
        tail = np.cumsum(counts[::-1])[::-1] / n
        assert tail[0] == 1 and np.all(np.diff(tail) <= 0)
        axs[0, 1].step(x, tail, where="pre", color=color, lw=1.7, label=name)
        ranked = np.sort(counts)[::-1]
        axs[1, 0].plot(
            np.arange(1, len(ranked) + 1),
            ranked,
            ".-",
            color=color,
            markersize=3,
            lw=1.2,
            label=name,
        )
        axs[1, 1].hist(values, bins=linear_bins, density=True, color=color, alpha=0.55, label=name)
        for p, style in ((95, ":"), (99, "--")):
            q = float(np.percentile(values, p))
            axs[1, 1].axvline(q, color=color, ls=style, lw=1.2)
        percentiles.append(
            f"{name}: P95={np.percentile(values,95):g}, " f"P99={np.percentile(values,99):g}"
        )
    axs[0, 0].set(
        xscale="log",
        yscale="log",
        xlabel=xlabel,
        ylabel="Density",
        title="(a) Log–log scale distribution",
    )
    axs[0, 1].set(
        xscale="log",
        yscale="log",
        xlabel=xlabel,
        ylabel=r"Tail probability $P(X \geq x)$",
        title="(b) Complementary CDF (log–log)",
    )
    axs[1, 0].set(
        xscale="log",
        yscale="log",
        xlabel="Rank of token-length frequency",
        ylabel="Molecules with that length",
        title="(c) Rank–frequency plot",
    )
    axs[1, 1].set(
        yscale="log",
        xlabel=xlabel,
        ylabel="Density (log scale)",
        title="(d) Linear–log with percentiles",
        xlim=(0, max_count + 1),
    )
    axs[1, 1].text(
        0.98,
        0.73,
        "\n".join(percentiles),
        transform=axs[1, 1].transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="none"),
    )
    for ax in axs.flat:
        ax.grid(alpha=0.20)
        ax.set_axisbelow(True)
        ax.legend(fontsize=8, framealpha=0.95)
    for ext in ("svg", "png", "pdf"):
        fig.savefig(output.with_suffix("." + ext), dpi=220, metadata={"Creator": "Matplotlib"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("graph", "brics", "demodiff", "text-tokenizer", "output-dir"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    keys, graph, brics = collect(args.graph, args.brics, args.demodiff)
    tokenizer = SAFETokenizer.load(str(args.text_tokenizer)).tokenizer
    tokenizer.no_padding()
    tokenizer.no_truncation()
    tokenizer.model.dropout = None
    fragment_series = [
        ("BRICS", [brics[k]["count"] for k in keys]),
        ("GraphBPE", [graph[k]["count"] for k in keys]),
    ]
    text_series = [
        (
            "BRICS+SAFE",
            [len(tokenizer.encode(brics[k]["safe"], add_special_tokens=False).ids) for k in keys],
        ),
        (
            "GraphBPE+SAFE",
            [len(tokenizer.encode(graph[k]["safe"], add_special_tokens=False).ids) for k in keys],
        ),
    ]
    report = dict(
        reference="https://arxiv.org/html/2510.08744v1",
        cohort="GraphBPE boundary-exact intersect BRICS successes and original-SMILES DemoDiff successes",
        matched_n=len(keys),
        cohort_key_sha256=hashlib.sha256("\n".join(map(str, keys)).encode()).hexdigest(),
        brics_changed_boundaries_in_cohort=sum(not brics[k]["exact"] for k in keys),
        brics_final_safe_counts=describe([brics[k]["safe_count"] for k in keys]),
        fragments={name: describe(vals) for name, vals in fragment_series},
        text={name: describe(vals) for name, vals in text_series},
        units=dict(
            fragments="One original BRICS/GraphBPE fragment per counting unit; no new BRICS ID vocabulary; connectivity excluded",
            text="Same existing SAFE text tokenizer; connectivity syntax included; no special tokens/padding/truncation/dropout",
        ),
        plotting=dict(
            log_histogram_bins=54,
            log_histogram_bin_edges=55,
            linear_histogram_bin_width=1,
            density="bin_count / (cohort_size * bin_width)",
            ccdf="P(X >= x)",
            rank="Rank of distinct length values by descending frequency, not motif-vocabulary frequency",
            percentiles="NumPy linear interpolation",
            tail="Full range retained; zero-height bins omitted by log axes",
        ),
        sources={
            role: fingerprint(path)
            for role, path in (
                ("graph", args.graph),
                ("brics", args.brics),
                ("demodiff", args.demodiff),
                ("text_tokenizer", args.text_tokenizer),
            )
        },
    )
    (args.output_dir / "token-distributions.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    for kind, series, unit in (
        ("fragments", fragment_series, "Number of fragment tokens"),
        ("text", text_series, "Number of SAFE text tokens"),
    ):
        plot_pair(
            series,
            f"BRICS vs GraphBPE — {len(keys):,} matched molecules\n"
            + (
                "Fragment counts (one fragment = one counting token)"
                if kind == "fragments"
                else "SAFE text lengths (same text tokenizer)"
            ),
            unit,
            args.output_dir / ("token-distributions-" + kind),
        )
    print(
        json.dumps(
            {
                "n": len(keys),
                "brics_boundary_changes": report["brics_changed_boundaries_in_cohort"],
                "stats": {
                    kind: {
                        name: {k: v for k, v in values.items() if k != "frequencies"}
                        for name, values in report[kind].items()
                    }
                    for kind in ("fragments", "text")
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
