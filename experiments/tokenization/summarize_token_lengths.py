"""Recompute full-cohort and matched-cohort statistics from saved experiment outputs.

The DemoDiff input MUST be export_demodiff_tokens.py output from original SMILES.
This script does not re-run GraphBPE or silently substitute the legacy SAFE run.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from safe.tokenizer import SAFETokenizer


def read_records(path):
    records = {}
    with Path(path).open(encoding="utf-8-sig") as stream:
        for line in stream:
            if not line.strip():
                continue
            row = json.loads(line)
            key = row["record_index"]
            if key in records:
                raise ValueError(f"Duplicate record_index {key} in {path}")
            records[key] = row
    return records


def same_molecule(left, right):
    if (left["chembl_id"], left["canonical_smiles"]) != (
        right["chembl_id"],
        right.get("canonical_smiles", right.get("smiles")),
    ):
        raise ValueError(f"Molecule mismatch for record {left['record_index']}")


def statistics(values):
    values = list(values)
    if not values:
        raise ValueError("Empty cohort")
    return dict(
        n=len(values),
        total=sum(values),
        mean=float(np.mean(values)),
        median=float(np.median(values)),
        p95=float(np.percentile(values, 95)),
        maximum=max(values),
    )


def fingerprint(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return dict(filename=path.name, bytes=path.stat().st_size, sha256=digest.hexdigest())


def summarize(graph_path, brics_path, demodiff_path, tokenizer_path):
    graph, brics, demo = map(read_records, (graph_path, brics_path, demodiff_path))
    successful = {}
    for key, row in demo.items():
        if "smiles" not in row or "safe" in row:
            raise ValueError("Expected original-SMILES DemoDiff exporter output")
        if row["status"] == "ok":
            if row["num_graph_tokens"] != len(row["node_ids"]):
                raise ValueError(f"DemoDiff token count mismatch: {key}")
            successful[key] = row
        elif row["status"] != "error":
            raise ValueError(f"Unknown DemoDiff status: {key}")
    exact = {}
    for key, row in graph.items():
        if type(row["bpe_boundary_exact"]) is not bool:
            raise ValueError(f"Invalid boundary flag: {key}")
        if not row["bpe_boundary_exact"]:
            continue
        if not (
            row["roundtrip_ok"] is True
            and len(row["bpe_node_ids"])
            == len(row["bpe_groups"])
            == len(row["safe_groups"])
            == len(row["safe"].split("."))
            == row["num_bpe_fragments"]
            == row["num_safe_fragments"]
            and sorted(map(sorted, row["bpe_groups"])) == sorted(map(sorted, row["safe_groups"]))
        ):
            raise ValueError(f"Invalid exact-boundary record: {key}")
        exact[key] = row
    for key in graph.keys() & successful.keys():
        same_molecule(graph[key], successful[key])
    for key in graph.keys() & brics.keys():
        same_molecule(graph[key], brics[key])
    paired = sorted(exact.keys() & brics.keys() & successful.keys())
    # SAFE JSON includes a custom pre-tokenizer: using Tokenizer.from_file
    # or removing its metadata would change the experiment's token boundaries.
    tokenizer = SAFETokenizer.load(str(tokenizer_path)).tokenizer
    tokenizer.no_padding()
    tokenizer.no_truncation()
    tokenizer.model.dropout = None

    def lengths(records):
        return {
            key: len(tokenizer.encode(row["safe"], add_special_tokens=False).ids)
            for key, row in records.items()
        }

    brics_lengths, graph_lengths = lengths(brics), lengths(graph)
    full = {
        "brics_safe_text": statistics(brics_lengths.values()),
        "demodiff_original_smiles_nodes": statistics(
            row["num_graph_tokens"] for row in successful.values()
        ),
        "graphbpe_safe_exact_motifs": statistics(
            len(row["bpe_node_ids"]) for row in exact.values()
        ),
        "graphbpe_safe_text_control": statistics(graph_lengths.values()),
    }
    matched = {
        "brics_safe_text": statistics(brics_lengths[key] for key in paired),
        "demodiff_original_smiles_nodes": statistics(
            successful[key]["num_graph_tokens"] for key in paired
        ),
        "graphbpe_safe_exact_motifs": statistics(len(exact[key]["bpe_node_ids"]) for key in paired),
        "graphbpe_safe_text_control": statistics(graph_lengths[key] for key in paired),
    }
    return dict(
        full=full,
        matched=matched,
        counts=dict(
            demodiff_processed=len(demo),
            demodiff_errors=len(demo) - len(successful),
            demodiff_rows_with_unknown_nodes=sum(
                row["unknown_node_count"] > 0 for row in successful.values()
            ),
            graph_records=len(graph),
            boundary_changed=len(graph) - len(exact),
            exact_records=len(exact),
            matched_records=len(paired),
        ),
        matched_length_disagreements=sum(
            successful[key]["num_graph_tokens"] != len(exact[key]["bpe_node_ids"]) for key in paired
        ),
        exact_excluded_from_matched=[
            dict(record_index=key, chembl_id=exact[key]["chembl_id"])
            for key in sorted(exact.keys() - set(paired))
        ],
        errors=[row for row in demo.values() if row["status"] == "error"],
        units=dict(
            brics_safe_text="Text IDs including connectivity syntax; no special tokens",
            demodiff_original_smiles_nodes="Node IDs only; edge/attachment arrays excluded",
            graphbpe_safe_exact_motifs="Recorded motif IDs; connectivity syntax excluded",
            graphbpe_safe_text_control="Same text tokenizer as BRICS",
        ),
        quantile_method="NumPy percentile, linear interpolation",
        sources={
            name: fingerprint(path)
            for name, path in zip(
                ("graph", "brics", "demodiff", "text_tokenizer"),
                (graph_path, brics_path, demodiff_path, tokenizer_path),
            )
        },
        note="Length agreement does not prove lossless reconstruction from motif IDs alone.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("graph", "brics", "demodiff", "text-tokenizer", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists; choose a new path")
    result = summarize(args.graph, args.brics, args.demodiff, args.text_tokenizer)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
    print(json.dumps(result["counts"], indent=2))


if __name__ == "__main__":
    main()
