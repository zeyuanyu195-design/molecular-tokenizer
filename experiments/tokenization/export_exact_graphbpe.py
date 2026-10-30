"""Export existing SAFE records whose boundaries exactly match GraphBPE.

Reuse the recorded DemoDiff-compatible node IDs; never re-run segmentation.
IDs retain bpe_groups order, NOT necessarily SAFE dot-fragment order.
"""

import argparse
import json
from pathlib import Path


def export_records(source, destination):
    source = Path(source).resolve()
    destination = Path(destination).resolve()
    # Refuse existing destinations so original results cannot be overwritten.
    destination.mkdir(parents=True, exist_ok=False)
    counts = dict(input_records=0, kept=0, excluded_boundary_changed=0, total_node_tokens=0)
    with source.open(encoding="utf-8-sig") as fin, (destination / "records.jsonl").open(
        "x", encoding="utf-8"
    ) as fout, (destination / "safe.txt").open("x", encoding="utf-8") as corpus:
        for line in fin:
            if not line.strip():
                continue
            row = json.loads(line)
            counts["input_records"] += 1
            if "bpe_boundary_exact" not in row:
                raise ValueError("Input must be the original GraphBPE records.jsonl")
            if row["bpe_boundary_exact"] is not True:
                counts["excluded_boundary_changed"] += 1
                continue
            ids, groups = row["bpe_node_ids"], row["bpe_groups"]
            safe_groups = row["safe_groups"]
            if not ids or not all(type(i) is int and i >= 0 for i in ids):
                raise ValueError("Invalid recorded node IDs")
            if not (
                len(ids)
                == len(groups)
                == len(safe_groups)
                == len(row["safe"].split("."))
                == row["num_bpe_fragments"]
                == row["num_safe_fragments"]
            ):
                raise ValueError("Fragment count mismatch")
            if sorted(map(sorted, groups)) != sorted(map(sorted, safe_groups)):
                raise ValueError("GraphBPE and SAFE partitions differ")
            if row.get("roundtrip_ok") is not True:
                raise ValueError("Input record has not passed the original round-trip check")
            row.update(
                node_ids=ids,
                num_node_tokens=len(ids),
                node_id_order="bpe_groups (not SAFE serialization order)",
                tokenization_source="recorded GraphBPE node-vocabulary lookup",
                ids_are_lossless=False,
            )
            fout.write(json.dumps(row, ensure_ascii=False) + "\n")
            corpus.write(row["safe"] + "\n")
            counts["kept"] += 1
            counts["total_node_tokens"] += len(ids)
    summary = dict(
        status="complete",
        input=str(source),
        **counts,
        mean_node_tokens=(counts["total_node_tokens"] / counts["kept"] if counts["kept"] else None),
        node_id_order="bpe_groups (not SAFE serialization order)",
        includes_connectivity_tokens=False
    )
    source_summary = source.parent / "summary.json"
    if source_summary.is_file():
        metadata = json.loads(source_summary.read_text(encoding="utf-8-sig"))
        summary["source_vocab_sha256"] = metadata.get("vocab_sha256")
    with (destination / "summary.json").open("x", encoding="utf-8") as out:
        json.dump(summary, out, ensure_ascii=False, indent=2)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(export_records(args.input, args.output_dir), ensure_ascii=False, indent=2))
