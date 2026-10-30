"""Export DemoDiff graph tokens from original SMILES (TSV or records.jsonl).

Uses the upstream API at commit 0acb85ef446c175a3e3a60bbad43de97d97ca76f.
No training; fixed vocabulary; one output record per input record, including errors.
"""

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from itertools import islice
from pathlib import Path


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        if path.suffix.lower() == ".jsonl":
            for line in f:
                if line.strip():
                    yield json.loads(line)
        else:
            yield from csv.DictReader(f, delimiter="\t")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--demodiff-root", type=Path, required=True)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--smiles-column", default="canonical_smiles")
    p.add_argument("--vocab-prefix", type=Path)
    p.add_argument("--limit", type=int)
    a = p.parse_args()
    if a.limit is not None and a.limit < 1:
        p.error("--limit must be positive")
    root = a.demodiff_root.resolve()
    source = root / "downstream" / "graphbpe.py"
    if not source.is_file():
        p.error(f"Missing upstream source: {source}")
    prefix = (a.vocab_prefix or root / "data/tokenizer/vocab3000ring300/pretrain-token").resolve()
    vocab_hashes = {}
    for suffix in (".node", ".edge", ".ring"):
        filename = Path(str(prefix) + suffix)
        vocab_hashes[suffix] = hashlib.sha256(filename.read_bytes()).hexdigest()
    sys.path.insert(0, str(root / "downstream"))
    from graphbpe import MolecularGraphTokenizer

    tokenizer = MolecularGraphTokenizer(kekulize=True, name="pretrain", simple_mode=False)
    tokenizer.load(str(prefix))
    counts = Counter()
    summary_path = Path(str(a.output) + ".summary.json")
    if a.output.exists() or summary_path.exists():
        p.error("Output already exists; choose a new output path")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("x", encoding="utf-8") as out:
        for index, row in enumerate(islice(rows(a.input), a.limit), 1):
            result = {
                "record_index": row.get("record_index", index),
                "chembl_id": row.get("chembl_id"),
                "smiles": row.get(a.smiles_column),
            }
            try:
                smiles = result["smiles"]
                if not isinstance(smiles, str) or not smiles.strip():
                    raise ValueError(f"Empty or missing {a.smiles_column}")
                ids, adj = tokenizer.encode(smiles, update_vocab_edge=False)
                bonds, positions = tokenizer.get_bond_position_by_vocab(adj)
                unknown = sum(i >= len(tokenizer.vocab_node) for i in ids)
                result.update(
                    status="ok",
                    node_ids=ids,
                    num_graph_tokens=len(ids),
                    adj_edge_ids=adj,
                    bond_adj=bonds,
                    position_adj=positions,
                    unknown_node_count=unknown,
                )
                counts["encoded"] += 1
                counts["graph_tokens"] += len(ids)
                counts["rows_with_unknown_nodes"] += bool(unknown)
            except Exception as exc:
                result.update(status="error", error=f"{type(exc).__name__}: {exc}")
                counts["errors"] += 1
            out.write(json.dumps(result, ensure_ascii=False) + "\n")
            counts["processed"] += 1
            if index % 1000 == 0:
                print(f"Processed {index}; errors={counts['errors']}", flush=True)
    summary = {
        "input": str(a.input.resolve()),
        "smiles_column": a.smiles_column,
        "vocab_prefix": str(prefix),
        "vocab_sha256": vocab_hashes,
        "graphbpe_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "counts": dict(counts),
        "mean_graph_tokens": (
            counts["graph_tokens"] / counts["encoded"] if counts["encoded"] else None
        ),
        "note": "Graph nodes, not SAFE text tokens. Encoding success does not establish lossless molecular reconstruction. Missing edges follow upstream frozen-vocabulary behavior; validate round trips before downstream use.",
    }
    with summary_path.open("x", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
