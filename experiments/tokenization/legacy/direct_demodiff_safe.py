"""Pass complete SAFE strings directly to the original DemoDiff encode() API.

Accepts records.jsonl (the 'safe' field) or plain text (one SAFE per line).
No SAFE decoding, fragment splitting, vocabulary fitting, or tokenizer edits.
"""

import argparse
import hashlib
import json
import sys
from collections import Counter
from itertools import islice
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demodiff-root", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be positive")
    root = args.demodiff_root.resolve()
    source = root / "downstream" / "graphbpe.py"
    prefix = root / "data" / "tokenizer" / "vocab3000ring300" / "pretrain-token"
    required = [source, root / "downstream/utils/token.py"]
    required += [Path(str(prefix) + suffix) for suffix in (".node", ".edge", ".ring")]
    for path in required:
        if not path.is_file():
            parser.error(f"Required DemoDiff file missing: {path}")
    summary_path = Path(str(args.output) + ".summary.json")
    if args.output.exists() or summary_path.exists():
        parser.error("Output exists; choose a new output filename")
    sys.path.insert(0, str(root / "downstream"))
    from graphbpe import MolecularGraphTokenizer

    tokenizer = MolecularGraphTokenizer(kekulize=True, name="pretrain", simple_mode=False)
    tokenizer.load(str(prefix))
    count = Counter(processed=0, success=0, failed=0, total_node_tokens=0, unknown_node_tokens=0)
    vocab_size = len(tokenizer.vocab_node)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(encoding="utf-8-sig") as fin, args.output.open(
        "x", encoding="utf-8"
    ) as fout:
        for line_number, line in enumerate(islice(fin, args.limit), 1):
            result = {"input_line": line_number}
            try:
                if args.input.suffix.lower() == ".jsonl":
                    row = json.loads(line)
                    safe_text = row["safe"]  # Deliberately do not read canonical_smiles.
                    result.update(
                        record_index=row.get("record_index"), chembl_id=row.get("chembl_id")
                    )
                else:
                    safe_text = line.rstrip("\r\n")
                result["safe"] = safe_text
                if not isinstance(safe_text, str) or not safe_text.strip():
                    raise ValueError("Missing or empty SAFE string")

                # Direct upstream call, matching the fixed-vocabulary downstream usage.
                node_ids, adj = tokenizer.encode(safe_text, update_vocab_edge=False)

                unknown_count = sum(i >= vocab_size for i in node_ids)
                result.update(
                    status="ok",
                    node_ids=node_ids,
                    adj_edge_ids=adj,
                    num_node_tokens=len(node_ids),
                    unknown_node_count=unknown_count,
                )
                count["success"] += 1
                count["total_node_tokens"] += len(node_ids)
                count["unknown_node_tokens"] += unknown_count
            except Exception as exc:
                result.update(status="error", error=f"{type(exc).__name__}: {exc}")
                count["failed"] += 1
            fout.write(json.dumps(result, ensure_ascii=False) + "\n")
            count["processed"] += 1
            if line_number % 1000 == 0:
                fout.flush()
                print(
                    f"Processed {line_number}; success={count['success']}; failed={count['failed']}",
                    flush=True,
                )
    summary = {
        "input": str(args.input.resolve()),
        "input_field": "safe" if args.input.suffix.lower() == ".jsonl" else "whole_line",
        "demodiff_root": str(root),
        "update_vocab_edge": False,
        "node_vocab_size": vocab_size,
        "edge_vocab_size": len(tokenizer.vocab_edge),
        "counts": dict(count),
        "mean_node_tokens": (
            count["total_node_tokens"] / count["success"] if count["success"] else None
        ),
        "file_sha256": {
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in required
        },
        "note": "Success means upstream encode() returned; no molecular round-trip check was performed. Upstream fixed-edge-vocabulary behavior is unchanged.",
    }
    with summary_path.open("x", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
