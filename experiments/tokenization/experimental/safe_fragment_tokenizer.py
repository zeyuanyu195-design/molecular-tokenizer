"""One dot-separated SAFE fragment per token; Python standard library only.

literal: IDs alone reconstruct SAFE. rings: IDs plus ring_labels reconstruct SAFE.
The rings mode normalizes ring-label spelling, not arbitrary molecular graph order.
"""

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

LEXEME = re.compile(r"\[[^\]]*\]|%\(\d+\)|%\d{2}|\d")
PLACEHOLDER = re.compile(r"\{(\d+)\}")


def split_safe(text):
    if not isinstance(text, str) or not text or text.strip() != text:
        raise ValueError("SAFE must be a nonempty string without outer whitespace")
    parts = text.split(".")
    if not all(parts):
        raise ValueError("Empty SAFE fragment")
    if "{" in text or "}" in text:
        raise ValueError("Input contains reserved template characters")
    return parts


def fragment_key(fragment, mode):
    if mode == "literal":
        return fragment, []
    label_index = {}

    def replace(match):
        label = match.group()
        if label.startswith("["):
            return label  # Preserve isotopes, atom maps, H counts, and charges.
        if label not in label_index:
            label_index[label] = len(label_index)
        return "{" + str(label_index[label]) + "}"

    return LEXEME.sub(replace, fragment), list(label_index)


def records(path):
    with Path(path).open(encoding="utf-8-sig") as f:
        for number, line in enumerate(f, 1):
            if not line.strip():
                raise ValueError(f"{path}:{number}: empty line")
            row = (
                json.loads(line)
                if Path(path).suffix.lower() == ".jsonl"
                else {"safe": line.rstrip("\r\n")}
            )
            split_safe(row["safe"])
            yield number, row


def encode(text, vocab, lookup):
    ids, labels = [], []
    for fragment in split_safe(text):
        key, side = fragment_key(fragment, vocab["mode"])
        if key not in lookup:
            raise ValueError(f"Unseen fragment: {fragment!r}; no UNK substitution is performed")
        ids.append(lookup[key])
        labels.append(side)
    result = {"input_ids": ids, "num_tokens": len(ids)}
    if vocab["mode"] == "rings":
        result["ring_labels"] = labels
    return result


def decode(encoded, vocab):
    ids = encoded["input_ids"]
    labels = encoded.get("ring_labels", [[] for _ in ids])
    if len(labels) != len(ids):
        raise ValueError("One ring_labels list is required per token")
    parts = []
    for token_id, side in zip(ids, labels):
        if not isinstance(token_id, int) or not 0 <= token_id < len(vocab["tokens"]):
            raise ValueError(f"Invalid token ID: {token_id}")
        part = vocab["tokens"][token_id]
        if vocab["mode"] == "rings":
            part = PLACEHOLDER.sub(lambda m: side[int(m.group(1))], part)
        parts.append(part)
    return ".".join(parts)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as f:
        json.dump(value, f, ensure_ascii=False, indent=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    fit = sub.add_parser("fit", help="Fit on training data only for ML evaluation")
    fit.add_argument("--input", nargs="+", required=True)
    fit.add_argument("--vocab", required=True)
    fit.add_argument("--mode", choices=("literal", "rings"), default="literal")
    enc = sub.add_parser("encode")
    enc.add_argument("--input", required=True)
    enc.add_argument("--vocab", required=True)
    enc.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "fit":
        if Path(args.vocab).exists():
            parser.error("Vocabulary exists; choose a new path")
        tokens, count = set(), 0
        for source in args.input:
            for _, row in records(source):
                tokens.update(fragment_key(f, args.mode)[0] for f in split_safe(row["safe"]))
                count += 1
        if not tokens:
            parser.error("No training fragments")
        vocab = {
            "format_version": 1,
            "mode": args.mode,
            "tokens": sorted(tokens),
            "special_tokens": [],
            "training_records": count,
            "training_sources": [str(Path(p).resolve()) for p in args.input],
        }
        write_json(args.vocab, vocab)
        print(json.dumps({"vocab_size": len(tokens), "training_records": count, "mode": args.mode}))
        return
    vocab_bytes = Path(args.vocab).read_bytes()
    vocab = json.loads(vocab_bytes)
    if vocab.get("format_version") != 1 or vocab.get("mode") not in ("literal", "rings"):
        parser.error("Unsupported vocabulary format")
    lookup = {token: i for i, token in enumerate(vocab["tokens"])}
    if len(lookup) != len(vocab["tokens"]):
        parser.error("Duplicate vocabulary entries")
    output = Path(args.output)
    summary_path = Path(str(output) + ".summary.json")
    if output.exists() or summary_path.exists():
        parser.error("Output exists; choose a new output path")
    # Preflight catches unknown fragments before creating any output.
    total = 0
    for number, row in records(args.input):
        try:
            item = encode(row["safe"], vocab, lookup)
            if decode(item, vocab) != row["safe"]:
                raise ValueError("Exact SAFE round trip failed")
        except Exception as exc:
            raise ValueError(f"Input line {number}: {exc}") from exc
        total += 1
    output.parent.mkdir(parents=True, exist_ok=True)
    sizes = Counter()
    with output.open("x", encoding="utf-8") as f:
        for number, row in records(args.input):
            item = encode(row["safe"], vocab, lookup)
            item.update(
                record_index=row.get("record_index", number), chembl_id=row.get("chembl_id")
            )
            f.write(json.dumps(item, ensure_ascii=False) + "\n")
            sizes[item["num_tokens"]] += 1
    summary = {
        "records": total,
        "vocab_size": len(lookup),
        "mode": vocab["mode"],
        "vocab_sha256": hashlib.sha256(vocab_bytes).hexdigest(),
        "exact_safe_roundtrip_checked": total,
        "mean_fragment_tokens": sum(k * v for k, v in sizes.items()) / total if total else None,
        "special_tokens_added": False,
        "padding": False,
        "truncation": False,
        "note": "Token count is fragment count; rings mode additionally requires ring_labels. This check reconstructs SAFE text, not original-SMILES stereochemistry.",
    }
    write_json(summary_path, summary)
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
