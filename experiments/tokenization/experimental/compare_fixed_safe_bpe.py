"""Compare existing BRICS/GraphBPE SAFE using one shared, dot-bounded text BPE.

No molecular graph parsing, DemoDiff calls, or changes to the SAFE strings.
Train on both representations of training molecules; evaluate held-out pairs.
"""

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

import tokenizers
from tokenizers import Tokenizer, Regex, models, pre_tokenizers, decoders, trainers


def read_records(path):
    result = {}
    with path.open(encoding="utf-8-sig") as f:
        for line_no, line in enumerate(f, 1):
            row = json.loads(line)
            key = (row["record_index"], row["chembl_id"])
            if key in result:
                raise ValueError(f"Duplicate record key at {path}:{line_no}: {key}")
            safe = row["safe"]
            smiles = row["canonical_smiles"]
            if not isinstance(safe, str) or not safe or not all(safe.split(".")):
                raise ValueError(f"Invalid or empty SAFE at {path}:{line_no}")
            if not isinstance(smiles, str) or not smiles:
                raise ValueError(f"Missing canonical_smiles at {path}:{line_no}")
            result[key] = {"safe": safe, "canonical_smiles": smiles}
    return result


def split_name(smiles, seed, fraction):
    digest = hashlib.sha256(f"{seed}\0{smiles}".encode("utf-8")).digest()
    value = int.from_bytes(digest[:8], "big") / 2**64
    return "train" if value < fraction else "test"


def make_tokenizer():
    tokenizer = Tokenizer(models.BPE(dropout=None))
    tokenizer.pre_tokenizer = pre_tokenizers.Sequence(
        [
            pre_tokenizers.Split(Regex(r"\."), behavior="isolated"),
            pre_tokenizers.ByteLevel(add_prefix_space=False, use_regex=False),
        ]
    )
    tokenizer.decoder = decoders.ByteLevel()
    tokenizer.no_padding()
    tokenizer.no_truncation()
    return tokenizer


def validate_encoding(text, enc, tokenizer):
    if tokenizer.decode(enc.ids, skip_special_tokens=False) != text:
        raise ValueError("SAFE text round trip failed")
    dot_count = 0
    for start, end in enc.offsets:
        segment = text[start:end]
        if "." in segment:
            if segment != ".":
                raise ValueError("A token crossed a SAFE fragment boundary")
            dot_count += 1
    if dot_count != text.count("."):
        raise ValueError("Fragment separators were not preserved individually")


def length_stats(values):
    values = sorted(values)
    return {
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "p95_nearest_rank": values[math.ceil(0.95 * len(values)) - 1],
        "min": values[0],
        "max": values[-1],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--graphbpe", type=Path, required=True)
    p.add_argument("--brics", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--vocab-size", type=int, default=8000)
    p.add_argument("--min-frequency", type=int, default=2)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--train-fraction", type=float, default=0.8)
    p.add_argument("--max-pairs", type=int, help="Limit aligned pairs for smoke tests only")
    p.add_argument("--save-ids", action="store_true")
    args = p.parse_args()
    if args.vocab_size < 256 or args.min_frequency < 1 or not 0 < args.train_fraction < 1:
        p.error("Require vocab-size >= 256, min-frequency >= 1, and 0 < train-fraction < 1")
    if args.max_pairs is not None and args.max_pairs < 1:
        p.error("max-pairs must be positive")
    if args.out_dir.exists():
        p.error("Output directory already exists; choose a new directory")
    graph = read_records(args.graphbpe)
    brics = read_records(args.brics)
    common_keys = sorted(graph.keys() & brics.keys())
    keys = common_keys[: args.max_pairs] if args.max_pairs is not None else common_keys
    splits = {"train": [], "test": []}
    for key in keys:
        if graph[key]["canonical_smiles"] != brics[key]["canonical_smiles"]:
            raise ValueError(f"Original molecule mismatch for {key}")
        name = split_name(graph[key]["canonical_smiles"], args.seed, args.train_fraction)
        splits[name].append(key)
    if not splits["train"] or not splits["test"]:
        p.error("Need nonempty training and test partitions; increase max-pairs or check inputs")
    print(
        f"Aligned {len(keys)} pairs; train={len(splits['train'])}; test={len(splits['test'])}",
        flush=True,
    )
    tokenizer = make_tokenizer()
    trainer = trainers.BpeTrainer(
        vocab_size=args.vocab_size,
        min_frequency=args.min_frequency,
        initial_alphabet=sorted(pre_tokenizers.ByteLevel.alphabet()),
        special_tokens=[],
        show_progress=True,
    )

    def corpus():
        for key in splits["train"]:
            yield graph[key]["safe"]
            yield brics[key]["safe"]

    tokenizer.train_from_iterator(corpus(), trainer=trainer, length=2 * len(splits["train"]))
    args.out_dir.mkdir(parents=True, exist_ok=False)
    tokenizer_path = args.out_dir / "tokenizer.json"
    tokenizer.save(str(tokenizer_path))
    # Verify the serialized tokenizer, which is what subsequent runs will load.
    tokenizer = Tokenizer.from_file(str(tokenizer_path))
    results = {}
    with (args.out_dir / "paired_results.jsonl").open("x", encoding="utf-8") as fout:
        for split, split_keys in splits.items():
            lengths = {"graphbpe": [], "brics": []}
            fragment_counts = {"graphbpe": [], "brics": []}
            for offset in range(0, len(split_keys), 128):
                batch_keys = split_keys[offset : offset + 128]
                texts = [data[k]["safe"] for k in batch_keys for data in (graph, brics)]
                encodings = tokenizer.encode_batch(texts, add_special_tokens=False)
                for index, key in enumerate(batch_keys):
                    item = {"record_index": key[0], "chembl_id": key[1], "split": split}
                    for j, method in enumerate(("graphbpe", "brics")):
                        text, enc = texts[2 * index + j], encodings[2 * index + j]
                        validate_encoding(text, enc, tokenizer)
                        n_fragments = text.count(".") + 1
                        lengths[method].append(len(enc.ids))
                        fragment_counts[method].append(n_fragments)
                        item[method + "_tokens"] = len(enc.ids)
                        item[method + "_fragments"] = n_fragments
                        item[method + "_separator_tokens"] = n_fragments - 1
                        if args.save_ids:
                            item[method + "_ids"] = enc.ids
                    fout.write(json.dumps(item, ensure_ascii=False) + "\n")
            g, b = lengths["graphbpe"], lengths["brics"]
            results[split] = {
                "paired_molecules": len(split_keys),
                "graphbpe": length_stats(g),
                "brics": length_stats(b),
                "mean_graphbpe_fragments": statistics.mean(fragment_counts["graphbpe"]),
                "mean_brics_fragments": statistics.mean(fragment_counts["brics"]),
                "mean_graphbpe_minus_brics": statistics.mean(x - y for x, y in zip(g, b)),
                "graphbpe_reduction_percent_ratio_of_means": 100 * (sum(b) - sum(g)) / sum(b),
                "graphbpe_shorter_fraction": sum(x < y for x, y in zip(g, b)) / len(g),
                "equal_length_fraction": sum(x == y for x, y in zip(g, b)) / len(g),
            }
            print(f"Verified {split}: {len(split_keys)} paired SAFE strings per method", flush=True)
    report = {
        "inputs": {"graphbpe": str(args.graphbpe.resolve()), "brics": str(args.brics.resolve())},
        "available_records": {
            "graphbpe": len(graph),
            "brics": len(brics),
            "common": len(common_keys),
        },
        "requested_vocab_size": args.vocab_size,
        "actual_vocab_size": tokenizer.get_vocab_size(),
        "min_frequency": args.min_frequency,
        "seed": args.seed,
        "train_fraction_target": args.train_fraction,
        "split_rule": "SHA256(seed + NUL + exact canonical_smiles); identical source SMILES stay together. Not a scaffold split.",
        "train_corpus": "Both SAFE representations for each training molecule; occurrence frequencies are retained.",
        "tokenizers_version": tokenizers.__version__,
        "tokenizer_sha256": hashlib.sha256(tokenizer_path.read_bytes()).hexdigest(),
        "max_pairs": args.max_pairs,
        "statistics": results,
        "counting": "Counts include every '.' separator; no special tokens, padding, or truncation. Merges cannot cross '.' boundaries.",
        "verification": "All evaluated SAFE strings round-trip exactly and all token offsets respect '.' boundaries. Molecular chemistry is not revalidated.",
        "status": "complete",
    }
    with (args.out_dir / "summary.json").open("x", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(json.dumps(report["statistics"]["test"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
