"""Frequency-ranked whole SAFE fragments with lossless base/byte fallback.

Train only on exact-boundary training molecules; evaluate held-out molecules.
Vocabulary IDs decode directly to SAFE bytes, with no original SAFE side table.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics

from rdkit import Chem, rdBase
from tokenizers import Tokenizer
from tokenizers.pre_tokenizers import PreTokenizer


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def split_group(key, seed):
    fraction = int.from_bytes(hashlib.sha256(f"{seed}\0{key}".encode()).digest()[:8], "big") / 2**64
    return "train" if fraction < .8 else "validation" if fraction < .9 else "test"


def molecule_identity(text, stereo=True):
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise ValueError("Invalid complete molecule")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=stereo)


def load_base(path, splitter_path):
    config = json.loads(Path(path).read_text())
    use_splitter = config.pop("custom_pre_tokenizer", False)
    config.pop("tokenizer_attrs", None)
    config.pop("tokenizer_type", None)
    tok = Tokenizer.from_str(json.dumps(config))
    if use_splitter:
        spec = importlib.util.spec_from_file_location("source_safe_splitter", splitter_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        tok.pre_tokenizer = PreTokenizer.custom(mod.SAFESplitter())
    tok.no_padding()
    tok.no_truncation()
    tok.model.dropout = None
    special = {x["id"] for x in config["added_tokens"] if x["special"]}
    payloads = [None if i in special else tok.decode([i], skip_special_tokens=False).encode()
                for i in range(tok.get_vocab_size())]
    return tok, payloads


def decode(ids, payloads):
    if any(not isinstance(i, int) or i < 0 or i >= len(payloads) for i in ids):
        raise ValueError("Invalid vocabulary ID")
    values = [payloads[i] for i in ids]
    if any(v is None for v in values):
        raise ValueError("Special model token encountered in molecular content")
    return b"".join(values).decode("utf-8")


def load_vocabulary(path):
    """Reload a persisted decoder without chemistry, input records or an edge table."""
    path = Path(path)
    config = json.loads(path.read_text(encoding="utf-8"))
    base = json.loads((path.parent / config["base_payloads_file"]).read_text(encoding="utf-8"))
    payloads = [None if p is None else bytes.fromhex(p) for p in base]
    if len(payloads) != config["M"]:
        raise ValueError("Incorrect base vocabulary size")
    payloads.extend(f.encode("utf-8") for f in config["added_fragments"])
    if len(payloads) != config["actual_size"] or len(payloads) > config["K"]:
        raise ValueError("Incorrect vocabulary budget")
    return payloads


class Codec:
    def __init__(self, tok, base_payloads):
        self.tok = tok
        self.base_size = len(base_payloads)
        self.payloads = list(base_payloads) + [bytes([i]) for i in range(256)]
        self.initial_size = len(self.payloads)
        self.base_by_text = {p.decode(): i for i, p in enumerate(base_payloads)
                             if p is not None and p}
        self.dot = self.base_by_text["."]
        self.cache = {}

    def fallback(self, fragment):
        if fragment in self.cache:
            return self.cache[fragment]
        byte_fallback = False
        try:
            ids = self.tok.encode(fragment, add_special_tokens=False).ids
            if decode(ids, self.payloads) != fragment:
                raise ValueError("Base token bytes are not an exact reconstruction")
        except Exception:
            ids = [self.base_size + b for b in fragment.encode("utf-8")]
            byte_fallback = True
        assert decode(ids, self.payloads) == fragment
        self.cache[fragment] = (ids, byte_fallback)
        return ids, byte_fallback

    def encode(self, safe, selected):
        fragments = safe.split(".")
        if not safe or not all(fragments):
            raise ValueError("Empty SAFE fragment")
        ids, added_hits, single_hits, byte_fragments = [], 0, 0, 0
        for n, fragment in enumerate(fragments):
            if n:
                ids.append(self.dot)
            if fragment in selected:
                part = [selected[fragment]]
                added_hits += 1
            elif fragment in self.base_by_text:
                part = [self.base_by_text[fragment]]
            else:
                part, used_bytes = self.fallback(fragment)
                byte_fragments += used_bytes
            single_hits += len(part) == 1
            ids.extend(part)
        return ids, added_hits, single_hits, byte_fragments


def length_stats(values):
    values = sorted(values)
    if not values:
        return {"n": 0}
    return {"n": len(values), "total": sum(values), "mean": statistics.mean(values),
            "median": statistics.median(values), "p95": values[math.ceil(.95 * len(values)) - 1],
            "min": values[0], "max": values[-1]}


def run(args):
    if args.output.exists():
        raise ValueError("Output exists; select a fresh directory")
    tok, payloads = load_base(args.base_tokenizer, args.splitter)
    codec = Codec(tok, payloads)
    rows, partitions, groups = [], defaultdict(list), defaultdict(set)
    for i, line in enumerate(args.input.open(encoding="utf-8-sig")):
        if args.limit and i >= args.limit:
            break
        r = json.loads(line)
        safe, original = r["safe"], r["canonical_smiles"]
        if molecule_identity(safe) != molecule_identity(original):
            raise ValueError(f"Molecular roundtrip failed at {r['record_index']}")
        group = molecule_identity(original, stereo=False)
        partition = split_group(group, args.seed)
        cohort = partition if r["bpe_boundary_exact"] else "boundary_changed_stress"
        groups[partition].add(group)
        row = dict(record_index=r["record_index"], chembl_id=r["chembl_id"], safe=safe,
                   canonical_smiles=original, group=group, split=partition, cohort=cohort,
                   num_fragments=len(safe.split(".")))
        rows.append(row)
        partitions[cohort].append(row)
        if (i + 1) % 10000 == 0:
            print(f"Validated molecular identity: {i + 1}", flush=True)
    assert not (groups["train"] & groups["test"] or groups["train"] & groups["validation"]
                or groups["test"] & groups["validation"])
    if any(not partitions[s] for s in ["train", "validation", "test"]):
        raise ValueError("Empty primary split")
    # Only exact-boundary TRAIN examples select new vocabulary entries.
    frequency = Counter(f for r in partitions["train"] for f in r["safe"].split("."))
    eligible = [(f, n) for f, n in frequency.items() if f not in codec.base_by_text]
    eligible.sort(key=lambda x: (-x[1], x[0]))
    args.output.mkdir(parents=True)
    dump(args.output / "base_payloads.json", [None if b is None else b.hex() for b in codec.payloads])
    (args.output / "base_tokenizer.json").write_bytes(args.base_tokenizer.read_bytes())
    (args.output / "source_splitter.py").write_bytes(args.splitter.read_bytes())
    with (args.output / "splits.jsonl").open("x", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps({k: v for k, v in row.items() if k != "safe"}) + "\n")
    with (args.output / "training_fragment_ranking.jsonl").open("x", encoding="utf-8") as f:
        for fragment, count in eligible:
            f.write(json.dumps({"fragment": fragment, "train_occurrences": count}) + "\n")
    for row in rows:
        ids, _, _, used_bytes = codec.encode(row["safe"], {})
        assert decode(ids, codec.payloads) == row["safe"]
        row["baseline_length"] = len(ids)
        row["baseline_byte_fragments"] = used_bytes
        try:
            base_ids = tok.encode(row["safe"], add_special_tokens=False).ids
            row["original_base_lossless"] = tok.decode(base_ids, skip_special_tokens=False) == row["safe"]
            row["original_base_length"] = len(base_ids)
        except Exception:
            row["original_base_lossless"] = False
            row["original_base_length"] = None
    summary = {"base_vocab_size": codec.base_size, "byte_fallback_entries": 256,
               "initial_vocab_M": codec.initial_size, "input_records": len(rows),
               "split_seed": args.seed, "split_policy": "sha256(nonisomeric canonical SMILES), 80/10/10",
               "split_note": "Identity-disjoint, not scaffold-disjoint; pre-existing base vocabulary training provenance unknown.",
               "training_scope": "Only boundary-exact TRAIN records; all changed-boundary rows are a separate stress set.",
               "train_unique_fragments": len(frequency), "eligible_new_fragments": len(eligible),
               "length_scope": "Includes explicit dots and all label/bond/stereo and fallback tokens; no BOS/EOS, padding or truncation.",
               "decoder": "Concatenate vocabulary byte payloads, decode UTF-8; no input-specific side information.",
               "ranking": "Training occurrence frequency descending, lexical tie break; no test-based selection.",
               "versions": {"rdkit": rdBase.rdkitVersion},
               "sources": {str(p): sha(p) for p in [args.input, args.base_tokenizer, args.splitter, Path(__file__)]},
               "cohorts": {}, "budgets": []}
    for name, entries in partitions.items():
        stress_groups = {r["group"] for r in entries}
        summary["cohorts"][name] = {
            "records": len(entries), "unique_nonisomeric_graphs": len(stress_groups),
            "fragment_occurrences": sum(r["num_fragments"] for r in entries),
            "baseline": length_stats([r["baseline_length"] for r in entries]),
            "original_base": length_stats([r["original_base_length"] for r in entries if r["original_base_lossless"]]),
            "original_base_lossless_records": sum(r["original_base_lossless"] for r in entries),
            "byte_fallback_fragment_occurrences": sum(r["baseline_byte_fragments"] for r in entries),
            "stress_overlap_with_primary_train_groups": len(stress_groups & {r["group"] for r in partitions["train"]}) if name=="boundary_changed_stress" else None,
        }
    for K in sorted({codec.initial_size, *args.budgets}):
        if K < codec.initial_size:
            raise ValueError(f"Budget {K} below M={codec.initial_size}")
        chosen = eligible[:K - codec.initial_size]
        selected = {text: codec.initial_size + i for i, (text, _) in enumerate(chosen)}
        vocab = codec.payloads + [text.encode() for text, _ in chosen]
        folder = args.output / f"K{K}"
        folder.mkdir()
        dump(folder / "vocab.json", {"K": K, "actual_size": len(vocab), "M": codec.initial_size,
                                   "base_payloads_file": "../base_payloads.json", "added_fragments": [f for f, _ in chosen]})
        stats, examples = {}, []
        with (folder / "lengths.jsonl").open("x", encoding="utf-8") as f, \
             (folder / "test_ids.jsonl").open("x", encoding="utf-8") as ids_file:
            for name, entries in partitions.items():
                lengths, add_hits, single_hits, byte_count, rt = [], 0, 0, 0, 0
                for row in entries:
                    ids, ah, sh, bc = codec.encode(row["safe"], selected)
                    restored = decode(ids, vocab)
                    if restored != row["safe"]:
                        raise ValueError(f"Exact string roundtrip failed: {K}/{row['record_index']}")
                    rt += 1
                    lengths.append(len(ids)); add_hits += ah; single_hits += sh; byte_count += bc
                    record = {"record_index": row["record_index"], "cohort": name,
                              "baseline_length": row["baseline_length"], "length": len(ids),
                              "new_fragment_hits": ah, "one_token_fragments": sh}
                    f.write(json.dumps(record) + "\n")
                    if name == "test":
                        ids_file.write(json.dumps({"record_index": row["record_index"], "ids": ids}) + "\n")
                        if len(examples) < 3:
                            examples.append(dict(record, safe=row["safe"], ids=ids))
                total_fragments = sum(r["num_fragments"] for r in entries)
                baseline_total = sum(r["baseline_length"] for r in entries)
                stats[name] = dict(length_stats(lengths),
                    length_reduction=1-sum(lengths)/baseline_total,
                    new_fragment_occurrence_coverage=add_hits/total_fragments,
                    one_token_fragment_coverage=single_hits/total_fragments,
                    byte_fallback_fragment_occurrences=byte_count,
                    exact_safe_roundtrips=rt,
                    molecular_identity_verified_by_exact_SAFE_and_input_check=rt)
        item = {"K": K, "actual_size": len(vocab), "new_entries": len(selected), "cohorts": stats}
        summary["budgets"].append(item)
        dump(folder / "metrics.json", item)
        dump(folder / "examples.json", examples)
        print(json.dumps({"K": K, "actual_size": len(vocab), "test": stats["test"]}), flush=True)
    dump(args.output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--base-tokenizer", type=Path, required=True)
    p.add_argument("--splitter", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--budgets", type=int, nargs="+", default=[4000, 8000, 16000, 32000, 64000])
    p.add_argument("--limit", type=int)
    run(p.parse_args())
