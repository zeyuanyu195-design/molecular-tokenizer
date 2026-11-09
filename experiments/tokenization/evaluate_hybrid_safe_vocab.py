"""Fixed-budget whole SAFE fragments plus within-fragment BPE, fit on TRAIN only.

BPE operates on reversible private-use symbols for existing base token IDs.
Every training iterator item is ONE fragment, so no merge crosses a dot.
Select positive whole/subword mixtures by validation total length before test.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import shutil

import tokenizers
from tokenizers import Tokenizer, models, trainers
from evaluate_bounded_safe_vocab import Codec, decode, dump, load_base, length_stats, sha

OFFSET = 0xE000


def jsonl(path):
    with Path(path).open(encoding="utf-8-sig") as stream:
        yield from map(json.loads, stream)


def symbols(ids):
    return "".join(chr(OFFSET+i) for i in ids)


def symbol_payload(text, base):
    return b"".join(base[ord(char)-OFFSET] for char in text)


def build_variant(base_codec, ranked, pool, K, whole_fraction):
    """Share duplicate payloads; budget counts MODEL IDs, not BPE internal IDs."""
    M = base_codec.initial_size
    whole_count = int((K-M)*whole_fraction)
    payloads = list(base_codec.payloads)
    by_payload = {}
    for i, payload in enumerate(payloads):
        if payload is not None:
            by_payload.setdefault(payload, i)
    whole = {}
    for fragment in ranked[:whole_count]:
        payload = fragment.encode("utf-8")
        if payload in by_payload:
            raise ValueError("Whole-fragment ranking contains an existing payload")
        whole[fragment] = len(payloads)
        by_payload[payload] = len(payloads)
        payloads.append(payload)
    internal_vocab, external_ids = {}, []
    for i, payload in enumerate(base_codec.payloads):
        if payload is not None:
            internal_vocab[chr(OFFSET+i)] = len(internal_vocab)
            external_ids.append(i)
    merges, substring_count = [], 0
    substring_budget = K-len(payloads)
    for left, right in pool:
        if substring_count == substring_budget:
            break
        combined = left+right
        assert left in internal_vocab and right in internal_vocab
        payload = symbol_payload(combined, base_codec.payloads)
        if b"." in payload:
            raise ValueError("Learned a substring crossing a dot")
        if payload not in by_payload:
            by_payload[payload] = len(payloads)
            payloads.append(payload)
            substring_count += 1
        if combined not in internal_vocab:
            internal_vocab[combined] = len(internal_vocab)
            external_ids.append(by_payload[payload])
        merges.append((left, right))
    if len(payloads) != K:
        raise ValueError(f"Insufficient learned substrings: {len(payloads)} != {K}")
    tokenizer = Tokenizer(models.BPE(vocab=internal_vocab, merges=merges, dropout=None))
    config = {"K": K, "actual_size": len(payloads), "M": M,
        "whole_fraction_of_additional_budget": whole_fraction, "whole_entries": len(whole),
        "substring_entries": substring_count, "internal_merge_count": len(merges),
        "whole_fragment_ids": whole, "internal_to_external_ids": external_ids,
        "payloads_hex": [None if p is None else p.hex() for p in payloads],
        "symbol_offset": OFFSET, "dot_id": base_codec.dot}
    return HybridCodec(base_codec, tokenizer, config)


class HybridCodec:
    def __init__(self, base, tokenizer, config):
        self.base = base
        self.tokenizer = tokenizer
        self.config = config
        self.payloads = [None if p is None else bytes.fromhex(p) for p in config["payloads_hex"]]
        self.whole = config["whole_fragment_ids"]
        self.external = config["internal_to_external_ids"]
        self.cache = {}

    def fragment(self, text):
        if text in self.cache:
            return self.cache[text]
        if text in self.whole:
            ids = [self.whole[text]]
        elif text in self.base.base_by_text:
            ids = [self.base.base_by_text[text]]
        else:
            base_ids = self.base.fallback(text)[0]
            encoded = self.tokenizer.encode(symbols(base_ids), add_special_tokens=False)
            ids = [self.external[i] for i in encoded.ids]
        if decode(ids, self.payloads) != text:
            raise ValueError("Fragment reconstruction failed")
        self.cache[text] = ids
        return ids

    def encode(self, text):
        parts = text.split(".")
        if not text or not all(parts):
            raise ValueError("Empty SAFE fragment")
        ids = []
        for i, part in enumerate(parts):
            if i:
                ids.append(self.base.dot)
            ids.extend(self.fragment(part))
        return ids

    def save(self, folder):
        folder.mkdir(parents=True)
        dump(folder / "vocab.json", self.config)
        self.tokenizer.save(str(folder / "substring_bpe.json"))


def load_hybrid(folder, base):
    config = json.loads((Path(folder) / "vocab.json").read_text(encoding="utf-8"))
    assert config["actual_size"] == len(config["payloads_hex"]) == config["K"]
    return HybridCodec(base, Tokenizer.from_file(str(Path(folder)/"substring_bpe.json")), config)


def evaluate(codec, rows, old_lengths, ids_path=None):
    stats, examples = {}, []
    stream = ids_path.open("x", encoding="utf-8") if ids_path else None
    try:
        for name, entries in rows.items():
            lengths, full_hits, one_token, fragments, unchanged_or_better = [], 0, 0, 0, 0
            for row in entries:
                ids = codec.encode(row["safe"])
                if decode(ids, codec.payloads) != row["safe"]:
                    raise ValueError(f"SAFE roundtrip failed: {row['record_index']}")
                lengths.append(len(ids))
                parts = row["safe"].split(".")
                fragments += len(parts)
                full_hits += sum(p in codec.whole for p in parts)
                one_token += sum(len(codec.fragment(p)) == 1 for p in parts)
                unchanged_or_better += len(ids) <= old_lengths[row["record_index"]]
                if stream:
                    stream.write(json.dumps({"record_index": row["record_index"], "cohort": name,
                        "length": len(ids), "whole_only_length": old_lengths[row["record_index"]], "ids": ids})+"\n")
                if row["record_index"] == 39224:
                    examples.append({"record_index": row["record_index"], "chembl_id": row["chembl_id"],
                        "safe": row["safe"], "length": len(ids), "ids": ids,
                        "fragments": [{"text": p, "whole_entry": p in codec.whole,
                            "ids": codec.fragment(p), "tokens": [codec.payloads[i].decode() for i in codec.fragment(p)]}
                            for p in parts]})
            old_total = sum(old_lengths[r["record_index"]] for r in entries)
            stats[name] = dict(length_stats(lengths),
                whole_only_total=old_total, reduction_vs_whole_only=1-sum(lengths)/old_total,
                exact_safe_roundtrips=len(entries), fragment_occurrences=fragments,
                explicit_whole_entry_hits=full_hits, one_token_fragment_coverage=one_token/fragments,
                molecules_no_longer_than_whole_only=unchanged_or_better)
    finally:
        if stream:
            stream.close()
    return stats, examples


def run(args):
    if args.output.exists():
        raise ValueError("Output exists; select a fresh directory")
    parent = args.experiment
    original = json.loads((parent/"summary.json").read_text(encoding="utf-8"))
    assert sha(args.input) == original["sources"][str(args.input)]
    tok, base_payloads = load_base(parent/"base_tokenizer.json", parent/"source_splitter.py")
    base = Codec(tok, base_payloads)
    splits = {r["record_index"]: r["cohort"] for r in jsonl(parent/"splits.jsonl")}
    rows = defaultdict(list)
    counts = Counter()
    for record in jsonl(args.input):
        cohort = splits[record["record_index"]]
        rows[cohort].append(record)
        if cohort == "train":
            counts.update(record["safe"].split("."))
    ranked = sorted(((f,n) for f,n in counts.items() if f not in base.base_by_text), key=lambda x:(-x[1],x[0]))
    saved_ranked = [(r["fragment"], r["train_occurrences"]) for r in jsonl(parent/"training_fragment_ranking.jsonl")]
    assert ranked == saved_ranked
    ranked = [f for f,n in ranked]
    args.output.mkdir(parents=True)
    for name in ["base_tokenizer.json", "source_splitter.py", "splits.jsonl"]:
        shutil.copyfile(parent/name, args.output/name)
    dump(args.output/"train_fragment_counts.json", dict(sorted(counts.items())))
    print(f"Training BPE on {sum(counts.values())} fragment occurrences from {len(rows['train'])} TRAIN records", flush=True)
    alphabet = [chr(OFFSET+i) for i,p in enumerate(base.payloads) if p is not None]
    trainer = trainers.BpeTrainer(vocab_size=len(alphabet)+2*(max(args.budgets)-base.initial_size),
        min_frequency=1, show_progress=False, initial_alphabet=alphabet)
    trained = Tokenizer(models.BPE())
    def corpus():
        for fragment, n in sorted(counts.items()):
            text = symbols(base.fallback(fragment)[0])
            for _ in range(n):
                yield text
    trained.train_from_iterator(corpus(), trainer=trainer, length=sum(counts.values()))
    trained.save(str(args.output/"trained_substring_pool.json"))
    pool = json.loads(trained.to_str())["model"]["merges"]
    print(f"Trained {len(pool)} merges; assessing allocations on TRAIN and VALIDATION only", flush=True)
    summary = {"M": base.initial_size, "budgets": [], "selected": {},
        "method": "TRAIN-only BPE over base-token ID symbols, no cross-fragment merges; exact whole entry takes priority",
        "pool_training": "All primary TRAIN fragment occurrences; same ordered merge pool for every allocation",
        "selection": "Minimum validation total tokens among 25%, 50%, 75% whole-entry allocations; ties prefer fewer whole entries",
        "test_note": "Existing test set was explored in previous diagnosis; not a pristine final holdout. No test data used to fit pool or select allocation in this run.",
        "base_provenance_note": original["split_note"], "versions": {"tokenizers": tokenizers.__version__},
        "sources": {str(p):sha(p) for p in [args.input, parent/"summary.json", parent/"splits.jsonl", Path(__file__)]}}
    chosen = {}
    for K in args.budgets:
        previous = {r["record_index"]:r["length"] for r in jsonl(parent/f"K{K}"/"lengths.jsonl")}
        candidates = []
        for ratio in [0.0, 0.25, 0.5, 0.75, 1.0]:
            name = f"K{K}_whole{int(ratio*100)}"
            folder = args.output/name
            codec = build_variant(base, ranked, pool, K, ratio)
            codec.save(folder)
            stats,_ = evaluate(codec, {n: rows[n] for n in ["train","validation"]}, previous)
            if ratio == 1:
                assert all(s["total"] == s["whole_only_total"] for s in stats.values())
            result = {"name": name, "K":K, "actual_size":len(codec.payloads), "whole_fraction":ratio,
                "whole_entries":len(codec.whole), "substring_entries":codec.config["substring_entries"],
                "internal_merge_count":codec.config["internal_merge_count"], "cohorts":stats}
            dump(folder/"metrics.json", result)
            summary["budgets"].append(result)
            if 0 < ratio < 1:
                candidates.append(result)
            print(json.dumps({"name":name,"validation_mean":stats['validation']['mean'],
                "train_mean":stats['train']['mean']}), flush=True)
        winner = min(candidates, key=lambda r:(r["cohorts"]["validation"]["total"],r["whole_entries"]))
        chosen[str(K)] = winner["name"]
    # Persist the selection BEFORE any test-set codec evaluation.
    dump(args.output/"selection_before_test.json", {"rule":summary["selection"], "selected":chosen})
    summary["selected"] = chosen
    for result in summary["budgets"]:
        K,name = result["K"],result["name"]
        if name != chosen[str(K)] and result["whole_fraction"] != 0:
            continue
        previous = {r["record_index"]:r["length"] for r in jsonl(parent/f"K{K}"/"lengths.jsonl")}
        folder = args.output/name
        codec = load_hybrid(folder,base)
        selected_rows = rows if name == chosen[str(K)] else {"test": rows["test"]}
        stats,examples = evaluate(codec, selected_rows, previous, folder/"encoded_ids.jsonl")
        result["cohorts"].update(stats)
        dump(folder/"metrics.json",result)
        dump(folder/"CHEMBL45183.json",examples)
        print(json.dumps({"evaluated":name,"test":stats['test']}),flush=True)
    dump(args.output/"summary.json",summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--experiment",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--budgets",type=int,nargs="+",default=[16000,32000])
    run(parser.parse_args())
