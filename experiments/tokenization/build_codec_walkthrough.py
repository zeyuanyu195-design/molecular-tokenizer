"""Attach real, frozen whole/hybrid vocabulary traces to verified graph examples.

No training or graph repartitioning happens here. The browser receives only the
used dictionary entries (with hashes of the full vocabularies), not invented IDs.
"""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem
import safe as safe_library
import safe.converter as safe_converter_module
from evaluate_bounded_safe_vocab import Codec, decode, load_base, load_vocabulary
from evaluate_hybrid_safe_vocab import load_hybrid


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(text):
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise ValueError("Invalid complete molecular representation")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def build(graphs, whole_dir, hybrid_dir):
    data = read(graphs)
    tok, payloads = load_base(whole_dir / "base_tokenizer.json", whole_dir / "source_splitter.py")
    base = Codec(tok, payloads)
    data["schema_version"] = 2
    data["codec_metadata"] = {
        "policy": "Frozen TRAIN vocabularies; whole-fragment lookup first; fallback never crosses dots",
        "graph_trace_sha256": digest(graphs),
        "base_tokenizer_sha256": digest(whole_dir / "base_tokenizer.json"),
        "splitter_sha256": digest(whole_dir / "source_splitter.py"),
        "safe_decoder_module_sha256": digest(safe_converter_module.__file__),
        "budgets": {},
    }
    for budget, name in [(16000, "K16000_whole50"), (32000, "K32000_whole25")]:
        whole_path = whole_dir / f"K{budget}" / "vocab.json"
        whole = read(whole_path)
        whole_payloads = load_vocabulary(whole_path)
        selected = {s: whole["M"] + i for i, s in enumerate(whole["added_fragments"])}
        hybrid = load_hybrid(hybrid_dir / name, base)
        dictionary = {"whole": {}, "hybrid": {}}
        data["codec_metadata"]["budgets"][str(budget)] = {
            "whole_sha256": digest(whole_path),
            "whole_base_payloads_sha256": digest(whole_dir / "base_payloads.json"),
            "hybrid_sha256": digest(hybrid_dir / name / "vocab.json"),
            "hybrid_bpe_sha256": digest(hybrid_dir / name / "substring_bpe.json"),
            "base_entries": base.initial_size,
            "whole_entries": len(selected),
            "hybrid_whole_entries": len(hybrid.whole),
            "hybrid_substring_entries": hybrid.config["substring_entries"],
            "dictionary": dictionary,
        }
        for example in data["examples"]:
            safe = example["safe"]
            assert identity(safe) == identity(example["smiles"])
            assert len(example["fragments"]) == len(safe.split("."))
            result = {}
            for scheme, table, decoder in [("whole", selected, whole_payloads),
                                            ("hybrid", hybrid.whole, hybrid.payloads)]:
                parts, full_ids, offset = [], [], 0
                for fragment_index, fragment in enumerate(safe.split(".")):
                    direct_id = table.get(fragment, base.base_by_text.get(fragment))
                    base_ids, byte_fallback = base.fallback(fragment)
                    ids = ([direct_id] if direct_id is not None else base_ids) if scheme == "whole" else hybrid.fragment(fragment)
                    mode = ("whole_hit" if fragment in table else "base_hit") if direct_id is not None else ("base_fallback" if scheme == "whole" else "substring_fallback")
                    assert decode(ids, decoder) == fragment
                    if direct_id is not None:
                        assert ids == [direct_id], "A complete hit must not be split"
                    tokens, cursor = [], offset
                    for token_id in ids:
                        payload = decoder[token_id]
                        text = payload.decode("utf-8")
                        assert "." not in text
                        kind = "base" if token_id < base.initial_size else ("whole" if token_id < base.initial_size + len(table) else "substring")
                        entry = dict(text=text, hex=payload.hex(), kind=kind)
                        assert str(token_id) not in dictionary[scheme] or dictionary[scheme][str(token_id)] == entry
                        dictionary[scheme][str(token_id)] = entry
                        tokens.append(dict(id=token_id, text=text, start=cursor, end=cursor + len(text), kind=kind))
                        cursor += len(text)
                    assert cursor == offset + len(fragment)
                    parts.append(dict(text=fragment, mode=mode, lookup_id=direct_id, ids=ids, tokens=tokens,
                                      base_ids=base_ids, base_tokens=[base.payloads[i].decode() for i in base_ids],
                                      byte_fallback=byte_fallback))
                    if fragment_index:
                        full_ids.append(base.dot)
                    full_ids.extend(ids)
                    offset += len(fragment) + 1
                dictionary[scheme][str(base.dot)] = dict(text=".", hex="2e", kind="separator")
                expected_ids = base.encode(safe, selected)[0] if scheme == "whole" else hybrid.encode(safe)
                assert full_ids == expected_ids
                restored = decode(full_ids, decoder)
                assert restored == safe
                # Exercise the real SAFE decoder on the ID-restored string only.
                decoded = safe_library.decode(restored, canonical=True, fix=False)
                assert identity(decoded) == identity(example["smiles"])
                result[scheme] = dict(fragments=parts, ids=full_ids, length=len(full_ids), dot_id=base.dot,
                                      content_length=len(full_ids) - len(parts) + 1,
                                      exact_safe_roundtrip=True, isomeric_roundtrip=True, decoded_smiles=decoded)
            example.setdefault("codecs", {})[str(budget)] = result
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--graphs", type=Path, required=True)
    p.add_argument("--whole-dir", type=Path, required=True)
    p.add_argument("--hybrid-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    data = build(args.graphs, args.whole_dir, args.hybrid_dir)
    args.output.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for e in data["examples"]:
        print(e["chembl_id"], {k: [v[s]["length"] for s in ("whole", "hybrid")] for k, v in e["codecs"].items()})


if __name__ == "__main__":
    main()
