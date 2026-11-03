"""Build a verified, offline SAFE connection walkthrough from two saved records.

The graph/atom metadata is for visualization and validation, not decoder input.
The decoder receives only the complete SAFE string. No merge history is inferred.
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def identity(mol):
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def partition(groups):
    return sorted(tuple(sorted(g)) for g in groups)


def edges(pairs):
    return {tuple(sorted(p)) for p in pairs}


def sha(path, normalize_text=False):
    content = Path(path).read_bytes()
    if normalize_text:
        content = content.replace(b"\r\n", b"\n")
    return hashlib.sha256(content).hexdigest()


def lex_safe(text, lexeme, atom_pattern):
    """Track SMILES lexical atoms and paired ring occurrences, including label reuse."""
    tokens, pairs, pending, stack = [], [], {}, []
    current, fragment, cursor, atom_count = None, 0, 0, 0
    for match in lexeme.finditer(text):
        assert match.start() == cursor, "Unrecognized syntax"
        value = match.group()
        token = dict(text=value, start=match.start(), end=match.end(), fragment=fragment)
        if atom_pattern.fullmatch(value):
            current = atom_count
            atom_count += 1
            token.update(kind="atom", atom=current)
        elif value == "(":
            stack.append(current)
            token["kind"] = "syntax"
        elif value == ")":
            current = stack.pop()
            token["kind"] = "syntax"
        elif value == ".":
            current = None
            token["kind"] = "separator"
            fragment += 1
        elif re.fullmatch(r"\d|%\d{2}|%\(\d+\)", value):
            assert current is not None
            token.update(kind="ring", atom=current)
            if value in pending:
                first = pending.pop(value)
                pair_id = len(pairs)
                tokens[first]["pair"] = pair_id
                token["pair"] = pair_id
                pairs.append(dict(id=pair_id, label=value, tokens=[first, len(tokens)]))
            else:
                pending[value] = len(tokens)
        else:
            token["kind"] = "syntax"
        tokens.append(token)
        cursor = match.end()
    assert cursor == len(text) and not pending and not stack
    return tokens, pairs, atom_count


def build(row, pipeline, mapper, lexeme, atom_pattern):
    # Replay the current pipeline; abort if historical boundaries or SAFE changed.
    replay = pipeline.convert(row["canonical_smiles"])
    for field in ("safe", "bpe_groups", "bpe_node_ids", "requested_cuts", "used_cuts"):
        assert (
            json.loads(json.dumps(replay[field])) == row[field]
        ), f"Replay differs: {row['chembl_id']} {field}"
    source = Chem.MolFromSmiles(row["indexed_source_smiles"])
    source = Chem.RenumberAtoms(
        source,
        sorted(range(source.GetNumAtoms()), key=lambda i: source.GetAtomWithIdx(i).GetAtomMapNum()),
    )
    for i, atom in enumerate(source.GetAtoms()):
        assert atom.GetAtomMapNum() == i + 1
        atom.SetAtomMapNum(0)
    safe_mol = Chem.MolFromSmiles(row["safe"])
    assert identity(source) == identity(safe_mol) == row["normalized_smiles"]
    tokens, pairs, atom_count = lex_safe(row["safe"], lexeme, atom_pattern)
    assert atom_count == safe_mol.GetNumAtoms() == source.GetNumAtoms()
    owners = [t["fragment"] for t in tokens if t["kind"] == "atom"]
    nfrags = len(row["safe"].split("."))
    mapping = None
    for match in source.GetSubstructMatches(
        safe_mol, uniquify=False, useChirality=True, maxMatches=10000
    ):
        groups = [[match[i] for i, f in enumerate(owners) if f == j] for j in range(nfrags)]
        if partition(groups) == partition(row["safe_groups"]):
            mapping = match
            break
    assert mapping is not None, "No stereochemistry/partition-compatible atom mapping"
    for token in tokens:
        if "atom" in token:
            token["atom"] = mapping[token["atom"]]
    # Both lexical occurrences must point to a real source bond.
    for pair in pairs:
        first, second = [tokens[i] for i in pair["tokens"]]
        u, v = first["atom"], second["atom"]
        bond = source.GetBondBetweenAtoms(u, v)
        assert bond is not None
        pair.update(
            atoms=[u, v],
            fragments=[first["fragment"], second["fragment"]],
            cross_fragment=first["fragment"] != second["fragment"],
            bond=bond.GetIdx(),
            bond_type=str(bond.GetBondType()),
        )
        if pair["cross_fragment"]:
            pair["dummy_isotope"] = next(
                i + 1 for i, uv in enumerate(row["used_cuts"]) if set(uv) == {u, v}
            )
    cross = [p for p in pairs if p["cross_fragment"]]
    assert edges(p["atoms"] for p in cross) == edges(row["used_cuts"])
    lookup = mapper.encode(row["safe"], strict=False)
    fragment_rows = []
    for j, entry in enumerate(lookup["fragments"]):
        group = sorted(mapping[i] for i, f in enumerate(owners) if f == j)
        motifs = [k for k, g in enumerate(row["bpe_groups"]) if set(g) <= set(group)]
        assert set().union(*(set(row["bpe_groups"][k]) for k in motifs)) == set(group)
        fragment_rows.append(dict(**entry, atoms=group, bpe_groups=motifs))
        if row["bpe_boundary_exact"]:
            assert len(motifs) == 1
            assert entry["node_id"] == row["bpe_node_ids"][motifs[0]]
    # No attachment table, atom mapping, or graph metadata is passed to decoder.
    decoded = pipeline.converter.decoder(row["safe"], canonical=True, remove_dummies=False)
    assert identity(Chem.MolFromSmiles(decoded)) == identity(source)
    rdDepictor.Compute2DCoords(source)
    drawer = rdMolDraw2D.MolDraw2DSVG(700, 430)
    options = drawer.drawOptions()
    options.addAtomIndices = True
    options.annotationFontScale = 0.8
    options.padding = 0.12
    options.useBWAtomPalette()
    drawer.DrawMolecule(source)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg[svg.index("<svg") :]
    atoms = []
    for atom in source.GetAtoms():
        i = atom.GetIdx()
        xy = drawer.GetDrawCoords(i)
        atoms.append(
            dict(
                id=i,
                element=atom.GetSymbol(),
                x=round(xy.x, 3),
                y=round(xy.y, 3),
                bpe=next(j for j, g in enumerate(row["bpe_groups"]) if i in g),
                fragment=next(j for j, f in enumerate(fragment_rows) if i in f["atoms"]),
            )
        )
    requested, used = edges(row["requested_cuts"]), edges(row["used_cuts"])
    bonds = []
    for b in source.GetBonds():
        uv = tuple(sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]))
        bonds.append(
            dict(
                id=b.GetIdx(),
                atoms=list(uv),
                type=str(b.GetBondType()),
                stereo=str(b.GetStereo()),
                requested=uv in requested,
                used=uv in used,
                protected=uv in requested - used,
            )
        )
    return dict(
        chembl_id=row["chembl_id"],
        record_index=row["record_index"],
        smiles=row["canonical_smiles"],
        safe=row["safe"],
        decoded=decoded,
        exact=row["bpe_boundary_exact"],
        bpe_groups=row["bpe_groups"],
        bpe_node_ids=row["bpe_node_ids"],
        fragments=fragment_rows,
        tokens=tokens,
        pairs=pairs,
        atoms=atoms,
        bonds=bonds,
        svg=svg,
        roundtrip_ok=True,
        lookup_fully_covered=lookup["fully_covered"],
        validation=dict(
            replayed_encoder=True,
            paired_labels_match_used_cuts=True,
            isomeric_roundtrip=True,
            atom_mapping="stereo- and partition-constrained graph isomorphism",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--vocab-prefix", type=Path, required=True)
    parser.add_argument(
        "--records", type=Path, default=Path(__file__).with_name("connection_demo_records.json")
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo_root.resolve()))
    sys.path.insert(0, str(args.repo_root.resolve() / "experiments/tokenization"))
    from safe.batch_chembl import Pipeline
    from demodiff_fixed_safe import ATOM, LEXEME, FixedSAFEVocabMapper

    pipeline = Pipeline(args.vocab_prefix, allow_coarsening=True)
    mapper = FixedSAFEVocabMapper(pipeline.slicer.tokenizer)
    records = json.loads(args.records.read_text(encoding="utf-8"))
    examples = [build(r, pipeline, mapper, LEXEME, ATOM) for r in records]
    report = dict(
        schema_version=1,
        rdkit_version=rdBase.rdkitVersion,
        source_commit="61e3a02f371e31f67ff1e04197b7462df6ca0972",
        fingerprint_policy="Source and record text normalized to LF; vocabulary files hashed as bytes",
        records_sha256=sha(args.records, normalize_text=True),
        source_sha256={
            p: sha(args.repo_root / p, normalize_text=True)
            for p in [
                "safe/converter.py",
                "safe/graphbpe_slicer.py",
                "safe/graphbpe.py",
                "experiments/tokenization/demodiff_fixed_safe.py",
            ]
        },
        vocab_sha256={s: sha(str(args.vocab_prefix) + s) for s in (".node", ".edge", ".ring")},
        examples=examples,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for e in examples:
        print(
            json.dumps(
                dict(
                    id=e["chembl_id"],
                    graphbpe=len(e["bpe_groups"]),
                    safe=len(e["fragments"]),
                    pairs=sum(p["cross_fragment"] for p in e["pairs"]),
                    protected=[b for b in e["bonds"] if b["protected"]],
                    ids=[f["node_id"] for f in e["fragments"]],
                    validation=e["validation"],
                )
            )
        )


if __name__ == "__main__":
    main()
