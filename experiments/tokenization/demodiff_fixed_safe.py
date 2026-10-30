"""Fixed SAFE fragment -> existing DemoDiff node-vocabulary lookup.

No GraphBPE merge, vocabulary update, or edge-table construction is performed.
This is a motif lookup adapter, NOT a lossless node-ID-only molecular codec.
The original SAFE string remains the authoritative connectivity representation.
"""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
import re
import sys

from rdkit import Chem

# Standard SMILES atoms; an entire bracket expression counts as one atom.
# Other recognized lexemes never create atoms. Extended ring labels are retained.
LEXEME = re.compile(r"\[[^\[\]\s]+\]|Cl|Br|[BCNOPSFIbcnops*]|%\(\d+\)|%\d{2}|[0-9().=#$:/\\~+@-]")
ATOM = re.compile(r"\[[^\[\]\s]+\]|Cl|Br|[BCNOPSFIbcnops*]")


class UnknownFragmentsError(ValueError):
    def __init__(self, result):
        self.result = result
        super().__init__(
            "Existing DemoDiff vocabulary does not cover fragments: "
            + repr(result["missing_fragments"])
        )


def _parse(text):
    params = Chem.SmilesParserParams()
    params.removeHs = False  # Preserve explicit atoms and their fragment membership.
    params.allowCXSMILES = False
    params.parseName = False
    mol = Chem.MolFromSmiles(text, params)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("Input must be a valid, complete, nonempty SAFE/SMILES string")
    if any(a.GetAtomMapNum() for a in mol.GetAtoms()):
        raise ValueError("Atom-mapped inputs are unsupported; remove maps explicitly first")
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("Open fragments/dummy atoms are unsupported")
    return mol


def _fragment_owners(text):
    owners = []
    cursor = 0
    fragment = 0
    for match in LEXEME.finditer(text):
        if match.start() != cursor:
            raise ValueError(f"Unsupported SAFE syntax at character {cursor}")
        token = match.group()
        if token == ".":
            fragment += 1
        elif ATOM.fullmatch(token):
            owners.append(fragment)
        cursor = match.end()
    if cursor != len(text):
        raise ValueError(f"Unsupported SAFE syntax at character {cursor}")
    return owners


def _reparse_with_owners(mol, owners, kekule=False):
    text = Chem.MolToSmiles(mol, canonical=True, kekuleSmiles=kekule)
    order = list(mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
    result = _parse(text)
    if len(order) != len(owners) or result.GetNumAtoms() != len(owners):
        raise ValueError("Atom count changed during canonicalization")
    return result, [owners[i] for i in order]


class FixedSAFEVocabMapper:
    """Reuse a loaded original MolecularGraphTokenizer without calling encode()."""

    def __init__(self, tokenizer):
        self.tokenizer = tokenizer
        module = importlib.import_module(type(tokenizer).__module__)
        self.get_sub_molecule = module.get_sub_molecule
        self.molecule_to_smiles = module.molecule_to_smiles
        self.node_to_id = {key: i for i, key in enumerate(tokenizer.vocab_node)}
        if len(self.node_to_id) != len(tokenizer.vocab_node):
            raise ValueError("Duplicate node-vocabulary keys")

    def encode(self, safe_text, strict=True):
        if not isinstance(safe_text, str) or not safe_text or any(c.isspace() for c in safe_text):
            raise ValueError("Pass one SAFE string without whitespace or a molecule name")
        fragments = safe_text.split(".")
        if any(not fragment for fragment in fragments):
            raise ValueError("Empty SAFE fragment")

        # Parse the COMPLETE string: individual fragments may have open labels.
        # RDKit atom order follows the atom-lexeme order in the input SMILES.
        mol = _parse(safe_text)
        owners = _fragment_owners(safe_text)
        if len(owners) != mol.GetNumAtoms():
            raise ValueError("SAFE atom lexemes do not match parsed atom count")
        kekulize = self.tokenizer.kekulize
        if kekulize:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        # Match DemoDiff's canonicalization while tracking original dot boundaries.
        mol, owners = _reparse_with_owners(mol, owners)
        if kekulize and mol.GetNumAtoms() > 1:
            mol, owners = _reparse_with_owners(mol, owners, kekule=True)
        if kekulize:
            Chem.Kekulize(mol, clearAromaticFlags=True)

        groups = [[] for _ in fragments]
        for atom_index, fragment_index in enumerate(owners):
            groups[fragment_index].append(atom_index)

        rows = []
        for index, (fragment, group) in enumerate(zip(fragments, groups)):
            submol = self.get_sub_molecule(mol, group, kekulize=kekulize)
            if submol.GetNumAtoms() != len(group):
                raise ValueError(f"Subgraph extraction lost atoms in fragment {index}")
            # Match the upstream get_nei_smis/get-rings call signature exactly.
            # Upstream uses its second positional argument as `canonical`.
            key = self.molecule_to_smiles(submol, kekulize)
            rows.append(
                {
                    "fragment_index": index,
                    "safe_fragment": fragment,
                    "vocab_key": key,
                    "node_id": self.node_to_id.get(key),
                }
            )

        missing = [row for row in rows if row["node_id"] is None]
        result = {
            "safe": safe_text,
            "fragments": rows,
            "node_ids": [row["node_id"] for row in rows],
            "num_fragments": len(rows),
            "matched_fragments": len(rows) - len(missing),
            "fully_covered": not missing,
            "num_node_tokens": len(rows) if not missing else None,
            "missing_fragments": missing,
            "ids_are_lossless": False,
        }
        if missing and strict:
            raise UnknownFragmentsError(result)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demodiff-root", type=Path)
    source.add_argument(
        "--graphbpe-dir", type=Path, help="Directory containing graphbpe.py and its utils package"
    )
    parser.add_argument("--vocab-prefix", type=Path)
    parser.add_argument("--safe", required=True)
    args = parser.parse_args()
    if args.demodiff_root:
        root = args.demodiff_root.resolve()
        module_dir = root / "downstream"
        prefix = args.vocab_prefix or root / "data/tokenizer/vocab3000ring300/pretrain-token"
    else:
        module_dir = args.graphbpe_dir.resolve()
        if args.vocab_prefix is None:
            parser.error("--graphbpe-dir requires --vocab-prefix")
        prefix = args.vocab_prefix
    if not (module_dir / "graphbpe.py").is_file():
        parser.error(f"Cannot find graphbpe.py in {module_dir}")
    sys.path.insert(0, str(module_dir))
    from graphbpe import MolecularGraphTokenizer

    tokenizer = MolecularGraphTokenizer(True, "pretrain", False)
    tokenizer.load(str(prefix))
    mapper = FixedSAFEVocabMapper(tokenizer)
    result = mapper.encode(args.safe, strict=False)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if not result["fully_covered"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
