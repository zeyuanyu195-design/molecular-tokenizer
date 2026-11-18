"""Draw the published SAFE examples and map actual codec outputs to graph atoms."""
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from demodiff_fixed_safe import ATOM, LEXEME

DOCS = Path(__file__).resolve().parents[2] / "docs"
DEST = DOCS / "example-graphs"


def lexical_atoms(text):
    words = LEXEME.findall(text)
    assert "".join(words) == text
    return [word for word in words if ATOM.fullmatch(word)]


def emitted_groups(model, piece):
    """Each ordinary ID or complete fallback frame covers one emitted group."""
    groups, frame = [], None
    for token in piece["ids"]:
        entry = model["dictionary"][str(token)]
        if entry["text"] == "BLOCK_BEGIN":
            assert frame is None
            frame = [token]
        elif entry["text"] == "BLOCK_END":
            assert frame is not None
            groups.append(frame + [token])
            frame = None
        elif frame is None:
            groups.append([token])
        else:
            frame.append(token)
    assert frame is None
    return groups


def token_locations(safe, mol):
    """Map SAFE character spans to RDKit atom indices or explicit bond endpoints."""
    rows, branches, rings = [], [], {}
    last, next_atom, pending_bond = None, 0, None
    def assign_bond(indices, a, b):
        assert mol.GetBondBetweenAtoms(a, b) is not None
        for index in indices:
            if index is not None:
                rows[index][3] = [a, b]

    for match in LEXEME.finditer(safe):
        word = match.group()
        index = len(rows)
        kind, atoms = "syntax", []
        rows.append([match.start(), match.end(), kind, atoms])
        if ATOM.fullmatch(word):
            kind, atoms = "atom", [next_atom]
            if pending_bond is not None:
                assign_bond([pending_bond], last, next_atom)
            last, next_atom, pending_bond = next_atom, next_atom + 1, None
        elif word == ".":
            kind, last = "dot", None
        elif word == "(":
            branches.append(last)
            kind, atoms = "branch", [last]
        elif word == ")":
            last = branches.pop()
            kind, atoms = "branch", [last]
        elif word[0].isdigit() or word.startswith("%"):
            rows[index][2] = "ring"
            if word in rings:
                other, opening, bond = rings.pop(word)
                assign_bond([index, opening, bond, pending_bond], other, last)
            else:
                rings[word] = (last, index, pending_bond)
            pending_bond = None
            continue
        else:
            kind, pending_bond = "bond", index
        rows[index][2:] = [kind, atoms]
    assert next_atom == mol.GetNumAtoms() and not rings and not branches and pending_bond is None
    assert all(atoms or kind == "dot" for _, _, kind, atoms in rows)
    return rows


def build(example):
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(example["safe"], params)
    original = Chem.MolFromSmiles(example["smiles"], params)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(original)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(720, 440)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.addAtomIndices = True
    options.padding = 0.075
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    points = [drawer.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
    bonds = [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()]
    fragments, offset = [], 0
    for text in example["fragments"]:
        stop = offset + len(lexical_atoms(text))
        fragments.append(list(range(offset, stop)))
        offset = stop
    assert offset == mol.GetNumAtoms()
    protected = [("ring", list(ring)) for ring in mol.GetRingInfo().AtomRings()]
    for name, pattern in [("amide", "[CX3](=[OX1])[NX3]"), ("acid/ester", "[CX3](=[OX1])[OX2]")]:
        protected.extend((name, list(group)) for group in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
    layouts = {}
    for budget, models in example["budgets"].items():
        layouts[budget] = {}
        for mode in ("A", "B"):
            model = models[mode]
            records = []
            for fi, piece in enumerate(model["fragments"]):
                groups = emitted_groups(model, piece)
                spans = ["".join(model["dictionary"][str(i)]["payload"] for i in group) for group in groups]
                assert "".join(spans) == "".join(piece["units"]) == example["fragments"][fi]
                views = {}
                for view, texts in [("encoded", spans), ("allowed", piece["units"])]:
                    start = fragments[fi][0]
                    rows = []
                    for j, text in enumerate(texts):
                        stop = start + len(lexical_atoms(text))
                        assert stop > start
                        rows.append(dict(text=text, atoms=list(range(start, stop)), ids=groups[j] if view=="encoded" else []))
                        start = stop
                    assert start == fragments[fi][-1]+1
                    views[view] = rows
                    if mode == "B":
                        for _, group in protected:
                            if set(group) <= set(fragments[fi]):
                                assert any(set(group) <= set(row["atoms"]) for row in rows)
                records.append(views)
            layouts[budget][mode] = records
    return dict(svg=drawer.GetDrawingText(), points=[[round(p.x,2),round(p.y,2)] for p in points],
                bonds=bonds, fragments=fragments, protected=protected, layouts=layouts,
                locations=token_locations(example["safe"], mol))


if __name__ == "__main__":
    DEST.mkdir(exist_ok=True)
    for example in json.loads((DOCS/"long-examples.json").read_text(encoding="utf-8"))["examples"]:
        data = build(example)
        (DEST/(example["chembl_id"]+".json")).write_text(json.dumps(data,separators=(",",":")),encoding="utf-8")
        print(example["chembl_id"], len(data["points"]), "atoms")
