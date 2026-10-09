"""Partition original atom indices; never reconstruct atoms from token labels."""
from rdkit import Chem
from rdkit.Chem import BRICS
from .types import TokenizerError
from ._vendor.graphbpe import MolecularSubgraphProcessor
from ._vendor.graphbpe_utils import get_sub_molecule, molecule_to_smiles


def validate_groups(mol, groups):
    indices = [i for group in groups for i in group]
    if (any(not group for group in groups) or len(indices) != mol.GetNumAtoms()
            or set(indices) != set(range(mol.GetNumAtoms()))):
        raise TokenizerError('Fragmentation must partition every original atom exactly once')
    return groups


def brics_groups(mol):
    cuts = {tuple(sorted(pair)) for pair, _ in BRICS.FindBRICSBonds(mol)}
    remaining = set(range(mol.GetNumAtoms()))
    groups = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        group, pending = [start], [start]
        while pending:
            atom = pending.pop()
            for neighbor in mol.GetAtomWithIdx(atom).GetNeighbors():
                other = neighbor.GetIdx()
                if other in remaining and tuple(sorted((atom, other))) not in cuts:
                    remaining.remove(other)
                    group.append(other)
                    pending.append(other)
        groups.append(tuple(sorted(group)))
    return validate_groups(mol, groups)


class NPEPartitioner:
    def __init__(self, engine):
        self.engine = engine

    def groups(self, mol):
        # Copy without a SMILES round trip: SAFE's atom indices must not change.
        work = Chem.Mol(mol)
        Chem.Kekulize(work, clearAromaticFlags=True)
        processor = MolecularSubgraphProcessor(
            Chem.MolToSmiles(work), work, kekulize=True,
            allowable_rings=self.engine.initial_rings)
        for _ in range(work.GetNumAtoms() + 1):
            choices = [s for s in processor.get_nei_smis() if s in self.engine.vocab_node_stats]
            if not choices:
                break
            selected = max(choices, key=lambda s: self.engine.vocab_node_stats[s][1])
            before = len(processor.subgraphs)
            processor.merge(selected)
            if len(processor.subgraphs) >= before:
                raise TokenizerError('NPE merge did not reduce the number of fragments')
        else:
            raise TokenizerError('NPE merge iteration limit exceeded')
        groups = [tuple(sorted(atoms)) for atoms in processor.subgraphs.values()]
        return validate_groups(mol, groups)

    def __call__(self, mol):
        groups = self.groups(mol)
        owner = {atom: group_id for group_id, group in enumerate(groups) for atom in group}
        # Keep rare ring systems intact even when NPE falls back to atom nodes.
        # Cutting aromatic ring bonds can leave non-ring aromatic fragments that
        # SAFE/RDKit cannot sanitize. This is an explicit SAFE cut constraint.
        return [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()
                if not bond.IsInRing()
                and owner[bond.GetBeginAtomIdx()] != owner[bond.GetEndAtomIdx()]]


def brics_processor(smiles, mol):
    """Use the DemoDiff motif/attachment codec with a fixed BRICS partition."""
    processor = MolecularSubgraphProcessor(smiles, mol, kekulize=True, simple_mode=True)
    groups = brics_groups(mol)
    processor.subgraphs = {i: {a: mol.GetAtomWithIdx(a).GetSymbol() for a in group}
                           for i, group in enumerate(groups)}
    processor.subgraphs_smiles = {
        i: molecule_to_smiles(get_sub_molecule(mol, list(group), True), kekulize=True)
        for i, group in enumerate(groups)}
    processor.inversed_index = {a: i for i, group in enumerate(groups) for a in group}
    return processor
