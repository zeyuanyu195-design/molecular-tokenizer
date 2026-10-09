"""Small feasibility probe: NPE atom partitions -> SAFE custom cut function.

No NPE decoding and no vocabulary retraining. Original model files are unchanged.
"""
import json
from pathlib import Path
from rdkit import Chem, RDLogger
from safe import SAFEConverter
from molecular_tokenizer import MolecularTokenizer
from molecular_tokenizer._vendor.graphbpe import MolecularSubgraphProcessor

ROOT=Path(__file__).resolve().parents[1]
RDLogger.DisableLog('rdApp.*')
npe=MolecularTokenizer.load(ROOT/'trained_models/chembl37_npe_vocab350ring300')
safe=MolecularTokenizer.load(ROOT/'trained_models/chembl37_safe_vocab3000')
engine=npe.backend.engine

def npe_cuts(mol):
    # Preserve the index correspondence to the molecule supplied by SAFE.
    work=Chem.Mol(mol)
    Chem.Kekulize(work,clearAromaticFlags=True)
    p=MolecularSubgraphProcessor(Chem.MolToSmiles(work),work,kekulize=True,
                                simple_mode=False,allowable_rings=engine.initial_rings)
    for _ in range(work.GetNumAtoms()+1):
        candidates=[s for s in p.get_nei_smis() if s in engine.vocab_node_stats]
        if not candidates: break
        chosen=max(candidates,key=lambda s:engine.vocab_node_stats[s][1])
        before=len(p.subgraphs)
        p.merge(chosen)
        if len(p.subgraphs)>=before: raise RuntimeError('Non-progressing merge')
    else: raise RuntimeError('Merge limit exceeded')
    owner={idx:pid for pid,atoms in p.subgraphs.items() for idx in atoms}
    assert len(owner)==mol.GetNumAtoms()
    assert sum(len(atoms) for atoms in p.subgraphs.values())==mol.GetNumAtoms()
    return [(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()
            if owner[b.GetBeginAtomIdx()]!=owner[b.GetEndAtomIdx()]]

converter=SAFEConverter(slicer=npe_cuts,ignore_stereo=False)
cases=['CCO','c1ccccc1','CC(=O)Oc1ccccc1C(=O)O','N[C@@H](C)C(=O)O',
       '[13CH3]CO','C/C=C\\C','O=[N+]([O-])c1ccccc1','CC(=O)[O-].[Na+]']
rows=[]
for smiles in cases:
    row={'input':smiles}
    try:
        encoded=converter.encoder(smiles,canonical=True,randomize=False,allow_empty=True)
        decoded=converter.decoder(encoded,canonical=True,fix=False,remove_dummies=False)
        canonical=lambda s:Chem.MolToSmiles(Chem.MolFromSmiles(s),canonical=True,isomericSmiles=True)
        old=safe.backend.representation(smiles)
        row.update(npe_safe=encoded,decoded=decoded,fidelity=canonical(decoded)==canonical(smiles),
                   brics_safe=old,npe_safe_blocks=len(encoded.split('.')),brics_safe_blocks=len(old.split('.')),
                   npe_safe_tokens_with_existing_lexical_bpe=len(safe.backend.tokenizer.encode(encoded).ids),
                   brics_safe_tokens_with_existing_lexical_bpe=len(safe.backend.tokenizer.encode(old).ids))
    except Exception as ex: row['error']=f'{type(ex).__name__}: {ex}'
    rows.append(row)
out={'scope':'8 diagnostic molecules, original graph retained, existing 350-node NPE partition; no retraining or broad evaluation',
     'warning':'SAFE may veto stereo-unsafe cuts, so final blocks can differ from the raw NPE partition.',
     'results':rows}
(ROOT/'analysis_units/npe_safe_feasibility.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=False,indent=2))
