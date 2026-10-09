"""Reproducible descriptive audit; not a controlled held-out benchmark."""
import collections
import csv
import json
import random
import re
import statistics
import time
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS
from molecular_tokenizer import MolecularTokenizer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis_units'
MODELS = ROOT / 'trained_models'
RDLogger.DisableLog('rdApp.*')

def stats(values):
    values = sorted(values)
    if not values:
        return {'n': 0}
    return dict(n=len(values), mean=statistics.mean(values), median=statistics.median(values),
                min=values[0], max=values[-1], p90=values[int(.9*(len(values)-1))])

def safe_kind(token):
    if token == '[UNK]': return 'special'
    if re.fullmatch(r'\[[^\]]+\]', token) or token in ('B','C','N','O','S','P','F','I','Br','Cl','b','c','n','o','s','p','*'):
        return 'complete_atom_lexeme'
    if re.fullmatch(r'[0-9]|%[0-9]{2,5}|%\([0-9]{1,5}\)', token): return 'ring_or_attachment_label'
    if token in ('(',')','.','=','#','-','+','\\','/',':','~','@','?','>','>>','$'): return 'syntax'
    return 'partial_lexeme_or_alphabet_fallback'

def mol_props(token):
    m = Chem.MolFromSmiles(token)
    if m is None and re.fullmatch(r'[A-Z][a-z]?', token):
        m = Chem.MolFromSmiles('['+token+']')
    if m is None: return {'parseable': False}
    return {'parseable': True, 'atoms': m.GetNumAtoms(),
            'heavy_atoms': m.GetNumHeavyAtoms(), 'ring_count': m.GetRingInfo().NumRings(),
            'dummy_atoms': sum(a.GetAtomicNum()==0 for a in m.GetAtoms())}

def main():
    OUT.mkdir(exist_ok=True)
    started = time.time()
    models = {k: MolecularTokenizer.load(MODELS / p) for k,p in
              [('safe','chembl37_safe_vocab3000'),('npe','chembl37_npe_vocab350ring300')]}
    sv = models['safe'].backend.tokenizer.get_vocab()
    engine = models['npe'].backend.engine
    nv = engine.vocab_node
    safe_details = [{'token':t, 'id':i, 'category':safe_kind(t)} for t,i in sorted(sv.items(), key=lambda x:x[1])]
    npe_details = []
    for i,t in enumerate(nv):
        kind = 'initial_ring' if t in engine.initial_rings else ('initial_element' if re.fullmatch(r'[A-Z][a-z]?',t) else 'learned_merge')
        npe_details.append(dict(token=t, id=i, category=kind, **mol_props(t),
                                stored_atom_count=engine.vocab_node_stats[t][0]))
    # Uniform sample of source row indices, independent of the NPE training seed.
    selected = set(random.Random(20261008).sample(range(1,2897819+1),2000))
    training = set()
    with (MODELS/'chembl37_npe_vocab350ring300/training_selection.tsv').open(encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'): training.add(int(r['source_row']))
    samples=[]
    with (ROOT.parent/'chembl_37_smiles.csv').open(encoding='utf-8',newline='') as f:
        for i,r in enumerate(csv.DictReader(f),1):
            if i in selected: samples.append(dict(source_row=i, **r))
    (OUT/'sample_selection.json').write_text(json.dumps(samples,ensure_ascii=False,indent=2),encoding='utf-8')
    counts = {k:collections.Counter() for k in ['safe','npe','brics']}
    doccounts = {k:collections.Counter() for k in counts}
    records=[]
    def evaluate(smiles):
        m=Chem.MolFromSmiles(smiles)
        original=Chem.MolToSmiles(m,True)
        row={'smiles':smiles, 'heavy_atoms':m.GetNumHeavyAtoms(),
             'has_stereo':('@' in smiles or '/' in smiles or '\\' in smiles),
             'has_charge':any(a.GetFormalCharge() for a in m.GetAtoms()),
             'has_isotope':any(a.GetIsotope() for a in m.GetAtoms())}
        for kind,model in models.items():
            result={}
            try:
                ids,edges=model.backend.encode(smiles)
                tokens=[model.backend.tokenizer.id_to_token(i) if kind=='safe' else nv[i] for i in ids]
                result.update(encoded=True, tokens=tokens, units=len(ids), connections=len(edges))
                if kind=='safe':
                    representation=model.backend.representation(smiles)
                    result['representation']=representation
                    result['safe_blocks']=len(representation.split('.'))
                try:
                    decoded=model.backend.decode(ids,edges)
                    result['decoded']=decoded
                    result['fidelity']=decoded==original
                    if not result['fidelity']:
                        dm=Chem.MolFromSmiles(decoded)
                        result['equal_without_stereo_and_isotopes']=Chem.MolToSmiles(dm,isomericSmiles=False)==Chem.MolToSmiles(m,isomericSmiles=False)
                except Exception as ex: result.update(fidelity=False,decode_error=f'{type(ex).__name__}: {ex}')
            except Exception as ex: result.update(encoded=False,fidelity=False,encode_error=f'{type(ex).__name__}: {ex}')
            row[kind]=result
        try:
            frags=Chem.GetMolFrags(BRICS.BreakBRICSBonds(m),asMols=True)
            ts=[Chem.MolToSmiles(x,True) for x in frags]
            row['brics']={'tokens':ts,'units':len(ts), 'heavy_atoms':[x.GetNumHeavyAtoms() for x in frags]}
        except Exception as ex: row['brics']={'error':str(ex)}
        return row
    for i,s in enumerate(samples):
        try: row=evaluate(s['canonical_smiles'])
        except Exception as ex:
            row={'smiles':s['canonical_smiles'],'input_error':str(ex)}
        row.update(source_row=s['source_row'],chembl_id=s['chembl_id'],in_npe_training=s['source_row'] in training)
        records.append(row)
        for k in counts:
            ts=row.get(k,{}).get('tokens',[])
            counts[k].update(ts); doccounts[k].update(set(ts))
        if (i+1)%200==0: print(f'{i+1}/{len(samples)} elapsed={time.time()-started:.1f}s',flush=True)
    examples={name:evaluate(s) for name,s in {
        'ethanol':'CCO','aspirin':'CC(=O)Oc1ccccc1C(=O)O','benzene':'c1ccccc1',
        'alanine':'N[C@@H](C)C(=O)O','isotope_ethanol':'[13CH3]CO',
        'tetramethylammonium':'C[N+](C)(C)C','trans_alkene':'C/C=C/C',
        'nitrobenzene':'O=[N+]([O-])c1ccccc1'}.items()}
    summary={'sample_n':len(samples),'sample_seed':20261008,
             'sample_is_heldout':False, 'npe_training_overlap':sum(r['in_npe_training'] for r in records),
             'sample_input_errors':sum('input_error' in r for r in records),
             'safe_vocab_categories':dict(collections.Counter(r['category'] for r in safe_details)),
             'npe_vocab_categories':dict(collections.Counter(r['category'] for r in npe_details)),
             'npe_vocab_atom_sizes':stats([r['heavy_atoms'] for r in npe_details if r['parseable']]),
             'npe_unparseable_tokens':[r['token'] for r in npe_details if not r['parseable']],
             'npe_stored_atom_count_disagreements':[r for r in npe_details if r['parseable'] and r['atoms']!=r['stored_atom_count']],
             'safe_vocab_size':len(sv),'npe_vocab_size':len(nv),'npe_edge_vocab_size':len(engine.vocab_edge)}
    for k in counts:
        total=sum(counts[k].values())
        summary[k]={'observed_unique_units':len(counts[k]), 'total_occurrences':total,
                    'top20':[{'token':t,'occurrences':n,'molecules':doccounts[k][t]} for t,n in counts[k].most_common(20)],
                    'top10_occurrence_fraction':sum(n for _,n in counts[k].most_common(10))/total,
                    'units_per_molecule':stats([r[k]['units'] for r in records if 'units' in r.get(k,{})])}
        if k in models:
            summary[k].update(encode_success=sum(r.get(k,{}).get('encoded',False) for r in records),
                              fidelity_success=sum(r.get(k,{}).get('fidelity',False) for r in records),
                              encode_errors=sum('encode_error' in r.get(k,{}) for r in records),
                              decode_errors=sum('decode_error' in r.get(k,{}) for r in records),
                              mismatch=sum('decoded' in r.get(k,{}) and not r[k]['fidelity'] for r in records))
            summary[k]['by_input_feature']={f:{'n':sum(bool(r.get(f)) for r in records),
                   'fidelity_success':sum(bool(r.get(f)) and r.get(k,{}).get('fidelity',False) for r in records)}
                   for f in ['has_stereo','has_charge','has_isotope']}
    summary['safe']['observed_category_occurrences']=dict(sum((collections.Counter({safe_kind(t):n}) for t,n in counts['safe'].items()),collections.Counter()))
    summary['safe']['blocks_per_molecule']=stats([r['safe']['safe_blocks'] for r in records if 'safe_blocks' in r.get('safe',{})])
    summary['npe']['connections_per_encoded_molecule']=stats([r['npe']['connections'] for r in records if 'connections' in r.get('npe',{})])
    both=[r for r in records if all(r.get(k,{}).get('fidelity') for k in models)]
    summary['joint_fidelity_subset']={'n':len(both),**{k:stats([r[k]['units'] for r in both]) for k in models}}
    summary['brics']['observed_vocab_heavy_atoms']=stats([mol_props(t)['heavy_atoms'] for t in counts['brics']])
    summary['brics']['singleton_types']=sum(n==1 for n in counts['brics'].values())
    summary['brics']['occurrence_weighted_heavy_atoms']=stats([n for r in records for n in r.get('brics',{}).get('heavy_atoms',[])])
    for name,data in [('summary',summary),('vocabulary_units',{'safe':safe_details,'npe':npe_details}),
                      ('molecule_results',records),('examples',examples),
                      ('unit_frequencies',{k:dict(v) for k,v in counts.items()})]:
        (OUT/f'{name}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__': main()
