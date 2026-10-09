"""Render figures and validate the descriptive audit's internal consistency."""
import collections
import hashlib
import importlib.metadata
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from molecular_tokenizer import MolecularTokenizer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis_units'
s=json.loads((OUT/'summary.json').read_text(encoding='utf-8'))
r=json.loads((OUT/'molecule_results.json').read_text(encoding='utf-8'))
v=json.loads((OUT/'vocabulary_units.json').read_text(encoding='utf-8'))
e=json.loads((OUT/'examples.json').read_text(encoding='utf-8'))
f=json.loads((OUT/'unit_frequencies.json').read_text(encoding='utf-8'))
assert len(r)==2000 and len({x['source_row'] for x in r})==2000
assert sum(s['safe_vocab_categories'].values())==1176
assert sum(s['npe_vocab_categories'].values())==350
for k in ['safe','npe','brics']:
    assert sum(f[k].values())==sum(len(x[k]['tokens']) for x in r)
for k in ['safe','npe']:
    assert s[k]['encode_success']==s[k]['fidelity_success']+s[k]['mismatch']+s[k]['decode_errors']
    model=MolecularTokenizer.load(ROOT/'trained_models'/('chembl37_safe_vocab3000' if k=='safe' else 'chembl37_npe_vocab350ring300'))
    for example in e.values():
        accepted=False
        try:
            enc=model.encode(example['smiles'])
            model.decode(enc)
            accepted=True
        except Exception:
            pass
        assert accepted==example[k]['fidelity']

extra={
    'safe_non_atom_occurrence_fraction':1-s['safe']['observed_category_occurrences']['complete_atom_lexeme']/s['safe']['total_occurrences'],
    'npe_category_occurrences':dict(sum((collections.Counter({x['category']:f['npe'].get(x['token'],0)}) for x in v['npe']),collections.Counter())),
    'npe_non_stereo_charge_isotope_mismatches':sum(not any(x[a] for a in ['has_stereo','has_charge','has_isotope']) and not x['npe']['fidelity'] for x in r),
    'npe_mismatch_equal_without_stereo_isotopes':sum(x['npe'].get('equal_without_stereo_and_isotopes',False) for x in r),
    'safe_brics_fragment_count_differences':sum(x['safe']['safe_blocks']!=x['brics']['units'] for x in r),
    'npe_stored_atom_count_disagreements_n':len(s['npe_stored_atom_count_disagreements']),
    'npe_learned_merges':[x['token'] for x in v['npe'] if x['category']=='learned_merge'],
    'dependencies':{p:importlib.metadata.version(p) for p in ['rdkit','safe-mol','tokenizers','matplotlib']},
    'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'examples/analyze_units.py',ROOT/'molecular_tokenizer/backends.py',ROOT/'molecular_tokenizer/_vendor/graphbpe.py',OUT/'sample_selection.json']},
    'checks':'Count conservation, unique source rows, category partitions, and strict facade acceptance on eight diagnostic molecules passed.',
}
(OUT/'audit_details.json').write_text(json.dumps(extra,ensure_ascii=False,indent=2),encoding='utf-8')
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
ax=axes[0,0]
ax.barh(['Atom lexemes','Ring / attachment labels','Syntax','Partial / fallback','Special'],[191,745,15,224,1],color='#377eb8')
ax.invert_yaxis(); ax.set_title('SAFE: 1,176 vocabulary entries'); ax.set_xlabel('Number of types')
ax=axes[0,1]
ax.barh(['Initial rings','Learned merges','Initial elements'],[300,24,26],color='#e18d38')
ax.invert_yaxis(); ax.set_title('NPE: 350 node entries (+31 edge types)'); ax.set_xlabel('Number of types')
ax=axes[1,0]
ax.hist([[x['heavy_atoms'] for x in v['npe']], [__import__('rdkit').Chem.MolFromSmiles(t).GetNumHeavyAtoms() for t in f['brics']]],bins=[0,1,2,3,5,7,10,15,25,50,150],label=['NPE vocabulary','Observed BRICS dictionary'],density=True,color=['#e18d38','#4daf4a'])
ax.set_xlabel('Heavy atoms per distinct unit (different dictionaries)'); ax.set_ylabel('Density'); ax.set_title('Chemical unit size; excludes BRICS dummy atoms'); ax.legend()
ax=axes[1,1]
ax.bar(['SAFE','NPE'],[100,57.7],color=['#377eb8','#e18d38'])
ax.set_ylim(0,110); ax.set_ylabel('Exact round-trip (%)'); ax.set_title('Same 2,000 molecules; descriptive audit')
for i,val in enumerate([100,57.7]): ax.text(i,val+2,f'{val:g}%',ha='center')
fig.suptitle('Existing model audit: unequal training sets / budgets; not a held-out benchmark',fontsize=12)
fig.savefig(OUT/'unit_comparison.png',dpi=180)
fig.savefig(OUT/'unit_comparison.pdf')
plt.close(fig)
print(json.dumps(extra,ensure_ascii=False,indent=2))
