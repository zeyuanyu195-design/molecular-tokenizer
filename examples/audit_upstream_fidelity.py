"""Check diagnostic molecules using the unmodified pinned DemoDiff code."""
import json
import sys
from pathlib import Path
from rdkit import Chem, RDLogger

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parent/'DemoDiff'
sys.path.insert(0,str(REPO/'downstream'))
from graphbpe import MolecularGraphTokenizer
RDLogger.DisableLog('rdApp.*')

def canonical(s):
    return Chem.MolToSmiles(Chem.MolFromSmiles(s),canonical=True,isomericSmiles=True)

cases=['CCO','[13CH3]CO','N[C@@H](C)C(=O)O','N[C@H](C)C(=O)O',
       'C/C=C/C','C/C=C\\C','O=[N+]([O-])c1ccccc1','O=[N+](O)c1ccccc1',
       'C[N+](C)(C)C','CC(=O)[O-].[Na+]']
output={'upstream_commit':'0acb85ef446c175a3e3a60bbad43de97d97ca76f',
        'rdkit_version':__import__('rdkit').__version__,'code_file':str(REPO/'downstream/graphbpe.py'),'models':{}}
for name,prefix in [('local_350',ROOT/'trained_models/chembl37_npe_vocab350ring300/graph'),
                    ('author_pretrained',REPO/'data/tokenizer/vocab3000ring300/pretrain-token')]:
    tokenizer=MolecularGraphTokenizer(kekulize=True,name='upstream-audit')
    tokenizer.load(str(prefix))
    rows=[]
    for s in cases:
        tokenizer.unknown.clear()
        row={'input':s,'canonical_input':canonical(s)}
        try:
            nodes,adj=tokenizer.encode(s,update_vocab_edge=False)
            bonds,positions=tokenizer.get_bond_position_by_vocab(adj)
            decoded,_,unknown=tokenizer.decode(nodes,bonds,positions,replace_unknown_with_random=False)
            row.update(nodes=nodes,adjacency=adj,decoded=decoded,unknown=unknown,
                       fidelity=decoded is not None and canonical(decoded)==canonical(s))
        except Exception as ex:
            row['error']=f'{type(ex).__name__}: {ex}'
        rows.append(row)
    output['models'][name]={'node_vocab_size':len(tokenizer.vocab_node),'results':rows}
(ROOT/'analysis_units/upstream_fidelity_audit.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
for name,data in output['models'].items():
    print(name,data['node_vocab_size'])
    for r in data['results']: print(r['input'],'=>',r.get('decoded',r.get('error')),r.get('fidelity'))
    for i,j in [(2,3),(4,5),(6,7)]:
        a,b=data['results'][i],data['results'][j]
        print('encoding collision',i,j,a.get('nodes')==b.get('nodes') and a.get('adjacency')==b.get('adjacency'))
