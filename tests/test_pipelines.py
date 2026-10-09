from dataclasses import replace
import json
import pytest
from rdkit import Chem
from molecular_tokenizer import MolecularTokenizer, TokenizerError, FidelityError
from molecular_tokenizer.fragmentation import NPEPartitioner

DATA = ['CCO','CCCO','CCN','CCCN','CC(=O)O','CCOC(=O)C','c1ccccc1',
        'Cc1ccccc1','Oc1ccccc1','c1ccncc1','CCCl','CCF','CCOC(=O)c1ccccc1','CCN(CC)CC']
DIAGNOSTICS = ['N[C@@H](C)C(=O)O','N[C@H](C)C(=O)O','C/C=C/C','C/C=C\\C',
               '[13CH3]CO','CC(=O)[O-].[Na+]','O=[N+]([O-])c1ccccc1',
               'C[N+](C)(C)C','C1CCC2(CC1)CCCC2','c1ccc2ccccc2c1']

def canonical(s):
    return Chem.MolToSmiles(Chem.MolFromSmiles(s),isomericSmiles=True)

@pytest.fixture(scope='module', params=['brics','npe'])
def sequence_model(request):
    t=MolecularTokenizer(fragmentation=request.param,representation='safe',bpe_scope='fragment')
    options={'motif_vocab_size':20,'ring_vocab_size':2,'num_workers':2} if request.param=='npe' else {}
    t.train(iter(DATA*3),256,**options)
    return t

@pytest.mark.parametrize('smiles',DATA+DIAGNOSTICS)
def test_hybrid_preserves_original_graph(sequence_model,smiles):
    result=sequence_model.encode(smiles)
    assert result.connections==() and result.fidelity_checked
    assert canonical(sequence_model.decode(result))==canonical(smiles)
    assert canonical(sequence_model.to_safe(smiles))==canonical(smiles)

def test_pipeline_save_load_and_no_source_smiles(sequence_model,tmp_path):
    before=sequence_model.tokenizer_id
    encoded=sequence_model.encode(DIAGNOSTICS[0])
    assert set(encoded.to_dict())=={'backend','tokenizer_id','token_ids','connections','fidelity_checked'}
    loaded=MolecularTokenizer.load(sequence_model.save(tmp_path/'model'))
    assert loaded.tokenizer_id==before
    assert loaded.fragmentation==sequence_model.fragmentation
    assert loaded.representation=='safe'
    for s in DATA+DIAGNOSTICS:
        assert loaded.encode(s)==sequence_model.encode(s)
    assert sequence_model._fingerprint()==before

def test_npe_partition_is_index_preserving(sequence_model):
    if sequence_model.fragmentation!='npe': return
    p=NPEPartitioner(sequence_model.backend.npe.engine)
    m=Chem.MolFromSmiles('[13CH3][C@@H](O)C(=O)[O-]')
    signature=Chem.MolToMolBlock(m)
    groups=p.groups(m)
    assert sorted(a for group in groups for a in group)==list(range(m.GetNumAtoms()))
    assert Chem.MolToMolBlock(m)==signature

def test_fragment_bpe_actually_merges_across_atoms():
    models=[]
    for scope in ['lexical','fragment']:
        t=MolecularTokenizer(fragmentation='brics',representation='safe',bpe_scope=scope)
        t.train(['CCO']*20,128)
        models.append(t)
    assert len(models[0].encode('CCO').token_ids)==3
    assert len(models[1].encode('CCO').token_ids)==1
    assert models[1].decode(models[1].encode('CCO'))=='CCO'

@pytest.mark.parametrize('fragmentation',['brics','npe'])
def test_demodiff_choice(fragmentation,tmp_path):
    t=MolecularTokenizer(fragmentation=fragmentation,representation='demodiff')
    options={'ring_vocab_size':2} if fragmentation=='npe' else {}
    t.train(DATA,32,**options)
    for s in DATA:
        assert canonical(t.decode(t.encode(s)))==canonical(s)
    loaded=MolecularTokenizer.load(t.save(tmp_path/'graph'))
    assert loaded.fragmentation==fragmentation
    assert loaded.representation=='demodiff'
    assert loaded.encode('CCO')==t.encode('CCO')
    with pytest.raises(TokenizerError,match='SAFE'):
        t.to_safe('CCO')

def test_brics_graph_oov_and_budget_explicit():
    t=MolecularTokenizer(fragmentation='brics',representation='demodiff')
    t.train(['CCO'],16)
    with pytest.raises(TokenizerError,match='unknown'):
        t.encode('CCN')
    before=t.tokenizer_id
    with pytest.raises(TokenizerError,match='increase vocab_size'):
        t.train(['CCO','CCN'],1)
    assert t.tokenizer_id==before

def test_hybrid_reuses_npe_model_and_is_self_contained(tmp_path):
    graph=MolecularTokenizer('npe')
    graph.train(DATA,20,ring_vocab_size=2)
    source=graph.save(tmp_path/'npe_source')
    t=MolecularTokenizer(fragmentation='npe',representation='safe')
    report=t.train(DATA,256,npe_model=source)
    assert report.motif_vocab_size==graph.vocab_size
    assert json.dumps(t.backend.npe.state(),sort_keys=True)==json.dumps(graph.backend.state(),sort_keys=True)
    saved=t.save(tmp_path/'hybrid')
    source.rename(tmp_path/'source_moved')
    loaded=MolecularTokenizer.load(saved)
    assert loaded.encode('[13CH3]CO')==t.encode('[13CH3]CO')
    (saved/'graph.ring').write_text('corruption')
    with pytest.raises(TokenizerError,match='checksum'):
        MolecularTokenizer.load(saved)

def test_rejection_indices_and_transactionality(sequence_model):
    t=MolecularTokenizer(fragmentation='npe',representation='safe')
    rejected=[]
    report=t.train(iter(['CCO','invalid','CCN','*CC*']),128,
                   motif_vocab_size=6,ring_vocab_size=0,skip_invalid=True,
                   on_reject_with_index=lambda s,r,i: rejected.append(i))
    assert sorted(rejected)==[2,4]
    assert report.molecules==2 and report.rejected_molecules==2
    before=sequence_model.tokenizer_id
    with pytest.raises(TokenizerError):
        sequence_model.train(['invalid'],128)
    assert sequence_model.tokenizer_id==before

def test_vocab_binding_includes_partition_model(sequence_model):
    encoded=sequence_model.encode('CCO')
    with pytest.raises(TokenizerError,match='different backend or vocabulary'):
        sequence_model.decode(replace(encoded,tokenizer_id='another-model'))

@pytest.mark.parametrize('options',[
    {'fragmentation':'bad','representation':'safe'},
    {'fragmentation':'npe','representation':'bad'},
    {'fragmentation':'npe','representation':'demodiff','bpe_scope':'fragment'},
    {'fragmentation':'npe','representation':'safe','bpe_scope':'bad'},
])
def test_invalid_options(options):
    with pytest.raises(TokenizerError): MolecularTokenizer(**options)
