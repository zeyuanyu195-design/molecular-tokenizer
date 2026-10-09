from dataclasses import replace
import json
from pathlib import Path
import shutil

import pytest
from rdkit import Chem
from molecular_tokenizer import (MolecularTokenizer, MoleculeEncoding, Connection,
                                 FidelityError, TokenizerError, NotTrainedError,
                                 UnsupportedMoleculeError)

DATA = ['CCO', 'CCCO', 'CCN', 'CCCN', 'CC(=O)O', 'CCOC(=O)C', 'c1ccccc1',
        'Cc1ccccc1', 'Oc1ccccc1', 'c1ccncc1', 'CCCl', 'CCF',
        'CCOC(=O)c1ccccc1', 'CCN(CC)CC']
SPECIAL = ['N[C@@H](C)C(=O)O', 'C/C=C/C', 'CC(=O)[O-].[Na+]', '[13CH3]CO']
PRETRAINED = Path(__file__).resolve().parent / 'fixtures/demodiff/pretrain-token'

def canonical(value):
    return Chem.MolToSmiles(Chem.MolFromSmiles(value), canonical=True, isomericSmiles=True)

@pytest.fixture(scope='module')
def safe():
    t = MolecularTokenizer('safe')
    t.train(iter(DATA), 256)
    return t

@pytest.fixture(scope='module')
def npe():
    return MolecularTokenizer.from_graphbpe(PRETRAINED)

@pytest.mark.parametrize('backend', ['safe', 'npe'])
@pytest.mark.parametrize('smiles', ['CCCl', 'CCOC(=O)c1ccccc1', 'c1ccncc1', 'CCN(CC)CC'])
def test_molecular_roundtrip(backend, smiles, safe, npe):
    model = safe if backend == 'safe' else npe
    encoded = model.encode(smiles)
    assert encoded.fidelity_checked
    assert canonical(model.decode(encoded)) == canonical(smiles)
    restored = MoleculeEncoding.from_dict(json.loads(json.dumps(encoded.to_dict())))
    assert restored == encoded
    assert model.decode(restored) == model.decode(encoded)

@pytest.mark.parametrize('smiles', SPECIAL)
def test_safe_preserves_unseen_stereo_charge_isotope_and_dummy(safe, smiles):
    assert canonical(safe.decode(safe.encode(smiles))) == canonical(smiles)

def test_safe_explicitly_rejects_open_dummy_attachments(safe):
    with pytest.raises(UnsupportedMoleculeError, match='dummy'):
        safe.encode('*CC*')

@pytest.mark.parametrize('smiles', ['CC(=O)[O-].[Na+]', '[13CH3]CO'])
def test_npe_rejects_upstream_information_loss(npe, smiles):
    with pytest.raises(FidelityError):
        npe.encode(smiles)

def test_frozen_graph_vocab_and_explicit_connections(npe):
    before = json.dumps(npe.backend.state(), sort_keys=True)
    encoded = npe.encode('CCCl')
    assert len(encoded.connections) == 1
    assert len(encoded.token_ids) == 2
    assert before == json.dumps(npe.backend.state(), sort_keys=True)
    assert npe.decode(encoded) == 'CCCl'

def test_unseen_edges_rejected_even_in_relaxed_mode(tmp_path):
    prefix = tmp_path / 'incomplete'
    for suffix in ['node', 'ring']:
        shutil.copyfile(str(PRETRAINED) + '.' + suffix, str(prefix) + '.' + suffix)
    Path(str(prefix) + '.edge').write_text('', encoding='utf-8')
    t = MolecularTokenizer.from_graphbpe(prefix, strict=False)
    with pytest.raises(UnsupportedMoleculeError, match='frozen edge'):
        t.encode('CCCl')

@pytest.mark.parametrize('backend', ['safe', 'npe'])
def test_save_load(backend, safe, npe, tmp_path):
    t = safe if backend == 'safe' else npe
    path = t.save(tmp_path / backend)
    saved = t.encode('CCCl')
    loaded = MolecularTokenizer.load(path)
    assert loaded.tokenizer_id == t.tokenizer_id
    assert loaded.encode('CCCl') == saved
    assert loaded.decode(saved) == t.decode(saved)
    with pytest.raises(FileExistsError):
        t.save(path)

def test_ring_free_npe_training_with_windows_multiprocessing(tmp_path):
    t = MolecularTokenizer('npe')
    report = t.train(['CCO', 'CCCO', 'CCN', 'CCCl'], 8, ring_vocab_size=0, num_workers=2)
    assert report.molecules == 4 and report.ring_vocab_size == 0
    for smiles in ['CCO', 'CCN', 'CCCl']:
        assert canonical(t.decode(t.encode(smiles))) == canonical(smiles)
    loaded = MolecularTokenizer.load(t.save(tmp_path / 'ringfree'))
    assert loaded.encode('CCCl') == t.encode('CCCl')

def test_empty_edge_vocabulary_persistence(tmp_path):
    t = MolecularTokenizer('npe')
    t.train(['C'], 1, ring_vocab_size=0)
    assert not t.backend.engine.vocab_edge
    loaded = MolecularTokenizer.load(t.save(tmp_path / 'atom'))
    assert loaded.decode(loaded.encode('C')) == 'C'

def test_checksum_detects_corruption(safe, tmp_path):
    path = safe.save(tmp_path / 'corrupt')
    (path / 'sequence.json').write_text('{}', encoding='utf-8')
    with pytest.raises(TokenizerError, match='checksum'):
        MolecularTokenizer.load(path)

def test_vocab_and_backend_binding(safe, npe):
    with pytest.raises(TokenizerError, match='different backend'):
        npe.decode(safe.encode('CCO'))
    other = MolecularTokenizer('safe')
    other.train(['CCCl', 'CCCl'], 256, min_frequency=1)
    with pytest.raises(TokenizerError, match='different backend or vocabulary'):
        other.decode(safe.encode('CCCl'))

@pytest.mark.parametrize('backend', ['safe', 'npe'])
def test_untrained_and_empty_input(backend):
    t = MolecularTokenizer(backend)
    with pytest.raises(NotTrainedError):
        t.encode('CCO')
    options = {} if backend == 'safe' else {'ring_vocab_size': 0}
    with pytest.raises(TokenizerError, match='empty'):
        t.train([], 256, **options)
    with pytest.raises(NotTrainedError):
        t.encode('CCO')

def test_failed_retraining_is_transactional(safe):
    previous = safe.tokenizer_id
    with pytest.raises(Exception):
        safe.train(['CCO', 'invalid'], 256)
    assert safe.tokenizer_id == previous
    assert safe.decode(safe.encode('CCO')) == 'CCO'

@pytest.mark.parametrize('smiles', ['', 'invalid', ' CCO'])
def test_invalid_smiles_raise(safe, smiles):
    with pytest.raises(TokenizerError):
        safe.encode(smiles)

def test_unknown_ids_and_bad_graph_positions(npe):
    encoded = npe.encode('CCCl')
    with pytest.raises(UnsupportedMoleculeError):
        npe.decode(replace(encoded, token_ids=(npe.vocab_size, *encoded.token_ids[1:])))
    edge = encoded.connections[0]
    with pytest.raises(TokenizerError, match='outside'):
        npe.decode(replace(encoded, connections=(replace(edge, source_attachment=10000),)))

def test_relaxed_mode_is_explicitly_unchecked():
    t = MolecularTokenizer.from_graphbpe(PRETRAINED, strict=False)
    encoded = t.encode('[13CH3]CO')
    assert not encoded.fidelity_checked
    assert t.decode(encoded) == 'CCO'

def test_training_rejects_bare_string(safe):
    with pytest.raises(TokenizerError):
        safe.train('CCO', 256)

def test_lazy_encode_many(safe):
    encoded = list(safe.encode_many(iter(['CCO', 'CCN'])))
    assert [safe.decode(x) for x in encoded] == ['CCO', 'CCN']

def test_invalid_serialized_payload(safe):
    encoded = safe.encode('CCO').to_dict()
    encoded['token_ids'] = [True]
    with pytest.raises(TokenizerError):
        MoleculeEncoding.from_dict(encoded)
