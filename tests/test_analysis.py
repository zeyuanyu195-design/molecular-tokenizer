"""Check audit denominators, paired comparisons, and identity handling."""
from examples.analyze_corpus import Aggregate, analyze_row, initialize
from molecular_tokenizer import MolecularTokenizer


def test_failures_stay_in_denominator_and_pairs_require_every_model():
    stats = Aggregate(['a', 'b'])
    exact = dict(status='exact', tokens=2, connections=1, payload_bytes=40, ids=[1, 2])
    stats.add(dict(row=1, input_valid=True, smiles_characters=8, atoms=6,
                   flags=['isotope', 'charge'], models={'a':exact, 'b':exact}))
    stats.add(dict(row=2, input_valid=True, smiles_characters=8, atoms=6,
                   flags=['charge'], models={'a':exact, 'b':dict(status='mismatch', reason='changed')}))
    stats.add(dict(row=3, input_valid=False, smiles_characters=0, models={}))
    result = stats.result()
    assert result['rows'] == 3 and result['valid_molecules'] == 2
    assert result['models']['a']['exact_fraction_all_rows'] == 2/3
    assert result['models']['b']['exact_fraction_valid_inputs'] == 1/2
    assert result['models']['a']['common_exact_subset']['n'] == 1
    assert result['models']['a']['token_usage'] == {1:2, 2:2}
    assert result['models']['b']['token_usage'] == {1:1, 2:1}
    assert result['models']['b']['groups']['charge'] == {'exact':1, 'mismatch':1}
    assert result['models']['a']['mean_atom_token_ratio_on_exact'] == 3


def test_audit_decodes_isotopes_and_identifies_invalid_source(tmp_path):
    model = MolecularTokenizer(fragmentation='brics', representation='safe')
    model.train(['CCO','CCN','[13CH3]CO'], 128)
    model.save(tmp_path/'model')
    initialize([('safe', tmp_path/'model')])
    valid = analyze_row((1, '[13CH3]CO'))
    assert valid['atoms'] == 3
    assert valid['flags'] == ['isotope']
    assert valid['models']['safe']['status'] == 'exact'
    assert valid['models']['safe']['connections'] == 0
    assert valid['models']['safe']['payload_bytes'] > 0
    invalid = analyze_row((2, 'not a smiles'))
    assert invalid['input_valid'] is False and invalid['models'] == {}
