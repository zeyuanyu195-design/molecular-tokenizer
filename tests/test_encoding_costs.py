import csv
import gzip
import json

from examples.summarize_encoding_costs import summarize


def test_mismatched_encodings_remain_in_length_statistics(tmp_path):
    names = ['npe_safe', 'brics_safe']
    (tmp_path / 'summary.json').write_text(json.dumps({'rows': 2, 'source_sha256': 'test',
        'models': {k: {'configuration': {'saved_model_bytes': 10}} for k in names}}))
    with gzip.open(tmp_path / 'molecules.csv.gz', 'wt', newline='') as stream:
        fields = ['source_row', 'model', 'status', 'tokens', 'payload_bytes', 'atoms', 'safe_blocks']
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for row, status, values in [(1, 'exact', [2, 4]), (2, 'mismatch', [6, 4])]:
            for name, tokens in zip(names, values):
                writer.writerow(dict(source_row=row, model=name, status=status, tokens=tokens,
                    payload_bytes=30, atoms=5, safe_blocks=2 if row == 1 else ''))
    result = summarize(tmp_path)
    assert result['models']['npe_safe']['encoded_molecules'] == 2
    assert result['models']['npe_safe']['mean_tokens'] == 4
    assert result['paired_length_outcomes'] == {'npe_shorter': 1, 'brics_shorter': 1}
    assert result['models']['npe_safe']['safe_blocks_recorded_molecules'] == 1
