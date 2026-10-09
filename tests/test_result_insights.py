import csv
import pytest

from examples.result_insights import histogram_quantile, summarize


def test_quantiles_weight_counts_and_keep_empty_explicit():
    assert histogram_quantile({'10': 2, '2': 3, '100': 1}, .5) == 2
    assert histogram_quantile({'10': 2, '2': 3, '100': 1}, .95) == 100
    assert histogram_quantile({}, .5) is None


def test_findings_use_common_subset_and_string_overlap(tmp_path):
    data = {'models': {}}
    for name, units, tokens in [('npe_safe', ['[PAD]', 'CC', 'N'], 10),
                                 ('brics_safe', ['[PAD]', 'C', 'N'], 20)]:
        filename = name + '.csv'
        with (tmp_path / filename).open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['unit', 'exact_reconstruction_usage'])
            writer.writerows(zip(units, [0, 4, 2]))
        data['models'][name] = {
            'token_histogram': {'5': 2, '10': 1}, 'token_usage': {'1': 4, '2': 2},
            'vocabulary_file': filename, 'status': {'exact': 3, 'mismatch': 1},
            'mean_tokens_on_exact': 20 / 3,
            'common_exact_subset': {'n': 2, 'tokens': tokens, 'payload_bytes': tokens * 4},
        }
    result = summarize(data, tmp_path)
    assert result['paired_comparison']['npe_token_reduction_percent'] == 50
    assert result['paired_comparison']['molecules'] == 2
    assert result['vocabulary_comparison']['shared_string_units'] == 2
    assert result['vocabulary_comparison']['jaccard'] == .5
    assert result['models']['npe_safe']['unused_vocabulary_entries'] == 1
    assert result['vocabulary_comparison']['exclusive_most_used']['npe_safe'][0]['unit'] == 'CC'
    data['models']['brics_safe']['common_exact_subset']['n'] = 3
    with pytest.raises(ValueError, match='denominators'):
        summarize(data, tmp_path)
