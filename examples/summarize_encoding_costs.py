"""Summarize every successful encoding, without filtering by reconstruction identity."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

try:
    from .result_insights import histogram_quantile
except ImportError:
    from result_insights import histogram_quantile


def summarize(audit):
    source = json.loads((audit / 'summary.json').read_text())
    models = {k: {'n': 0, 'tokens': Counter(), 'bytes': 0, 'blocks': 0,
                   'blocks_n': 0, 'atoms': 0, 'ratios': Counter(), 'joint': Counter()} for k in source['models']}
    paired = Counter()
    pending = {}
    with gzip.open(audit / 'molecules.csv.gz', 'rt', encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream):
            name = row['model']
            # Metrics recorded before decoding remain eligible even when identity differs.
            if row['tokens']:
                stats = models[name]
                tokens = int(row['tokens'])
                stats['n'] += 1
                stats['tokens'][tokens] += 1
                stats['bytes'] += int(row['payload_bytes'])
                stats['atoms'] += int(row['atoms'])
                stats['ratios'][int(int(row['atoms']) / tokens * 10)] += 1
                stats['joint'][(int(row['atoms']) // 5, tokens // 5)] += 1
                if row['safe_blocks']:
                    stats['blocks'] += int(row['safe_blocks'])
                    stats['blocks_n'] += 1
            pending[name] = row
            if len(pending) == len(models):
                if len({r['source_row'] for r in pending.values()}) != 1:
                    raise ValueError('Per-molecule rows are not aligned')
                if {'npe_safe', 'brics_safe'} <= pending.keys() and all(r['tokens'] for r in pending.values()):
                    n = int(pending['npe_safe']['tokens'])
                    b = int(pending['brics_safe']['tokens'])
                    paired['npe_shorter' if n < b else 'brics_shorter' if n > b else 'equal'] += 1
                pending = {}
    if pending:
        raise ValueError('Incomplete final molecule')
    result = {'scope': 'All successful encodings, regardless of reconstruction identity',
              'source_rows': source['rows'], 'source_sha256': source['source_sha256'],
              'paired_length_outcomes': dict(paired), 'models': {}}
    for name, stats in models.items():
        hist = stats['tokens']; total = sum(k * v for k, v in hist.items()); n = stats['n']
        result['models'][name] = {
            'encoded_molecules': n, 'mean_tokens': total / n,
            'median_tokens': histogram_quantile(hist, .5),
            'p95_tokens': histogram_quantile(hist, .95),
            'p99_tokens': histogram_quantile(hist, .99), 'maximum_tokens': max(hist),
            'mean_compact_json_bytes': stats['bytes'] / n,
            'ratio_of_total_atoms_to_total_tokens': stats['atoms'] / total,
            'token_histogram': dict(hist),
            'atom_token_ratio_histogram_tenths': dict(stats['ratios']),
            'atom_token_joint_bins5': [[a, t, count] for (a, t), count in stats['joint'].items()],
            'mean_safe_blocks_where_recorded': stats['blocks'] / stats['blocks_n'],
            'safe_blocks_recorded_molecules': stats['blocks_n'],
            'saved_model_bytes': source['models'][name]['configuration']['saved_model_bytes'],
        }
    if {'npe_safe', 'brics_safe'} <= models.keys():
        a, b = result['models']['npe_safe'], result['models']['brics_safe']
        if a['encoded_molecules'] != b['encoded_molecules'] or sum(paired.values()) != a['encoded_molecules']:
            raise ValueError('All-encoding comparison requires a common encoded set')
        result['npe_mean_token_reduction_percent'] = 100 * (1 - a['mean_tokens'] / b['mean_tokens'])
        result['npe_mean_json_byte_reduction_percent'] = 100 * (1 - a['mean_compact_json_bytes'] / b['mean_compact_json_bytes'])
    result['definitions'] = {
        'json_bytes': 'Compact JSON encoding of IDs and external edges, excluding model storage and metadata; not an optimal binary format.',
        'safe_blocks': 'Dot-separated SAFE blocks, including disconnected components. Reported only where the original audit recorded this metric; denominator is explicit.',
        'length_scope': 'No identity-based filtering; invalid inputs and failed encodings have no token sequence.',
    }
    (audit / 'encoding_costs.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.audit)
    print(json.dumps({k: {a: b for a, b in v.items() if a != 'token_histogram'} for k, v in result['models'].items()}, indent=2))
