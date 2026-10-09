"""Descriptive statistics for a completed audit; never infer held-out quality."""
import csv
import math
from pathlib import Path


def histogram_quantile(histogram, q):
    """Nearest-rank quantile of an exact integer-frequency histogram."""
    total = sum(histogram.values())
    if not total:
        return None
    target = max(1, math.ceil(q * total))
    cumulative = 0
    for value, count in sorted((int(k), v) for k, v in histogram.items()):
        cumulative += count
        if cumulative >= target:
            return value


def summarize(data, directory):
    results = {'scope': 'In-corpus descriptive analysis; not held-out evaluation',
               'models': {}, 'vocabulary_comparison': None, 'paired_comparison': None}
    vocabularies = {}
    for name, model in data['models'].items():
        histogram = model['token_histogram']
        counts = sorted(model['token_usage'].values(), reverse=True)
        total = sum(counts)
        filename = Path(directory) / model['vocabulary_file']
        with filename.open(encoding='utf-8') as stream:
            units = list(csv.DictReader(stream))
        vocabularies[name] = {r['unit']: int(r['exact_reconstruction_usage']) for r in units}
        results['models'][name] = {
            'exact': model['status'].get('exact', 0),
            'mismatch': model['status'].get('mismatch', 0),
            'error': model['status'].get('error', 0),
            'invalid_input': model['status'].get('invalid_input', 0),
            'mean_tokens': model['mean_tokens_on_exact'],
            'median_tokens': histogram_quantile(histogram, .5),
            'p95_tokens': histogram_quantile(histogram, .95),
            'p99_tokens': histogram_quantile(histogram, .99),
            'maximum_tokens': max(map(int, histogram), default=None),
            'used_vocabulary_entries': len(counts),
            'unused_vocabulary_entries': len(units) - len(counts),
            'top100_occurrence_fraction': sum(counts[:100]) / total if total else None,
            'vocabulary_entries': len(units),
        }
    if {'npe_safe', 'brics_safe'} <= data['models'].keys():
        n = data['models']['npe_safe']['common_exact_subset']
        b = data['models']['brics_safe']['common_exact_subset']
        if n['n'] != b['n']:
            raise ValueError('Paired model statistics must have identical denominators')
        if n['n']:
            results['paired_comparison'] = {
                'molecules': n['n'],
                'npe_mean_tokens': n['tokens'] / n['n'],
                'brics_mean_tokens': b['tokens'] / b['n'],
                'npe_token_reduction_percent': 100 * (1 - n['tokens'] / b['tokens']),
                'npe_json_byte_reduction_percent': 100 * (1 - n['payload_bytes'] / b['payload_bytes']),
                'definition': 'Positive reduction means NPE is shorter; ratio of aggregate totals on the common exact subset.',
            }
        a, b = vocabularies['npe_safe'], vocabularies['brics_safe']
        shared = a.keys() & b.keys()
        union = a.keys() | b.keys()
        exclusive = {}
        for name, own, other in [('npe_safe', a, b), ('brics_safe', b, a)]:
            entries = sorted(((unit, count) for unit, count in own.items() if unit not in other),
                             key=lambda item: (-item[1], item[0]))
            exclusive[name] = [{'unit': unit, 'occurrences': count} for unit, count in entries[:12]]
        results['vocabulary_comparison'] = {
            'shared_string_units': len(shared), 'union_string_units': len(union),
            'jaccard': len(shared) / len(union) if union else None,
            'exclusive_most_used': exclusive,
            'definition': 'Exact string equality, including special tokens. Shared text need not have the same token ID or usage context. Usage counts include exact reconstructions only.',
        }
    return results
