"""Compare the trained SAFE and NPE vocabularies for ChEMBL 37."""
import json
from pathlib import Path

from molecular_tokenizer import MolecularTokenizer

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / 'trained_models'
SAFE = MODELS / 'chembl37_safe_vocab3000'
NPE = MODELS / 'chembl37_npe_vocab350ring300'
OUTPUT = ROOT / 'training_runs' / 'vocabulary_comparison.json'

safe_model = MolecularTokenizer.load(SAFE)
npe_model = MolecularTokenizer.load(NPE)
safe_vocab = safe_model.backend.tokenizer.get_vocab()
npe_tokens = npe_model.backend.engine.vocab_node
safe_set, npe_set = set(safe_vocab), set(npe_tokens)
intersection = safe_set & npe_set
union = safe_set | npe_set

summary = {
    'source_dataset': 'ChEMBL 37',
    'safe_tokenizer_id': safe_model.tokenizer_id,
    'npe_tokenizer_id': npe_model.tokenizer_id,
    'safe_vocabulary_size': len(safe_set),
    'npe_vocabulary_size': len(npe_set),
    'shared_token_strings': len(intersection),
    'jaccard_similarity': len(intersection) / len(union) if union else 1.0,
    'safe_only_examples': sorted(safe_set - npe_set)[:20],
    'npe_only_examples': sorted(npe_set - safe_set)[:20],
    'vocabularies_differ': safe_set != npe_set,
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
