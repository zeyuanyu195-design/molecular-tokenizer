"""Fresh, matched SAFE-pipeline training on a reproducible existing audit sample.

Requires the local analysis_units/sample_selection.json from the ChEMBL audit.
Models and reports are local outputs, excluded from source control.
"""
from dataclasses import asdict
import json
from pathlib import Path
import random
import statistics
import time
from molecular_tokenizer import MolecularTokenizer

ROOT=Path(__file__).resolve().parents[1]

def main():
    rows=json.loads((ROOT/'analysis_units/sample_selection.json').read_text(encoding='utf-8'))
    # Split by canonical SMILES identity, never duplicate molecules across splits.
    from molecular_tokenizer.backends import canonical_molecule
    unique={canonical_molecule(r['canonical_smiles']):r for r in rows}
    molecules=sorted(unique)
    random.Random(20261009).shuffle(molecules)
    train,test=molecules[:1600],molecules[1600:]
    output=ROOT/'validation/composable_v02'
    if output.exists(): raise FileExistsError(output)
    output.mkdir(parents=True)
    (output/'split.json').write_text(json.dumps({'seed':20261009,'train':train,'test':test},indent=2),encoding='utf-8')
    summaries={}
    for fragmentation in ['brics','npe']:
        model=MolecularTokenizer(fragmentation=fragmentation,representation='safe',bpe_scope='fragment')
        options={'motif_vocab_size':350,'ring_vocab_size':100,'num_workers':4} if fragmentation=='npe' else {}
        start=time.perf_counter()
        report=model.train(train,1024,**options)
        elapsed=time.perf_counter()-start
        directory=model.save(output/f'{fragmentation}_safe')
        loaded=MolecularTokenizer.load(directory)
        results=[]
        for smiles in test:
            row={'smiles':smiles}
            try:
                encoded=loaded.encode(smiles)
                row.update(ok=True,tokens=len(encoded.token_ids),connections=len(encoded.connections),
                           blocks=len(loaded.to_safe(smiles).split('.')))
            except Exception as error:
                row.update(ok=False,error=f'{type(error).__name__}: {error}')
            results.append(row)
        passed=[r for r in results if r['ok']]
        summaries[fragmentation]={'training_report':asdict(report),'training_seconds':elapsed,
                                 'test_molecules':len(test),'roundtrip_success':len(passed),
                                 'mean_tokens':statistics.mean(r['tokens'] for r in passed) if passed else None,
                                 'median_tokens':statistics.median(r['tokens'] for r in passed) if passed else None,
                                 'tokenizer_id':loaded.tokenizer_id}
        (output/f'{fragmentation}_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
        print(fragmentation,json.dumps(summaries[fragmentation]),flush=True)
    summaries['scope']='Fresh training; same 1600 training molecules; remaining unique molecules held out; sequence vocab target 1024. Additional NPE motif budget 350, ring budget 100. Small descriptive validation, not downstream quality evidence.'
    (output/'summary.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')

if __name__=='__main__': main()
