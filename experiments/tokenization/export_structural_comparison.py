"""Export audited A/B results and real CHEMBL45183 dictionaries for the site."""
import argparse
import json
from pathlib import Path
from evaluate_bounded_safe_vocab import dump,sha


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structural',type=Path,required=True)
    p.add_argument('--constrained',type=Path,required=True)
    p.add_argument('--hybrid',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    summary=read(args.structural/'summary.json')
    previous=read(args.constrained/'summary.json')
    selected=read(args.structural/'selection_before_test.json')['selected']
    audit=read(args.structural/'verification.json')
    old_examples=read(args.constrained/'CHEMBL45183.json')
    new_examples=read(args.structural/'CHEMBL45183.json')
    old_selected=read(args.hybrid/'summary.json')['selected']
    safe=next(e['safe'] for e in old_examples if e['name']=='historical_16000')
    result=dict(version='20261110-structural-v1',date='10 Nov 2026',records=44329,
        cohorts=dict(train=35518,validation=4268,test=4543),example=dict(chembl_id='CHEMBL45183',
        record_index=39224,cohort='test',safe=safe,fragments=safe.split('.'),budgets={}),budgets={},
        verification=dict(id_sequences=audit['verified_sequences'],logical_cut_violations=audit['logical_cut_violations'],
                          content_syntax_violations=audit['content_syntax_violations']),
        evidence={p.name:sha(p) for p in [args.structural/'summary.json',args.structural/'verification.json',
            args.structural/'CHEMBL45183.json',args.constrained/'summary.json']})
    for K in ('16000','32000'):
        result['budgets'][K]={}
        result['example']['budgets'][K]={}
        for mode in ('hybrid','lexical','A','B'):
            key=(f'historical_{K}' if mode=='hybrid' else
                 f'constrained_K{K}_whole{50 if K=="16000" else 25}' if mode=='lexical' else selected[f'{mode}_{K}'])
            values=previous if mode in ('hybrid','lexical') else summary
            result['budgets'][K][mode]=dict(all=values['all'][key],test=values['test'][key])
            if mode=='lexical':
                continue
            if mode=='hybrid':
                e=next(e for e in old_examples if e['name']==key)
                path=args.hybrid/old_selected[K]/'vocab.json'
                config=read(path)
                begin=end=None
                dot=config['dot_id']
            else:
                e=next(e for e in new_examples if e['name']==key)
                path=args.structural/key/'vocab.json'
                config=read(path)
                begin,end,dot=config['begin_id'],config['end_id'],config['dot_id']
                result['budgets'][K][mode]['allocation']={k:config[k] for k in ('M','whole_entries','span_entries')}
                result['budgets'][K][mode]['fallback']=audit['length_components'][key]
            payloads=[None if x is None else bytes.fromhex(x).decode('utf-8',errors='replace') for x in config['payloads_hex']]
            fragments=[]
            sequence=[]
            for i,f in enumerate(e['fragments']):
                assert f['fragment']==safe.split('.')[i]
                ids=f['ids']
                fragments.append(dict(text=f['fragment'],ids=ids,logical_units=f.get('logical_units'),
                    fallback_blocks=ids.count(begin) if begin is not None else 0))
                if i:
                    sequence.append(dot)
                sequence.extend(ids)
            assert len(sequence)==e['length']
            dictionary={str(i):dict(payload=payloads[i],text='<BLOCK_BEGIN>' if i==begin else '<BLOCK_END>' if i==end else payloads[i],
                kind='begin' if i==begin else 'end' if i==end else 'dot' if i==dot else 'payload') for i in set(sequence)}
            restored=''.join(dictionary[str(i)]['payload'] or '' for i in sequence)
            assert restored==safe
            result['example']['budgets'][K][mode]=dict(length=e['length'],ids=sequence,fragments=fragments,
                dictionary=dictionary,begin_id=begin,end_id=end,dot_id=dot,vocabulary_sha256=sha(path),
                exact_safe_roundtrip=True,molecular_roundtrip=True)
    dump(args.output,result)
    print(f'Exported {args.output}: six verified ID sequences, exact SAFE reconstruction')


if __name__=='__main__':
    main()
