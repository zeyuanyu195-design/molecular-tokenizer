"""Freeze whole-fragment K16k IDs; append train-only A/B spans at K20k/K24k."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
from types import SimpleNamespace
import time

from evaluate_structural_span_vocab import (
    SpanCodec, make_foundation, partition_record, collect_candidates, lex,
    MAX_UNITS, MAX_ATOMS, MAX_LEXEMES, PROTECTION_SMARTS,
)
from evaluate_bounded_safe_vocab import Codec, load_base, load_vocabulary, dump, sha, length_stats
from report_structural_span_vocab import logical_outputs, check_content
from rdkit import Chem, rdBase
from safe.converter import SAFEConverter
from safe.batch_chembl import identity


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def append_codec(payloads, lookup, train_lexemes, ranked, extra, mode, source_sha):
    foundation=make_foundation(SimpleNamespace(payloads=payloads,base_by_text=lookup),train_lexemes)
    values,content,primitives,begin,extras=foundation
    K=len(payloads)+extra
    spans=[]
    for text in ranked:
        if len(values)>=K:
            break
        if text not in content:
            content[text]=len(values)
            values.append(text.encode())
            spans.append(text)
    assert len(values)==K, 'Insufficient eligible candidates for exact total budget'
    assert values[:len(payloads)]==payloads
    assert all(content[t]==i for t,i in lookup.items())
    return SpanCodec(dict(K=K,mode=mode,ranking='frequency_on_train_whole_misses',
        frozen_size=len(payloads),append_budget=extra,original_M=2136,
        whole_entries=13864,extra_atomic_lexemes=len(extras),span_entries=len(spans),
        begin_id=begin,end_id=begin+1,dot_id=primitives['.'],primitive_ids=primitives,
        content_ids=content,payloads_hex=[None if x is None else x.hex() for x in values],
        frozen_vocabulary_sha256=source_sha,
        limits=dict(max_units=MAX_UNITS,max_atoms=MAX_ATOMS,max_lexemes=MAX_LEXEMES)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--records',type=Path,required=True)
    p.add_argument('--whole',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():
        raise ValueError('Use a fresh output folder; prior results are preserved')
    rows=[json.loads(x) for x in args.records.read_text(encoding='utf-8').splitlines()]
    assert len(rows)==len({r['record_index'] for r in rows})==44329
    assert all(r['boundary_exact'] for r in rows)
    assert Counter(r['cohort'] for r in rows)==Counter(train=35518,validation=4268,test=4543)
    source=args.whole/'K16000/vocab.json'
    frozen=load_vocabulary(source)
    config=read(source)
    assert len(frozen)==config['actual_size']==16000 and config['M']==2136
    tok,base_payloads=load_base(args.whole/'base_tokenizer.json',args.whole/'source_splitter.py')
    baseline=Codec(tok,base_payloads)
    selected={text:config['M']+i for i,text in enumerate(config['added_fragments'])}
    assert len(selected)==13864
    lookup={**baseline.base_by_text,**selected}
    baseline.payloads=frozen
    train_lexemes={'.'}|{t for r in rows if r['cohort']=='train' for t in lex(r['safe'])}
    foundation=make_foundation(SimpleNamespace(payloads=frozen,base_by_text=lookup),train_lexemes)
    args.output.mkdir(parents=True)
    sources=[args.records,source,args.whole/'base_payloads.json',args.whole/'base_tokenizer.json',
        args.whole/'source_splitter.py',Path(__file__),Path(__file__).with_name('evaluate_structural_span_vocab.py'),
        Path(__file__).with_name('evaluate_bounded_safe_vocab.py'),Path(__file__).parents[2]/'safe/converter.py']
    dump(args.output/'provenance.json',dict(sources={str(x.resolve()):sha(x) for x in sources},
        rdkit=rdBase.rdkitVersion,protection_smarts=PROTECTION_SMARTS,
        training='Only TRAIN fragments missing the frozen whole-fragment lookup; occurrence-frequency ranking fixed before TEST',
        scope='44329 boundary-exact records, same pre-existing splits; test split previously explored',
        budget='Keep all 16000 original IDs, including unused legacy entries. Added lexical units and frame markers count toward 4000/8000.',
        semantics='Logical blocks preserve A/B boundaries; fallback syntax IDs are not standalone motifs. No byte fallback.'))
    dump(args.output/'foundation.json',dict(frozen_size=16000,extra_lexemes=foundation[4],
        begin_id=foundation[3],end_id=foundation[3]+1))
    parts={}; corpora={'A':Counter(),'B':Counter()}; protection=Counter()
    start=time.monotonic()
    with (args.output/'partition_cache.jsonl').open('w',encoding='utf-8') as out:
        for i,r in enumerate(rows,1):
            part,stats=partition_record(r['safe'])
            parts[r['record_index']]=part
            protection.update(stats)
            assert all(t in foundation[2] for t in lex(r['safe'])), 'Unsupported lexeme; no held-out vocabulary additions'
            if r['cohort']=='train':
                for mode in corpora:
                    corpora[mode].update(tuple(u) for u in part[mode] if ''.join(u) not in lookup)
            out.write(json.dumps(dict(record_index=r['record_index'],partitions=part))+'\n')
            if i%5000==0:
                print(f'Partitions {i}/44329, {time.monotonic()-start:.1f}s',flush=True)
    dump(args.output/'partition_audit.json',dict(records=44329,protection_counts=dict(protection)))
    models={}
    for mode in ('A','B'):
        counts,_=collect_candidates(corpora[mode],foundation[1])
        ranked=sorted(counts,key=lambda t:(-counts[t],t))
        with (args.output/f'{mode}_training_candidates.jsonl').open('w',encoding='utf-8') as out:
            for text in ranked:
                out.write(json.dumps(dict(text=text,train_occurrences=counts[text]))+'\n')
        for extra in (4000,8000):
            name=f'{mode}_plus{extra}'
            c=append_codec(frozen,lookup,train_lexemes,ranked,extra,mode,sha(source))
            folder=args.output/name; folder.mkdir()
            dump(folder/'vocab.json',c.config)
            models[name]=SpanCodec(read(folder/'vocab.json'))
            print(f'{name}: K={c.config["K"]}, appended lexemes={c.config["extra_atomic_lexemes"]}, markers=2, spans={c.config["span_entries"]}',flush=True)
        assert models[f'{mode}_plus8000'].payloads[:20000]==models[f'{mode}_plus4000'].payloads
    dump(args.output/'frozen_before_evaluation.json',dict(models={n:sha(args.output/n/'vocab.json') for n in models},
        selection='No test-driven selection; fixed TRAIN frequency ranking for all four configurations'))
    lengths={n:defaultdict(list) for n in ['whole_16000',*models]}
    coverage={n:Counter() for n in models}
    comparisons={n:Counter() for n in models}
    examples=[]; decoder=SAFEConverter()
    with (args.output/'per_molecule.jsonl').open('w',encoding='utf-8') as out:
        for i,r in enumerate(rows,1):
            old_ids,*_=baseline.encode(r['safe'],selected)
            assert b''.join(frozen[j] for j in old_ids).decode()==r['safe']
            lengths['whole_16000'][r['cohort']].append(len(old_ids))
            result={k:r[k] for k in ('record_index','chembl_id','cohort')}
            result.update(baseline_length=len(old_ids),methods={})
            for name,c in models.items():
                part=parts[r['record_index']][c.config['mode']]
                ids=c.encode(part)
                restored=c.decode(ids)
                assert restored==r['safe']
                assert all(0<=j<c.config['K'] for j in ids)
                frames=ids.count(c.begin)
                coverage[name].update(records=1,fragments=len(part),fallback_blocks=frames,
                    molecules_without_fallback=int(frames==0))
                for units in part:
                    text=''.join(units); fids=c.encode_fragment(units)
                    if text in lookup:
                        assert fids==[lookup[text]], 'Frozen complete hit was changed or split'
                        coverage[name]['frozen_full_fragment_hits']+=1
                    coverage[name]['single_token_fragments']+=len(fids)==1
                    coverage[name]['fragments_without_fallback']+=c.begin not in fids
                lengths[name][r['cohort']].append(len(ids))
                comparisons[name]['shorter' if len(ids)<len(old_ids) else 'longer' if len(ids)>len(old_ids) else 'equal']+=1
                result['methods'][name]=dict(length=len(ids),ids=ids,fallback_blocks=frames,text_roundtrip=True)
                if r['record_index']==39224:
                    examples.append(dict(name=name,length=len(ids),safe=r['safe'],fragments=[dict(
                        fragment=''.join(u),logical_units=u,ids=c.encode_fragment(u),
                        tokens=['<BLOCK_BEGIN>' if j==c.begin else '<BLOCK_END>' if j==c.end else c.payloads[j].decode() for j in c.encode_fragment(u)]) for u in part]))
            # Each of the four reconstructed strings equals this exact SAFE.
            # Decode that shared string once; no repair/dummy removal is allowed.
            molecule=decoder.decoder(restored,canonical=True,fix=False,remove_dummies=False)
            assert identity(Chem.MolFromSmiles(molecule))==identity(Chem.MolFromSmiles(r['canonical_smiles']))
            result['molecular_roundtrip']=True
            for mode in ('A','B'):
                assert result['methods'][f'{mode}_plus8000']['length']<=result['methods'][f'{mode}_plus4000']['length']
            out.write(json.dumps(result)+'\n')
            if i%5000==0:
                print(f'Evaluated {i}/44329, {time.monotonic()-start:.1f}s',flush=True)
    summary={cohort:{n:length_stats([v for group in values.values() for v in group] if cohort=='all' else values[cohort]) for n,values in lengths.items()} for cohort in ('all','train','validation','test')}
    assert summary['all']['whole_16000']['total']==1555440, 'Baseline or cohort drift'
    dump(args.output/'summary.json',summary)
    dump(args.output/'coverage.json',coverage)
    dump(args.output/'paired_comparison.json',comparisons)
    dump(args.output/'CHEMBL45183.json',examples)
    # Independent persisted-ID audit: reload dictionaries and saved output, then
    # check every logical cut against the independently built A/B partitions.
    models={n:SpanCodec(read(args.output/n/'vocab.json')) for n in models}
    index={r['record_index']:r for r in rows}; seen=set(); totals=Counter(); components={n:Counter() for n in models}
    for line in (args.output/'per_molecule.jsonl').open(encoding='utf-8'):
        r=json.loads(line); idx=r['record_index']; assert idx not in seen; seen.add(idx)
        for name,value in r['methods'].items():
            c=models[name]; assert c.payloads[:16000]==frozen
            assert c.decode(value['ids'])==index[idx]['safe']
            assert len(value['ids'])==value['length']
            allowed=[]
            for units in parts[idx][c.config['mode']]:
                pos=0; boundaries={0}
                for u in units: pos+=len(u); boundaries.add(pos)
                allowed.append(boundaries)
            fi=pos=0
            for text,frame in logical_outputs(c,value['ids']):
                if text=='.':
                    assert pos==max(allowed[fi]); fi+=1; pos=0; components[name]['dots']+=1
                    continue
                check_content(text); pos+=len(text); assert pos in allowed[fi]
                if frame:
                    components[name]['frame_markers']+=2
                    components[name]['primitive_ids']+=len(lex(text))
                else: components[name]['literal_ids']+=1
            assert fi==len(allowed)-1 and pos==max(allowed[fi])
            totals[name]+=value['length']
    assert seen==set(index)
    for n,c in models.items():
        assert totals[n]==sum(components[n].values())==summary['all'][n]['total']
        for text,i in c.lookup.items():
            if i>c.end: check_content(text)
    for path,expected in read(args.output/'provenance.json')['sources'].items(): assert sha(path)==expected
    dump(args.output/'verification.json',dict(records=44329,verified_id_sequences=177316,
        exact_safe_roundtrips=177316,molecular_roundtrip_records=44329,
        molecular_check='Each of four ID decodings equals the same source SAFE; SAFE decoder and canonical isomeric identity checked once per record.',
        unchanged_first_16000_payloads=True,unchanged_frozen_fragment_ids=True,
        logical_cut_violations=0,content_syntax_violations=0,unseen_lexemes=0,
        nested_vocabularies=True,length_components=components,per_molecule_sha256=sha(args.output/'per_molecule.jsonl')))
    print(json.dumps(dict(all=summary['all'],test=summary['test'],verification='passed'),indent=2),flush=True)


if __name__=='__main__':
    main()
