"""A: complete syntactic blocks; B: graph-checked protected structural blocks.

Keep original SAFE unchanged. Full-fragment lookup precedes shortest-path span
encoding. OOV blocks use framed atomic lexemes, not arbitrary substring BPE.
Whole motifs, learned spans and framed fallback blocks are logical units;
syntax tokens inside a fallback frame are NOT claimed to be chemical motifs.
"""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from rdkit import Chem, rdBase
from safe.converter import SAFEConverter
from safe.batch_chembl import identity
from demodiff_fixed_safe import LEXEME,ATOM
from evaluate_bounded_safe_vocab import Codec,load_base,dump,sha,length_stats

MAX_UNITS=16
MAX_ATOMS=32
MAX_LEXEMES=128
PROTECTION_SMARTS={
    'amide':'[CX3](=[OX1])[NX3]',
    'carboxylic_acid_or_ester':'[CX3](=[OX1])[OX2]',
}
PATTERNS={name:Chem.MolFromSmarts(s) for name,s in PROTECTION_SMARTS.items()}
BONDS=set('=#$:/\\~-')


@lru_cache(maxsize=300000)
def lex(text):
    parts=tuple(m.group() for m in LEXEME.finditer(text))
    if ''.join(parts)!=text:
        raise ValueError(f'Unsupported SAFE lexeme: {text}')
    return parts


def syntactic_units(fragment,first_atom=0):
    tokens=lex(fragment)
    depth=0
    starts=[]
    atoms={}
    atom=first_atom
    for i,t in enumerate(tokens):
        if ATOM.fullmatch(t):
            atoms[i]=atom
            atom+=1
            if depth==0:
                start=i-1 if i and tokens[i-1] in BONDS else i
                starts.append(start)
        elif t=='(':
            depth+=1
        elif t==')':
            depth-=1
            if depth<0:
                raise ValueError('Unbalanced fragment branch')
    if depth or not starts or starts[0]!=0:
        raise ValueError('Invalid fragment root or unbalanced branches')
    starts.append(len(tokens))
    units=[]
    for a,b in zip(starts,starts[1:]):
        ids=[atoms[i] for i in range(a,b) if i in atoms]
        assert ids
        units.append(dict(text=''.join(tokens[a:b]),atoms=ids))
    assert ''.join(u['text'] for u in units)==fragment
    return units,atom


def connected(ids,adjacency):
    allowed=set(ids)
    reached=set()
    todo=[ids[0]]
    while todo:
        i=todo.pop()
        if i not in reached:
            reached.add(i)
            todo.extend((adjacency[i]&allowed)-reached)
    return reached==allowed


def partition_record(safe):
    params=Chem.SmilesParserParams()
    params.removeHs=False
    params.allowCXSMILES=False
    params.parseName=False
    mol=Chem.MolFromSmiles(safe,params)
    if mol is None:
        raise ValueError('Invalid complete SAFE molecule')
    fragments=[]
    offset=0
    for fragment in safe.split('.'):
        units,offset=syntactic_units(fragment,offset)
        fragments.append(units)
    assert offset==mol.GetNumAtoms()
    adjacency=[set(a.GetIdx() for a in atom.GetNeighbors()) for atom in mol.GetAtoms()]
    owner={a:(fi,ui) for fi,units in enumerate(fragments) for ui,u in enumerate(units) for a in u['atoms']}
    groups=[('ring',list(g)) for g in mol.GetRingInfo().AtomRings()]
    for name,pattern in PATTERNS.items():
        groups.extend((name,list(g)) for g in mol.GetSubstructMatches(pattern))
    blocked=[set() for _ in fragments]
    stats=Counter()
    for name,group in groups:
        fragment_ids={owner[a][0] for a in group}
        if len(fragment_ids)>1:
            stats['cross_fragment_'+name]+=1
            continue
        fi=next(iter(fragment_ids))
        chunk_ids=[owner[a][1] for a in group]
        blocked[fi].update(range(min(chunk_ids)+1,max(chunk_ids)+1))
        stats['protected_'+name]+=1
    result={'A':[],'B':[]}
    for fi,units in enumerate(fragments):
        for u in units:
            if not connected(u['atoms'],adjacency):
                raise ValueError('Syntactic block does not induce a connected subgraph')
        for left,right in zip(units,units[1:]):
            if not any(adjacency[a]&set(right['atoms']) for a in left['atoms']):
                raise ValueError('Adjacent blocks have no graph edge')
        # Connected blocks joined by at least one edge imply every consecutive
        # span is connected. No need to repeat BFS for every vocabulary span.
        result['A'].append([u['text'] for u in units])
        merged=[]
        for ui,u in enumerate(units):
            if ui in blocked[fi]:
                merged[-1]['text']+=u['text']
                merged[-1]['atoms']+=u['atoms']
            else:
                merged.append(dict(text=u['text'],atoms=list(u['atoms'])))
        for u in merged:
            assert connected(u['atoms'],adjacency)
        b_owner={a:i for i,u in enumerate(merged) for a in u['atoms']}
        for name,group in groups:
            if all(owner[a][0]==fi for a in group):
                assert len({b_owner[a] for a in group})==1
        result['B'].append([u['text'] for u in merged])
        stats['A_blocks']+=len(units)
        stats['B_blocks']+=len(merged)
    return result,dict(stats)


@lru_cache(maxsize=300000)
def size(text):
    tokens=lex(text)
    return len(tokens),sum(bool(ATOM.fullmatch(t)) for t in tokens)


def spans(units,start):
    text=''
    nlex=0
    natoms=0
    for j in range(start,min(len(units),start+MAX_UNITS)):
        text+=units[j]
        a,b=size(units[j])
        nlex+=a
        natoms+=b
        if nlex>MAX_LEXEMES or natoms>MAX_ATOMS:
            break
        yield j+1,text,nlex


class SpanCodec:
    def __init__(self,config):
        self.config=config
        self.payloads=[None if s is None else bytes.fromhex(s) for s in config['payloads_hex']]
        self.lookup=config['content_ids']
        self.primitives=config['primitive_ids']
        self.begin=config['begin_id']
        self.end=config['end_id']
        self.dot=config['dot_id']
        self.cache={}
        assert len(self.payloads)==config['K']

    def encode_fragment(self,units):
        key=tuple(units)
        if key in self.cache:
            return self.cache[key]
        text=''.join(units)
        if text in self.lookup:
            ids=[self.lookup[text]]
        else:
            n=len(units)
            costs=[0]*(n+1)
            choices=[None]*n
            for i in range(n-1,-1,-1):
                primitive=[self.primitives[t] for t in lex(units[i])]
                fallback=[self.begin,*primitive,self.end]
                if units[i] in self.lookup:
                    fallback=[self.lookup[units[i]]]
                costs[i]=len(fallback)+costs[i+1]
                choices[i]=(i+1,fallback)
                for j,payload,_ in spans(units,i):
                    idx=self.lookup.get(payload)
                    if idx is not None and 1+costs[j]<=costs[i]:
                        # Stable tie break: prefer the longer block span.
                        costs[i]=1+costs[j]
                        choices[i]=(j,[idx])
            ids=[]
            i=0
            while i<n:
                i,part=choices[i]
                ids.extend(part)
            assert len(ids)==costs[0]
        assert self.decode(ids)==text
        self.cache[key]=ids
        return ids

    def encode(self,fragments):
        ids=[]
        for units in fragments:
            if ids:
                ids.append(self.dot)
            ids.extend(self.encode_fragment(units))
        return ids

    def decode(self,ids):
        output=[]
        inside=False
        frame=[]
        primitive_ids=set(self.primitives.values())
        for idx in ids:
            if not isinstance(idx,int) or not 0<=idx<len(self.payloads):
                raise ValueError('Invalid ID')
            if idx==self.begin:
                if inside:
                    raise ValueError('Nested block frame')
                inside=True
                frame=[]
            elif idx==self.end:
                if not inside or not frame:
                    raise ValueError('Unmatched/empty block frame')
                tokens=[self.payloads[i].decode() for i in frame]
                text=''.join(tokens)
                if list(lex(text))!=tokens:
                    raise ValueError('Broken atomic lexeme boundary')
                # The leading bond, if present, links this block to its neighbor.
                depth=0
                for t in tokens:
                    depth+=(t=='(')-(t==')')
                    if depth<0:
                        raise ValueError('Unbalanced block branch')
                if depth or not any(ATOM.fullmatch(t) for t in tokens):
                    raise ValueError('Invalid structural block')
                output.append(text)
                inside=False
            elif inside:
                if idx not in primitive_ids or idx==self.dot:
                    raise ValueError('Non-primitive token in structural frame')
                frame.append(idx)
            else:
                p=self.payloads[idx]
                if p is None:
                    raise ValueError('Special/reserved token in content')
                output.append(p.decode())
        if inside:
            raise ValueError('Unclosed block frame')
        return ''.join(output)


def make_foundation(base,train_lexemes):
    payloads=list(base.payloads)
    lookup=dict(base.base_by_text)
    extras=sorted(t for t in train_lexemes if t not in lookup)
    for text in extras:
        lookup[text]=len(payloads)
        payloads.append(text.encode())
    # Preserve ALL pre-existing complete lexical units, including rare atoms
    # absent from this TRAIN split. This is reuse, not held-out vocabulary fit.
    primitive_texts=set(train_lexemes)|{t for t in base.base_by_text if LEXEME.fullmatch(t)}
    primitives={t:lookup[t] for t in primitive_texts}
    begin=len(payloads)
    payloads.extend([None,None])
    return payloads,lookup,primitives,begin,extras


def collect_candidates(corpus,base_lookup):
    counts=Counter()
    gains=Counter()
    for number,(units,weight) in enumerate(corpus.items(),1):
        unit_costs=[1 if u in base_lookup else len(lex(u))+2 for u in units]
        for i in range(len(units)):
            cumulative=0
            for j,text,_ in spans(units,i):
                cumulative+=unit_costs[j-1]
                if text not in base_lookup:
                    counts[text]+=weight
                    gains[text]+=weight*max(0,cumulative-1)
        if number%20000==0:
            print(f'Candidate signatures {number}/{len(corpus)}, distinct={len(counts)}',flush=True)
    return counts,gains


def build_codec(K,mode,ranking,foundation,whole_ranked,candidates):
    base_payloads,base_lookup,primitives,begin,extras=foundation
    payloads=list(base_payloads)
    lookup=dict(base_lookup)
    whole_count={16000:6932,32000:7466}[K]
    for text in whole_ranked[:whole_count]:
        assert text not in lookup
        lookup[text]=len(payloads)
        payloads.append(text.encode())
    selected=[]
    for text in candidates:
        if len(payloads)==K:
            break
        if text not in lookup:
            lookup[text]=len(payloads)
            payloads.append(text.encode())
            selected.append(text)
    if len(payloads)!=K:
        raise ValueError(f'Insufficient eligible span vocabulary {mode}/{K}: {len(payloads)}')
    return SpanCodec(dict(K=K,mode=mode,ranking=ranking,
        M=len(base_payloads),original_M=2136,extra_atomic_lexemes=len(extras),
        whole_entries=whole_count,span_entries=len(selected),begin_id=begin,end_id=begin+1,
        content_ids=lookup,primitive_ids=primitives,dot_id=primitives['.'],
        payloads_hex=[None if p is None else p.hex() for p in payloads],
        limits=dict(max_units=MAX_UNITS,max_atoms=MAX_ATOMS,max_lexemes=MAX_LEXEMES)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records',type=Path,required=True)
    parser.add_argument('--base',type=Path,required=True)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise ValueError('Use a fresh output folder')
    rows=[json.loads(s) for s in args.records.read_text().splitlines()]
    assert len(rows)==44329 and len({r['record_index'] for r in rows})==44329
    assert all(r['boundary_exact'] for r in rows)
    assert Counter(r['cohort'] for r in rows)==Counter(train=35518,validation=4268,test=4543)
    tok,payloads=load_base(args.base/'base_tokenizer.json',args.base/'source_splitter.py')
    base=Codec(tok,payloads)
    whole_counts=Counter(f for r in rows if r['cohort']=='train' for f in r['safe'].split('.'))
    whole_ranked=[f for f,n in sorted(whole_counts.items(),key=lambda x:(-x[1],x[0])) if f not in base.base_by_text]
    assert whole_ranked==[json.loads(s)['fragment'] for s in (args.base/'training_fragment_ranking.jsonl').read_text().splitlines()]
    train_lexemes={t for f in whole_counts for t in lex(f)}|{'.'}
    foundation=make_foundation(base,train_lexemes)
    args.output.mkdir(parents=True)
    sources=[Path(__file__),args.records,args.base/'base_tokenizer.json',args.base/'source_splitter.py',
        args.previous/'summary.json',Path(__file__).with_name('demodiff_fixed_safe.py'),
        Path(__file__).with_name('evaluate_bounded_safe_vocab.py'),ROOT/'safe/converter.py']
    dump(args.output/'provenance.json',dict(sources={str(p.resolve()):sha(p) for p in sources},
        rdkit=rdBase.rdkitVersion,protection_smarts=PROTECTION_SMARTS,
        protection_scope='ring atom sets and named SMARTS wholly inside an original SAFE fragment only',
        A='top-level atom plus attached complete branches/ring labels; incoming bond kept with following atom',
        B='A blocks coalesced across any protected group; graph connectivity checked',
        lexical_fallback='complete atomic lexemes inside explicit BEGIN/END block; no byte/text fallback',
        graph_checks='each A block connected and consecutive blocks joined by edge; hence every contiguous span connected',
        test_note='Previously explored test split, not a pristine holdout; train-only vocabulary, validation-only selection',
        limits=dict(max_units=MAX_UNITS,max_atoms=MAX_ATOMS,max_lexemes=MAX_LEXEMES)))
    dump(args.output/'foundation.json',dict(M=len(foundation[0]),extra_lexemes=foundation[4],
        begin_id=foundation[3],end_id=foundation[3]+1,train_lexemes=sorted(train_lexemes)))
    partitions={}
    corpora={'A':Counter(),'B':Counter()}
    protection=Counter()
    start=time.monotonic()
    with (args.output/'partition_cache.jsonl').open('w',encoding='utf-8') as out:
        for row in rows:
            parts,stats=partition_record(row['safe'])
            for token in lex(row['safe']):
                if token not in foundation[2]:
                    raise ValueError(f'Unseen atomic lexeme {token!r}; do not add held-out entries')
            partitions[row['record_index']]=parts
            protection.update(stats)
            if row['cohort']=='train':
                for mode in corpora:
                    corpora[mode].update(tuple(units) for units in parts[mode])
            out.write(json.dumps(dict(record_index=row['record_index'],cohort=row['cohort'],partitions=parts,stats=stats))+'\n')
            if len(partitions)%5000==0:
                print(f'Partitions {len(partitions)}, seconds={time.monotonic()-start:.1f}',flush=True)
    dump(args.output/'partition_audit.json',dict(records=len(rows),counts=dict(protection),
        supported_atomic_lexemes=len(foundation[2]),unseen_atomic_lexemes=0))
    grid=[]
    selected={}
    all_names=[]
    for mode in ('A','B'):
        counts,gains=collect_candidates(corpora[mode],foundation[1])
        dump(args.output/f'{mode}_candidate_stats.json',dict(unique=len(counts),occurrences=sum(counts.values()),
            training_fragment_signatures=len(corpora[mode])))
        ranking_lists={
            'frequency':sorted(counts,key=lambda t:(-counts[t],t)),
            'gain':sorted(counts,key=lambda t:(-gains[t],-counts[t],t)),
        }
        with (args.output/f'{mode}_training_candidates.jsonl').open('w',encoding='utf-8') as f:
            for t in ranking_lists['frequency']:
                f.write(json.dumps(dict(text=t,train_occurrences=counts[t],estimated_gain=gains[t]))+'\n')
        for K in (16000,32000):
            options=[]
            for ranking,words in ranking_lists.items():
                name=f'{mode}_K{K}_{ranking}'
                codec=build_codec(K,mode,ranking,foundation,whole_ranked,words)
                folder=args.output/name
                folder.mkdir()
                dump(folder/'vocab.json',codec.config)
                lengths=[len(codec.encode(partitions[r['record_index']][mode])) for r in rows if r['cohort']=='validation']
                stats=length_stats(lengths)
                item=dict(name=name,mode=mode,K=K,ranking=ranking,validation=stats,
                    whole_entries=codec.config['whole_entries'],span_entries=codec.config['span_entries'],M=codec.config['M'])
                options.append(item)
                grid.append(item)
                print(json.dumps(item),flush=True)
            winner=min(options,key=lambda r:(r['validation']['total'],r['ranking']))
            selected[f'{mode}_{K}']=winner['name']
    dump(args.output/'validation_grid.json',grid)
    dump(args.output/'selection_before_test.json',dict(selected=selected,
        rule='minimum validation total length; lexical ranking-name tie break',
        candidates='frequency vs estimated primitive-block saving; same fixed whole-entry counts'))
    # Reload persisted selected vocabularies before held-out evaluation.
    codecs={name:SpanCodec(json.loads((args.output/name/'vocab.json').read_text())) for name in selected.values()}
    previous={r['record_index']:r for r in map(json.loads,(args.previous/'per_molecule.jsonl').read_text().splitlines())}
    lengths={name:defaultdict(list) for name in codecs}
    fallback=Counter()
    examples=[]
    converter=SAFEConverter()
    started=time.monotonic()
    with (args.output/'per_molecule.jsonl').open('w',encoding='utf-8') as out:
        for i,row in enumerate(rows,1):
            result={k:row[k] for k in ('record_index','chembl_id','cohort')}
            result['methods']={}
            restored=converter.decoder(row['safe'],canonical=True,fix=False,remove_dummies=False)
            assert identity(Chem.MolFromSmiles(restored))==identity(Chem.MolFromSmiles(row['canonical_smiles']))
            for name,codec in codecs.items():
                parts=partitions[row['record_index']][codec.config['mode']]
                ids=codec.encode(parts)
                assert codec.decode(ids)==row['safe']
                assert all(0<=k<codec.config['K'] for k in ids)
                frames=ids.count(codec.begin)
                assert frames==ids.count(codec.end)
                fallback[name]+=frames
                lengths[name][row['cohort']].append(len(ids))
                result['methods'][name]=dict(length=len(ids),ids=ids,fallback_blocks=frames,
                    text_roundtrip=True,molecular_roundtrip=True)
                if row['record_index']==39224:
                    detail=[]
                    for units in parts:
                        part_ids=codec.encode_fragment(units)
                        display=['<BLOCK_BEGIN>' if k==codec.begin else '<BLOCK_END>' if k==codec.end else codec.payloads[k].decode() for k in part_ids]
                        detail.append(dict(fragment=''.join(units),logical_units=units,ids=part_ids,tokens=display))
                    examples.append(dict(name=name,length=len(ids),fragments=detail))
            result['previous_lengths']={k:previous[row['record_index']]['methods'][v]['length'] for k,v in {
                'hybrid_16000':'historical_16000','hybrid_32000':'historical_32000',
                'lexical_16000':'constrained_K16000_whole50','lexical_32000':'constrained_K32000_whole25'}.items()}
            out.write(json.dumps(result)+'\n')
            if i%5000==0:
                print(f'Evaluated {i}, seconds={time.monotonic()-started:.1f}',flush=True)
    summary={}
    for cohort in ('all','train','validation','test'):
        summary[cohort]={}
        for name,values in lengths.items():
            vals=[n for group in values.values() for n in group] if cohort=='all' else values[cohort]
            summary[cohort][name]=dict(length_stats(vals),exact_safe_roundtrips=len(vals),molecular_roundtrips=len(vals))
    dump(args.output/'summary.json',summary)
    dump(args.output/'fallback_counts.json',dict(fallback))
    dump(args.output/'CHEMBL45183.json',examples)
    print(json.dumps(dict(selected=selected,all=summary['all'],test=summary['test']),indent=2),flush=True)


if __name__=='__main__':
    main()
