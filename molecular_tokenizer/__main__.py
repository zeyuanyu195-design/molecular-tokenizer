"""Portable CLI: python -m molecular_tokenizer train|encode|evaluate."""
import argparse
from contextlib import redirect_stdout
from dataclasses import asdict
import hashlib
import itertools
import json
from pathlib import Path
import sys

from . import MolecularTokenizer, TokenizerError


def samples(path, limit=None):
    # One unmodified SMILES per line; blank/invalid rows are handled explicitly.
    with Path(path).open(encoding='utf-8') as stream:
        for line in itertools.islice(stream, limit):
            yield line.rstrip('\r\n')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    train=commands.add_parser('train')
    train.add_argument('--input',required=True,type=Path)
    train.add_argument('--output',required=True,type=Path)
    train.add_argument('--fragmentation',choices=['brics','npe'],default='npe')
    train.add_argument('--representation',choices=['safe','demodiff'],default='safe')
    train.add_argument('--bpe-scope',choices=['lexical','fragment'])
    train.add_argument('--vocab-size',type=int,default=3000)
    train.add_argument('--motif-vocab-size',type=int,default=350)
    train.add_argument('--ring-vocab-size',type=int,default=300)
    train.add_argument('--npe-model',type=Path)
    train.add_argument('--num-workers',type=int,default=1)
    train.add_argument('--max-molecules',type=int)
    train.add_argument('--skip-invalid',action='store_true')
    encode=commands.add_parser('encode')
    encode.add_argument('--model',required=True,type=Path)
    encode.add_argument('--smiles',required=True)
    evaluate=commands.add_parser('evaluate')
    evaluate.add_argument('--model',required=True,type=Path)
    evaluate.add_argument('--input',required=True,type=Path)
    evaluate.add_argument('--output',required=True,type=Path)
    evaluate.add_argument('--max-molecules',type=int)
    args=parser.parse_args(argv)
    if getattr(args,'max_molecules',None) is not None and args.max_molecules<1:
        parser.error('--max-molecules must be positive')
    if args.command=='train':
        if args.output.exists(): parser.error('Output directory already exists')
        if args.npe_model is not None and (args.fragmentation,args.representation)!=('npe','safe'):
            parser.error('--npe-model requires --fragmentation npe --representation safe')
        model=MolecularTokenizer(fragmentation=args.fragmentation,representation=args.representation,
                                  bpe_scope=args.bpe_scope)
        rejected=[]
        options=dict(num_workers=args.num_workers,skip_invalid=args.skip_invalid,
                     on_reject_with_index=lambda smi,reason,index:rejected.append(
                         {'source_row':index,'smiles':smi,'reason':reason}))
        if args.fragmentation=='npe':
            options['ring_vocab_size']=args.ring_vocab_size
            if args.representation=='safe':
                options.update(motif_vocab_size=args.motif_vocab_size,npe_model=args.npe_model)
        with redirect_stdout(sys.stderr):
            report=model.train(samples(args.input,args.max_molecules),args.vocab_size,**options)
        model.save(args.output)
        with args.input.open('rb') as stream:
            checksum=hashlib.file_digest(stream,'sha256').hexdigest()
        metadata={'report':asdict(report),'source_file':args.input.name,'source_sha256':checksum,
                  'selection':'source order; first max_molecules rows if limited',
                  'max_molecules':args.max_molecules,'rejected':rejected,
                  'tokenizer_id':model.tokenizer_id}
        (args.output/'training_report.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
        print(json.dumps(metadata,indent=2))
    elif args.command=='encode':
        model=MolecularTokenizer.load(args.model)
        result=model.encode(args.smiles)
        output={'encoding':result.to_dict(),'decoded':model.decode(result)}
        if model.representation=='safe': output['safe_string']=model.to_safe(args.smiles)
        print(json.dumps(output,indent=2))
    else:
        if args.output.exists(): parser.error('Output file already exists')
        model=MolecularTokenizer.load(args.model)
        model.strict=True
        count=success=tokens=connections=0
        failures=[]
        for row,smi in enumerate(samples(args.input,args.max_molecules),1):
            count+=1
            try:
                encoded=model.encode(smi)
                success+=1
                tokens+=len(encoded.token_ids)
                connections+=len(encoded.connections)
            except Exception as error:
                failures.append({'source_row':row,'smiles':smi,'error':f'{type(error).__name__}: {error}'})
        if not count: raise TokenizerError('Evaluation corpus is empty')
        output={'molecules':count,'strict_roundtrip_success':success,
                'mean_tokens_on_successes':tokens/success if success else None,
                'mean_connections_on_successes':connections/success if success else None,
                'tokenizer_id':model.tokenizer_id,'failures':failures}
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x',encoding='utf-8') as stream: json.dump(output,stream,indent=2)
        print(json.dumps({k:v for k,v in output.items() if k!='failures'},indent=2))
    return 0


if __name__=='__main__':
    try: sys.exit(main())
    except (TokenizerError,ValueError,OSError) as error:
        print(f'{type(error).__name__}: {error}',file=sys.stderr)
        sys.exit(1)
