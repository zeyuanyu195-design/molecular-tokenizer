"""Run full-corpus training/audit with a bounded-memory NPE checkpoint.

The orchestrator survives as a normal local process while this command is running.
If interrupted, rerun with --resume. NPE graph merges resume from committed batches;
sequence BPE and audit restart if incomplete. No machine sleep is disabled.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from profile_full_training import supervise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--metadata',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--max-rss-gib',type=float,default=12)
    parser.add_argument('--reserve-gib',type=float,default=4)
    parser.add_argument('--min-free-disk-gib',type=float,default=20)
    parser.add_argument('--timeout',type=float,default=7*24*3600)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--train-only',action='store_true',help='Finish both full-corpus models without audit, figures, or publication')
    parser.add_argument('--publish',action='store_true',help='Commit/push only the generated docs after all scientific stages succeed')
    args=parser.parse_args()
    if args.train_only and args.publish:
        parser.error('--train-only cannot be combined with --publish')
    root=Path(__file__).resolve().parents[1]
    args.input=args.input.resolve(); args.metadata=args.metadata.resolve(); args.output=args.output.resolve()
    if args.output.exists() and not args.resume:
        parser.error('Output exists; use --resume')
    args.output.mkdir(parents=True,exist_ok=True)
    expected=json.loads(args.metadata.read_text(encoding='utf-8'))
    with args.input.open('rb') as stream:
        checksum=hashlib.file_digest(stream,'sha256').hexdigest()
    if checksum!=expected['outputs']['chembl_37_smiles.smi']['sha256']:
        raise ValueError('Input differs from verified ChEMBL 37 manifest')
    state=dict(status='running',source_rows=expected['exported_rows'],source_sha256=checksum,
               started_utc=datetime.now(timezone.utc).isoformat(),steps={})
    if args.resume and (args.output/'run_state.json').exists():
        state=json.loads((args.output/'run_state.json').read_text())
        if state['source_sha256']!=checksum:
            raise ValueError('Resume input does not match run')
        state['status']='running'
    state['scope']='training_only' if args.train_only else 'training_and_analysis'
    state.pop('error',None)
    def save():
        state['updated_utc']=datetime.now(timezone.utc).isoformat()
        (args.output/'run_state.json').write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8')
    save()
    models={name:args.output/name/'model' for name in ('npe_safe','brics_safe')}
    try:
        for name,frag in [('npe_safe','npe'),('brics_safe','brics')]:
            directory=args.output/name
            directory.mkdir(exist_ok=True)
            if models[name].exists():
                from molecular_tokenizer import MolecularTokenizer
                restored=MolecularTokenizer.load(models[name])
                manifest=json.loads((models[name]/'training_report.json').read_text())
                if manifest['source_sha256']!=checksum or manifest['max_molecules'] is not None or restored.fragmentation!=frag:
                    raise ValueError(f'{name} does not contain the requested full-corpus model')
                continue
            state['current_stage']=name; save()
            command=[sys.executable,'-m','molecular_tokenizer','train','--input',str(args.input),
                '--output',str(models[name]),'--fragmentation',frag,'--representation','safe',
                '--vocab-size','3000','--num-workers',str(args.workers),'--skip-invalid']
            if frag=='npe':
                checkpoint=directory/'checkpoint'
                command+=['--motif-vocab-size','3000','--ring-vocab-size','300',
                          '--npe-storage','sqlite','--work-dir',str(checkpoint),'--batch-size','256']
                if args.resume and (checkpoint/'npe.sqlite3').exists(): command+=['--resume']
            print(f'Full-corpus stage: {name}',flush=True)
            state['steps'][name]=supervise(command,directory,args); save()
            if state['steps'][name]['exit_code']!=0:
                raise RuntimeError(f"{name} stopped: {state['steps'][name]}")
        if args.train_only:
            state['status']='complete'; state['current_stage']='complete'; save()
            print(json.dumps(state,indent=2),flush=True)
            return
        state['current_stage']='audit'; save()
        audit=args.output/'audit'
        if not (audit/'summary.json').exists():
            # Preserve partial audits by choosing a new output directory on resume.
            attempt=1
            while audit.exists():
                attempt+=1; audit=args.output/f'audit_attempt{attempt}'
            logdir=args.output/'audit_logs'; logdir.mkdir(exist_ok=True)
            command=[sys.executable,str(root/'examples/analyze_corpus.py'),'--input',str(args.input),
                     '--metadata',str(args.metadata),'--output',str(audit),'--workers',str(args.workers),
                     '--scope','training_corpus']
            for name,path in models.items(): command+=['--model',f'{name}={path}']
            state['steps']['audit']=supervise(command,logdir,args); save()
            if state['steps']['audit']['exit_code']!=0:
                raise RuntimeError('Full audit failed or exceeded a resource limit')
        state['audit_directory']=str(audit)
        state['current_stage']='render'; save()
        subprocess.run([sys.executable,str(root/'examples/render_full_results.py'),
                        '--summary',str(audit/'summary.json'),'--output',str(root/'docs')],check=True)
        if args.publish:
            state['current_stage']='publish'; save()
            staged=subprocess.check_output(['git','diff','--cached','--name-only'],cwd=root,text=True)
            if staged.strip():
                raise RuntimeError('Other staged work exists; preserve it and publish results manually')
            subprocess.run(['git','add','--','docs'],cwd=root,check=True)
            if subprocess.run(['git','diff','--cached','--quiet'],cwd=root).returncode:
                subprocess.run(['git','commit','-m','Publish completed full ChEMBL tokenizer analysis'],cwd=root,check=True)
                subprocess.run(['git','push','origin','main'],cwd=root,check=True)
        state['status']='complete'; state['current_stage']='complete'; save()
        print(json.dumps(state,indent=2),flush=True)
    except BaseException as error:
        state['status']='stopped'; state['error']=str(error); save()
        raise


if __name__=='__main__': main()
