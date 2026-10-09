"""Render completed all-row tokenizer audits; refuses partial source selections."""
import argparse
import csv
import html
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
try:
    from .result_insights import summarize
except ImportError:
    from result_insights import summarize


def distribution(histogram):
    pairs=sorted((int(k),v) for k,v in histogram.items())
    return np.array([k for k,v in pairs]),np.array([v for k,v in pairs])


def render(args):
    data=json.loads(args.summary.read_text(encoding='utf-8'))
    if data['selection']!='all rows' or data['rows']!=data['source_rows'] or not data['models']:
        raise ValueError('Final full-corpus page requires an all-row audit with trained models')
    for name,model in data['models'].items():
        provenance=model.get('training_provenance',{})
        if provenance.get('source_sha256')!=data['source_sha256'] or provenance.get('max_molecules') is not None:
            raise ValueError(f'Missing full-corpus training provenance for {name}')
        report=provenance['report']
        if report['molecules']+report['rejected_molecules']!=data['source_rows']:
            raise ValueError(f'Training row coverage mismatch for {name}')
    costs_path=args.summary.parent/'encoding_costs.json'
    costs=json.loads(costs_path.read_text()) if costs_path.exists() else None
    output=args.output
    (output/'figures').mkdir(parents=True,exist_ok=True)
    (output/'data').mkdir(exist_ok=True)
    if (output/'index.html').exists() and not (output/'preflight.html').exists():
        shutil.copyfile(output/'index.html',output/'preflight.html')
    public={k:v for k,v in data.items() if k not in ('models','definitions','chemistry_groups')}
    public['models']={}
    for name,model in data['models'].items():
        public['models'][name]={k:v for k,v in model.items() if k not in ('status','reasons','failure_examples','groups','training_provenance','exact_fraction_all_rows','exact_fraction_valid_inputs')}
        public['models'][name]['recorded_usage_molecules']=sum(model['token_histogram'].values())
    (output/'data/full_analysis.json').write_text(json.dumps(public,indent=2)+'\n',encoding='utf-8')
    labels=list(data['models'])
    insights=summarize(data,args.summary.parent)
    public_insights=json.loads(json.dumps(insights))
    for model in public_insights['models'].values():
        for key in ('exact','mismatch','error','invalid_input'): model.pop(key,None)
    (output/'data/key_findings.json').write_text(json.dumps(public_insights,indent=2)+'\n',encoding='utf-8')
    names={k:k.replace('_',' + ').upper() for k in labels}
    colors=dict(zip(labels,['#087f8c','#bd5b32','#754fa1','#559047']))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
    def save(fig,name):
        fig.savefig(output/f'figures/{name}.png',dpi=180)
        fig.savefig(output/f'figures/{name}.svg')
        svg=output/f'figures/{name}.svg'
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
        plt.close(fig)
    costs_html=''
    if costs:
        shutil.copyfile(costs_path,output/'data/encoding_costs.json')
        rows=''
        for name in labels:
            c=costs['models'][name]
            rows+=f"<tr><td>{html.escape(names[name])}</td><td>{c['encoded_molecules']:,}</td><td>{c['mean_tokens']:.3f}</td><td>{c['median_tokens']}</td><td>{c['p95_tokens']}</td><td>{c['p99_tokens']}</td><td>{c['mean_compact_json_bytes']:.2f}</td></tr>"
        fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
        for ax,key,title in zip(axs,['mean_tokens','p95_tokens'],['Mean sequence length','95th-percentile sequence length']):
            bars=ax.bar([names[k] for k in labels],[costs['models'][k][key] for k in labels],color=[colors[k] for k in labels])
            ax.bar_label(bars,fmt='%.2f',padding=4); ax.margins(y=.15)
            ax.set(title=title,ylabel='Sequence tokens')
        fig.suptitle('All successful encodings · no reconstruction-identity filter',fontsize=12)
        save(fig,'full_encoding_costs')
        block_text=''.join(f"<p>{names[k]}: mean {costs['models'][k]['mean_safe_blocks_where_recorded']:.3f} SAFE blocks across {costs['models'][k]['safe_blocks_recorded_molecules']:,} recorded molecule counts.</p>" for k in labels)
        costs_html=f'<section><h2>Encoding efficiency: all encoded molecules</h2><p>The current focus is sequence length, vocabulary use, and storage. This table includes every successfully encoded molecule. </p><div class="scroll"><table><tr><th>Pipeline</th><th>Encoded molecules</th><th>Mean tokens</th><th>Median</th><th>p95</th><th>p99</th><th>Mean JSON bytes</th></tr>{rows}</table></div><p>NPE mean sequence reduction: {costs["npe_mean_token_reduction_percent"]:.2f}%; mean compact JSON reduction: {costs["npe_mean_json_byte_reduction_percent"]:.2f}%. JSON size excludes model storage and is not a binary compression benchmark.</p><figure><img src="figures/full_encoding_costs.png" alt="Mean and p95 sequence lengths for all encoded molecules"><figcaption><a href="figures/full_encoding_costs.svg">SVG</a> · <a href="data/encoding_costs.json">All-encoding statistics, SAFE block counts and definitions</a></figcaption></figure><h3>SAFE block counts</h3>{block_text}<p>Blocks are dot-separated SAFE pieces, including disconnected components. The available block-count denominator is reported separately from the all-encoding length statistics.</p></section>'
    fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    sets=[('Atom nodes (all valid inputs)',data['atom_histogram'],'#70808b')]+[
        (names[k],(costs['models'][k]['token_histogram'] if costs else data['models'][k]['token_histogram']),colors[k]) for k in labels]
    for name,hist,color in sets:
        x,y=distribution(hist)
        if not len(x): continue
        axs[0].plot(x,y,lw=1.2,color=color,label=name)
        axs[1].step(x,np.cumsum(y[::-1])[::-1]/sum(y),where='post',color=color,label=name)
    axs[0].set(xscale='log',yscale='log',xlabel='Atom nodes / sequence tokens / graph nodes',ylabel='Molecules',title='Representation length distribution')
    axs[1].set(xscale='log',yscale='log',xlabel='Representation length',ylabel='Fraction with at least this length',title='Complementary cumulative distribution')
    for ax in axs: ax.grid(alpha=.15); ax.legend(frameon=False,fontsize=8)
    save(fig,'full_length_distributions')
    fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained')
    for k in labels:
        hist=(costs['models'][k] if costs else data['models'][k])['atom_token_ratio_histogram_tenths']
        filled={i:hist.get(str(i),0) for i in range(max(map(int,hist),default=0)+1)}
        x,y=distribution(filled)
        if len(x) and sum(y): ax.step(x/10,y/sum(y),where='post',color=colors[k],label=names[k])
    ax.axvline(1,color='#788',ls=':',label='One atom per representation unit')
    ax.set(yscale='log',xlabel='Atom count / representation-unit count (bins of 0.1)',ylabel='Fraction of recorded molecules',title='Length ratio, not byte compression')
    ax.legend(frameon=False); ax.grid(alpha=.15)
    save(fig,'full_length_ratios')
    fig,axs=plt.subplots(1,len(labels),figsize=(6*len(labels),4.8),layout='constrained',squeeze=False)
    all_bins=[b for k in labels for b in (costs['models'][k] if costs else data['models'][k])['atom_token_joint_bins5']]
    color_max=max(1,np.ceil(np.log10(max((b[2] for b in all_bins),default=1))))
    for ax,k in zip(axs.flat,labels):
        pairs=np.array((costs['models'][k] if costs else data['models'][k])['atom_token_joint_bins5'])
        if pairs.size:
            dots=ax.scatter(pairs[:,0]*5+2.5,pairs[:,1]*5+2.5,c=np.log10(pairs[:,2]),s=12,cmap='viridis',vmin=0,vmax=color_max,rasterized=True)
            fig.colorbar(dots,ax=ax,label='log10(molecules per 5 × 5 bin)')
        ax.set(xscale='log',yscale='log',xlabel='Atom nodes (bin centers)',ylabel='Representation units (bin centers)',title=names[k])
    save(fig,'full_size_vs_length')
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for k in labels:
        counts=np.array(sorted(data['models'][k]['token_usage'].values(),reverse=True))
        if not len(counts): continue
        ranks=np.arange(1,len(counts)+1)
        axs[0].plot(ranks,counts,color=colors[k],label=names[k])
        axs[1].plot(ranks,np.cumsum(counts)/sum(counts),color=colors[k],label=names[k])
    axs[0].set(xscale='log',yscale='log',xlabel='Observed token rank',ylabel='Token occurrences',title='Vocabulary rank–frequency')
    axs[1].set(xlabel='Most frequent observed tokens',ylabel='Cumulative token occurrence coverage',title='Vocabulary usage concentration',ylim=(0,1.02))
    for ax in axs: ax.grid(alpha=.15); ax.legend(frameon=False,fontsize=9)
    save(fig,'full_vocabulary_usage')
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    pairs=[data['models'][k]['common_exact_subset'] for k in labels]
    for ax,key,title in zip(axs,['tokens','payload_bytes'],['Mean representation units','Mean compact JSON bytes']):
        means=[p[key]/p['n'] if p['n'] else float('nan') for p in pairs]
        bars=ax.bar([names[k] for k in labels],means,color=[colors[k] for k in labels])
        ax.bar_label(bars,fmt='%.2f',padding=4)
        ax.margins(y=.15)
        ax.set(title=title,ylabel='Per molecule, common exact subset')
    fig.suptitle(f"Paired comparison | {pairs[0]['n']:,} molecules in the recorded comparison subset",fontsize=12)
    save(fig,'full_paired_costs')
    fig,ax=plt.subplots(figsize=(10,4.6),layout='constrained')
    bottom=np.zeros(len(labels))
    for status,color in [('exact','#087f8c'),('mismatch','#bd5b32'),('error','#ad83ac'),('invalid_input','#b4bec4')]:
        values=np.array([data['models'][k]['status'].get(status,0)/data['rows']*100 for k in labels])
        ax.bar([names[k] for k in labels],values,bottom=bottom,label=status,color=color)
        bottom+=values
    ax.set(ylabel='Percent of ALL source rows',title='Reconstruction outcome, with failures in the denominator',ylim=(0,103))
    ax.legend(frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(.5,-.10))
    plt.close(fig)
    table=[]; vocab_blocks=[]
    for k in labels:
        model=data['models'][k]; p=model['common_exact_subset']; cfg=model['configuration']
        table.append(f"<tr><td>{html.escape(names[k])}</td><td>{p['n']:,}</td><td>{cfg['vocab_size']:,}</td><td>{cfg.get('partition_vocab_size') or '—'}</td><td>{(p['tokens']/p['n']) if p['n'] else float('nan'):.3f}</td><td>{(p['payload_bytes']/p['n']) if p['n'] else float('nan'):.1f}</td><td>{cfg['saved_model_bytes']:,}</td></tr>")
        filename=model['vocabulary_file']
        text=(args.summary.parent/filename).read_text(encoding='utf-8')
        (output/'data'/filename).write_text(text.replace('exact_reconstruction_usage','recorded_occurrences'),encoding='utf-8')
        with (args.summary.parent/filename).open(encoding='utf-8') as stream:
            units=sorted(csv.DictReader(stream),key=lambda r:int(r['exact_reconstruction_usage']),reverse=True)[:12]
        items=''.join(f"<tr><td><code>{html.escape(r['unit'])}</code></td><td>{int(r['exact_reconstruction_usage']):,}</td></tr>" for r in units)
        vocab_blocks.append(f'<h3>{html.escape(names[k])}</h3><table><tr><th>Unit</th><th>Occurrences on recorded molecules</th></tr>{items}</table><p><a href="data/{filename}">Complete vocabulary CSV</a></p>')
    figure_specs=[
        ('full_length_distributions','Inspired by DemoDiff Fig. 3(b)/9. Atom counts use all valid inputs; each tokenizer curve uses its own recorded molecules. String tokens and graph nodes are different units.'),
        ('full_length_ratios','Inspired by DemoDiff Fig. 10. Ratios below one mean more representation units than atom nodes. This ratio is not byte compression.'),
        ('full_size_vs_length','Inspired by DemoDiff Fig. 11. Every exact reconstruction contributes to a 5-by-5 density bin, with logarithmic axes.'),
        ('full_paired_costs','Both methods are measured on the SAME recorded subset. Payload JSON includes IDs and five fields per graph connection; it excludes model storage and metadata. SAFE connection information is in its string.'),
        ('full_vocabulary_usage','Observed counts on recorded molecules. Unused vocabulary entries are absent from rank plots but retained in the CSV files.')]
    figures=''.join(f'<figure><img src="figures/{name}.png" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)} <a href="figures/{name}.svg">SVG</a></figcaption></figure>' for name,caption in figure_specs if not costs or name!='full_fidelity')
    groups=''
    for k in labels:
        rows=''.join(f"<tr><td>{html.escape(g)}</td><td>{counts.get('exact',0):,}</td><td>{sum(counts.values()):,}</td><td>{counts.get('mismatch',0):,}</td><td>{counts.get('error',0):,}</td></tr>" for g,counts in data['models'][k]['groups'].items())
        groups+=f'<h3>{html.escape(names[k])}</h3><table><tr><th>Overlapping subset</th><th>Exact</th><th>Valid inputs</th><th>Changed</th><th>Error</th></tr>{rows}</table>'
    findings=[]
    paired=insights['paired_comparison']
    if paired:
        delta=paired['npe_token_reduction_percent']
        direction='fewer' if delta>=0 else 'more'
        findings.append(f"On {paired['molecules']:,} molecules in the recorded comparison subset, NPE + SAFE uses {paired['npe_mean_tokens']:.3f} tokens per molecule versus {paired['brics_mean_tokens']:.3f} for BRICS + SAFE: {abs(delta):.2f}% {direction} tokens. This comparison includes the connection syntax inside SAFE strings.")
        findings.append('This result measures the complete partition-to-SAFE-to-BPE pipeline. It does not isolate the effect of partitioning alone, because the two sequence vocabularies are learned separately.')
        n=insights['models']['npe_safe']; b=insights['models']['brics_safe']
        if n['p95_tokens']>b['p95_tokens']:
            findings.append(f"The average advantage is not uniform: the NPE + SAFE p95 length is {n['p95_tokens']} tokens versus {b['p95_tokens']} for BRICS + SAFE, and p99 is {n['p99_tokens']} versus {b['p99_tokens']}. NPE has a longer high-length tail in this experiment.")
    for k, m in insights['models'].items():
        findings.append(f"{names[k]} median / p95 / p99 token lengths on the recorded subset: {m['median_tokens']} / {m['p95_tokens']} / {m['p99_tokens']}.")
    findings.append('Shorter sequences can reduce the number of positions processed by a downstream model. No downstream training-speed, generation-quality, or property-prediction improvement has been measured here.')
    verification_path=args.summary.parent/'verification.json'
    if verification_path.exists():
        verified=json.loads(verification_path.read_text())
        outcomes=verified['paired_length_outcomes']; total=sum(outcomes.values())
        findings.append(f"A separate pass over the per-molecule audit CSV verified all {verified['verified_source_rows']:,} source rows and the aggregate counts. On the common exact subset, NPE is shorter for {100*outcomes.get('npe_shorter',0)/total:.2f}% of molecules, BRICS is shorter for {100*outcomes.get('brics_shorter',0)/total:.2f}%, and {100*outcomes.get('equal',0)/total:.2f}% tie.")
    findings_html=''.join(f'<p>{html.escape(x)}</p>' for x in findings)
    vocab_detail=''
    overlap=insights['vocabulary_comparison']
    if overlap:
        vocab_detail+=f"<p>The dictionaries share {overlap['shared_string_units']:,} exact string units; union size {overlap['union_string_units']:,}, Jaccard overlap {overlap['jaccard']:.3f}. This includes special tokens and compares text, not token IDs or chemical equivalence.</p>"
        for k, entries in overlap['exclusive_most_used'].items():
            items=''.join(f"<tr><td><code>{html.escape(r['unit'])}</code></td><td>{r['occurrences']:,}</td></tr>" for r in entries)
            vocab_detail+=f'<h3>Frequent units exclusive to {html.escape(names[k])}</h3><table><tr><th>String unit</th><th>Occurrences</th></tr>{items}</table>'
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>Full ChEMBL tokenizer results</title><style>body{{margin:0;color:#18323f;background:#f4f7f8;font:16px/1.65 system-ui,sans-serif}}header{{background:#142f3c;color:white;padding:48px max(24px,calc((100vw - 1060px)/2))}}h1{{font-size:38px;line-height:1.2}}main{{max-width:1060px;margin:auto;padding:24px}}section{{background:white;border:1px solid #dce5e8;border-radius:12px;padding:26px;margin:24px 0}}a{{color:#087f8c}}header a{{color:#86dfdb}}img{{width:100%;height:auto}}figure{{margin:28px 0}}figcaption{{font-size:14px;color:#536975}}table{{width:100%;border-collapse:collapse;font-size:14px}}td,th{{padding:10px;border-bottom:1px solid #dce5e8;text-align:left}}.scroll{{overflow:auto}}code{{overflow-wrap:anywhere}}.notice{{background:#e6f3f2;padding:20px;border-radius:8px}}@media(max-width:600px){{main{{padding:12px}}section{{padding:16px}}h1{{font-size:28px}}}}</style></head><body>
    <header><p>MOLECULAR TOKENIZER · FULL-CORPUS ENCODING STATISTICS</p><h1>BRICS and NPE on all ChEMBL 37 records</h1>
    <p>{data['rows']:,} input records · {data['valid_molecules']:,} parsed molecules · {data['invalid_molecules']:,} parser failures</p>
    <a href="https://github.com/zeyuanyu195-design/molecular-tokenizer">Code and reproduction instructions</a></header><main>
    <div class="notice">This is an in-corpus descriptive experiment. Every available source SMILES was submitted to training. The same database is used for these descriptive measurements. These results do not establish held-out generalization or molecule-generation quality.</div>
    {costs_html}
    <section><h2>What the experiment shows</h2>{findings_html}<p><a href="data/key_findings.json">Download descriptive statistics</a></p></section>
    <section><h2>Measured results</h2><div class="scroll"><table><tr><th>Pipeline</th><th>Recorded comparison molecules</th><th>Sequence vocabulary</th><th>Additional NPE motifs</th><th>Paired mean units</th><th>Paired mean JSON bytes</th><th>Model bytes</th></tr>{''.join(table)}</table></div>
    <p>NPE has an additional learned motif dictionary. A matched sequence vocabulary is not matched total model capacity.
    Means are computed over the recorded common subset; coverage is reported over all source records.</p></section>
    <section><h2>Length distributions and representation cost</h2>{figures}</section>
    <section><h2>Vocabulary units and their differences</h2><p>SAFE sequence units are BPE string pieces. They need not be independently valid molecules or complete chemical fragments. BRICS uses fixed chemistry rules for partitioning; NPE learns frequent adjacent subgraphs. Both are then serialized as SAFE and receive fragment-scoped sequence BPE. SAFE stores connections inside its string, so an empty external edge table does not mean connections cost nothing.</p>{vocab_detail}<h3>Most-used units overall</h3>{''.join(vocab_blocks)}</section>
    <section><h2>Reproducibility and scope</h2><p><a href="data/full_analysis.json">Encoding statistics and model fingerprints (JSON)</a></p>
    <p>Source SHA-256: <code>{html.escape(data['source_sha256'])}</code>. Original ChEMBL 37 records are retained, including duplicates. Dataset: CC BY-SA 3.0, ChEMBL/EMBL-EBI.
    Figure design references <a href="https://arxiv.org/html/2510.08744v1">DemoDiff, Section 3.1 / Appendix B.2</a>; plots show our measurements, not copied paper results.
    No vocabulary-budget sweep or downstream neural-model experiment is claimed. Audit timestamp: {html.escape(data['created_utc'])}.</p></section></main></body></html>'''
    (output/'index.html').write_text(page,encoding='utf-8')
    (output/'.nojekyll').touch()
    print(f'Wrote completed full-corpus report: {output}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=Path('docs'))
    render(parser.parse_args())
