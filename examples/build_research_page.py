"""Build a static, dependency-free research page from measured JSON artifacts.

This preflight page deliberately does not substitute probe models for requested
full-corpus experiments. Figures are original matplotlib plots of our data.
"""
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


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def build(args):
    resources = json.loads(args.resources.read_text(encoding='utf-8'))
    dataset = json.loads(args.dataset.read_text(encoding='utf-8'))
    census = json.loads(args.census.read_text(encoding='utf-8')) if args.census else None
    disk = json.loads(args.disk_resources.read_text(encoding='utf-8')) if args.disk_resources else None
    parity = json.loads(args.disk_parity.read_text(encoding='utf-8')) if args.disk_parity else None
    running = json.loads(args.run_state.read_text(encoding='utf-8')) if args.run_state else None
    out = args.output
    (out/'figures').mkdir(parents=True, exist_ok=True)
    (out/'data').mkdir(exist_ok=True)
    (out/'.nojekyll').touch()
    cases = resources['cases']
    complete = [c for c in cases if c['resources']['exit_code']==0 and 'training' in c]
    npe = sorted([c for c in complete if c['fragmentation']=='npe' and c['motifs']==350],key=lambda c:c['rows'])
    if len(npe) < 2:
        raise ValueError('Need at least two successful NPE probe sizes')
    n = np.array([c['rows'] for c in npe],dtype=float)
    ram = np.array([c['resources']['peak_tree_rss_gib'] for c in npe])
    elapsed = np.array([c['training']['train_seconds'] for c in npe])
    slopes = [(ram[j]-ram[i])/(n[j]-n[i]) for i in range(len(n)) for j in range(i+1,len(n))]
    full = dataset['exported_rows']
    projections = [float(ram[-1]+s*(full-n[-1])) for s in slopes]
    slope, intercept = np.polyfit(n,ram,1)
    estimate = float(intercept+slope*full)
    scenarios = []
    for case in complete:
        scenarios.append(dict(fragmentation=case['fragmentation'], motifs=case['motifs'],
            probe_rows=case['rows'], full_hours_linear=case['training']['train_seconds']*full/case['rows']/3600))
    forecast = dict(scope='Extrapolation only; not a measured full training run or confidence interval',
        full_rows=full, memory_fit_gib=estimate, memory_pairwise_slope_range_gib=[min(projections),max(projections)],
        extrapolation_factor_from_largest_probe=full/n[-1],
        memory_budget_measured='350 NPE motifs, 300 rings, 3000 SAFE sequence tokens',
        hardware=resources['hardware'], linear_time_scenarios=scenarios,
        caveats=['3000-motif training has only a small resource probe, not a full-scale memory measurement.',
                 'Peak RSS sums child processes and can count shared pages more than once.',
                 'RAM scaling, rare-motif diversity, disk paging and CPU contention can change the extrapolation.',
                 'This is CPU tokenizer training; no GPU model training is included.',
                 'No bootstrap or inferential confidence intervals are claimed from single-run timings.'])
    dump(out/'data/resource_forecast.json',forecast)
    dump(out/'data/resource_measurements.json',resources)
    dump(out/'data/dataset_manifest.json',dataset)
    if census:
        shutil.copyfile(args.census,out/'data/corpus_census.json')
    disk_html=''
    if disk:
        dump(out/'data/disk_probe.json',disk)
        if parity: dump(out/'data/disk_parity.json',parity)
        case=disk['cases'][0]
        disk_html=f'''<section><h2>Disk-backed NPE: implemented and validated</h2>
        <p>Graph states, candidate frequencies and the inverted index now live in SQLite. Each batch is committed
        atomically; interrupted initialization or either phase of a merge can resume. Global frequency selection and tie-breaking retain the in-memory implementation's ordering.</p>
        <p>On the same 5,000 ChEMBL records, with 3,000 NPE motifs and 300 rings, the complete model fingerprint matches the original implementation:
        <strong>{'verified identical' if parity and parity['complete_model_fingerprint_identical'] else 'not yet verified'}</strong>.
        This probe took {case['training']['train_seconds']:.1f} seconds and peaked at {case['resources']['peak_tree_rss_gib']:.2f} GiB.
        SAFE conversion also uses workers with a frozen copy of the learned NPE engine, so its time is not directly comparable to the earlier serial-conversion probes.</p>
        <p><a href="data/disk_probe.json">Disk probe measurements</a> · <a href="data/disk_parity.json">Model parity evidence</a></p></section>'''
    if running:
        snapshot={k:v for k,v in running.items() if k in ['status','current_stage','source_rows','source_sha256','updated_utc']}
        progress=args.run_state.parent/'npe_safe/checkpoint/progress.jsonl'
        if progress.exists():
            with progress.open(encoding='utf-8') as stream:
                for line in stream:
                    try: snapshot['last_npe_event']=json.loads(line)
                    except json.JSONDecodeError: pass
        dump(out/'data/full_run_snapshot.json',snapshot)
    with (out/'data/resource_measurements.csv').open('w',encoding='utf-8',newline='') as sink:
        writer=csv.writer(sink)
        writer.writerow(['fragmentation','rows','motifs','sequence_vocab','train_seconds','peak_tree_rss_gib','status'])
        for c in cases:
            writer.writerow([c['fragmentation'],c['rows'],c['motifs'],
                c.get('training',{}).get('report',{}).get('vocab_size'),
                c.get('training',{}).get('train_seconds'),c['resources']['peak_tree_rss_gib'],
                'complete' if c['resources']['exit_code']==0 else c['resources'].get('stop_reason','failed')])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    ax[0].plot(n/1000,ram,'o-',color='#087f8c',lw=2)
    ax[0].set(xlabel='Resource-probe molecules (thousands)',ylabel='Peak process-tree RSS (GiB)',title='Measured memory')
    ax[1].plot(n/1000,elapsed/60,'o-',color='#b55431',lw=2)
    ax[1].set(xlabel='Resource-probe molecules (thousands)',ylabel='Training time (minutes)',title='Measured NPE + SAFE training')
    for axis in ax: axis.grid(alpha=.18)
    fig.suptitle('Resource probes only | NPE motifs 350, rings 300, sequence vocabulary 3,000',fontsize=12)
    fig.savefig(out/'figures/resource_scaling.png',dpi=180)
    fig.savefig(out/'figures/resource_scaling.svg')
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
    xmax=full/1e6
    xs=np.array([n[-1],full])
    ax.plot(xs/1e6,intercept+slope*xs,'--',color='#087f8c',label='Linear extrapolation (unvalidated)')
    ax.scatter(n/1e6,ram,color='#087f8c',s=45,zorder=3,label='Measured probes')
    ax.axhline(resources['hardware']['total_ram_gib'],color='#b55431',label='Installed RAM')
    ax.axhline(resources['hardware']['available_ram_gib'],color='#b55431',ls=':',label='Available RAM at probe start')
    ax.set(xlabel='Training molecules (millions)',ylabel='Process-tree RSS (GiB)',
           title='Can the current in-memory NPE trainer fit the full corpus?',xlim=(-.05,xmax*1.03),ylim=(0,max(estimate*1.2,40)))
    ax.grid(alpha=.18)
    ax.legend(loc='upper left',frameon=False)
    ax.annotate(f'{estimate:.0f} GiB estimated\n350 motifs only',xy=(xmax,estimate),
                xytext=(-150,20),textcoords='offset points',arrowprops={'arrowstyle':'->','color':'#455'})
    fig.savefig(out/'figures/full_memory_projection.png',dpi=180)
    fig.savefig(out/'figures/full_memory_projection.svg')
    plt.close(fig)
    if census:
        hist=sorted((int(k),v) for k,v in census['atom_histogram'].items())
        x=np.array([p[0] for p in hist]); y=np.array([p[1] for p in hist])
        fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        ax[0].plot(x,y,color='#087f8c',lw=1.2)
        ax[0].set(xscale='log',yscale='log',xlabel='Atom nodes per molecule',ylabel='Molecules',title='Full-corpus atom-count distribution')
        ax[1].step(x,np.cumsum(y[::-1])[::-1]/sum(y),where='post',color='#087f8c')
        ax[1].set(xscale='log',yscale='log',xlabel='Atom nodes per molecule',ylabel='Fraction with at least this many atoms',title='Full-corpus complementary CDF')
        for axis in ax: axis.grid(alpha=.18)
        fig.suptitle(f"ChEMBL 37 | {census['valid_molecules']:,} RDKit-readable records | no tokenizer compression results yet",fontsize=12)
        fig.savefig(out/'figures/corpus_atom_distribution.png',dpi=180)
        fig.savefig(out/'figures/corpus_atom_distribution.svg')
        plt.close(fig)
    rows=''.join(f"<tr><td>{c['fragmentation'].upper()} + SAFE</td><td>{c['rows']:,}</td><td>{c['motifs'] or 'Fixed rules'}</td><td>{c['training']['train_seconds']:.1f}</td><td>{c['resources']['peak_tree_rss_gib']:.2f}</td></tr>" for c in complete)
    scope='''<div class="notice"><strong>Preflight report — full-corpus tokenizer results are pending.</strong>
      <p>The requested full NPE training has not been completed. Subset runs below measure resource requirements only.
      They are not substitutes for full training, compression benchmarks, or held-out evaluation.</p></div>'''
    if running:
        scope=f'''<div class="notice"><strong>Full ChEMBL training launched with disk-backed NPE. Final tokenizer results are pending.</strong>
        <p>Run status at the published snapshot: {html.escape(running['status'])}; stage: {html.escape(running.get('current_stage','unknown'))}.
        Corpus and resource figures below are measured; full tokenizer-comparison figures will only appear after training and the full audit finish.
        <a href="data/full_run_snapshot.json">Timestamped run snapshot</a> (not a live monitor).</p></div>'''
    census_html=''
    if census:
        groups=''.join(f'<tr><td>{html.escape(k)}</td><td>{v:,}</td></tr>' for k,v in census['chemistry_groups'].items())
        histogram=sorted((int(k),v) for k,v in census['atom_histogram'].items())
        n_valid=census['valid_molecules']
        def quantile(q):
            cumulative=0
            for size,count in histogram:
                cumulative+=count
                if cumulative>=q*n_valid: return size
        mean_atoms=sum(size*count for size,count in histogram)/n_valid
        long_count=sum(count for size,count in histogram if size>100)
        census_html=f'''<section id="corpus"><h2>Full-corpus chemistry census</h2>
        <p>Every one of the {census['rows']:,} input records was visited. RDKit parsed {census['valid_molecules']:,};
        {census['invalid_molecules']:,} failed parsing. No tokenizer model is required for this census.
        Explicit atom nodes are counted; implicit hydrogens are excluded. Records are not deduplicated.</p>
        <p><strong>Observed molecule sizes:</strong> mean {mean_atoms:.2f} atom nodes, median {quantile(.5)},
        95th percentile {quantile(.95)}, 99th percentile {quantile(.99)}, maximum {histogram[-1][0]}.
        There are {long_count:,} molecules with more than 100 atoms. This long tail motivates reporting both typical lengths and tail behavior.</p>
        <figure><img src="figures/corpus_atom_distribution.png" alt="Full ChEMBL atom-count frequency and complementary cumulative distribution on logarithmic axes">
        <figcaption>Original corpus-level plots inspired by the distribution analyses in DemoDiff Fig. 9.
        These show atom counts only; motif/token curves will be added after full training.</figcaption></figure>
        <table><thead><tr><th>Overlapping molecular subset</th><th>Records</th></tr></thead><tbody>{groups}</tbody></table>
        <p>Subsets overlap. “Stereo” means specified atom chirality or bond stereochemistry;
        “charge” means at least one formally charged atom, including net-neutral zwitterions.</p>
        <p>{census['chemistry_groups'].get('stereo',0)/n_valid:.2%} of readable molecules contain specified stereochemistry.
        Fidelity on these molecules is therefore a substantial part of the evaluation, rather than just a few edge cases.
        These input statistics do not establish a compression advantage for either tokenizer.</p>
        <a href="data/corpus_census.json">Download the complete census aggregates</a></section>'''
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>Molecular Tokenizer · ChEMBL research report</title><style>
    :root{{--ink:#172c36;--muted:#526873;--teal:#087f8c;--line:#dce5e8;--paper:#f4f7f8}}
    *{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.65 system-ui,sans-serif}}
    header{{background:#142f3c;color:white;padding:54px max(24px,calc((100vw - 1050px)/2)) 40px}}
    header p{{max-width:760px;color:#c8dce2}}h1{{font-size:clamp(28px,4vw,46px);line-height:1.15;margin:12px 0 20px;letter-spacing:-.02em}}
    .eyebrow{{font-size:12px;letter-spacing:.15em;text-transform:uppercase;color:#7bddda}}nav{{display:flex;flex-wrap:wrap;gap:20px}}nav a{{color:#b6f0ed}}
    main{{max-width:1100px;margin:auto;padding:28px 24px 60px}}section{{padding:28px;margin:24px 0;background:white;border:1px solid var(--line);border-radius:12px}}
    h2{{font-size:25px;line-height:1.25;margin:0 0 18px}}h3{{margin:24px 0 8px}}a{{color:#096b79}}.notice{{background:#fff4db;border-left:5px solid #b87818;padding:22px;border-radius:6px}}
    .notice p{{margin-bottom:0}}.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:24px 0}}.card{{background:white;border:1px solid var(--line);border-radius:10px;padding:20px}}.number{{display:block;font-size:29px;font-weight:700;color:var(--teal)}}.small{{font-size:14px;color:var(--muted)}}
    table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{text-align:left;padding:12px;border-bottom:1px solid var(--line)}}th{{background:#eef4f5}}
    .scroll{{overflow:auto}}figure{{margin:24px 0}}img{{width:100%;height:auto;display:block}}figcaption{{font-size:14px;color:var(--muted);margin-top:10px}}
    pre{{overflow:auto;background:#142f3c;color:#e0f1f5;padding:20px;border-radius:8px;font-size:13px}}code{{font-family:ui-monospace,monospace;overflow-wrap:anywhere}}li{{margin:8px 0}}footer{{font-size:13px;color:var(--muted)}}
    @media(max-width:640px){{.cards{{grid-template-columns:1fr}}section{{padding:20px}}main{{padding:20px 14px}}th,td{{padding:8px}}}}
    </style></head><body><header><span class="eyebrow">Molecular Tokenizer / Research notebook / 09 October 2026</span>
    <h1>ChEMBL, from molecular graphs<br>to reusable token vocabularies.</h1>
    <p>BRICS and NPE fragmentation, paired with SAFE strings or DemoDiff graph encoding.
    Reproducible measurements, explicit reconstruction checks, and transparent resource limits.</p>
    <nav><a href="#status">Experiment status</a><a href="#resources">Resource assessment</a><a href="#design">Analysis design</a><a href="https://github.com/zeyuanyu195-design/molecular-tokenizer">Source code</a></nav></header>
    <main><div id="status">{scope}</div>
    <div class="cards"><div class="card"><span class="number">{full:,}</span>ChEMBL 37 source records</div>
    <div class="card"><span class="number">{resources['hardware']['total_ram_gib']:.1f} GiB</span>Installed system memory</div>
    <div class="card"><span class="number">{estimate:.0f} GiB*</span>Original in-memory NPE estimate<br><span class="small">*Not the new SQLite trainer</span></div></div>
    {disk_html}
    <section><h2>What has actually been run?</h2><p>The original ChEMBL 37 SMILES file has been verified against its SHA-256 manifest.
    The following models were trained on nested uniform random subsets selected from the entire source, seed {resources['seed']}.
    Each run uses {resources['workers']} workers, 300 ring types, a 3,000-token SAFE vocabulary target, and the requested NPE motif budget.
    BRICS uses fixed rules and has no motif-training budget.</p><div class="scroll"><table><thead><tr><th>Pipeline</th><th>Probe records</th><th>NPE motifs</th><th>Training seconds</th><th>Peak RSS / GiB</th></tr></thead><tbody>{rows}</tbody></table></div>
    <p class="small">Single-run wall times; end-to-end training includes SAFE conversion and sequence BPE. Sampling, imports, model saving, and smoke checks are excluded from training seconds.
    Peak RSS includes the training process and descendants. It is sampled every 0.5 seconds and can double-count shared pages.
    The first 100 training rows are used only as a smoke check, never as held-out accuracy evidence.</p></section>
    {census_html}
    <section id="resources"><h2>Full training: resource feasibility</h2>
    <figure><img src="figures/resource_scaling.png" alt="Measured NPE SAFE peak memory and training time against resource probe size"><figcaption>Measured resource scaling at 350 motifs. These are completed subset runs, not full-corpus results.</figcaption></figure>
    <figure><img src="figures/full_memory_projection.png" alt="Linear memory extrapolation to the full ChEMBL corpus against installed and available RAM"><figcaption>Dashed line: extrapolation over {full/n[-1]:.0f}× beyond the largest probe.
    The pairwise-slope scenarios give {min(projections):.0f}–{max(projections):.0f} GiB; this range is not a confidence interval.
    A 3,000-motif full run can require more memory. The current graph corpus and inverted indexes are kept in RAM.</figcaption></figure>
    <p>The original in-memory implementation exceeds this machine's available RAM at the extrapolated scale.
    The disk-backed replacement avoids retaining the complete graph corpus in Python memory and has separate measurements above when available.
    A smaller motif vocabulary or subset would change the requested experiment and is not silently substituted.</p>
    <p>Linear time scenarios for the original implementation are provided in <a href="data/resource_forecast.json">the forecast JSON</a>. They are planning estimates, not completion promises or estimates for the SQLite replacement.
    A GPU does not directly accelerate the current Python/RDKit tokenizer implementation.</p></section>
    <section id="design"><h2>Analysis design after full training</h2><div class="scroll"><table><thead><tr><th>Question</th><th>Measurements</th><th>Reference / boundary</th></tr></thead><tbody>
    <tr><td>How much shorter?</td><td>Atom nodes, graph motifs, SAFE sequence lengths; mean, median, maximum and CCDF.</td><td>DemoDiff Fig. 3(b), 9. Graph nodes and string tokens are distinct units.</td></tr>
    <tr><td>Does a larger vocabulary help?</td><td>Vocabulary-budget sweep and sequence/graph length curves.</td><td>DemoDiff Fig. 8. Each budget requires a specified training configuration.</td></tr>
    <tr><td>Which molecules compress?</td><td>Per-molecule atom-to-unit ratio and atom-count versus unit-count density.</td><td>DemoDiff Fig. 10–11. Not a byte compression claim.</td></tr>
    <tr><td>What is the connection cost?</td><td>Edge records, SAFE block counts, compact JSON payload bytes.</td><td>SAFE stores connectivity in its string; no claim of zero connection cost.</td></tr>
    <tr><td>Is the molecule preserved?</td><td>Canonical isomeric-SMILES identity; all errors and mismatches; isotope, charge, stereo subsets.</td><td>Failures remain in the denominator. Common-success subsets support paired length comparisons.</td></tr>
    <tr><td>What is learned?</td><td>Vocabulary usage, rank-frequency, entropy and token examples.</td><td>A SAFE BPE token need not be an independently valid molecular fragment.</td></tr></tbody></table></div>
    <p>With all ChEMBL records used for training, evaluation on the same database is an <strong>in-corpus descriptive analysis</strong>, not a held-out generalization experiment.
    Neural generation quality and downstream task performance are outside this tokenizer-only experiment.</p></section>
    <section><h2>Reproduce the resource assessment</h2><p>Python 3.11, Git, and the separately downloaded ChEMBL corpus are required.
    Commands below use a POSIX shell; on the existing Windows environment use <code>.\\.venv\\Scripts\\python.exe</code> instead of <code>python</code>.</p>
    <pre>python -m pip install -e ".[test,analysis]"
python examples/profile_full_training.py \\
  --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json \\
  --output training_runs/resource_probe --workers 4

# Census visits every source record; no tokenizer model is required.
python examples/analyze_corpus.py \\
  --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json \\
  --output analysis_runs/census --workers 4</pre>
    <p>Full-training and analysis commands, resource assumptions, and publication instructions are in
    <a href="https://github.com/zeyuanyu195-design/molecular-tokenizer/blob/main/RESEARCH.md">RESEARCH.md</a>.</p>
    <p><a href="data/resource_measurements.csv">Resource table (CSV)</a> · <a href="data/resource_measurements.json">Measured runs (JSON)</a> ·
    <a href="data/resource_forecast.json">Forecast assumptions (JSON)</a> · <a href="data/dataset_manifest.json">Dataset provenance (JSON)</a> ·
    <a href="figures/resource_scaling.svg">Resource figure (SVG)</a></p></section>
    <section><h2>Sources and attribution</h2><p>Liu et al., <a href="https://arxiv.org/html/2510.08744v1">Graph Diffusion Transformers are In-Context Molecular Designers</a>,
    Section 3.1 and Appendix B.2, arXiv v1 (2025). The requested <a href="https://openreview.net/pdf?id=lJ87GN5zJc">OpenReview PDF</a>
    was access-gated during review; the accessible author preprint was used. This is an adapted analysis plan, not a reproduction of their reported dataset or results.</p>
    <p>ChEMBL 37 chemical representations: <a href="{html.escape(dataset['source_url'])}">official EBI release</a>,
    CC BY-SA 3.0. Source records are preserved, including duplicates. No raw SMILES corpus or trained model is embedded in this page.
    All figures on this page were generated from the measured local data.</p></section>
    <footer>Preflight status is explicit. Full training, full tokenizer comparison, and final performance conclusions remain pending.
    Source hash: <code>{html.escape(dataset['outputs']['chembl_37_smiles.smi']['sha256'])}</code>.</footer></main></body></html>'''
    (out/'index.html').write_text(page,encoding='utf-8')
    for svg in (out/'figures').glob('*.svg'):
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(forecast,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resources',type=Path,required=True)
    parser.add_argument('--dataset',type=Path,required=True)
    parser.add_argument('--census',type=Path)
    parser.add_argument('--disk-resources',type=Path)
    parser.add_argument('--disk-parity',type=Path)
    parser.add_argument('--run-state',type=Path)
    parser.add_argument('--output',type=Path,default=Path('docs'))
    build(parser.parse_args())
