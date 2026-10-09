# ChEMBL 37 full-corpus experiment

Current scope (2026-10-09, resumed after training): **full-corpus audit, analysis,
figures, and GitHub Pages publication**. Both NPE + SAFE and BRICS + SAFE training
finished at 22:28 Asia/Shanghai; the models are reused without retraining.
The current audit command is the command below with `--train-only` omitted.
Publication follows inspection of the generated report and plots. For a
training-only run that deliberately skips analysis:

```bash
python examples/run_full_experiment.py --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json --output training_runs/chembl37_full_sqlite_20261009 --workers 8 --max-rss-gib 12 --reserve-gib 4 --resume --train-only
```

For a new training-only run, choose a fresh output directory and omit `--resume`.
`--train-only` cannot be combined with `--publish`. It records
`scope: training_only` in the run state and finishes after saving both models.

The requested experiment uses **all 2,897,819 available ChEMBL 37 SMILES records**.
Subset runs in the resource report are feasibility probes only. They are not the
final tokenizer experiment. No completed full NPE model or compression result is claimed until
its trained model, input manifest, and evaluation outputs have been produced.

The static research page is in `docs/index.html`. It can be opened locally without
a web server. Data behind each published figure are under `docs/data/`.

## Inputs and provenance

- Original export: [ChEMBL 37 chemical representations](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_chemreps.txt.gz).
- Dataset license: CC BY-SA 3.0. Credit ChEMBL/EMBL-EBI when redistributing derived data.
- Required local files: `../chembl_37_smiles.smi` and `../dataset_metadata.json`.
- SMILES SHA-256: `87d10f64ed147c6371ff6d385f46762ef77df6368cc68768de58df134a05cbfb`.
- Keep source order and duplicates. The source has 2,897,657 unique SMILES strings;
  row-level counts describe the released records, not unique molecular identities.
- Invalid/unconvertible rows must be reported, including source row numbers.
  A complete source pass does not imply every molecule is accepted by every codec.

## Environment

Use Python 3.11 and Git. In the repository root:

```bash
python -m pip install -e ".[test,analysis]"
python -m pytest -q
```

On the prepared Windows machine, replace `python` below with
`.\.venv\Scripts\python.exe`. Commands use fresh output directories and do not
overwrite previous models. Downloaded data, model files, and per-molecule results
remain outside version control.

## 1. Resource assessment

```bash
python examples/profile_full_training.py --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json --output training_runs/resource_probe --workers 4
```

This verifies the entire input checksum, draws nested uniform subsets (seed
20261009), then restores their source order. The default cases are:

| Fragmentation | Records | NPE motifs | Ring budget | SAFE vocabulary |
|---|---:|---:|---:|---:|
| NPE | 5,000 | 350 | 300 | 3,000 |
| NPE | 20,000 | 350 | 300 | 3,000 |
| NPE | 50,000 | 350 | 300 | 3,000 |
| NPE | 5,000 | 3,000 | 300 | 3,000 |
| BRICS | 20,000 | n/a | n/a | 3,000 |

Each case runs in a fresh process. The supervisor measures process-tree resident
memory every 0.5 seconds and stops its own child process tree at 10 GiB RSS, less
than 5 GiB free system memory, or 1,200 seconds. These are resource-probe limits,
not suitable default limits for full training. Peak summed RSS can double-count
shared pages; the Windows virtual-environment launcher has little memory of its
own, so parent-only RSS must not be used as the training footprint.

The timing includes the real tokenizer training code, SAFE conversion, and BPE.
It excludes sampling, imports, model saving, and the 100-training-row smoke check.
Probe timings are single measurements, not replicated speed benchmarks. Resource
extrapolations are scenarios, not confidence intervals or guaranteed completion
times. The 350-motif probes cannot establish an exact 3,000-motif full-run memory
requirement. No neural-network training or GPU work is included.

## 2. Full-corpus chemistry census

```bash
python examples/analyze_corpus.py --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json --output analysis_runs/census --workers 4
```

No tokenizer is loaded in this mode. Every row is parsed with RDKit. The output
contains atom-count and string-length distributions and overlapping specified
stereochemistry, formal-charge, isotope and disconnected-component subsets.
Explicit atom nodes are counted; implicit hydrogens are excluded. Parsing success
alone does not establish SAFE or DemoDiff reconstruction fidelity.

## 3. Full training on a sufficiently provisioned machine

The primary comparison is **BRICS + SAFE** against **NPE + SAFE**, both with the
same 3,000-token sequence budget. NPE additionally has a 3,000-motif budget and
300 initial rings, matching the motif configuration discussed in DemoDiff.
This gives NPE an additional learned dictionary; it is not an equal-total-model-
capacity comparison. Publish motif counts and model sizes alongside sequence
vocabulary sizes.

Review the resource report before executing these commands. The original memory
mode holds every graph and inverted index in RAM. The new SQLite mode stores them
on disk, limits graph batches to 256, and checkpoints initialization and both phases
of every merge. Global frequencies and frequency-tie ordering match the in-memory
algorithm. The complete 3,000-motif model fingerprint matched on the same 5,000
ChEMBL records; tests also exercise different batch sizes and mid-merge recovery.

SQLite checkpoints are trusted local artifacts, not portable model exchange files.
Use `model.save/load` for model exchange. On restart, use exactly the same input
and configuration with `--resume`; the input checksum and RDKit version are checked.
Sequence BPE and incomplete corpus audits restart rather than resuming mid-pass.

```bash
python -m molecular_tokenizer train --input ../chembl_37_smiles.smi --output trained_models/chembl37_full_brics_safe_v3000 --fragmentation brics --representation safe --bpe-scope fragment --vocab-size 3000 --num-workers 4 --skip-invalid

python -m molecular_tokenizer train --input ../chembl_37_smiles.smi --output trained_models/chembl37_full_npe_safe_v3000_m3000_r300 --fragmentation npe --representation safe --bpe-scope fragment --vocab-size 3000 --motif-vocab-size 3000 --ring-vocab-size 300 --num-workers 8 --skip-invalid --npe-storage sqlite --work-dir training_runs/npe_checkpoint
```

These commands pass every source row to training. Review each saved
`training_report.json` for the source checksum, accepted/rejected counts and
rejection reasons. Do not describe a probe model as a full-corpus model.

The automated full run trains NPE+SAFE, then BRICS+SAFE, audits both on every source
row, and renders the completed report. It monitors the process tree and stops at
12 GiB RSS, 4 GiB remaining system memory, or 20 GiB remaining disk space. It does
not disable sleep or change system settings. Keep the computer awake for execution.

```bash
python examples/run_full_experiment.py --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json --output training_runs/chembl37_full_sqlite_20261009 --workers 8
# If interrupted, rerun the same command with --resume.
```

Use `--publish` only when committing/pushing the generated `docs/` results to the
configured GitHub repository is intended. The running task was launched with this
option under the user's request to synchronize the research results. Final figures
are produced only after both full models and the audit succeed. This is a local
process, not a hosted job or a recurring scheduled task.

## 4. Full-corpus tokenizer audit

After both full models exist:

```bash
python examples/analyze_corpus.py --input ../chembl_37_smiles.smi --metadata ../dataset_metadata.json --model brics_safe=trained_models/chembl37_full_brics_safe_v3000 --model npe_safe=trained_models/chembl37_full_npe_safe_v3000_m3000_r300 --output analysis_runs/full_comparison --workers 4 --scope training_corpus
```

Outputs:

- `summary.json`: distributions, reconstruction rates, failure groups, usage
  frequencies, entropy, model fingerprints, and common-success comparisons.
- `molecules.csv.gz`: source-row metrics for every model, including failures.
- `vocabulary_*.csv`: vocabulary units and observed usage on exact reconstructions.

The audit explicitly decodes each encoding and compares canonical isomeric
SMILES, even though the internal strict flag is disabled to record diagnostic
information for mismatches. Lossy reconstructions are never counted as successful.
All source rows, including parse failures, remain in the overall success
denominator. Model-specific length statistics are conditional on exact recovery;
the common-success subset avoids comparing different surviving molecule sets.

Graph node counts and SAFE sequence token counts are different quantities. The
atom-to-unit ratio describes representation length, not storage compression.
Compact JSON payload size includes node/sequence IDs and all five integer fields
per graph connection. It excludes model files and metadata, and is not an optimal
binary codec or Transformer memory measurement. SAFE connection information is
inside its string even when the extra edge table is empty.

Because every ChEMBL record is used for training, this is **in-corpus descriptive
analysis**, not a held-out generalization result. A later external dataset or
separate train/test experiment is needed for that claim. Likewise, tokenizer
compression does not establish improved molecule generation or property prediction.

## 5. Figure plan and paper reference

Reference: Liu et al., [Graph Diffusion Transformers are In-Context Molecular
Designers](https://arxiv.org/html/2510.08744v1), Section 3.1 and Appendix B.2.
The requested [OpenReview PDF](https://openreview.net/pdf?id=lJ87GN5zJc) was
access-gated, so the accessible author preprint (v1) was used. Do not assert that
its experiments or manuscript version are identical to the current OpenReview PDF.

| Analysis | Paper reference | Our adaptation |
|---|---|---|
| Vocabulary size versus node count | Fig. 8 | Separate motif budget and SAFE sequence budget; include dictionary size |
| Node-count frequency and CCDF | Fig. 3(b), 9 | Show atoms, graph motifs and string tokens with distinct labels |
| Compression-ratio distribution | Fig. 10 | State denominator and report failures; avoid byte-compression claims |
| Molecule size versus representation length | Fig. 11 | Aggregate all molecules into density bins, avoiding cherry-picked examples |
| Fidelity and connection costs | Additional audit | Isotope/charge/stereo groups, failures, edges and compact payload size |
| Vocabulary units | Additional audit | Rank-frequency, coverage and representative units with usage counts |

The preflight page has measured resource plots and, when complete, a full-corpus
atom-count baseline. It intentionally does not display unrun motif/sequence
compression curves or claim a completed budget sweep. `build_research_page.py`
builds this preflight report; final comparison figures must be added from the
actual full-model audit outputs.

## 6. Build and publish the preflight page

```bash
python examples/build_research_page.py --resources training_runs/resource_probe/manifest.json --dataset ../dataset_metadata.json --census analysis_runs/census/summary.json --output docs
```

Use GitHub Pages with branch `main` and folder `/docs` (Settings → Pages → Deploy
from a branch). No Jekyll build is needed; `.nojekyll` is included. Only aggregate
data, original figures, and static HTML are published. Do not upload the full
SMILES file, virtual environment, or raw experiment directories.

GitHub Free supports Pages for public repositories; private repositories require
an eligible paid plan. A private source repository does not by itself make the
published Pages site private. See [GitHub Pages documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages).
Changing repository visibility or buying a plan is a separate user decision.
