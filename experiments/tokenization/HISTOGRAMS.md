# BRICS / GraphBPE length distributions

The four-panel layout follows Appendix B.2, Figure 9 of [DemoDiff (arXiv:2510.08744v1)](https://arxiv.org/html/2510.08744v1). These figures use this project's saved experimental records, not the paper's data. The paper compares atom/motif representations; this comparison uses BRICS/GraphBPE.

Both comparisons use the same 44,328 molecules: GraphBPE boundary-exact records intersected with successful BRICS records and original-SMILES DemoDiff records. Record index, ChEMBL ID and original SMILES must agree. The GraphBPE counts must match the saved DemoDiff counts.

| Counted unit | BRICS mean | GraphBPE mean | BRICS median | GraphBPE median |
|---|---:|---:|---:|---:|
| Original fragmentation (one fragment = one counting token) | 6.102915 | 8.412832 | 5 | 5 |
| Final SAFE strings, same existing text tokenizer | 60.222320 | 68.847861 | 52 | 56 |

The fragment comparison excludes connectivity overhead and does not construct a BRICS ID vocabulary. It uses original `num_brics_fragments` / `num_bpe_fragments`, checked against saved groups. In 304 matched molecules, SAFE changes the BRICS boundaries; final BRICS SAFE fragment counts average 6.095809 instead. The text comparison uses final SAFE strings and includes connectivity syntax, with no special tokens, padding, truncation or dropout.

Panels: (a) common 54 log-spaced bins (55 edges), normalized by cohort size and bin width, with median lines; (b) empirical P(X >= x), including ties; (c) distinct length values ranked by their molecule frequencies, not motif-vocabulary frequencies; (d) unit-width density histogram with P95/P99 lines. All observations and the full range are retained. Zero bins cannot appear on log axes. Percentiles use NumPy's linear interpolation.

From the repository root, in the experiment environment described in [README](README.md), install Matplotlib if necessary and run:

```sh
python experiments/tokenization/plot_token_distributions.py --graph safe/data/output/records.jsonl --brics safe/data/output_brics/records.jsonl --demodiff runs/demodiff_original_smiles/smiles_50000_tokens.jsonl --text-tokenizer expts/tokenizer/tokenizer-custom.json --output-dir runs/token-distributions
```

Use a new output directory. The script reads saved records, does not retrain vocabularies or rerun GraphBPE, and writes PNG, vector PDF, SVG and a JSON frequency/statistics report with input hashes. The website publishes PNG/PDF figures and aggregated JSON; full molecular records are not added by this update.

- [Fragment-count PNG](../../docs/figures/token-distributions-fragments.png)
- [Fragment-count PDF](../../docs/figures/token-distributions-fragments.pdf)
- [SAFE text-count PNG](../../docs/figures/token-distributions-text.png)
- [SAFE text-count PDF](../../docs/figures/token-distributions-text.pdf)
- [Statistics and exact length frequencies](../../docs/figures/token-distributions.json)

Matplotlib version used: 3.11.2. Recomputed SAFE text totals (2,669,535 and 3,051,888) and GraphBPE motif total (372,924) match the published matched-cohort table.
