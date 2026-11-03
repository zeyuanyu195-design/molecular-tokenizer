# Table-aligned counted-length histogram

The current interactive page shows BRICS+SAFE text-token counts against GraphBPE+SAFE fixed-motif ID counts for the same 44,328 matched molecules. It has no counting-unit selector or figure-download controls.

- BRICS+SAFE: 2,669,535 text tokens; mean 60.22231997834326.
- GraphBPE+SAFE: 372,924 motif IDs; mean 8.412831618841365.
- Text counts include connectivity syntax. Motif counts exclude connectivity and do not form a standalone reversible encoding.

The plot therefore shows a difference in counted lengths, not a complete-encoding compression ratio. Two overlaid histograms share 54 logarithmic bins across the full data range; both axes are logarithmic and dashed lines mark means. Density is bin count / (cohort size × bin width).

Regenerate from the saved aggregate frequency report without rerunning tokenization:

```sh
python experiments/tokenization/plot_counted_length_comparison.py --input-report docs/figures/token-distributions.json --output docs/figures/counted-length-comparison.png
```

The script reads `text.BRICS+SAFE` and `fragments.GraphBPE`, verifies sample counts, totals and means, and produces the figure. The previous same-unit BRICS/GraphBPE comparison script and aggregate data remain available for reproducibility.
