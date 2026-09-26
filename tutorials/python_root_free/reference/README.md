# Published-curve reference points

The CSVs contain coordinates extracted from the vector graphics of
[arXiv:1707.00913v2](https://arxiv.org/pdf/1707.00913v2), not numerical
tables supplied by the authors. Each curve has 15 points. Experimental
points and uncertainty bands are not digitized here.

The PDF used has SHA-256:

`003e9ce2f68787bea780b8f74e041282d2f3d8e74040a53d66403fde78c0431b`

The extraction used Poppler 26.07.0 and visually checked axis ticks.
Figure 4 is on physical PDF page 9; Figure 5 is on page 10. Figure 4 uses
a logarithmic horizontal axis and linear vertical axes; Figure 5 uses a
linear horizontal axis and logarithmic vertical axis. `extract.py` records
the exact path indices and tick transformations. Tiny zero-error marker
segments supply central values for curves without a convenient continuous
path. The transformations already account for the inverted SVG y axis.

To regenerate, place that PDF in your own workspace and, from the tutorial
folder, run:

```sh
mkdir -p results/reference_check
pdftocairo -svg -f 9 -l 9 /path/to/1707.00913v2.pdf results/reference_check/page9.svg
pdftocairo -svg -f 10 -l 10 /path/to/1707.00913v2.pdf results/reference_check/page10.svg
python3 reference/extract.py --page9 results/reference_check/page9.svg \
  --page10 results/reference_check/page10.svg --output results/reference_check
```

The extractor is specific to that PDF and SVG layout; a different Poppler
version can change path ordering. Its assertions detect some such changes,
but visual inspection is still necessary. It is not a general PDF digitizer.
The PDF and full-page renders are not duplicated in this tutorial. Small
vector rounding and any fluctuations already present in the published
curves limit the precision of comparisons.
