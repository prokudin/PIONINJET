# Points to check with Alexei

The ROOT-free tutorial uses the uploaded native kernels and grids at
`3945b91086dd153b4d940531e2f09f8f70c0496f`. The Python integration agrees
with a separate four-dimensional integration of the original C++ NLO
routine within its integration errors. Numerical refinement changes the
asymmetries by less than 0.2%, while the Figure 4 central curves remain
different. No fit parameters or normalizations have been adjusted.

For `sqrt(s)=500 GeV`, `pt=31 GeV`, `0<eta<1`, `jt=0.1 GeV`:

| z | Pion | Tutorial A | Figure 4 central A (vector extraction) |
|---:|:---:|---:|---:|
| 0.13 | + | 0.00135923 | 0.00098552 |
| 0.13 | − | -0.00245117 | -0.00171671 |
| 0.37 | + | 0.00594409 | 0.00429761 |
| 0.37 | − | -0.00856644 | -0.00589064 |

Useful things to confirm:

1. Which steering file, PDF/FF order, and central parameter files generated
   Figures 4 and 5? Are the two supplied KANG2015 files the exact ones used?
2. Is the root `PINLO.GRID` the intended evolved FF input? It differs from
   the one in `FORTRAN/frag/`, and contains three repeated blocks.
3. Was the finite delta matching term `-2 CF alpha_s/pi` included in the
   figures? The code includes it; the printed Eqs. (8),(10) do not.
4. Should `Spert` in Eq. (14) use `dmu^2/mu^2` before the factor `1/2` in
   Eq. (6)? The uploaded antiderivative follows that convention.
5. Was outgoing gluon fragmentation included in the Gaussian Figure 5
   curve? Including it in the legacy driver gives agreement within 0.52%;
   excluding it leaves about 37.5% of the plotted result.
6. Does Figure 5 use the driver's `2pt*2jt` measure, despite the Cartesian
   differential label? Both quantities are stored in the tutorial outputs.
7. Are numerical tables for the central curves, and the missing
   `transmaxlo_new.grid` / `transmaxnlo_new.grid`, available?

These are questions about reconstructing the historical configuration,
not a conclusion that the paper's calculation is incorrect. The
[audit](AUDIT.md), [Figure 4 comparison](example_results/comparison/figure4_comparison.png),
and [Figure 5 comparison](example_results/figure5/figure5_comparison.png)
contain the evidence and exact settings.
