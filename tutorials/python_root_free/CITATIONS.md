# Attribution and scientific references

The scientific kernels, fit parameters, data grids and custom HOPPET code
remain the work of their original authors. This tutorial adds a Python
interface, a small C++ bridge, build preparation, numerical checks and
documentation. It does not relabel the original scientific work or assign
a new license to the upstream files. Preserve the notices in each retained
source. HOPPET's supplied source identifies its GNU public license and
requests citation of its paper.

For scientific use, cite the original pion-in-jet calculation, the fit and
all PDF/FF/evolution ingredients actually used:

- Z.-B. Kang, A. Prokudin, F. Ringer and F. Yuan,
  [Collins azimuthal asymmetries of hadron production inside jets](https://arxiv.org/abs/1707.00913).
- Z.-B. Kang, A. Prokudin, P. Sun and F. Yuan,
  [Extraction of Quark Transversity Distribution and Collins Fragmentation Functions with QCD Evolution](https://arxiv.org/abs/1505.05589).
- G. P. Salam and J. Rojo,
  [A Higher Order Perturbative Parton Evolution Toolkit (HOPPET)](https://arxiv.org/abs/0804.3755).
- A. Buckley et al.,
  [LHAPDF6: parton density access in the LHC precision era](https://arxiv.org/abs/1412.7420).

Also use the CT10/CTEQ, DSS fragmentation and DSSV references identified
by the corresponding supplied data and source files. The ambiguity between
the two `PINLO.GRID` files is documented in [AUDIT.md](AUDIT.md); confirm
the historical FF-table identity before assigning a dataset-specific claim.

When sharing results, identify the source revision, tutorial revision,
input hashes and explicit settings. Describe the present output as an
audited implementation of the uploaded source with an unresolved figure
discrepancy, not an exact reproduction of the published central curves.
