# Central PDF inputs

This directory bundles the exact central PDF data used in the tutorial:

| Set | Member | Data version | Reference |
|---|---:|---:|---|
| CT10nlo | 0 | 4 | [arXiv:1007.2241](https://arxiv.org/abs/1007.2241) |
| cteq6l1 | 0 | 4 | [hep-ph/0201195](https://arxiv.org/abs/hep-ph/0201195) |

The original `.info` metadata and central `.dat` files are unchanged.
The `.info` for CT10nlo describes the full 53-member set, but only member
zero is included here. Error members and PDF uncertainty evaluation are
not part of this central-fit tutorial.

These are CTEQ PDF datasets distributed in the LHAPDF6 format, not new fits
made by this project. Their authors and references remain in the original
metadata. Cite both the applicable PDF paper and LHAPDF in scientific use.
The canonical distribution is the
[LHAPDF PDF-set archive](https://lhapdfsets.web.cern.ch/lhapdfsets/current/).
`SOURCES.json` records the dataset archive URLs and the checksums of the
specific files used here; hashes identify this snapshot even if the remote
archives change.

`Backend` checks these hashes and adds this directory to the native library's
search path before initialization. It always requests member zero. The
shared LHAPDF library and its C++ headers remain a build dependency.
