# pTC-MRSF-CIS manuscript

`ptc_mrsf_cis.tex` — theory manuscript for projective transcorrelated MRSF-CIS.

## Status / how this was produced

The requested Overleaf project
(`https://git.overleaf.com/6a337be49b4f9baacf6d936c`) could **not** be reached
from the remote execution environment: the network egress policy does not
allowlist `git.overleaf.com` (HTTP 403, "Host not in allowlist"). The manuscript
is therefore committed here instead.

To sync it to Overleaf, either:

1. Add `git.overleaf.com` to the environment's network egress allowlist, then I
   can `git clone` the Overleaf project and push `ptc_mrsf_cis.tex` into it; or
2. Copy `ptc_mrsf_cis.tex` directly into the Overleaf project (upload, or paste
   into the main `.tex`).

## Compilation

Written in REVTeX 4-2 (`aps`/`prb`); a one-line fallback to the standard
`article` class is documented in the file header. It depends only on
`amsmath`, `amssymb`, `braket`, `graphicx`, `booktabs` — all standard on
Overleaf. No figures are required to compile.

## Please verify before submission

- **References.** Bibliographic details (volumes/pages) were drafted from memory
  and must be checked against the primary sources before any submission. The
  central pTC reference is S. Ten-no, *J. Chem. Phys.* **159**, 171103 (2023).
- **Author/affiliation** are placeholders.
- The numerical values quoted (CH2 CASCI match, <S^2>, the ~7e-15 eigensolver
  gate, 100% Hubbard correlation recovery) come from the validated prototypes in
  `../prototype/` and can be regenerated with `python3 *.py` there.
