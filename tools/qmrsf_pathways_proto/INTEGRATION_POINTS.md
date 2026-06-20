# QMRSF dual-pathway integration points

Copy-pasteable wire-in for the two stub modules created on `feat/qmrsf-dual-pathways`:

- `source/modules/tdhf_qmrsf_icpt2.F90`  (Pathway I, C entry `tdhf_qmrsf_icpt2`)
- `source/modules/tdhf_qmrsf_dk.F90`     (Pathway II, C entry `tdhf_qmrsf_dk`)

**These are DOCUMENTATION / diffs only.** Per the branch coordination rule
(DESIGN_QMRSF_DUAL_PATHWAYS.md §0), the existing files are **not** edited yet, so
this branch stays free of merge conflicts with the concurrent backbone /
spin-projection work and the build cannot break. Apply the edits below once the
backbone's S^2 projection lands. Every change is **additive**: with
`qmrsf_pathway=none` (the default) the code path is bit-identical to today.

The two stub modules themselves are convention-matched but **unbuilt** (a full
OpenQP build is heavy and was not run); they mirror `tdhf_mrsf_energy.F90` /
`tdhf_sf_energy.F90` closely enough to compile as drop-ins.

---

## (a) CMake registration -- NO EDIT REQUIRED

`source/modules/CMakeLists.txt` already globs every module source:

```cmake
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB SOURCES_modules CONFIGURE_DEPENDS "*.F90" "*.f90")

list(APPEND SOURCES ${SOURCES_modules})
set(SOURCES ${SOURCES} PARENT_SCOPE)
```

Because the glob uses `CONFIGURE_DEPENDS` and matches `*.F90`, the two new files
`tdhf_qmrsf_icpt2.F90` and `tdhf_qmrsf_dk.F90` are picked up automatically the
next time CMake configures. **No CMakeLists edit is needed.**

If you prefer an explicit (non-glob) listing instead, replace the `file(GLOB ...)`
line with an explicit set that includes the two new files, e.g.:

```cmake
# (optional explicit form -- only if the glob is ever removed)
list(APPEND SOURCES_modules
  "${dir}/tdhf_qmrsf_icpt2.F90"
  "${dir}/tdhf_qmrsf_dk.F90")
```

---

## (b) New input flag `qmrsf_pathway` (additive; default `none` = no-op)

Values: `none | icpt2 | dk`  (default `none`). Companion knobs from the design
note -- `qmrsf_0os_diag`, `qmrsf_icpt2_h0`, `qmrsf_dk_gamma` -- follow the same
pattern and are listed at the end of this section.

The flag is carried as a small integer code in the C-interop `tddft_parameters`
struct, parsed from the `[tdhf]` input section, and consumed by the Fortran
dispatch in (c). Four files change.

### b.1 -- `source/types.F90` : add a field to `tddft_parameters`

The derived type is `bind(c)`, so the new field must be a C-interop kind.
Add it as the **last** field (append-only keeps the C struct layout backward
compatible). Surrounding lines (currently end of the type):

```fortran
    integer(c_int64_t) :: z_solver = 0     !< z-vector solver: 0 (CG), 1 (GMRES legacy), 2 (MINRES), 3 (AUTO)
    integer(c_int64_t) :: gmres_dim = 50   !< The Restart dimension of GMRES
    logical(c_bool) :: umrsf = .false.     !< UMRSF branch calculations switch in td_mrsf_energy module
  end type tddft_parameters
```

Change to (append the new fields before `end type`):

```fortran
    integer(c_int64_t) :: z_solver = 0     !< z-vector solver: 0 (CG), 1 (GMRES legacy), 2 (MINRES), 3 (AUTO)
    integer(c_int64_t) :: gmres_dim = 50   !< The Restart dimension of GMRES
    logical(c_bool) :: umrsf = .false.     !< UMRSF branch calculations switch in td_mrsf_energy module
    integer(c_int64_t) :: qmrsf_pathway = 0  !< QMRSF dynamic-correlation pathway: 0=none, 1=icpt2, 2=dk
    integer(c_int64_t) :: qmrsf_0os_diag = 0 !< 0OS diagonal source: 0=backbone, 1=seq, 2=hve
    integer(c_int64_t) :: qmrsf_icpt2_h0 = 0 !< icPT2 zeroth-order H: 0=dyall, 1=fink
    real(c_double)    :: qmrsf_dk_gamma = 0.0_dp !< DK dressed-kernel strength / frequency parameter
  end type tddft_parameters
```

### b.2 -- `include/oqp.h` : declare the two C entry points

So pyoqp's cffi parser (`pyoqp/oqp/__init__.py` reads `include/oqp.h` via
`ffi.cdef`) exposes them as `oqp.tdhf_qmrsf_icpt2` / `oqp.tdhf_qmrsf_dk`.
(No struct change needed in the header: `tddft_parameters` is mirrored on the
Python side as the cffi struct, and append-only new fields are read by the
existing layout-following code; see b.3.) Surrounding lines:

```c
void tdhf_mrsf_energy(struct oqp_handle_t *inf);
void tdhf_umrsf_energy(struct oqp_handle_t *inf);
void tdhf_mrsf_ekt_ip(struct oqp_handle_t *inf);
```

Add the two prototypes (anywhere in that block; here right after the MRSF ones):

```c
void tdhf_mrsf_energy(struct oqp_handle_t *inf);
void tdhf_umrsf_energy(struct oqp_handle_t *inf);
void tdhf_qmrsf_icpt2(struct oqp_handle_t *inf);
void tdhf_qmrsf_dk(struct oqp_handle_t *inf);
void tdhf_mrsf_ekt_ip(struct oqp_handle_t *inf);
```

Because `pyoqp/oqp/__init__.py` binds **every** callable symbol in the library
into the `oqp` namespace automatically, no further Python binding code is needed
to make `oqp.tdhf_qmrsf_icpt2(mol)` / `oqp.tdhf_qmrsf_dk(mol)` callable.

### b.3 -- `pyoqp/oqp/molecule/oqpdata.py` : declare + parse + set the keyword

Three edits in this file.

(i) Keyword declaration -- in the `'tdhf'` config block (currently ends):

```python
        'z_solver': {'type': int, 'default': '0'},  # 0: CG, 1: GMRES (legacy), 2: MINRES, 3: AUTO
        'gmres_dim': {'type': int, 'default': '50'},  # Dimension for GMRES during Z-vector
    },
```

Change to (append the new keywords; string-valued so users type `none|icpt2|dk`):

```python
        'z_solver': {'type': int, 'default': '0'},  # 0: CG, 1: GMRES (legacy), 2: MINRES, 3: AUTO
        'gmres_dim': {'type': int, 'default': '50'},  # Dimension for GMRES during Z-vector
        'qmrsf_pathway': {'type': string, 'default': 'none'},   # none | icpt2 | dk
        'qmrsf_0os_diag': {'type': string, 'default': 'backbone'},  # backbone | seq | hve
        'qmrsf_icpt2_h0': {'type': string, 'default': 'dyall'}, # dyall | fink
        'qmrsf_dk_gamma': {'type': float, 'default': '0.0'},
    },
```

(ii) Setter dispatch -- in the `"tdhf"` entry of the setter map (currently ends):

```python
            "z_solver": "set_tdhf_z_solver",
            "gmres_dim": "set_tdhf_gmres_dim",
        },
    }
```

Change to:

```python
            "z_solver": "set_tdhf_z_solver",
            "gmres_dim": "set_tdhf_gmres_dim",
            "qmrsf_pathway": "set_tdhf_qmrsf_pathway",
            "qmrsf_0os_diag": "set_tdhf_qmrsf_0os_diag",
            "qmrsf_icpt2_h0": "set_tdhf_qmrsf_icpt2_h0",
            "qmrsf_dk_gamma": "set_tdhf_qmrsf_dk_gamma",
        },
    }
```

(iii) Setter methods -- add next to the existing `set_tdhf_*` methods (e.g. right
after `set_tdhf_gmres_dim`). These map the user string to the integer code that
matches `source/types.F90` (b.1):

```python
    def set_tdhf_qmrsf_pathway(self, qmrsf_pathway):
        """Select QMRSF dynamic-correlation pathway: none | icpt2 | dk (default none = bare backbone)."""
        codes = {'none': 0, 'icpt2': 1, 'dk': 2}
        key = str(qmrsf_pathway).strip().lower()
        if key not in codes:
            raise ValueError(f"qmrsf_pathway must be none, icpt2, or dk; got {qmrsf_pathway}")
        self._data.tddft.qmrsf_pathway = codes[key]

    def set_tdhf_qmrsf_0os_diag(self, qmrsf_0os_diag):
        """0OS diagonal source for QMRSF: backbone | seq | hve (default backbone)."""
        codes = {'backbone': 0, 'seq': 1, 'hve': 2}
        key = str(qmrsf_0os_diag).strip().lower()
        if key not in codes:
            raise ValueError(f"qmrsf_0os_diag must be backbone, seq, or hve; got {qmrsf_0os_diag}")
        self._data.tddft.qmrsf_0os_diag = codes[key]

    def set_tdhf_qmrsf_icpt2_h0(self, qmrsf_icpt2_h0):
        """Zeroth-order H for the icPT2 external-Q downfold: dyall | fink (default dyall)."""
        codes = {'dyall': 0, 'fink': 1}
        key = str(qmrsf_icpt2_h0).strip().lower()
        if key not in codes:
            raise ValueError(f"qmrsf_icpt2_h0 must be dyall or fink; got {qmrsf_icpt2_h0}")
        self._data.tddft.qmrsf_icpt2_h0 = codes[key]

    def set_tdhf_qmrsf_dk_gamma(self, qmrsf_dk_gamma):
        """DK dressed-kernel strength / frequency parameter."""
        self._data.tddft.qmrsf_dk_gamma = qmrsf_dk_gamma
```

### b.4 -- `pyoqp/oqp/utils/input_checker.py` (optional but recommended)

Add `qmrsf_pathway` to the validated `[tdhf]` keyword set / value checks here so a
typo is rejected early (mirror the existing `td_type` value guard). Allowed
values `none|icpt2|dk`; this is a no-op for correctness but improves UX.

---

## (c) Dispatch after the backbone solve

Two equivalent options; **option 1 (Python dispatch in `single_point.py`) is
recommended** because it keeps the backbone Fortran (the file the concurrent
chat is editing) untouched.

### Option 1 (recommended) -- `pyoqp/oqp/library/single_point.py`

(i) Register the two C entries in the `energy_func` map (currently):

```python
        self.energy_func = {
            'hf': oqp.hf_energy,
            'rpa': oqp.tdhf_energy,
            'tda': oqp.tdhf_energy,
            'sf': oqp.tdhf_sf_energy,
            'mrsf': oqp.tdhf_mrsf_energy,
            'umrsf': oqp.tdhf_umrsf_energy,
            'mrsf_ekt_ip': oqp.tdhf_mrsf_ekt_ip,
            'mrsf_ekt_ea': oqp.tdhf_mrsf_ekt_ea,
        }
```

Change to:

```python
        self.energy_func = {
            'hf': oqp.hf_energy,
            'rpa': oqp.tdhf_energy,
            'tda': oqp.tdhf_energy,
            'sf': oqp.tdhf_sf_energy,
            'mrsf': oqp.tdhf_mrsf_energy,
            'umrsf': oqp.tdhf_umrsf_energy,
            'mrsf_ekt_ip': oqp.tdhf_mrsf_ekt_ip,
            'mrsf_ekt_ea': oqp.tdhf_mrsf_ekt_ea,
            'qmrsf_icpt2': oqp.tdhf_qmrsf_icpt2,
            'qmrsf_dk': oqp.tdhf_qmrsf_dk,
        }
```

(ii) Dispatch right AFTER the backbone solve. The backbone runs at the end of the
`tdhf_energy` method (currently):

```python
        # do TDDFT
        dump_log(self.mol, title='PyOQP: TDDFT steps', section='tdhf')
        self.energy_func[self.td](self.mol)
```

Change to (the pathway layers on top of the converged MRSF backbone; `none`
leaves behavior identical):

```python
        # do TDDFT
        dump_log(self.mol, title='PyOQP: TDDFT steps', section='tdhf')
        self.energy_func[self.td](self.mol)

        # QMRSF dynamic-correlation pathway (additive; runs only on an MRSF backbone)
        qmrsf_pathway = str(self.mol.config['tdhf'].get('qmrsf_pathway', 'none')).strip().lower()
        if qmrsf_pathway != 'none':
            if self.td not in ('mrsf', 'umrsf'):
                raise ValueError('qmrsf_pathway requires an MRSF backbone ([tdhf] type=mrsf)')
            dump_log(self.mol, title='PyOQP: QMRSF pathway (%s)' % qmrsf_pathway, section='tdhf')
            self.energy_func['qmrsf_%s' % qmrsf_pathway](self.mol)
```

This consumes the backbone's `OQP::td_bvec_mo` / `OQP::td_energies` (already in
`infos%dat`) exactly as documented in each module's tagarray contract.

### Option 2 (alternative) -- inside `tdhf_mrsf_energy` (Fortran)

If a single-call Fortran dispatch is preferred over a second Python call, add the
following just before the final `call int2_driver%clean()` in
`source/modules/tdhf_mrsf_energy.F90` (this DOES edit the backbone file, so defer
it until the concurrent work merges):

```fortran
    ! --- QMRSF dynamic-correlation pathway (additive; none = no-op) ---
    select case (infos%tddft%qmrsf_pathway)
    case (1)  ! icpt2
      call tdhf_qmrsf_icpt2(infos)
    case (2)  ! dk
      call tdhf_qmrsf_dk(infos)
    case default
      ! 0 = none: bare backbone, unchanged
    end select
```

with, at the top of the module's `use` list:

```fortran
    use tdhf_qmrsf_icpt2_mod, only: tdhf_qmrsf_icpt2
    use tdhf_qmrsf_dk_mod, only: tdhf_qmrsf_dk
```

---

## (d) New producer tagarray symbols (for when the physics lands)

The stubs are no-ops and reserve nothing. When the real downfold/kernel is
implemented, add these tag constants to `source/tagarray_driver.F90` next to the
existing `OQP_td_*` definitions, and add them to the module's public list:

```fortran
  character(len=*), parameter, public :: OQP_qmrsf_icpt2_energies = OQP_prefix // "qmrsf_icpt2_energies"
  character(len=*), parameter, public :: OQP_qmrsf_icpt2_vec      = OQP_prefix // "qmrsf_icpt2_vec"
  character(len=*), parameter, public :: OQP_qmrsf_dk_energies    = OQP_prefix // "qmrsf_dk_energies"
  character(len=*), parameter, public :: OQP_qmrsf_icpt2_energies_comment = "QMRSF-icPT2 dressed state energies"
  character(len=*), parameter, public :: OQP_qmrsf_dk_energies_comment    = "QMRSF-DK dressed state energies"
```

reserve them in the inner driver via `infos%dat%reserve_data(...)` (mirroring the
`OQP_td_energies` reservation in `tdhf_mrsf_energy`), and write the corrected
energies into them.

---

## Summary checklist

| Step | File | Edit |
|------|------|------|
| (a)  | `source/modules/CMakeLists.txt` | none (glob auto-discovers) |
| b.1  | `source/types.F90` | append 4 fields to `tddft_parameters` |
| b.2  | `include/oqp.h` | declare `tdhf_qmrsf_icpt2`, `tdhf_qmrsf_dk` |
| b.3  | `pyoqp/oqp/molecule/oqpdata.py` | keyword decl + setter map + 4 setter methods |
| b.4  | `pyoqp/oqp/utils/input_checker.py` | optional value validation |
| (c)  | `pyoqp/oqp/library/single_point.py` | `energy_func` entries + post-backbone dispatch |
| (d)  | `source/tagarray_driver.F90` | producer tags (deferred until physics lands) |

With every step applied and `[tdhf] qmrsf_pathway=none` (default), behavior is
bit-identical to the current backbone.
