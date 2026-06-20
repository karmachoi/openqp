"""QMRSF-icPT2 results post-processing.

The live Fortran routine ``tdhf_qmrsf_icpt2`` writes a fixed-name validation
dump ``qmrsf_icpt2_full_live.dat`` into the run's working directory.  This module
turns that raw dump into clean, consumable results: a parser, a results builder,
a JSON writer, and a formatted log table.

Dump file format (free-format ``es24.16`` text, written by
``source/modules/tdhf_qmrsf_icpt2.F90``)::

    line 1            : "<norb_w> <nPd>"        two ints (nPd = number of dressed roots, 36)
    next norb_w lines : h_win rows               (norb_w values each)
    next norb_w**3    : eri_win blocks (p,q,r,:) (norb_w values each)
    1 line            : ecore                    (nuc. repulsion + frozen-core const, Hartree)
    1 line            : eps_win                  (norb_w ROHF MO energies)
    1 line            : eP                        (nPd bare CAS-root ELECTRONIC energies)
    1 line            : edr_en                    (nPd icPT2 Epstein-Nesbet ELECTRONIC energies)
    1 line            : edr_dy                    (nPd icPT2 Dyall ELECTRONIC energies)

All ``e*`` arrays are electronic energies; add ``ecore`` for totals.
"""

import json

# CODATA Hartree -> eV (matches the value used elsewhere in the QMRSF stack).
HARTREE_TO_EV = 27.211386245988


def _read_floats(line):
    """Parse a whitespace-separated line of Fortran ``es24.16`` floats."""
    return [float(tok) for tok in line.split()]


def parse_qmrsf_icpt2_dump(path):
    """Read the scalar results from a ``qmrsf_icpt2_full_live.dat`` dump.

    Only the scalar spectra are returned (norb_w, nPd, ecore, eps_win, eP,
    edr_en, edr_dy).  The large ``h_win`` / ``eri_win`` integral blocks are
    skipped, but the exact number of their lines is consumed so the file
    pointer reaches the trailing scalar records.

    Parameters
    ----------
    path : str
        Filesystem path to the dump file.

    Returns
    -------
    dict
        ``{"norb_w", "nPd", "ecore", "eps_win", "eP", "edr_en", "edr_dy"}``.

    Raises
    ------
    ValueError
        If the file is truncated or the header is malformed.
    """
    with open(path, 'r') as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError('QMRSF-icPT2 dump is empty: %s' % path)

    header = lines[0].split()
    if len(header) < 2:
        raise ValueError('QMRSF-icPT2 dump header malformed (expected "<norb_w> <nPd>"): %r' % lines[0])
    norb_w = int(header[0])
    nPd = int(header[1])
    if norb_w <= 0 or nPd <= 0:
        raise ValueError('QMRSF-icPT2 dump header has non-positive sizes: norb_w=%d nPd=%d' % (norb_w, nPd))

    # Skip the integral blocks: norb_w rows of h_win + norb_w**3 rows of eri_win.
    skip = norb_w + norb_w ** 3
    idx = 1 + skip  # first trailing scalar record (ecore)

    # 5 trailing scalar records: ecore, eps_win, eP, edr_en, edr_dy.
    needed = idx + 5
    if len(lines) < needed:
        raise ValueError(
            'QMRSF-icPT2 dump is truncated: have %d lines, need at least %d '
            '(norb_w=%d, nPd=%d)' % (len(lines), needed, norb_w, nPd)
        )

    ecore = _read_floats(lines[idx])[0]
    eps_win = _read_floats(lines[idx + 1])
    eP = _read_floats(lines[idx + 2])
    edr_en = _read_floats(lines[idx + 3])
    edr_dy = _read_floats(lines[idx + 4])

    for name, arr, n in (
        ('eps_win', eps_win, norb_w),
        ('eP', eP, nPd),
        ('edr_en', edr_en, nPd),
        ('edr_dy', edr_dy, nPd),
    ):
        if len(arr) != n:
            raise ValueError(
                'QMRSF-icPT2 dump record %s has %d values, expected %d' % (name, len(arr), n)
            )

    return {
        'norb_w': norb_w,
        'nPd': nPd,
        'ecore': ecore,
        'eps_win': eps_win,
        'eP': eP,
        'edr_en': edr_en,
        'edr_dy': edr_dy,
    }


def build_qmrsf_results(dump, ref_energy):
    """Assemble a clean results dict from a parsed dump and the reference energy.

    Total energies are ``E_electronic + ecore``.  Excitation energies are
    measured relative to each method's own ground state (state 0) and converted
    Hartree -> eV.

    Parameters
    ----------
    dump : dict
        Output of :func:`parse_qmrsf_icpt2_dump`.
    ref_energy : float
        The converged quintet ROHF reference energy (Hartree).

    Returns
    -------
    dict
        Results dict with per-state totals (Hartree) and excitation energies (eV).
    """
    ecore = dump['ecore']
    nPd = dump['nPd']
    eP = dump['eP']
    edr_en = dump['edr_en']
    edr_dy = dump['edr_dy']

    # Each method's ground-state total (state 0) anchors its excitation energies.
    cas0 = eP[0] + ecore
    en0 = edr_en[0] + ecore
    dy0 = edr_dy[0] + ecore

    states = []
    for i in range(nPd):
        e_cas = eP[i] + ecore
        e_en = edr_en[i] + ecore
        e_dy = edr_dy[i] + ecore
        states.append({
            'index': i,
            'E_CAS': e_cas,
            'E_icPT2_EN': e_en,
            'E_icPT2_Dyall': e_dy,
            'exc_CAS_eV': (e_cas - cas0) * HARTREE_TO_EV,
            'exc_EN_eV': (e_en - en0) * HARTREE_TO_EV,
            'exc_Dyall_eV': (e_dy - dy0) * HARTREE_TO_EV,
        })

    return {
        'method': 'QMRSF-icPT2',
        'reference_energy': float(ref_energy),
        'ecore': ecore,
        'n_window_orbitals': dump['norb_w'],
        'n_states': nPd,
        'states': states,
    }


def write_qmrsf_json(results, json_path):
    """Write the results dict to ``json_path`` as pretty-printed JSON."""
    with open(json_path, 'w') as fh:
        json.dump(results, fh, indent=2)
    return json_path


def format_qmrsf_log_table(results, max_states=10):
    """Render an aligned text table of the lowest states for the log.

    Columns: state index, the three total energies (Hartree), and the
    Epstein-Nesbet icPT2 excitation energy (eV, relative to state 0).
    """
    states = results.get('states', [])
    n_show = min(max_states, len(states))

    header_lines = [
        'QMRSF-icPT2 results',
        'reference (quintet ROHF) = %18.10f Hartree' % results['reference_energy'],
        'E_core (nuc + frozen)    = %18.10f Hartree' % results['ecore'],
        'window orbitals = %d   dressed roots = %d   (showing lowest %d)'
        % (results['n_window_orbitals'], results['n_states'], n_show),
        '',
    ]

    col = '%5s %18s %18s %18s %14s'
    rule = '-' * 76
    table = [
        col % ('state', 'E_CAS(Eh)', 'E_icPT2-EN(Eh)', 'E_icPT2-Dyall(Eh)', 'exc-EN(eV)'),
        rule,
    ]
    rowfmt = '%5d %18.10f %18.10f %18.10f %14.4f'
    for st in states[:n_show]:
        table.append(rowfmt % (
            st['index'], st['E_CAS'], st['E_icPT2_EN'], st['E_icPT2_Dyall'], st['exc_EN_eV'],
        ))

    return '\n'.join(header_lines + table)
