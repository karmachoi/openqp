!> @file tdhf_qmrsf_icpt2.F90
!> @brief QMRSF-icPT2 dynamic-correlation pathway (Pathway I) -- STUB.
!>
!> @details
!>   QMRSF = Quintet Mixed-Reference Spin-Flip.  This module implements the
!>   *wavefunction-picture* dynamic-correlation layer that is applied AFTER the
!>   determinant-union MRSF backbone (`tdhf_mrsf_ensemble_sigma` in
!>   source/modules/tdhf_mrsf_energy.F90) has converged.
!>
!>   Physics target (see DESIGN_QMRSF_DUAL_PATHWAYS.md, Pathway I):
!>     bare spin-pure DSF/RAS-SF(2) backbone (the P space) dressed by an
!>     internally-contracted external-Q second-order self-energy downfold,
!>
!>         H_eff(E) = H_PP + Sigma(E),
!>         Sigma(E) = H_PQ (E - H_QQ)^{-1} H_QP,
!>
!>     with a Dyall/Fink zeroth-order H0, a Hermitized (des Cloizeaux /
!>     NEVPT2-style partitioned-denominator) effective Hamiltonian, then a small
!>     final diagonalization.  One-shot, O(N^5).
!>
!>   STATUS: STUB.  This routine is a compile-time no-op.  It prints a banner and
!>   returns without touching the backbone results, so that wiring it in behind
!>   `qmrsf_pathway=icpt2` is safe and `qmrsf_pathway=none` stays bit-identical to
!>   the current backbone.  No physics is implemented yet.
!>
!>   CONVENTIONS MIRRORED FROM:
!>     - source/modules/tdhf_mrsf_energy.F90  (C-binding pattern, `information`
!>       handle via c_interop, oqp_tagarray_driver usage, print_module_info,
!>       module_name parameter, log-file open/close discipline).
!>     - source/modules/tdhf_sf_energy.F90    (minimal module skeleton).
!>
!>   INTEGRATION: see tools/qmrsf_pathways_proto/INTEGRATION_POINTS.md for the
!>   exact (documentation-only) edits to register this module, add the
!>   `qmrsf_pathway` input flag, and dispatch to it after the backbone solve.
module tdhf_qmrsf_icpt2_mod

  implicit none

  private
  public :: tdhf_qmrsf_icpt2_C
  public :: tdhf_qmrsf_icpt2

  character(len=*), parameter :: module_name = "tdhf_qmrsf_icpt2_mod"

contains

  !> @brief C-bound entry point (matches the `void f(struct oqp_handle_t*)` ABI
  !>        declared in include/oqp.h and parsed by pyoqp via cffi).
  !> @details Mirrors tdhf_mrsf_energy_C: dereference the opaque C handle to the
  !>          Fortran `information` object, then call the inner driver.
  subroutine tdhf_qmrsf_icpt2_C(c_handle) bind(C, name="tdhf_qmrsf_icpt2")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_qmrsf_icpt2(inf)
  end subroutine tdhf_qmrsf_icpt2_C

  !> @brief Inner driver for the QMRSF-icPT2 downfold (STUB / no-op).
  !>
  !> @param[inout] infos  OpenQP run container.  Consumed (read-only) inputs and
  !>                       produced outputs are exchanged through `infos%dat`
  !>                       (the tagarray container); see the tagarray contract
  !>                       block below.
  !>
  !> @details TAGARRAY CONTRACT (documentation; not yet enforced in the stub)
  !>
  !>   CONSUMES (must already be present, produced by the backbone solve):
  !>     OQP_td_bvec_mo   - converged P-space response eigenvectors
  !>                        {C_I^(state)} (the trial/Ritz vectors in MO basis).
  !>     OQP_td_energies  - converged backbone (bare P-space) excitation energies;
  !>                        the reference energies E that enter Sigma(E).
  !>     OQP_VEC_MO_A / OQP_VEC_MO_B  - alpha/beta MO coefficients (define the
  !>                        active O1..O4 space and the external-Q C/O/V classes).
  !>     OQP_E_MO_A / OQP_E_MO_B      - MO energies (Dyall/Fink H0 denominators).
  !>     OQP_FOCK_A / OQP_FOCK_B, OQP_DM_A / OQP_DM_B, OQP_SM - reference Fock,
  !>                        density and overlap (one-body H0 pieces + metric).
  !>     (Active-space 1-/2-RDMs over O1..O4 are cheap -- 4 electrons -- and are
  !>      rebuilt here from the converged P vectors; no separate tag required.)
  !>
  !>   PRODUCES (new producer tags, to be reserved when the physics lands; the
  !>   stub does NOT reserve them so it remains a pure no-op):
  !>     OQP_qmrsf_icpt2_energies - dressed (icPT2-corrected) state energies.
  !>     OQP_qmrsf_icpt2_vec      - dressed P-space eigenvectors (optional).
  !>   See INTEGRATION_POINTS.md for the tagarray_driver.F90 symbol additions.
  subroutine tdhf_qmrsf_icpt2(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use precision, only: dp
    use printing, only: print_module_info
    use qmrsf_ao2mo_mod, only: qmrsf_active_integrals
    use qmrsf_cas_mod, only: qmrsf_cas_solve, QMRSF_NACT, QMRSF_NDET
    use qmrsf_icpt2_engine_mod, only: qmrsf_icpt2_dress_contracted, &
                                      qmrsf_icpt2_count_perturbers

    implicit none

    integer(8), parameter :: MAXQ = 4000000_8      ! contracted perturber-count guard

    type(information), target, intent(inout) :: infos

    integer :: nact, ncore, nbf, norb_w, i, p, q, r, s
    integer :: nQ
    integer(8) :: nQ_est
    integer, allocatable :: act_w(:)
    real(dp), allocatable :: h_win(:,:), eri_win(:,:,:,:), eps_win(:)
    real(dp) :: h4(QMRSF_NACT,QMRSF_NACT)
    real(dp) :: eri4(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp) :: eP(QMRSF_NDET), edr_en(QMRSF_NDET), edr_dy(QMRSF_NDET), evals(QMRSF_NDET)
    real(dp) :: ecore, herm
    logical  :: do_downfold

    ! --- Open the main log file (append), matching the backbone discipline. ---
    open(unit=iw, file=infos%log_filename, position="append")

    call print_module_info('QMRSF_icPT2', &
         'Internally-contracted external-Q PT2 self-energy downfold')

    ! ---- Active space from the quintet (S=2) reference -------------------
    !   nelec_A = ncore + 4 , nelec_B = ncore. The four SOMOs are MOs
    !   ncore+1..ncore+4; the icPT2 window is the frozen-core-dressed set of
    !   active+virtual orbitals MOs ncore+1..nbf (4 electrons, M_s=0).
    nbf   = int(infos%basis%nbf)
    nact  = QMRSF_NACT
    ncore = int(infos%mol_prop%nelec_B)
    if (int(infos%mol_prop%nelec_A) - int(infos%mol_prop%nelec_B) /= nact) then
      write(iw,'(/,5x,a)') 'QMRSF-icPT2 ERROR: reference is not a quintet (S=2) '// &
           'high-spin ROHF (need nelec_A - nelec_B = 4). Aborting pathway.'
      call flush(iw); close(iw); return
    end if
    norb_w = nbf - ncore
    allocate(act_w(norb_w))
    do i = 1, norb_w
      act_w(i) = ncore + i
    end do

    write(iw,'(/,5x,a,i0)') 'QMRSF-icPT2: basis functions          = ', nbf
    write(iw,'(5x,a,i0)')   'QMRSF-icPT2: frozen-core MOs          = ', ncore
    write(iw,'(5x,a,i0)')   'QMRSF-icPT2: active+virtual window    = ', norb_w
    write(iw,'(5x,a,i0,a)') 'QMRSF-icPT2: CAS active (SOMOs)       = ', nact, ' (MOs ncore+1..ncore+4)'

    ! ---- STAGE 0: window MO integrals via int2 reuse (frozen-core dressed) --
    allocate(h_win(norb_w,norb_w), eri_win(norb_w,norb_w,norb_w,norb_w))
    call qmrsf_active_integrals(infos, norb_w, act_w, ncore, h_win, eri_win, ecore)
    h4  = h_win(1:nact,1:nact)
    eri4 = eri_win(1:nact,1:nact,1:nact,1:nact)

    ! window ROHF MO energies (Dyall H0 denominators)
    allocate(eps_win(norb_w))
    block
      real(dp), contiguous, pointer :: emo(:)
      call tagarray_get_data(infos%dat, OQP_E_MO_A, emo)
      do i = 1, norb_w
        eps_win(i) = emo(act_w(i))
      end do
    end block

    ! ---- decide contracted-engine feasibility (perturber count ~ nvirt^4) ----
    nQ_est = qmrsf_icpt2_count_perturbers(norb_w, 2, 2, nact)
    do_downfold = (norb_w > nact) .and. (nQ_est <= MAXQ)

    if (do_downfold) then
      ! ---- backbone + CONTRACTED external-Q downfold (EN and Dyall) ----
      call qmrsf_icpt2_dress_contracted(h_win, eri_win, eps_win, norb_w, 2, 2, nact, &
                                        QMRSF_NDET, eP, edr_en, edr_dy, nQ, herm)
      write(iw,'(/,5x,a,f18.10)') 'QMRSF-icPT2: E_core (nuc + frozen core) = ', ecore
      write(iw,'(5x,a,i0)')       'QMRSF-icPT2: CAS P=36 ; external-Q perturbers nQ = ', nQ
      write(iw,'(5x,a,es10.2)')   'QMRSF-icPT2: H_eff Hermiticity          = ', herm
      write(iw,'(5x,a)')          'QMRSF-icPT2: state   E_CAS(total)     E_icPT2-EN(total)   E_icPT2-Dyall'
      do i = 1, min(8, QMRSF_NDET)
        write(iw,'(7x,i3,3f18.10)') i-1, eP(i)+ecore, edr_en(i)+ecore, edr_dy(i)+ecore
      end do
      evals = eP
    else
      ! ---- no virtuals (window==CAS) or perturber count over guard: CAS only ----
      call qmrsf_cas_solve(h4, eri4, evals, herm=herm)
      eP = evals; edr_en = evals; edr_dy = evals
      write(iw,'(/,5x,a,f18.10)') 'QMRSF-icPT2: E_core (nuc + frozen core) = ', ecore
      write(iw,'(5x,a,es10.2)')   'QMRSF-icPT2: CAS Hamiltonian |H-H^T|    = ', herm
      if (norb_w <= nact) then
        write(iw,'(5x,a)')        'QMRSF-icPT2: no virtual orbitals (window=CAS); CAS = FCI in-space.'
      else
        write(iw,'(5x,a,i0,a)')   'QMRSF-icPT2: external-Q downfold SKIPPED (perturber count ', &
             int(nQ_est), ' exceeds guard).'
      end if
      write(iw,'(5x,a)')          'QMRSF-icPT2: lowest CAS state totals:'
      do i = 1, min(8, QMRSF_NDET)
        write(iw,'(7x,a,i3,a,f18.10)') 'state ', i-1, '  E = ', evals(i)+ecore
      end do
    end if

    ! ---- Validation dumps ------------------------------------------------
    !  (a) qmrsf_cact_live.dat + qmrsf_icpt2_live.dat: the 4-active block gated
    !      against the closed-form AO oracle (route_a_oracle.py).
    !  (b) qmrsf_icpt2_full_live.dat: the full window integrals + spectra, gated
    !      against the NumPy det-CI+downfold oracle.
    block
      use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A
      real(dp), contiguous, pointer :: mo_a(:,:)
      integer :: mu
      call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
      open(unit=98, file='qmrsf_cact_live.dat', status='replace', action='write')
      write(98,'(i0,1x,i0)') nbf, nact
      do mu = 1, nbf
        write(98,'(*(es24.16))') (mo_a(mu, ncore+i), i=1,nact)
      end do
      close(98)
    end block
    open(unit=97, file='qmrsf_icpt2_live.dat', status='replace', action='write')
    write(97,'(i0)') nact
    do p = 1, nact; write(97,'(*(es24.16))') (h4(p,q), q=1,nact); end do
    do p = 1, nact; do q = 1, nact; do r = 1, nact
      write(97,'(*(es24.16))') (eri4(p,q,r,s), s=1,nact)
    end do; end do; end do
    write(97,'(es24.16)') ecore
    write(97,'(i0)') QMRSF_NDET
    write(97,'(*(es24.16))') (evals(i), i=1,QMRSF_NDET)
    close(97)

    open(unit=96, file='qmrsf_icpt2_full_live.dat', status='replace', action='write')
    write(96,'(i0,1x,i0)') norb_w, QMRSF_NDET
    do p = 1, norb_w; write(96,'(*(es24.16))') (h_win(p,q), q=1,norb_w); end do
    do p = 1, norb_w; do q = 1, norb_w; do r = 1, norb_w
      write(96,'(*(es24.16))') (eri_win(p,q,r,s), s=1,norb_w)
    end do; end do; end do
    write(96,'(es24.16)') ecore
    write(96,'(*(es24.16))') (eps_win(i), i=1,norb_w)
    write(96,'(*(es24.16))') (eP(i), i=1,QMRSF_NDET)
    write(96,'(*(es24.16))') (edr_en(i), i=1,QMRSF_NDET)
    write(96,'(*(es24.16))') (edr_dy(i), i=1,QMRSF_NDET)
    close(96)
    write(iw,'(/,5x,a)') 'QMRSF-icPT2: wrote validation dumps (qmrsf_icpt2_{live,full_live}.dat).'

    call flush(iw)
    close(iw)

  end subroutine tdhf_qmrsf_icpt2

end module tdhf_qmrsf_icpt2_mod
