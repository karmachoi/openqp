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

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_qmrsf_icpt2"

    type(information), target, intent(inout) :: infos

    ! --- Open the main log file (append), matching the backbone discipline. ---
    open(unit=iw, file=infos%log_filename, position="append")

    call print_module_info('QMRSF_icPT2', &
         'Internally-contracted external-Q PT2 self-energy downfold (STUB)')

    write(iw,'(/,5x,a)') 'QMRSF-icPT2 stub: dynamic-correlation downfold not yet implemented.'
    write(iw,'(5x,a)')   'Backbone (P-space) results are passed through unchanged.'

    ! =====================================================================
    ! REAL STAGES TO FILL IN (all post-backbone; operate on the converged
    ! P-space vectors + active-space RDMs consumed from infos%dat):
    !
    !   STAGE 1 -- Build the internally-contracted external-Q perturbers.
    !     Form the first-order interacting space by contracting the C->O, O->V,
    !     C->V and core/virtual double-excitation classes against the P-space
    !     reference RDMs (1- and 2-RDM over the O1..O4 active space).
    !
    !   STAGE 2 -- Construct the Dyall/Fink zeroth-order Hamiltonian H0.
    !     Active two-body part retained; inactive/virtual reduced to one-body.
    !     Orthonormalize the perturbers against their overlap metric S_Q and
    !     drop near-singular directions (intruder-state guard).
    !
    !   STAGE 3 -- Assemble the self-energy / downfold Sigma(E).
    !     Sigma(E) = H_PQ (E - H_QQ)^{-1} H_QP on the orthonormalized perturbers
    !     (partitioned resolvent), giving the dressed P-space matrix
    !     H_eff = H_PP + Sigma(E).
    !
    !   STAGE 4 -- Hermitize.
    !     State-independent effective Hamiltonian via des Cloizeaux /
    !     NEVPT2-style symmetrization of H_eff.
    !
    !   STAGE 5 -- Diagonalize.
    !     Diagonalize the small Hermitian dressed P matrix -> corrected
    !     energies/vectors; write OQP_qmrsf_icpt2_energies (and _vec).
    !
    !   H0 SELECTION: honor infos%tddft%qmrsf_icpt2_h0 (dyall|fink) once added.
    ! =====================================================================

    call flush(iw)
    close(iw)

  end subroutine tdhf_qmrsf_icpt2

end module tdhf_qmrsf_icpt2_mod
