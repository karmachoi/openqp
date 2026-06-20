!> @file tdhf_qmrsf_dk.F90
!> @brief QMRSF-DK dressed-kernel dynamic-correlation pathway (Pathway II) -- STUB.
!>
!> @details
!>   QMRSF = Quintet Mixed-Reference Spin-Flip.  This module implements the
!>   *density-functional-picture* dynamic-correlation layer that is applied on
!>   top of the determinant-union MRSF backbone
!>   (`tdhf_mrsf_ensemble_sigma` in source/modules/tdhf_mrsf_energy.F90).
!>
!>   Physics target (see DESIGN_QMRSF_DUAL_PATHWAYS.md, Pathway II):
!>     KS orbitals + a dressed / frequency-dependent quadratic exchange-
!>     correlation kernel g_xc(omega) added to the SIX closed-shell (0OS)
!>     double-spin-flip *diagonal* elements only.  The frequency dependence
!>     restores the double-excitation pole an adiabatic kernel cannot produce.
!>     Off-diagonals stay bare; spin-adaptation is untouched.  Correlation comes
!>     from the functional -- there is NO external-Q here (avoids double
!>     counting with v_xc/f_xc already present in the backbone diagonals).
!>     Cost ~ SCF + diagonalization.
!>
!>   STATUS: STUB.  This routine is a compile-time no-op.  It prints a banner and
!>   returns without modifying the 0OS diagonals, so wiring it in behind
!>   `qmrsf_pathway=dk` is safe and `qmrsf_pathway=none` stays bit-identical to
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
module tdhf_qmrsf_dk_mod

  implicit none

  private
  public :: tdhf_qmrsf_dk_C
  public :: tdhf_qmrsf_dk

  character(len=*), parameter :: module_name = "tdhf_qmrsf_dk_mod"

contains

  !> @brief C-bound entry point (matches the `void f(struct oqp_handle_t*)` ABI
  !>        declared in include/oqp.h and parsed by pyoqp via cffi).
  !> @details Mirrors tdhf_mrsf_energy_C: dereference the opaque C handle to the
  !>          Fortran `information` object, then call the inner driver.
  subroutine tdhf_qmrsf_dk_C(c_handle) bind(C, name="tdhf_qmrsf_dk")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_qmrsf_dk(inf)
  end subroutine tdhf_qmrsf_dk_C

  !> @brief Inner driver for the QMRSF dressed-kernel pathway (STUB / no-op).
  !>
  !> @param[inout] infos  OpenQP run container.  Consumed (read-only) inputs and
  !>                       produced outputs are exchanged through `infos%dat`
  !>                       (the tagarray container); see the tagarray contract
  !>                       block below.
  !>
  !> @details TAGARRAY CONTRACT (documentation; not yet enforced in the stub)
  !>
  !>   CONSUMES (must already be present, produced by the backbone solve):
  !>     OQP_td_bvec_mo   - converged P-space response eigenvectors (the 0OS
  !>                        double-spin-flip configurations live here).
  !>     OQP_td_energies  - converged backbone excitation energies; supply the
  !>                        frequency omega at which g_xc(omega) is evaluated.
  !>     OQP_VEC_MO_A / OQP_VEC_MO_B  - KS (alpha/beta) MO coefficients.
  !>     OQP_E_MO_A / OQP_E_MO_B      - KS MO energies (pole structure of g_xc).
  !>     OQP_DM_A / OQP_DM_B          - reference densities (kernel evaluation).
  !>     OQP_SM                       - AO overlap metric.
  !>
  !>   PRODUCES (new producer tags, to be reserved when the physics lands; the
  !>   stub does NOT reserve them so it remains a pure no-op):
  !>     OQP_qmrsf_dk_energies - dressed (DK-corrected) state energies.
  !>   See INTEGRATION_POINTS.md for the tagarray_driver.F90 symbol additions.
  subroutine tdhf_qmrsf_dk(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use precision, only: dp
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_qmrsf_dk"

    type(information), target, intent(inout) :: infos

    ! --- Open the main log file (append), matching the backbone discipline. ---
    open(unit=iw, file=infos%log_filename, position="append")

    call print_module_info('QMRSF_DK', &
         'Frequency-dependent dressed xc kernel on the 0OS diagonals (STUB)')

    write(iw,'(/,5x,a)') 'QMRSF-DK stub: dressed-kernel g_xc(omega) term not yet implemented.'
    write(iw,'(5x,a)')   'Backbone (0OS diagonals) results are passed through unchanged.'

    ! =====================================================================
    ! REAL STAGES TO FILL IN (post-backbone; modifies ONLY the six 0OS
    ! double-spin-flip diagonal elements -- off-diagonals and spin-adaptation
    ! stay bare):
    !
    !   STAGE 1 -- Identify the six closed-shell 0OS double-spin-flip diagonal
    !     configurations within the converged P space.
    !
    !   STAGE 2 -- Evaluate the dressed / frequency-dependent quadratic kernel
    !     g_xc(omega) for those diagonals on the KS orbitals.  Use the backbone
    !     excitation energy as the frequency omega (the term that restores the
    !     double-excitation pole; an adiabatic kernel gives exactly zero here).
    !     Generalize the Maitra one-single/one-double dressing to the coupled
    !     20-singlet block; document the chosen prescription and verify it does
    !     NOT split degenerate multiplets.
    !
    !   STAGE 3 -- Add g_xc(omega) to the 0OS diagonal ONLY.
    !     Assert no double counting with v_xc / f_xc already in the backbone
    !     diagonals (single correlation source).  Honor a 0OS-diagonal source
    !     selector (infos%tddft%qmrsf_0os_diag, backbone|seq|hve) and the kernel
    !     strength / frequency parameter (infos%tddft%qmrsf_dk_gamma) once added.
    !
    !   STAGE 4 -- Re-solve / re-evaluate the dressed P-space eigenproblem and
    !     write OQP_qmrsf_dk_energies.
    ! =====================================================================

    call flush(iw)
    close(iw)

  end subroutine tdhf_qmrsf_dk

end module tdhf_qmrsf_dk_mod
