module hf_cphf_response_mod

  use precision, only: dp

  implicit none

  character(len=*), parameter :: module_name = "hf_cphf_response_mod"

  private
  public :: hf_cphf_hessian_response

contains

!###############################################################################

  !> @brief Guarded scaffold for HF CPHF orbital-response Hessian terms.
  !> @details hf_cphf_hessian_response_scaffold partial_kernel:
  !>  A complete HF analytic Hessian needs the occupied-virtual orbital-response
  !>  amplitudes for each nuclear perturbation.  The validated solver must build
  !>  the Fock derivative RHS, divide/precondition by the orbital-energy denominator,
  !>  include two-electron response coupling, and pass
  !>  finite-difference gradient validation before this guard can be removed.
  subroutine hf_cphf_hessian_response(infos, hessian)
    use messages, only: show_message, WITH_ABORT
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(inout) :: hessian(:,:)

    integer :: nbasis, nocc, nvir, npert

    nbasis = int(infos%mol_prop%nbf)
    nocc = infos%mol_prop%nocc
    nvir = nbasis - nocc
    npert = 3 * int(infos%mol_prop%natom)

    hessian = hessian

    if (nocc <= 0 .or. nvir <= 0) then
      call show_message('Invalid HF CPHF occupied-virtual orbital space.', WITH_ABORT)
    end if
    if (size(hessian,1) /= npert .or. size(hessian,2) /= npert) then
      call show_message('Invalid HF CPHF Hessian response matrix shape.', WITH_ABORT)
    end if

    call show_message(&
      'Native HF CPHF Hessian partial_kernel reached: orbital-response equations are not complete yet.', &
      WITH_ABORT)
  end subroutine hf_cphf_hessian_response

end module hf_cphf_response_mod
