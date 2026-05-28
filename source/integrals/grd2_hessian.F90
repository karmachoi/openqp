module grd2_hessian

  use precision, only: dp
  use basis_tools, only: basis_set
  use grd2, only: grd2_compute_data_t, grd2_driver
  use grd2_rys, only: grd2_int_data_t

  implicit none

  character(len=*), parameter :: module_name = "grd2_hessian"

  private
  public :: grd2_hessian_driver

contains

!###############################################################################

  !> @brief Guarded scaffold for two-electron second-derivative Hessian work.
  !> @details two_electron_hessian_der2_scaffold partial_kernel:
  !>  The existing two-electron gradient path (`grd2_driver`) already owns the
  !>  first-derivative ERI contraction machinery.  This driver reserves a native
  !>  Hessian entry point and verifies that the Rys integral workspace can be
  !>  allocated for two-electron second derivatives (`nder=2`).  Production use
  !>  must remain guarded until the second-derivative contraction and finite-
  !>  difference validation against `grd2_driver` are implemented.
  subroutine grd2_hessian_driver(infos, basis, hessian, gcomp)
    use messages, only: show_message, WITH_ABORT
    use types, only: information

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(inout) :: hessian(:,:)
    class(grd2_compute_data_t), intent(inout) :: gcomp

    type(grd2_int_data_t) :: gdat
    integer :: maxang, iok
    real(kind=dp) :: dtol, dabcut

    ! Keep references to the existing first-derivative driver and compute-data
    ! object at this boundary so future finite-difference validation compares
    ! the new two-electron second derivatives against grd2_driver, not against
    ! an unrelated integral implementation.
    hessian = hessian
    if (gcomp%attenuated) hessian = hessian

    maxang = basis%mxam
    dtol = 0.0_dp
    dabcut = 0.0_dp
    call gdat%init(maxang, 2, dtol, dabcut, iok)
    if (iok /= 0) call show_message('Cannot allocate two-electron Hessian Rys workspace.', WITH_ABORT)
    call gdat%clean()

    call show_message(&
      'Native HF/DFT two-electron Hessian partial_kernel reached: two-electron second derivatives are not complete yet.', &
      WITH_ABORT)
  end subroutine grd2_hessian_driver

end module grd2_hessian
