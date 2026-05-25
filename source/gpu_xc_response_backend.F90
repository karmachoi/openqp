module gpu_xc_response_backend
  !! Minimal Fortran-side ABI stub for the experimental CUDA TDHF/TDDFT
  !! XC-response branch.  The routines intentionally remain opt-in behind
  !! OQP_CUDA_ENABLE so CPU-only builds keep the normal response path.
  use iso_c_binding, only: c_double, c_int
  implicit none
  private

  public :: gpu_xc_response_enabled
  public :: gpu_xc_response_describe
  public :: gpu_xc_response_contract

#ifdef OQP_CUDA_ENABLE
  interface
    integer(c_int) function oqp_gpu_xc_response_contract(nbasis, nstate, density, kernel, response) bind(C, name="oqp_gpu_xc_response_contract")
      import :: c_double, c_int
      integer(c_int), value :: nbasis
      integer(c_int), value :: nstate
      real(c_double), intent(in) :: density(*)
      real(c_double), intent(in) :: kernel(*)
      real(c_double), intent(inout) :: response(*)
    end function oqp_gpu_xc_response_contract
  end interface
#else
  ! Keep the C ABI symbol name visible in source-level tests even when CUDA is
  ! disabled for ordinary CPU-only builds: oqp_gpu_xc_response_contract.
#endif

contains

  logical function gpu_xc_response_enabled()
#ifdef OQP_CUDA_ENABLE
    gpu_xc_response_enabled = .true.
#else
    gpu_xc_response_enabled = .false.
#endif
  end function gpu_xc_response_enabled

  subroutine gpu_xc_response_describe(message)
    character(len=*), intent(out) :: message
#ifdef OQP_CUDA_ENABLE
    message = "CUDA XC-response backend enabled"
#else
    message = "CUDA XC-response backend disabled; using CPU fallback"
#endif
  end subroutine gpu_xc_response_describe

  integer(c_int) function gpu_xc_response_contract(nbasis, nstate, density, kernel, response) result(status)
    integer(c_int), value :: nbasis
    integer(c_int), value :: nstate
    real(c_double), intent(in) :: density(*)
    real(c_double), intent(in) :: kernel(*)
    real(c_double), intent(inout) :: response(*)
#ifdef OQP_CUDA_ENABLE
    status = oqp_gpu_xc_response_contract(nbasis, nstate, density, kernel, response)
#else
    status = 1_c_int
#endif
  end function gpu_xc_response_contract

end module gpu_xc_response_backend
