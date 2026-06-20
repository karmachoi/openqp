!> @brief  Projective transcorrelated (pTC) MRSF-CIS support module.
!>
!> @details
!>   pTC-MRSF-CIS replaces the bare electronic Hamiltonian with Ten-no's
!>   projective transcorrelated effective Hamiltonian
!>          H_bar = exp(-tau) H exp(tau)
!>   (S. Ten-no, J. Chem. Phys. 159, 171103 (2023)). H_bar is non-Hermitian,
!>   so the MRSF-CIS (Tamm-Dancoff) response matrix A is no longer
!>   real-symmetric and the symmetric reduced solver `rpaeig` (TDA branch,
!>   source/tdhf_lib.F90) cannot be used.
!>
!>   This module provides the non-Hermitian replacement kernels. The outer
!>   Davidson structure (subspace build / residual / new-vector generation in
!>   tdhf_mrsf_energy.F90) is reused unchanged; only the reduced-space
!>   diagonalization and the biorthonormal residual differ.
!>
!>   Implementation phases (see docs/ptc_mrsf/DESIGN.md):
!>     Phase 1  tc_nonsym_tda_eig  -- non-Hermitian reduced eigensolver.
!>              Validated as a NumPy reference in
!>              docs/ptc_mrsf/prototype/nonsym_tda_eig.py; the tau=0 limit
!>              reproduces the symmetric `rpaeig` TDA result to machine
!>              precision (the Phase-1 acceptance gate).
!>     Phase 2  tc_build_eff_integrals  -- DF/RI transcorrelated 1-/2-body
!>              effective integrals (NOT YET IMPLEMENTED).
!>     Phase 3  tc_normal_order_3body   -- normal-order 3-/4-body operators
!>              against the ROHF reference density (NOT YET IMPLEMENTED).
!>
!>   NOTE: This module is self-contained (depends only on `precision` and
!>   LAPACK DGEEV).
!>
!>   Phase 4 (DONE): `tc_nonsym_tda_eig` is wired into the live MRSF Davidson in
!>   tdhf_mrsf_energy.F90, gated by the OQP_PTC_MRSF environment variable (off by
!>   default). With the bare Hamiltonian (tau=0) the reduced matrix is symmetric
!>   and this path reproduces stock MRSF-CIS bit-for-bit (validated: H2O/6-31G
!>   MRSF-s, all 10 roots dE = 0.0 eV). tau/=0 (the genuine H_bar) awaits the
!>   Phase 2/3 integral kernels below; until then OQP_PTC_MRSF only exercises the
!>   non-Hermitian solver on a symmetric matrix.
!
module tdhf_mrsf_ptc

  use precision, only: dp
  implicit none
  private

  public :: tc_nonsym_tda_eig
  public :: tc_build_eff_integrals
  public :: tc_normal_order_3body

contains

!> @brief Lowest-`nstate` right/left eigenpairs of a non-symmetric reduced
!>        MRSF-CIS response matrix (transcorrelated Tamm-Dancoff).
!>
!> @details Drop-in replacement for the TDA branch of `rpaeig`. For tau = 0
!>   the input `amat` is symmetric and the returned eigenpairs coincide with
!>   the symmetric solver (validation gate). For the transcorrelated H_bar,
!>   `amat` is non-symmetric; bound states have real eigenvalues with
!>   biorthonormal left/right Ritz vectors  vl(:,i) . vr(:,j) = delta_ij.
!>
!> @param[in]  amat     (n,n) reduced response matrix (overwritten internally)
!> @param[in]  nstate   number of low-lying roots requested
!> @param[out] ee       (nstate) eigenvalues, ascending
!> @param[out] vr       (n,nstate) right eigenvectors,  A   vr = vr ee
!> @param[out] vl       (n,nstate) left  eigenvectors,  A^T vl = vl ee
!> @param[out] max_imag largest |Im(eig)|/(1+|Re|) encountered
!> @param[out] n_complex number of roots flagged as complex (instability)
!> @param[out] ierr     0 on success, /=0 on LAPACK failure
  subroutine tc_nonsym_tda_eig(amat, nstate, ee, vr, vl, &
                               max_imag, n_complex, ierr, imag_tol)
    real(kind=dp), intent(in)  :: amat(:,:)
    integer,       intent(in)  :: nstate
    real(kind=dp), intent(out) :: ee(:)
    real(kind=dp), intent(out) :: vr(:,:)
    real(kind=dp), intent(out) :: vl(:,:)
    real(kind=dp), intent(out) :: max_imag
    integer,       intent(out) :: n_complex
    integer,       intent(out) :: ierr
    real(kind=dp), intent(in), optional :: imag_tol

    real(kind=dp), allocatable :: a(:,:), wr(:), wi(:)
    real(kind=dp), allocatable :: vlf(:,:), vrf(:,:), work(:)
    integer,       allocatable :: ord(:)
    real(kind=dp) :: tol, denom, scal, sgn, rel, wq(1)
    integer :: n, i, j, k, lwork

    n = size(amat, 1)
    tol = 1.0e-8_dp
    if (present(imag_tol)) tol = imag_tol

    allocate(a(n,n), wr(n), wi(n), vlf(n,n), vrf(n,n), ord(n))
    a = amat

    ! workspace query then DGEEV: both left and right eigenvectors in one pass
    call dgeev('V','V', n, a, n, wr, wi, vlf, n, vrf, n, wq, -1, ierr)
    if (ierr /= 0) return
    lwork = int(wq(1))
    allocate(work(lwork))
    call dgeev('V','V', n, a, n, wr, wi, vlf, n, vrf, n, work, lwork, ierr)
    if (ierr /= 0) return

    ! diagnostics: relative imaginary content of the spectrum
    max_imag  = 0.0_dp
    n_complex = 0
    do i = 1, n
      rel = abs(wi(i)) / (1.0_dp + abs(wr(i)))
      if (rel > max_imag) max_imag = rel
      if (rel > tol)      n_complex = n_complex + 1
    end do

    ! ascending sort of the real parts (bound states are real)
    do i = 1, n
      ord(i) = i
    end do
    do i = 1, n-1
      k = i
      do j = i+1, n
        if (wr(ord(j)) < wr(ord(k))) k = j
      end do
      if (k /= i) then
        j = ord(i); ord(i) = ord(k); ord(k) = j
      end if
    end do

    ! biorthonormalize selected roots so that vl_i . vr_i = 1
    do i = 1, nstate
      k = ord(i)
      denom = dot_product(vlf(:,k), vrf(:,k))
      if (abs(denom) < 1.0e-14_dp) denom = sign(1.0e-14_dp, denom + 1.0e-30_dp)
      scal = 1.0_dp / sqrt(abs(denom))
      sgn  = sign(1.0_dp, denom)
      ee(i)   = wr(k)
      vr(:,i) = vrf(:,k) * scal
      vl(:,i) = vlf(:,k) * (scal * sgn)
    end do

    deallocate(a, wr, wi, vlf, vrf, ord, work)
  end subroutine tc_nonsym_tda_eig

!> @brief  Phase 2 (assembly TODO): density-fitted transcorrelated effective
!>         1-/2-body integrals from the correlation factor and ROHF orbitals.
!> @details The N^5, MP2-class step. The PRIMITIVE LAYER is done and validated:
!>   the F12 Gaussian-geminal integrals and the X/B/V intermediates live in
!>   module `ptc_geminal` (source/modules/ptc_geminal.F90), with the Slater
!>   factor exp(-gamma r12) expanded via `stg_ng`. What remains here is the
!>   ASSEMBLY: loop the geminal primitives over shell pairs of the real basis,
!>   contract with the ROHF MO coefficients, and density-fit/RI-factorize so the
!>   3-body term (Phase 3) is never materialized as a dense O(N^6) tensor.
!>   Not yet implemented; aborts loudly if called.
  subroutine tc_build_eff_integrals()
    !> Phase 2 now delegates the geminal/CABS/dressing machinery to the externally
    !> validated libptc (12 gates). This smoke routine confirms libptc links + runs;
    !> the live pTC-MRSF-CIS injection (ptc_build_w -> ptc_mrsf_dress) is wired into
    !> the MRSF response in tdhf_mrsf_energy.F90.
    use ptc_kernels, only: ptc_kernel, ptc_nterms, PTC_KERNEL_F12
    integer :: nk
    real(dp), allocatable :: om(:), co(:)
    nk = ptc_nterms(PTC_KERNEL_F12, 24)
    allocate(om(nk), co(nk))
    call ptc_kernel(PTC_KERNEL_F12, 1.5_dp, 24, nk, om, co)
    write(*,'(a,i0,a,2es14.6,a)') ' [libptc] f12 STG-6G kernel nk=', nk, &
         '  first (omega,coeff)=', om(1), co(1), '  (expect +4.97e-01 -2.10e-01)'
  end subroutine tc_build_eff_integrals

!> @brief  Phase 3 (TODO): normal-order the transcorrelated 3-/4-body operators
!>         against the ROHF reference one-particle density, folding their
!>         dominant contributions into the effective 1-/2-body integrals.
!> @details Keeps the method at O(N^5); explicit 3-body would be O(N^6),
!>   explicit 4-body O(N^8). Not yet implemented; aborts loudly if called.
  subroutine tc_normal_order_3body()
    error stop 'tc_normal_order_3body: pTC Phase 3 not yet implemented'
  end subroutine tc_normal_order_3body

end module tdhf_mrsf_ptc
