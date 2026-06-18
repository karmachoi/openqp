!> Standalone validation of tdhf_mrsf_ptc :: tc_nonsym_tda_eig.
!>
!> Mirrors the NumPy reference tests/ptc_mrsf/prototype/nonsym_tda_eig.py and
!> exercises the ACTUAL Fortran kernel:
!>   (1) tau=0 gate     : symmetric A -> eigenvalues match LAPACK DSYEV;
!>   (2) non-symmetric  : real spectrum, right-eigenpair residual ~0,
!>                        biorthonormality vl^T vr = I.
!>
!> Build:
!>   gfortran -J/tmp source/precision.F90 source/modules/tdhf_mrsf_ptc.F90 \
!>            tests/ptc_mrsf/prototype/tc_nonsym_eig_test.F90 -llapack -lblas -o /tmp/tctest
!>   /tmp/tctest
program tc_nonsym_eig_test
  use precision, only: dp
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none

  integer, parameter :: n = 12, nst = 4
  real(dp) :: a(n,n), asym(n,n), askew(n,n)
  real(dp) :: ee(nst), vr(n,nst), vl(n,nst)
  real(dp) :: ref(n), work(3*n), gram(nst,nst), res(n)
  real(dp) :: max_imag, perturb, biorth, maxres, maxgram, maxde
  integer  :: n_complex, ierr, i, j, k
  logical  :: ok

  call random_seed_fixed()

  ! ---- (1) tau = 0 gate: symmetric matrix ---------------------------------
  call random_number(a)
  asym = 0.5_dp*(a + transpose(a))
  do i = 1, n
    asym(i,i) = asym(i,i) + 2.0_dp*real(i,dp)   ! well-separated spectrum
  end do

  ref = pack_eigvals(asym, n, work)              ! LAPACK DSYEV reference
  call tc_nonsym_tda_eig(asym, nst, ee, vr, vl, max_imag, n_complex, ierr)
  if (ierr /= 0) stop 'DGEEV failed (sym case)'

  maxde = 0.0_dp
  do i = 1, nst
    maxde = max(maxde, abs(ee(i) - ref(i)))
  end do
  ok = (maxde < 1.0e-10_dp) .and. (n_complex == 0)
  write(*,'(a,es10.2,a,i0,a,l1)') 'tau=0 gate     max|dE|=', maxde, &
        '  n_complex=', n_complex, '   PASS=', ok
  if (.not. ok) stop 'FAIL: tau=0 gate'

  ! ---- (2) non-symmetric (transcorrelated) matrix -------------------------
  call random_number(askew)
  askew = askew - transpose(askew)               ! antisymmetric part
  perturb = 0.05_dp
  a = asym + perturb*askew

  call tc_nonsym_tda_eig(a, nst, ee, vr, vl, max_imag, n_complex, ierr)
  if (ierr /= 0) stop 'DGEEV failed (nonsym case)'

  ! right-eigenpair residual  A vr - vr ee
  maxres = 0.0_dp
  do j = 1, nst
    res = matmul(a, vr(:,j)) - ee(j)*vr(:,j)
    maxres = max(maxres, maxval(abs(res)))
  end do
  ! biorthonormality  vl^T vr = I
  gram = matmul(transpose(vl), vr)
  maxgram = 0.0_dp
  do i = 1, nst
    do j = 1, nst
      maxgram = max(maxgram, abs(gram(i,j) - merge(1.0_dp,0.0_dp,i==j)))
    end do
  end do
  ok = (maxres < 1.0e-9_dp) .and. (maxgram < 1.0e-9_dp) .and. (n_complex == 0)
  write(*,'(a,es10.2,a,es10.2,a,es10.2,a,l1)') 'non-symmetric  res=', &
        maxres, '  biorth=', maxgram, '  Im=', max_imag, '   PASS=', ok
  if (.not. ok) stop 'FAIL: non-symmetric case'

  write(*,'(a)') 'ALL FORTRAN KERNEL TESTS PASSED'

contains

  subroutine random_seed_fixed()
    integer :: sz
    integer, allocatable :: seed(:)
    call random_seed(size=sz)
    allocate(seed(sz)); seed = 20260618
    call random_seed(put=seed)
  end subroutine random_seed_fixed

  function pack_eigvals(m, nn, wk) result(w)
    integer, intent(in) :: nn
    real(dp), intent(in) :: m(nn,nn)
    real(dp), intent(inout) :: wk(:)
    real(dp) :: w(nn), mm(nn,nn)
    integer :: info
    mm = m
    call dsyev('N','U', nn, mm, nn, w, wk, size(wk), info)
    if (info /= 0) stop 'DSYEV failed'
  end function pack_eigvals

end program tc_nonsym_eig_test
