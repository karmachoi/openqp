! QMRSF-DK dressed-kernel pole search (standalone, build-testable) -- Stage A (DFT dressing).
! Demonstrates the DK mechanism on a single-spin-flip "single" coupled to closed-shell "doubles":
!   adiabatic kernel  : drop the frequency dependence -> only the single state, the 0OS doubles are MISSED;
!   dressed kernel    : solve the frequency-dependent secular equation
!                         omega = A + sum_d |V_d|^2 / (omega - omega_d)
!                       whose roots (bracketed bisection) recover the FULL single+double spectrum,
!                       i.e. the exact Feshbach downfold of the augmented Hamiltonian.
! Validated here against the exact diagonalization of that augmented matrix (the reference), matching
! the NumPy proto (qmrsf_dk_proto.py: dressed == exact, adiabatic misses the doubles).
module dkmod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: gA
  real(dp), allocatable :: gW(:), gV(:)
contains
  real(dp) function fsec(omega)            ! secular function: omega - A - sum V_d^2/(omega-omega_d)
    real(dp), intent(in) :: omega
    integer :: d
    fsec = omega - gA
    do d = 1, size(gW); fsec = fsec - gV(d)**2 / (omega - gW(d)); end do
  end function
  real(dp) function root_bisect(lo, hi)    ! root of fsec in (lo,hi), assumes sign change
    real(dp), intent(in) :: lo, hi
    real(dp) :: a, b, m, fa, fm
    integer :: it
    a=lo; b=hi; fa=fsec(a)
    do it=1,200
      m=0.5_dp*(a+b); fm=fsec(m)
      if (fa*fm <= 0.0_dp) then; b=m; else; a=m; fa=fm; end if
      if (b-a < 1.0d-14) exit
    end do
    root_bisect = 0.5_dp*(a+b)
  end function
end module

program qmrsf_dk
  use dkmod
  implicit none
  integer, parameter :: nd = 3
  integer :: n, i, j, info, lwork, nr
  real(dp) :: Haug(1+nd,1+nd), evx(1+nd), work(256), dressed(1+nd), wsort(nd), delta
  real(dp) :: adiab_err, dmax, lo, hi

  ! ---- model: one single (A) coupled to nd well-separated doubles (omega_d, V_d) ----
  gA = 5.0_dp
  allocate(gW(nd), gV(nd))
  gW = (/ 6.0_dp, 7.5_dp, 9.0_dp /)
  gV = (/ 0.40_dp, 0.30_dp, 0.50_dp /)

  ! ---- exact: diagonalize the augmented [[A, V^T],[V, diag(w)]] ----
  n = 1+nd; Haug = 0.0_dp
  Haug(1,1) = gA
  do i=1,nd; Haug(1,1+i)=gV(i); Haug(1+i,1)=gV(i); Haug(1+i,1+i)=gW(i); end do
  lwork=256; call dsyev('N','U', n, Haug, n, evx, work, lwork, info)   ! evx = exact spectrum

  ! ---- adiabatic kernel: single-only -> one state at A, the nd doubles are missed ----
  adiab_err = abs(gA - evx(1)); do i=1,n; adiab_err = min(adiab_err, abs(gA-evx(i))); end do

  ! ---- dressed kernel: bracket the nd+1 roots of fsec between consecutive poles ----
  wsort = gW; call dsort(wsort, nd)
  delta = 1.0d-9; nr = 0
  ! interval (-inf, w1)
  nr=nr+1; dressed(nr) = root_bisect(wsort(1)-1.0d6, wsort(1)-delta)
  do i=1,nd-1                                ! between consecutive poles
    nr=nr+1; dressed(nr) = root_bisect(wsort(i)+delta, wsort(i+1)-delta)
  end do
  ! interval (w_nd, +inf)
  nr=nr+1; dressed(nr) = root_bisect(wsort(nd)+delta, wsort(nd)+1.0d6)
  call dsort(dressed, n)

  dmax = 0.0_dp; do i=1,n; dmax = max(dmax, abs(dressed(i)-evx(i))); end do

  print '(a)',       "==== QMRSF-DK dressed-kernel pole search (Fortran) ===="
  print '(a,i0,a)',  "  model: 1 single + ", nd, " doubles"
  print '(a)',       "  exact spectrum   : "
  print '(6f12.6)',  (evx(i), i=1,n)
  print '(a)',       "  dressed (pole search): "
  print '(6f12.6)',  (dressed(i), i=1,n)
  print '(a,es12.3)',"  max|dressed - exact| = ", dmax
  print '(a,f10.5,a,i0,a,f10.5)', "  adiabatic single at A=",gA," misses ",nd," double(s); nearest-exact err=", adiab_err
  if (dmax < 1.0d-9 .and. adiab_err > 0.1_dp) then
     print '(a)', "  RESULT: PASS  (dressed kernel == exact; adiabatic misses the doubles)"
  else
     print '(a)', "  RESULT: FAIL"
  end if
contains
  subroutine dsort(a,m)
    real(dp),intent(inout)::a(*); integer,intent(in)::m
    integer::i,j; real(dp)::t
    do i=1,m-1; do j=1,m-i; if(a(j)>a(j+1))then; t=a(j);a(j)=a(j+1);a(j+1)=t; end if; end do; end do
  end subroutine
end program
