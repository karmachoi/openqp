!> End-to-end native validation of the pTC transcorrelation mechanism on the
!> exactly-solvable half-filled Hubbard dimer, driving the PRODUCTION solver
!> tc_nonsym_tda_eig (source/modules/tdhf_mrsf_ptc.F90). Pyscf-free port of
!> tests/ptc_mrsf/prototype/tc_hubbard_demo.py.
!>
!> Demonstrates the complete chain that pTC-MRSF-CIS exploits:
!>   H_bar = J^{-1} H J (non-Hermitian, Gutzwiller J=g^docc) evaluated in a single
!>   mean-field determinant recovers ~100% of the correlation energy, and the
!>   non-Hermitian solver returns the EXACT spectrum.
!>
!> Oracles (all closed-form / exact, no pyscf):
!>   E_exact = 1/2 (U - sqrt(U^2+16 t^2)),  E_RHF = U/2 - 2t.
!>
!> Build:
!>   gfortran -O2 source/precision.F90 source/modules/tdhf_mrsf_ptc.F90 \
!>     tests/ptc_mrsf/prototype/tc_hubbard_test.F90 -llapack -lblas -o /tmp/ptc_hub && /tmp/ptc_hub
program tc_hubbard_test
  use precision, only: dp
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none

  integer, parameter :: nd = 4
  integer  :: bas(nd), docc(nd)
  real(dp) :: H(nd,nd), Hbar(nd,nd), g2(nd), jd(nd), jpsi(nd)
  real(dp) :: ev(nd), ee(nd), vr(nd,nd), vl(nd,nd), biow(nd,nd)
  real(dp) :: t, U, e_exact, e_exact_cf, e_rhf, e_rhf_cf, e_corr
  real(dp) :: g, gopt, eg, eg_min, recov, e_proj, resid, maxim, bime, spece
  real(dp), allocatable :: wrk(:)
  integer  :: i, j, k, lwk, ierr, ncx
  integer  :: nfail

  nfail = 0
  t = 1.0_dp; U = 4.0_dp

  ! Ms=0 determinant basis (bit masks over spin-orbitals 0=1up,1=1dn,2=2up,3=2dn):
  !   |1u1d>=3, |2u2d>=12, |1u2d>=9, |2u1d>=6
  bas  = [3, 12, 9, 6]
  do i = 1, nd
    docc(i) = 0
    if (btest(bas(i),0) .and. btest(bas(i),1)) docc(i) = docc(i) + 1   ! site 1
    if (btest(bas(i),2) .and. btest(bas(i),3)) docc(i) = docc(i) + 1   ! site 2
  end do

  call build_hubbard(t, U, bas, docc, H)

  ! exact spectrum via LAPACK DSYEV
  ev = 0.0_dp
  call sym_eig(H, nd, ev)
  e_exact = ev(1)
  e_exact_cf = 0.5_dp*(U - sqrt(U*U + 16.0_dp*t*t))
  call chk('E_exact vs closed form ', e_exact, e_exact_cf, 1.0e-10_dp, nfail)

  ! mean-field (RHF) determinant |g^2> = (c1+c2)up (c1+c2)dn / 2 in this basis
  g2 = [0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]
  g2 = g2 / sqrt(dot_product(g2,g2))
  e_rhf = dot_product(g2, matmul(H, g2))
  e_rhf_cf = 0.5_dp*U - 2.0_dp*t
  call chk('E_RHF   vs closed form ', e_rhf, e_rhf_cf, 1.0e-10_dp, nfail)

  e_corr = e_exact - e_rhf
  write(*,'(a)') ''
  write(*,'(a,2f12.6)') 'Hubbard dimer (t,U)            = ', t, U
  write(*,'(a,f14.8)')  'exact ground energy           = ', e_exact
  write(*,'(a,f14.8)')  'mean-field (RHF) energy       = ', e_rhf
  write(*,'(a,f14.8)')  'correlation energy            = ', e_corr
  write(*,'(a)') ''

  ! variational Gutzwiller scan: E_G(g) = <g2|J H J|g2>/<g2|J^2|g2>
  eg_min = 1.0e30_dp; gopt = 0.0_dp
  do k = 1, 1000
    g = 0.05_dp + (1.0_dp - 0.05_dp)*real(k-1,dp)/999.0_dp
    do i = 1, nd
      jd(i) = g**docc(i)
    end do
    jpsi = jd*g2
    eg = dot_product(jpsi, matmul(H, jpsi)) / dot_product(jpsi, jpsi)
    if (eg < eg_min) then
      eg_min = eg; gopt = g
    end if
  end do
  recov = (eg_min - e_rhf)/e_corr*100.0_dp
  write(*,'(a,f10.4)') 'Gutzwiller optimal g          = ', gopt
  write(*,'(a,f14.8)') 'E_G(g_opt)                    = ', eg_min
  write(*,'(a,f8.2,a)') 'correlation recovered         = ', recov, ' %'
  if (eg_min < e_exact - 1.0e-9_dp) nfail = nfail + 1     ! variational bound
  if (recov <= 99.0_dp) nfail = nfail + 1
  write(*,'(a)') ''

  ! transcorrelated H_bar = J^{-1} H J at g_opt (non-Hermitian)
  do i = 1, nd
    jd(i) = gopt**docc(i)
  end do
  do i = 1, nd
    do j = 1, nd
      Hbar(i,j) = H(i,j)*jd(j)/jd(i)
    end do
  end do
  e_proj = dot_product(g2, matmul(Hbar, g2)) / dot_product(g2, g2)
  resid  = sqrt(sum((matmul(Hbar, g2) - e_proj*g2)**2))
  call chk('projective <g2|Hbar|g2> = E_exact', e_proj, e_exact, 5.0e-3_dp, nfail)
  write(*,'(a,es10.2)') 'right-eigenvector residual |Hbar g2 - E g2| = ', resid
  if (resid > 5.0e-3_dp) nfail = nfail + 1
  if (maxval(abs(Hbar - transpose(Hbar))) < 1.0e-10_dp) then
    write(*,'(a)') 'ERROR: H_bar is symmetric (should be non-Hermitian)'
    nfail = nfail + 1
  end if
  write(*,'(a)') ''

  ! production non-Hermitian solver on H_bar: must return the exact spectrum
  call tc_nonsym_tda_eig(Hbar, nd, ee, vr, vl, maxim, ncx, ierr)
  if (ierr /= 0) then
    write(*,'(a,i0)') 'tc_nonsym_tda_eig ierr = ', ierr
    nfail = nfail + 1
  end if
  call sort_asc(ee, nd)
  spece = maxval(abs(ee - ev))
  biow = matmul(transpose(vl), vr)
  do i = 1, nd
    biow(i,i) = biow(i,i) - 1.0_dp
  end do
  bime = maxval(abs(biow))
  write(*,'(a,es10.2)') 'non-Herm solver: max Im(eig)/(1+|Re|)      = ', maxim
  call chk('Hbar spectrum vs exact ', spece, 0.0_dp, 1.0e-9_dp, nfail)
  write(*,'(a,es10.2)') 'biorthonormality |vl^T vr - I|             = ', bime
  if (maxim > 1.0e-8_dp) nfail = nfail + 1
  if (bime > 1.0e-8_dp) nfail = nfail + 1
  write(*,'(a)') ''

  if (nfail == 0) then
    write(*,'(a)') 'ALL PASS: native end-to-end pTC transcorrelation on the Hubbard'
    write(*,'(a)') 'dimer -- H_bar + tc_nonsym_tda_eig recover the EXACT spectrum,'
    write(*,'(a)') 'and a mean-field determinant recovers ~100% of the correlation.'
  else
    write(*,'(a,i0,a)') 'FAILURES: ', nfail, ' check(s) failed.'
    error stop 1
  end if

contains

  subroutine chk(name, a, b, tol, nf)
    character(*), intent(in)    :: name
    real(dp),     intent(in)    :: a, b, tol
    integer,      intent(inout) :: nf
    real(dp) :: d
    d = abs(a - b)
    write(*,'(a,2f16.8,es10.2,a)') name, a, b, d, &
         merge('  PASS', '  FAIL', d <= tol)
    if (d > tol) nf = nf + 1
  end subroutine chk

  subroutine build_hubbard(t, U, bas, docc, H)
    real(dp), intent(in)  :: t, U
    integer,  intent(in)  :: bas(:), docc(:)
    real(dp), intent(out) :: H(:,:)
    integer :: hops(2,4), i, h_, sgn, ndet, p, q, jnew, jj
    H = 0.0_dp
    ndet = size(bas)
    hops = reshape([0,2, 2,0, 1,3, 3,1], [2,4])
    do i = 1, ndet
      H(i,i) = H(i,i) + U*real(docc(i),dp)
      do h_ = 1, 4
        p = hops(1,h_); q = hops(2,h_)
        sgn = adag_a(bas(i), p, q, jnew)
        if (sgn /= 0) then
          do jj = 1, ndet
            if (bas(jj) == jnew) H(jj,i) = H(jj,i) - t*real(sgn,dp)
          end do
        end if
      end do
    end do
  end subroutine build_hubbard

  !> c^dag_p c_q on a bitmask determinant; returns sign (0 if killed) and newdet.
  integer function adag_a(det, p, q, newdet) result(sgn)
    integer, intent(in)  :: det, p, q
    integer, intent(out) :: newdet
    integer :: m
    sgn = 0; newdet = 0
    if (.not. btest(det, q)) return
    sgn = parity_below(det, q)
    m = ibclr(det, q)
    if (btest(m, p)) then
      sgn = 0; return
    end if
    sgn = sgn * parity_below(m, p)
    newdet = ibset(m, p)
  end function adag_a

  !> (-1)^(number of set bits in positions 0..x-1 of d)
  integer function parity_below(d, x) result(s)
    integer, intent(in) :: d, x
    integer :: i, c
    c = 0
    do i = 0, x-1
      if (btest(d, i)) c = c + 1
    end do
    s = 1 - 2*mod(c, 2)
  end function parity_below

  subroutine sym_eig(A, n, w)
    real(dp), intent(in)  :: A(:,:)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: w(:)
    real(dp) :: Acpy(n,n), wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    Acpy = A(1:n,1:n)
    call dsyev('N','U', n, Acpy, n, w, wq, -1, info)
    lw = int(wq(1)); allocate(wk(lw))
    call dsyev('N','U', n, Acpy, n, w, wk, lw, info)
    deallocate(wk)
  end subroutine sym_eig

  subroutine sort_asc(a, n)
    real(dp), intent(inout) :: a(:)
    integer,  intent(in)    :: n
    integer :: i, j
    real(dp) :: tmp
    do i = 1, n-1
      do j = i+1, n
        if (a(j) < a(i)) then
          tmp = a(i); a(i) = a(j); a(j) = tmp
        end if
      end do
    end do
  end subroutine sort_asc

end program tc_hubbard_test
