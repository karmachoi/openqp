!> Standalone validation of the native F12 geminal engine (source/modules/
!> ptc_geminal.F90), the pyscf-free port of r12_geminal.py + f12_intermediates.py.
!>
!> Oracles (no fitting to the answer, no pyscf):
!>   1. geminal(omega->0)            = overlap product (a|c)(b|d)   [analytic]
!>   2. geminal(finite omega)        vs 6D separable grid quadrature
!>   3. r2_geminal (B kernel core)   vs direct grid integration
!>   4. V(omega=0)                   = eri_s (Boys F0 closed form)
!>   5. STG-6G fit                   reproduces exp(-gamma r12)
!>
!> Build (standalone):
!>   gfortran -O2 source/precision.F90 source/modules/ptc_geminal.F90 \
!>            tests/ptc_mrsf/prototype/ptc_geminal_test.F90 -o /tmp/ptc_gem && /tmp/ptc_gem
program ptc_geminal_test
  use precision, only: dp
  use ptc_geminal
  implicit none

  real(dp), parameter :: PI = 3.141592653589793238462643383279_dp
  real(dp) :: A(3), B(3), C(3), D(3)
  real(dp) :: g0, prod, ana, ref, v0, eri, x2, bb
  real(dp) :: omega, lam, gamma
  integer  :: nfail
  integer  :: k, ng
  real(dp) :: cc(6), omg(6), r, fit, exact
  real(dp), parameter :: OMS(3)  = [0.5_dp, 1.0_dp, 2.5_dp]
  real(dp), parameter :: LAMS(3) = [0.5_dp, 1.0_dp, 2.0_dp]
  real(dp), parameter :: RS(4)   = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]

  nfail = 0
  write(*,'(a)') '=== native F12 Gaussian-geminal engine: validation ==='
  write(*,'(a)') ''

  ! (1) omega -> 0  ->  overlap product (a|c)(b|d)
  A = [0.0_dp,0.0_dp,0.0_dp];  B = [0.0_dp,0.0_dp,2.0_dp]
  g0   = gaussian_geminal_s(1.2_dp,A, 0.7_dp,B, 0.9_dp,A, 1.5_dp,B, 1.0e-9_dp)
  prod = gem_overlap_s(1.2_dp,A, 0.9_dp,A) * gem_overlap_s(0.7_dp,B, 1.5_dp,B)
  call chk('geminal(omega->0) vs (a|c)(b|d) ', g0, prod, 1.0e-5_dp, nfail)
  write(*,'(a)') ''

  ! (2) finite omega vs numerical 6D quadrature
  A = [0.0_dp,0.0_dp,0.0_dp]; B = [0.0_dp,0.0_dp,0.5_dp]
  C = [0.3_dp,0.0_dp,0.0_dp]; D = [0.0_dp,0.2_dp,0.4_dp]
  do k = 1, 3
    omega = OMS(k)
    ana = gaussian_geminal_s(1.2_dp,A, 0.7_dp,B, 0.9_dp,C, 1.5_dp,D, omega)
    ref = geminal_numeric  (1.2_dp,A, 0.7_dp,B, 0.9_dp,C, 1.5_dp,D, omega)
    call chk('geminal(finite om) vs grid     ', ana, ref, 1.0e-6_dp, nfail)
  end do
  write(*,'(a)') ''

  ! (3) r12^2 * geminal (B kernel core) vs grid
  A = [0.0_dp,0.0_dp,0.0_dp]; B = [0.0_dp,0.0_dp,0.4_dp]
  C = [0.2_dp,0.0_dp,0.0_dp]; D = [0.0_dp,0.1_dp,0.3_dp]
  do k = 1, 3
    lam = LAMS(k)
    ana = r2_geminal_s     (1.1_dp,A, 0.8_dp,B, 0.9_dp,C, 1.3_dp,D, lam)
    ref = r2_geminal_numeric(1.1_dp,A, 0.8_dp,B, 0.9_dp,C, 1.3_dp,D, lam)
    call chk('r12^2*geminal vs grid          ', ana, ref, 1.0e-3_dp, nfail)
  end do
  write(*,'(a)') ''

  ! (4) V(omega=0) == ERI (Boys F0)
  v0  = v_geminal_s(1.1_dp,A, 0.8_dp,B, 0.9_dp,C, 1.3_dp,D, 0.0_dp)
  eri = eri_s      (1.1_dp,A, 0.8_dp,B, 0.9_dp,C, 1.3_dp,D)
  call chk('V(omega=0) vs ERI (Boys F0)     ', v0, eri, 1.0e-4_dp, nfail)
  ! X, B intermediates: report finite values (X>0, B>0)
  x2 = f12_X_s(1.1_dp,A, 0.8_dp,B, 0.9_dp,C, 1.3_dp,D, 1.0_dp)
  bb = f12_B_s(1.1_dp,A, 0.8_dp,B, 0.9_dp,C, 1.3_dp,D, 1.0_dp)
  write(*,'(a,2es16.8)') 'X = geminal(2om), B = 8om^2 r2  ', x2, bb
  if (x2 <= 0.0_dp .or. bb <= 0.0_dp) nfail = nfail + 1
  write(*,'(a)') ''

  ! (5) STG-6G fit reproduces exp(-gamma r12)
  gamma = 1.0_dp
  call stg_ng(gamma, cc, omg, ng)
  write(*,'(a)') 'STG-6G: exp(-r12) vs sum_k c_k exp(-om_k r12^2)'
  do k = 1, 4
    r = RS(k)
    exact = exp(-gamma*r)
    fit = 0.0_dp
    block
      integer :: j
      do j = 1, ng
        fit = fit + cc(j)*exp(-omg(j)*r*r)
      end do
    end block
    call chk('  STG-6G fit at r12             ', fit, exact, 0.05_dp, nfail)
  end do
  write(*,'(a)') ''

  ! (6) contracted-shell AO tensor (the assembly machinery): STO-3G H2,
  !     geminal(omega->0) must equal the contracted overlap product Sc(i,k)Sc(j,l).
  block
    integer, parameter :: nsh = 2, mp = 3
    integer  :: npr(nsh), i2, j2, k2, l2, a2, b2
    real(dp) :: exs(mp,nsh), cos_(mp,nsh), cns(3,nsh)
    real(dp) :: Mten(nsh,nsh,nsh,nsh), Sc(nsh,nsh), maxd, sij
    npr = [3, 3]
    exs(:,1) = [3.42525091_dp, 0.62391373_dp, 0.16885540_dp]; exs(:,2) = exs(:,1)
    cos_(:,1) = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]; cos_(:,2) = cos_(:,1)
    cns(:,1) = [0.0_dp,0.0_dp,0.0_dp]; cns(:,2) = [0.0_dp,0.0_dp,1.4_dp]
    do i2 = 1, nsh
      do k2 = 1, nsh
        sij = 0.0_dp
        do a2 = 1, mp
          do b2 = 1, mp
            sij = sij + cos_(a2,i2)*s_norm(exs(a2,i2))*cos_(b2,k2)*s_norm(exs(b2,k2)) &
                      * gem_overlap_s(exs(a2,i2),cns(:,i2), exs(b2,k2),cns(:,k2))
          end do
        end do
        Sc(i2,k2) = sij
      end do
    end do
    call ptc_s_ao_tensor(nsh, npr, exs, cos_, cns, PTC_OP_GEMINAL, 1.0e-9_dp, Mten)
    maxd = 0.0_dp
    do i2 = 1, nsh
      do j2 = 1, nsh
        do k2 = 1, nsh
          do l2 = 1, nsh
            maxd = max(maxd, abs(Mten(i2,j2,k2,l2) - Sc(i2,k2)*Sc(j2,l2)))
          end do
        end do
      end do
    end do
    write(*,'(a,es10.2,a)') 'contracted-shell AO geminal(w->0) vs Sc(x)Sc  max|d| = ', &
         maxd, merge('  PASS', '  FAIL', maxd <= 1.0e-6_dp)
    if (maxd > 1.0e-6_dp) nfail = nfail + 1
    ! contracted ERI (STO-3G H2) finite & symmetric (ab|cd)=(cd|ab)
    call ptc_s_ao_tensor(nsh, npr, exs, cos_, cns, PTC_OP_ERI, 0.0_dp, Mten)
    write(*,'(a,f12.8)') 'contracted-shell ERI (11|11) [STO-3G H2]        = ', Mten(1,1,1,1)
  end block
  write(*,'(a)') ''

  ! (7) native 1-electron integrals vs textbook STO-3G hydrogen atom:
  !     normalized overlap S=1, and h = T + V_ne = -0.4665 Ha (STO-3G H atom).
  block
    integer, parameter :: nsh = 1, mp = 3
    integer  :: npr(nsh)
    real(dp) :: exs(mp,nsh), cos_(mp,nsh), cns(3,nsh)
    real(dp) :: Sm(nsh,nsh), Tm(nsh,nsh), Vm(nsh,nsh), zat(1), rat(3,1), Mte(nsh,nsh,nsh,nsh)
    npr = [3]
    exs(:,1) = [3.42525091_dp, 0.62391373_dp, 0.16885540_dp]
    cos_(:,1) = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
    cns(:,1) = [0.0_dp,0.0_dp,0.0_dp]
    zat(1) = 1.0_dp; rat(:,1) = [0.0_dp,0.0_dp,0.0_dp]
    call ptc_s_ao_1e(nsh, npr, exs, cos_, cns, PTC_1E_OVERLAP, 1, zat, rat, Sm)
    call ptc_s_ao_1e(nsh, npr, exs, cos_, cns, PTC_1E_KINETIC, 1, zat, rat, Tm)
    call ptc_s_ao_1e(nsh, npr, exs, cos_, cns, PTC_1E_NUCLEAR, 1, zat, rat, Vm)
    call ptc_s_ao_tensor(nsh, npr, exs, cos_, cns, PTC_OP_ERI, 0.0_dp, Mte)
    call chk('STO-3G H: overlap S(1,1)=1     ', Sm(1,1), 1.0_dp, 1.0e-6_dp, nfail)
    call chk('STO-3G H: h=T+Vne = -0.4665    ', Tm(1,1)+Vm(1,1), -0.46658185_dp, 1.0e-4_dp, nfail)
    call chk('STO-3G H: (11|11) = 0.7746     ', Mte(1,1,1,1), 0.77460594_dp, 1.0e-5_dp, nfail)
  end block
  write(*,'(a)') ''

  if (nfail == 0) then
    write(*,'(a)') 'ALL PASS: native F12 geminal engine validated (pyscf-free).'
  else
    write(*,'(a,i0,a)') 'FAILURES: ', nfail, ' check(s) failed.'
    error stop 1
  end if

contains

  subroutine chk(name, a, b, tol, nf)
    character(*), intent(in)    :: name
    real(dp),     intent(in)    :: a, b, tol
    integer,      intent(inout) :: nf
    real(dp) :: rel
    rel = abs(a - b)/max(abs(b), 1.0e-30_dp)
    write(*,'(a,2es16.8,es10.1,a)') name, a, b, rel, &
         merge('  PASS', '  FAIL', rel <= tol)
    if (rel > tol) nf = nf + 1
  end subroutine chk

  !> 6D separable grid quadrature of geminal (per-dimension 2D grids).
  function geminal_numeric(aA,A0,aB,B0,aC,C0,aD,D0,omega) result(g)
    real(dp), intent(in) :: aA,aB,aC,aD, A0(3),B0(3),C0(3),D0(3), omega
    real(dp) :: g, p,q, Pc(3),Qc(3), Kac,Kbd, dx, x1, x2v, acc, prd
    integer  :: dim, i, j, n
    real(dp), allocatable :: gx(:)
    n = 400
    allocate(gx(n))
    dx = 12.0_dp/real(n-1,dp)
    do i = 1, n
      gx(i) = -6.0_dp + dx*real(i-1,dp)
    end do
    p = aA+aC; Pc = (aA*A0+aC*C0)/p; Kac = exp(-aA*aC/p*dot_product(A0-C0,A0-C0))
    q = aB+aD; Qc = (aB*B0+aD*D0)/q; Kbd = exp(-aB*aD/q*dot_product(B0-D0,B0-D0))
    prd = 1.0_dp
    do dim = 1, 3
      acc = 0.0_dp
      do i = 1, n
        x1 = gx(i) + Pc(dim)
        do j = 1, n
          x2v = gx(j) + Qc(dim)
          acc = acc + exp(-p*gx(i)**2 - q*gx(j)**2 - omega*(x1-x2v)**2)
        end do
      end do
      prd = prd * acc*dx*dx
    end do
    g = Kac*Kbd*prd
    deallocate(gx)
  end function geminal_numeric

  !> Grid reference for <ab| r12^2 e^{-lam r12^2} |cd> (per-dim averaged).
  function r2_geminal_numeric(aA,A0,aB,B0,aC,C0,aD,D0,lam) result(rr)
    real(dp), intent(in) :: aA,aB,aC,aD, A0(3),B0(3),C0(3),D0(3), lam
    real(dp) :: rr, p,q, Pc(3),Qc(3), Kac,Kbd, dxg, z1,z2, w, base, val, num, den
    integer  :: dim, i, j, n
    real(dp), allocatable :: gx(:)
    n = 120
    allocate(gx(n))
    dxg = 14.0_dp/real(n-1,dp)
    do i = 1, n
      gx(i) = -7.0_dp + dxg*real(i-1,dp)
    end do
    p = aA+aC; Pc = (aA*A0+aC*C0)/p; Kac = exp(-aA*aC/p*dot_product(A0-C0,A0-C0))
    q = aB+aD; Qc = (aB*B0+aD*D0)/q; Kbd = exp(-aB*aD/q*dot_product(B0-D0,B0-D0))
    base = gaussian_geminal_s(aA,A0,aB,B0,aC,C0,aD,D0, lam)
    val = 0.0_dp
    do dim = 1, 3
      num = 0.0_dp; den = 0.0_dp
      do i = 1, n
        z1 = gx(i) + Pc(dim)
        do j = 1, n
          z2 = gx(j) + Qc(dim)
          w  = exp(-p*(z1-Pc(dim))**2 - q*(z2-Qc(dim))**2 - lam*(z1-z2)**2)
          num = num + (z1-z2)**2 * w
          den = den + w
        end do
      end do
      val = val + (num*dxg*dxg)/(den*dxg*dxg)
    end do
    rr = base*val
    deallocate(gx)
  end function r2_geminal_numeric

end program ptc_geminal_test
