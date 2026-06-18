!> Validate the McMurchie-Davidson general-L integrals (ptc_md): for s functions
!> they must reproduce the s-only engine (ptc_geminal); for p functions they must
!> match independent numerical quadrature. ERI checked vs eri_s (ssss) + symmetry.
program tc_md_test
  use precision, only: dp
  use ptc_md
  use ptc_geminal, only: gem_overlap_s, kinetic_s, nuclear_s, eri_s, gaussian_geminal_s
  implicit none
  real(dp) :: A(3), B(3), C(3), D(3), ze, zf, zg, zh, zc
  integer  :: s0(3), px(3), py(3)
  integer  :: nfail
  real(dp) :: v1, v2

  nfail = 0
  s0 = [0,0,0]; px = [1,0,0]; py = [0,1,0]
  A = [0.0_dp,0.0_dp,0.0_dp]; B = [0.0_dp,0.0_dp,1.4_dp]
  C = [0.3_dp,0.0_dp,0.2_dp]; D = [0.0_dp,0.5_dp,0.9_dp]
  ze = 1.2_dp; zf = 0.8_dp; zg = 0.9_dp; zh = 1.3_dp; zc = 1.0_dp

  write(*,'(a)') '=== McMurchie-Davidson general-L integrals: validation ==='
  ! --- s functions vs s-only engine ---
  call chk('overlap ss vs s-engine ', overlap_cart(s0,A,ze,s0,B,zf), gem_overlap_s(ze,A,zf,B), 1e-12_dp, nfail)
  call chk('kinetic ss vs s-engine ', kinetic_cart(s0,A,ze,s0,B,zf), kinetic_s(ze,A,zf,B), 1e-12_dp, nfail)
  call chk('nuclear ss vs s-engine ', nuclear_cart(s0,A,ze,s0,B,zf,zc,C), nuclear_s(ze,A,zf,B,zc,C), 1e-12_dp, nfail)
  call chk('eri ssss vs s-engine   ', eri_cart(s0,A,ze,s0,B,zf,s0,C,zg,s0,D,zh), eri_s(ze,A,zg,C,zf,B,zh,D), 1e-10_dp, nfail)
  write(*,'(a)') ''
  ! --- p functions vs numerical quadrature ---
  v1 = overlap_cart(px,A,ze,px,B,zf); v2 = num_overlap(px,A,ze,px,B,zf)
  call chk('overlap px-px vs grid  ', v1, v2, 1e-6_dp, nfail)
  v1 = overlap_cart(px,A,ze,py,B,zf); v2 = num_overlap(px,A,ze,py,B,zf)
  if (abs(v1) > 1e-10_dp .or. abs(v2) > 1e-10_dp) nfail=nfail+1
  write(*,'(a,2es17.9,a)') 'overlap px-py (both ~0)  ', v1, v2, '  PASS'
  v1 = kinetic_cart(px,A,ze,px,B,zf); v2 = num_kinetic(px,A,ze,px,B,zf)
  call chk('kinetic px-px vs grid  ', v1, v2, 1e-5_dp, nfail)
  v1 = nuclear_cart(px,A,ze,px,B,zf,zc,C); v2 = num_nuclear(px,A,ze,px,B,zf,zc,C)
  call chk('nuclear px-px vs grid  ', v1, v2, 3e-3_dp, nfail)
  write(*,'(a)') ''
  ! --- ERI permutational symmetry with p functions ---
  v1 = eri_cart(px,A,ze,s0,B,zf,s0,C,zg,s0,D,zh)
  v2 = eri_cart(s0,B,zf,px,A,ze,s0,C,zg,s0,D,zh)
  call chk('eri (pa,s|s,s)=(s,pa|s,s)', v1, v2, 1e-12_dp, nfail)
  v2 = eri_cart(s0,C,zg,s0,D,zh,px,A,ze,s0,B,zf)
  call chk('eri bra<->ket symmetry  ', v1, v2, 1e-12_dp, nfail)
  v1 = eri_cart(px,A,ze,px,B,zf,py,C,zg,py,D,zh)
  v2 = num_eri(px,A,ze,px,B,zf,py,C,zg,py,D,zh)
  call chk('eri (px px|py py) vs grid', v1, v2, 6e-2_dp, nfail)
  write(*,'(a)') ''

  ! --- Gaussian geminal at general L ---
  v1 = geminal_cart(s0,A,ze,s0,B,zf,s0,C,zg,s0,D,zh, 0.8_dp)
  v2 = gaussian_geminal_s(ze,A,zg,C,zf,B,zh,D, 0.8_dp)
  call chk('geminal ssss vs s-engine', v1, v2, 1e-10_dp, nfail)
  v1 = geminal_cart(px,A,ze,px,B,zf,py,C,zg,py,D,zh, 0.8_dp)
  v2 = num_geminal(px,A,ze,px,B,zf,py,C,zg,py,D,zh, 0.8_dp)
  call chk('geminal pppp vs grid    ', v1, v2, 6e-2_dp, nfail)
  write(*,'(a)') ''

  if (nfail == 0) then
    write(*,'(a)') 'ALL PASS: McMurchie-Davidson overlap/kinetic/nuclear/ERI validated.'
  else
    write(*,'(a,i0,a)') 'FAILURES: ', nfail, ' check(s).'
    error stop 1
  end if

contains

  subroutine chk(name, x, y, tol, nf)
    character(*), intent(in) :: name
    real(dp), intent(in) :: x, y, tol
    integer, intent(inout) :: nf
    real(dp) :: rel
    rel = abs(x-y)/max(abs(y),1e-30_dp)
    write(*,'(a,2es17.9,es10.1,a)') name, x, y, rel, merge('  PASS','  FAIL', rel<=tol)
    if (rel > tol) nf = nf + 1
  end subroutine chk

  pure function gval(l, R0, e, r) result(g)
    integer, intent(in) :: l(3)
    real(dp), intent(in) :: R0(3), e, r(3)
    real(dp) :: g, dr(3)
    dr = r - R0
    g = dr(1)**l(1) * dr(2)**l(2) * dr(3)**l(3) * exp(-e*dot_product(dr,dr))
  end function gval

  function num_overlap(la,A,ze,lb,B,zf) result(s)
    integer, intent(in) :: la(3), lb(3)
    real(dp), intent(in) :: A(3), B(3), ze, zf
    real(dp) :: s, r(3), h, lo, hi
    integer :: i, j, k, n
    n = 60; lo = -6.0_dp; hi = 6.0_dp; h = (hi-lo)/n
    s = 0.0_dp
    do i=0,n; do j=0,n; do k=0,n
      r = [lo+i*h, lo+j*h, lo+k*h]
      s = s + gval(la,A,ze,r)*gval(lb,B,zf,r)
    end do; end do; end do
    s = s * h**3
  end function num_overlap

  function num_kinetic(la,A,ze,lb,B,zf) result(t)
    integer, intent(in) :: la(3), lb(3)
    real(dp), intent(in) :: A(3), B(3), ze, zf
    real(dp) :: t, r(3), h, lo, hi, lap, eps
    integer :: i, j, k, n, dd
    n = 70; lo = -6.0_dp; hi = 6.0_dp; h = (hi-lo)/n; eps = 1e-4_dp
    t = 0.0_dp
    do i=0,n; do j=0,n; do k=0,n
      r = [lo+i*h, lo+j*h, lo+k*h]
      lap = 0.0_dp
      do dd=1,3
        block
          real(dp) :: rp(3), rm(3)
          rp=r; rm=r; rp(dd)=r(dd)+eps; rm(dd)=r(dd)-eps
          lap = lap + (gval(lb,B,zf,rp)-2*gval(lb,B,zf,r)+gval(lb,B,zf,rm))/eps**2
        end block
      end do
      t = t + gval(la,A,ze,r)*(-0.5_dp*lap)
    end do; end do; end do
    t = t * h**3
  end function num_kinetic

  function num_nuclear(la,A,ze,lb,B,zf,zc,C) result(v)
    integer, intent(in) :: la(3), lb(3)
    real(dp), intent(in) :: A(3), B(3), ze, zf, zc, C(3)
    real(dp) :: v, r(3), h, lo, hi, rc
    integer :: i, j, k, n
    n = 80; lo = -6.0_dp; hi = 6.0_dp; h = (hi-lo)/n
    v = 0.0_dp
    do i=0,n; do j=0,n; do k=0,n
      r = [lo+i*h, lo+j*h, lo+k*h]
      rc = sqrt(dot_product(r-C,r-C))
      if (rc < 1e-6_dp) cycle
      v = v - zc*gval(la,A,ze,r)*gval(lb,B,zf,r)/rc
    end do; end do; end do
    v = v * h**3
  end function num_nuclear

  function num_eri(la,A,ze,lb,B,zf,lc,C,zg,ld,D,zh) result(g)
    integer, intent(in) :: la(3),lb(3),lc(3),ld(3)
    real(dp), intent(in) :: A(3),B(3),C(3),D(3), ze,zf,zg,zh
    real(dp) :: g, r1(3), r2(3), h, lo, hi, r12, rho1
    integer :: i1,j1,k1,i2,j2,k2, n
    n = 18; lo = -4.0_dp; hi = 4.0_dp; h = (hi-lo)/n
    g = 0.0_dp
    do i1=0,n; do j1=0,n; do k1=0,n
      r1 = [lo+i1*h, lo+j1*h, lo+k1*h]
      rho1 = gval(la,A,ze,r1)*gval(lb,B,zf,r1)
      if (abs(rho1) < 1e-14_dp) cycle
      do i2=0,n; do j2=0,n; do k2=0,n
        r2 = [lo+i2*h, lo+j2*h, lo+k2*h]
        r12 = sqrt(dot_product(r1-r2,r1-r2))
        if (r12 < 1e-6_dp) cycle
        g = g + rho1*gval(lc,C,zg,r2)*gval(ld,D,zh,r2)/r12
      end do; end do; end do
    end do; end do; end do
    g = g * h**6
  end function num_eri

  function num_geminal(la,A,ze,lb,B,zf,lc,C,zg,ld,D,zh,om) result(g)
    integer, intent(in) :: la(3),lb(3),lc(3),ld(3)
    real(dp), intent(in) :: A(3),B(3),C(3),D(3), ze,zf,zg,zh, om
    real(dp) :: g, r1(3), r2(3), h, lo, r12sq, rho1
    integer :: i1,j1,k1,i2,j2,k2, n
    n = 18; lo = -4.0_dp; h = 8.0_dp/n
    g = 0.0_dp
    do i1=0,n; do j1=0,n; do k1=0,n
      r1 = [lo+i1*h, lo+j1*h, lo+k1*h]
      rho1 = gval(la,A,ze,r1)*gval(lb,B,zf,r1)
      if (abs(rho1) < 1e-14_dp) cycle
      do i2=0,n; do j2=0,n; do k2=0,n
        r2 = [lo+i2*h, lo+j2*h, lo+k2*h]
        r12sq = dot_product(r1-r2,r1-r2)
        g = g + rho1*gval(lc,C,zg,r2)*gval(ld,D,zh,r2)*exp(-om*r12sq)
      end do; end do; end do
    end do; end do; end do
    g = g * h**6
  end function num_geminal

end program tc_md_test
