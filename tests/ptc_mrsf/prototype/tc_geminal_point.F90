!> Single-geometry validation of the genuine-geminal pTC-MRSF-CIS on H2/6-311G
!> at a stretched bond (multireference). Checks: ROHF-triplet reference; tau=0
!> gate (scale=0 -> bare); the geminal LOWERS S0 (recovers dynamic correlation);
!> and -- the central point -- the geminal also DRESSES THE TRIPLET (pTC != bare
!> for T1), unlike the singlet-only MP2-T2 proxy.
program tc_geminal_point
  use precision, only: dp
  use ptc_geminal
  use tc_geminal_engine
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none

  integer, parameter :: NS = 6, MP = 3
  integer  :: npr(NS), i, j, k, l
  real(dp) :: exs(MP,NS), cos_(MP,NS), cns(3,NS), rat(3,2)
  real(dp) :: Cmo(NS,NS), eps(NS), e_scf, enuc, h1ao(NS,NS), eri_c(NS,NS,NS,NS)
  real(dp) :: h1mo(NS,NS), eri_mo(NS,NS,NS,NS), Gmo(NS,NS,NS,NS)
  integer, allocatable :: dets(:), cas(:)
  integer  :: dim, hfidx, nc, iact(2), iext(4)
  real(dp), allocatable :: H(:,:), T2op(:,:), Hbar(:,:), Em(:,:), Ep(:,:), S2(:,:)
  real(dp), allocatable :: Hc(:,:), Hbc(:,:), S2c(:,:)
  real(dp) :: R, gamma, efci(3), ebare(3), eptc(3), eptc0(3)
  integer  :: nfail

  nfail = 0
  R = 2.4_dp        ! stretched H2 (multireference)
  gamma = 1.0_dp    ! Slater geminal exponent
  call set_h_basis(npr, exs, cos_, cns, R)
  rat(:,1) = [0.0_dp,0.0_dp,0.0_dp]; rat(:,2) = [0.0_dp,0.0_dp,R]

  ! ROHF triplet reference (2 alpha, 0 beta -> sigma_g sigma_u)
  call rohf_highspin(NS, 2, 0, npr, exs, cos_, cns, 2, rat, Cmo, eps, e_scf, enuc, h1ao, eri_c)
  write(*,'(a,f12.6)') 'ROHF triplet (3Su+) energy   = ', e_scf
  if (e_scf > -0.6_dp .or. e_scf < -1.2_dp) nfail = nfail + 1   ! sane range

  call ao2mo_1e(h1ao, Cmo, NS, h1mo)
  call ao2mo_2e(eri_c, Cmo, NS, eri_mo)
  call geminal_mo(NS, npr, exs, cos_, cns, gamma, Cmo, Gmo)

  ! Ms=0 target determinants (1 alpha, 1 beta) in the triplet orbitals
  call build_dets(NS, 1, 1, dets, dim, hfidx)
  allocate(H(dim,dim), T2op(dim,dim), Hbar(dim,dim), Em(dim,dim), Ep(dim,dim), S2(dim,dim))
  call build_fci_H(h1mo, eri_mo, enuc, NS, dets, dim, H)
  call build_s2(NS, dets, dim, S2)

  iact = [1, 2]                 ! frontier sigma_g, sigma_u (MRSF active)
  iext = [3, 4, 5, 6]           ! external (dynamic correlation)
  call cas22_compact(dets, dim, NS, iact, 2, cas, nc)
  allocate(Hc(nc,nc), Hbc(nc,nc), S2c(nc,nc))

  ! FCI and bare MRSF (ROHF, (2,2))
  call states3(H, S2, dim, .false., efci)
  do i=1,nc; do j=1,nc
    Hc(i,j) = H(cas(i),cas(j)); S2c(i,j) = S2(cas(i),cas(j))
  end do; end do
  call states3(Hc, S2c, nc, .false., ebare)

  ! genuine geminal pTC (scale=1)
  call build_geminal_T2(Gmo, NS, iact, 2, iext, 4, 1.0_dp, gamma, dets, dim, T2op)
  call expm_nilpotent(-1.0_dp, T2op, dim, Em)
  call expm_nilpotent( 1.0_dp, T2op, dim, Ep)
  Hbar = matmul(Em, matmul(H, Ep))
  do i=1,nc; do j=1,nc
    Hbc(i,j) = Hbar(cas(i),cas(j))
  end do; end do
  call states3(Hbc, S2c, nc, .true., eptc)

  ! tau=0 gate (scale=0 -> H_bar = H)
  call build_geminal_T2(Gmo, NS, iact, 2, iext, 4, 0.0_dp, gamma, dets, dim, T2op)
  call expm_nilpotent(-1.0_dp, T2op, dim, Em)
  call expm_nilpotent( 1.0_dp, T2op, dim, Ep)
  Hbar = matmul(Em, matmul(H, Ep))
  do i=1,nc; do j=1,nc
    Hbc(i,j) = Hbar(cas(i),cas(j))
  end do; end do
  call states3(Hbc, S2c, nc, .true., eptc0)

  write(*,'(a)') ''
  write(*,'(a)') '             S0           T1           S1'
  write(*,'(a,3f13.6)') 'FCI       ', efci
  write(*,'(a,3f13.6)') 'bare MRSF ', ebare
  write(*,'(a,3f13.6)') 'pTC(gem)  ', eptc
  write(*,'(a,3f13.6)') 'pTC(t=0)  ', eptc0
  write(*,'(a)') ''
  call chk('tau=0 gate S0 == bare        ', eptc0(1), ebare(1), 1.0e-7_dp, nfail)
  call chk('tau=0 gate T1 == bare        ', eptc0(2), ebare(2), 1.0e-7_dp, nfail)
  ! genuine geminal must LOWER S0 (recover dynamic correlation toward FCI)
  if (.not. (eptc(1) < ebare(1) - 1.0e-4_dp)) then
    write(*,'(a)') 'FAIL: geminal did not lower S0'; nfail = nfail + 1
  else
    write(*,'(a,f6.1,a)') 'PASS: geminal lowers S0 by ', (ebare(1)-eptc(1))*1000, ' mHa'
  end if
  ! THE POINT: the geminal also dresses the TRIPLET (unlike MP2-T2)
  if (abs(eptc(2) - ebare(2)) < 1.0e-4_dp) then
    write(*,'(a)') 'NOTE: triplet NOT dressed by the geminal (check spin channel)'
    nfail = nfail + 1
  else
    write(*,'(a,f6.1,a)') 'PASS: geminal DRESSES the triplet T1 by ', &
      (ebare(2)-eptc(2))*1000, ' mHa (the spin-resolved cusp at work)'
  end if

  write(*,'(a)') ''
  if (nfail == 0) then
    write(*,'(a)') 'ALL PASS: genuine geminal pTC-MRSF-CIS (ROHF ref) validated at one point.'
  else
    write(*,'(a,i0,a)') 'FAILURES: ', nfail, ' check(s).'
  end if

contains

  subroutine set_h_basis(npr, exs, cos_, cns, R)
    integer,  intent(out) :: npr(NS)
    real(dp), intent(out) :: exs(MP,NS), cos_(MP,NS), cns(3,NS)
    real(dp), intent(in)  :: R
    integer :: s
    do s = 1, NS
      npr(s) = 1
    end do
    do s = 1, NS, 3
      npr(s) = 3
      exs(:,s)  = [33.8650_dp, 5.094790_dp, 1.158790_dp]
      cos_(:,s) = [0.0254938_dp, 0.190373_dp, 0.852161_dp]
    end do
    do s = 2, NS, 3
      exs(1,s) = 0.325840_dp; cos_(1,s) = 1.0_dp
    end do
    do s = 3, NS, 3
      exs(1,s) = 0.102741_dp; cos_(1,s) = 1.0_dp
    end do
    cns(:,1:3) = spread([0.0_dp,0.0_dp,0.0_dp], 2, 3)
    cns(:,4:6) = spread([0.0_dp,0.0_dp,R], 2, 3)
  end subroutine set_h_basis

  subroutine states3(Hm, S2m, m, nonherm, out)
    integer, intent(in) :: m
    real(dp), intent(in) :: Hm(m,m), S2m(m,m)
    logical, intent(in) :: nonherm
    real(dp), intent(out) :: out(3)
    real(dp) :: w(m), Vv(m,m), vrr(m,m), vll(m,m), ss(m), mi, den
    integer :: kk, ncx, ierr, ns0, nt0, ns1, ord(m), i2, j2, tmp
    if (nonherm) then
      call tc_nonsym_tda_eig(Hm, m, w, vrr, vll, mi, ncx, ierr)
      do kk=1,m
        den = dot_product(vll(:,kk), vrr(:,kk))
        ss(kk) = dot_product(vll(:,kk), matmul(S2m, vrr(:,kk)))/den
      end do
    else
      Vv = Hm
      call sym_eig_vec(Vv, m, w)
      do kk=1,m
        ss(kk) = dot_product(Vv(:,kk), matmul(S2m, Vv(:,kk)))
      end do
    end if
    do i2=1,m
      ord(i2)=i2
    end do
    do i2=1,m-1; do j2=i2+1,m
      if (w(ord(j2))<w(ord(i2))) then
        tmp=ord(i2); ord(i2)=ord(j2); ord(j2)=tmp
      end if
    end do; end do
    out = huge(1.0_dp); ns0=0; nt0=0; ns1=0
    do i2=1,m
      kk=ord(i2)
      if (ss(kk) < 1.0_dp) then
        if (ns0==0) then
          out(1)=w(kk); ns0=1
        else if (ns1==0) then
          out(3)=w(kk); ns1=1
        end if
      else
        if (nt0==0) then
          out(2)=w(kk); nt0=1
        end if
      end if
    end do
  end subroutine states3

  subroutine chk(name, a, b, tol, nf)
    character(*), intent(in) :: name
    real(dp), intent(in) :: a, b, tol
    integer, intent(inout) :: nf
    real(dp) :: d
    d = abs(a-b)
    write(*,'(a,2f14.7,es10.2,a)') name, a, b, d, merge('  PASS','  FAIL', d<=tol)
    if (d>tol) nf = nf + 1
  end subroutine chk

end program tc_geminal_point
