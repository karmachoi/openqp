!> Native genuine-geminal pTC-MRSF-CIS H2 dissociation curves (pyscf-free):
!> ROHF-triplet reference, the cusp-fixed Slater geminal F12 correlator (spin-
!> resolved 1/2 opposite-spin, 1/4 same-spin), downfolded into the (2,2) frontier
!> MRSF space, vs bare MRSF-CIS (ROHF) and FCI, in 6-311G. Unlike the MP2-T2
!> proxy, the geminal dresses BOTH the singlet and the triplet states.
!> Writes pes_h2_geminal.dat: R  S0/T1/S1 (fci | bare | pTC-geminal).
!>
!> Build (OpenMP):
!>   gfortran -O2 -fopenmp source/precision.F90 source/modules/ptc_geminal.F90 \
!>     source/modules/tdhf_mrsf_ptc.F90 tests/ptc_mrsf/prototype/tc_geminal_engine.F90 \
!>     tests/ptc_mrsf/prototype/tc_h2_pes_geminal.F90 -llapack -lblas -o /tmp/pesg && /tmp/pesg
program tc_h2_pes_geminal
  use precision, only: dp
  use ptc_geminal
  use tc_geminal_engine
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none
  integer, parameter :: NS = 6, MP = 3
  integer  :: npr(NS), ir, nR, u
  real(dp) :: exs(MP,NS), cos_(MP,NS), cns(3,NS), rat(3,2)
  real(dp) :: R, gamma, efci(3), ebare(3), eptc(3)
  real(dp), allocatable :: Rs(:)
  gamma = 1.0_dp
  nR = 28
  allocate(Rs(nR))
  do ir = 1, nR
    Rs(ir) = 1.0_dp + 7.0_dp*real(ir-1,dp)/real(nR-1,dp)
  end do
  open(newunit=u, file='pes_h2_geminal.dat', status='replace')
  write(u,'(a)') '# R   S0_fci T1_fci S1_fci   S0_bare T1_bare S1_bare   S0_ptc T1_ptc S1_ptc'
  do ir = 1, nR
    R = Rs(ir)
    call one_point(R, gamma, efci, ebare, eptc)
    write(u,'(f8.3, 3f12.6, 3x, 3f12.6, 3x, 3f12.6)') R, efci, ebare, eptc
    write(*,'(a,f6.3,a,3f10.5,a,3f10.5,a,3f10.5)') 'R=',R, &
      ' FCI=', efci, ' bare=', ebare, ' pTCgem=', eptc
  end do
  close(u)
  write(*,'(a)') 'wrote pes_h2_geminal.dat'

contains

  subroutine one_point(R, gamma, efci, ebare, eptc)
    real(dp), intent(in)  :: R, gamma
    real(dp), intent(out) :: efci(3), ebare(3), eptc(3)
    real(dp) :: Cmo(NS,NS), eps(NS), e_scf, enuc, h1ao(NS,NS), eri_c(NS,NS,NS,NS)
    real(dp) :: h1mo(NS,NS), eri_mo(NS,NS,NS,NS), Gmo(NS,NS,NS,NS)
    integer, allocatable :: dets(:), cas(:)
    integer  :: dim, hfidx, nc, iact(2), iext(4), i, j
    real(dp), allocatable :: H(:,:), T2op(:,:), Hbar(:,:), Em(:,:), Ep(:,:), S2(:,:)
    real(dp), allocatable :: Hc(:,:), Hbc(:,:), S2c(:,:)
    call set_h_basis(npr, exs, cos_, cns, R)
    rat(:,1) = [0.0_dp,0.0_dp,0.0_dp]; rat(:,2) = [0.0_dp,0.0_dp,R]
    call rohf_highspin(NS, 2, 0, npr, exs, cos_, cns, 2, rat, Cmo, eps, e_scf, enuc, h1ao, eri_c)
    call ao2mo_1e(h1ao, Cmo, NS, h1mo)
    call ao2mo_2e(eri_c, Cmo, NS, eri_mo)
    call geminal_mo(NS, npr, exs, cos_, cns, gamma, Cmo, Gmo)
    call build_dets(NS, 1, 1, dets, dim, hfidx)
    allocate(H(dim,dim), T2op(dim,dim), Hbar(dim,dim), Em(dim,dim), Ep(dim,dim), S2(dim,dim))
    call build_fci_H(h1mo, eri_mo, enuc, NS, dets, dim, H)
    call build_s2(NS, dets, dim, S2)
    iact = [1,2]; iext = [3,4,5,6]
    call cas22_compact(dets, dim, NS, iact, 2, cas, nc)
    allocate(Hc(nc,nc), Hbc(nc,nc), S2c(nc,nc))
    call states3(H, S2, dim, .false., efci)
    do i=1,nc; do j=1,nc
      Hc(i,j)=H(cas(i),cas(j)); S2c(i,j)=S2(cas(i),cas(j))
    end do; end do
    call states3(Hc, S2c, nc, .false., ebare)
    call build_geminal_T2(Gmo, NS, iact, 2, iext, 4, 1.0_dp, gamma, dets, dim, T2op)
    call expm_nilpotent(-1.0_dp, T2op, dim, Em)
    call expm_nilpotent( 1.0_dp, T2op, dim, Ep)
    Hbar = matmul(Em, matmul(H, Ep))
    do i=1,nc; do j=1,nc
      Hbc(i,j)=Hbar(cas(i),cas(j))
    end do; end do
    call states3(Hbc, S2c, nc, .true., eptc)
    deallocate(dets,cas,H,T2op,Hbar,Em,Ep,S2,Hc,Hbc,S2c)
  end subroutine one_point

  subroutine set_h_basis(npr, exs, cos_, cns, R)
    integer,  intent(out) :: npr(NS)
    real(dp), intent(out) :: exs(MP,NS), cos_(MP,NS), cns(3,NS)
    real(dp), intent(in)  :: R
    integer :: s
    do s=1,NS
      npr(s)=1
    end do
    do s=1,NS,3
      npr(s)=3
      exs(:,s)=[33.8650_dp,5.094790_dp,1.158790_dp]
      cos_(:,s)=[0.0254938_dp,0.190373_dp,0.852161_dp]
    end do
    do s=2,NS,3
      exs(1,s)=0.325840_dp; cos_(1,s)=1.0_dp
    end do
    do s=3,NS,3
      exs(1,s)=0.102741_dp; cos_(1,s)=1.0_dp
    end do
    cns(:,1:3)=spread([0.0_dp,0.0_dp,0.0_dp],2,3)
    cns(:,4:6)=spread([0.0_dp,0.0_dp,R],2,3)
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
        den=dot_product(vll(:,kk),vrr(:,kk))
        ss(kk)=dot_product(vll(:,kk),matmul(S2m,vrr(:,kk)))/den
      end do
    else
      Vv=Hm; call sym_eig_vec(Vv,m,w)
      do kk=1,m
        ss(kk)=dot_product(Vv(:,kk),matmul(S2m,Vv(:,kk)))
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
    out=huge(1.0_dp); ns0=0; nt0=0; ns1=0
    do i2=1,m
      kk=ord(i2)
      if (ss(kk)<1.0_dp) then
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

end program tc_h2_pes_geminal
