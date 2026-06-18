!> How much of the pTC-FCI gap is just an un-optimized geminal exponent gamma?
!> Scan gamma for H2/6-311G at two bond lengths; report S0: bare, pTC(gamma), FCI,
!> and the residual pTC-FCI gap (mHa). gamma is the only knob (fixed amplitudes).
program tc_gamma_scan
  use precision, only: dp
  use ptc_geminal
  use tc_geminal_engine
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none
  integer, parameter :: NS = 6, MP = 3
  integer  :: npr(NS), ig
  real(dp) :: exs(MP,NS), cos_(MP,NS), cns(3,NS), rat(3,2)
  real(dp) :: Rlist(2), gam, s0b, s0p, s0f, R
  integer  :: ir
  Rlist = [1.5_dp, 2.5_dp]
  do ir = 1, 2
    R = Rlist(ir)
    write(*,'(a,f5.2,a)') '=== H2/6-311G  R = ', R, ' bohr ==='
    write(*,'(a)') ' gamma   bare S0     pTC S0      FCI S0    pTC-FCI(mHa)  recov%'
    do ig = 1, 11
      gam = 0.4_dp + 0.2_dp*real(ig-1,dp)     ! 0.4 .. 2.4
      call s0_at(R, gam, s0b, s0p, s0f)
      write(*,'(f6.2, 3f12.6, f12.2, f9.0)') gam, s0b, s0p, s0f, &
        (s0p-s0f)*1000.0_dp, (s0b-s0p)/(s0b-s0f)*100.0_dp
    end do
    write(*,'(a)') ''
  end do

contains

  subroutine s0_at(R, gamma, e_bare, e_ptc, e_fci)
    real(dp), intent(in)  :: R, gamma
    real(dp), intent(out) :: e_bare, e_ptc, e_fci
    real(dp) :: Cmo(NS,NS), eps(NS), e_scf, enuc, h1ao(NS,NS), eri_c(NS,NS,NS,NS)
    real(dp) :: h1mo(NS,NS), eri_mo(NS,NS,NS,NS), Gmo(NS,NS,NS,NS)
    integer, allocatable :: dets(:), cas(:)
    integer  :: dim, hfidx, nc, iact(2), iext(4), i, j, s
    real(dp), allocatable :: H(:,:), T2op(:,:), Hbar(:,:), Em(:,:), Ep(:,:), S2(:,:)
    real(dp), allocatable :: Hc(:,:), Hbc(:,:), S2c(:,:), w(:), vr(:,:), vl(:,:)
    real(dp) :: mi, eb(3), wf(NS*NS)
    integer  :: ncx, ierr
    do s=1,NS
      npr(s)=1
    end do
    do s=1,NS,3
      npr(s)=3; exs(:,s)=[33.8650_dp,5.094790_dp,1.158790_dp]
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
    rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
    call rohf_highspin(NS,2,0,npr,exs,cos_,cns,2,rat,Cmo,eps,e_scf,enuc,h1ao,eri_c)
    call ao2mo_1e(h1ao,Cmo,NS,h1mo)
    call ao2mo_2e(eri_c,Cmo,NS,eri_mo)
    call geminal_mo(NS,npr,exs,cos_,cns,gamma,Cmo,Gmo)
    call build_dets(NS,1,1,dets,dim,hfidx)
    allocate(H(dim,dim),T2op(dim,dim),Hbar(dim,dim),Em(dim,dim),Ep(dim,dim),S2(dim,dim))
    call build_fci_H(h1mo,eri_mo,enuc,NS,dets,dim,H)
    call build_s2(NS,dets,dim,S2)
    iact=[1,2]; iext=[3,4,5,6]
    call cas22_compact(dets,dim,NS,iact,2,cas,nc)
    allocate(Hc(nc,nc),Hbc(nc,nc),S2c(nc,nc),w(nc),vr(nc,nc),vl(nc,nc))
    ! FCI S0
    block
      real(dp) :: Hcp(dim,dim), wful(dim)
      Hcp=H; call sym_eig_vec(Hcp,dim,wful); e_fci=wful(1)
    end block
    do i=1,nc; do j=1,nc
      Hc(i,j)=H(cas(i),cas(j)); S2c(i,j)=S2(cas(i),cas(j))
    end do; end do
    block
      real(dp) :: Hcc(nc,nc), wb(nc)
      Hcc=Hc; call sym_eig_vec(Hcc,nc,wb); e_bare=wb(1)
    end block
    call build_geminal_T2(Gmo,NS,iact,2,iext,4,1.0_dp,gamma,dets,dim,T2op)
    call expm_nilpotent(-1.0_dp,T2op,dim,Em)
    call expm_nilpotent( 1.0_dp,T2op,dim,Ep)
    Hbar=matmul(Em,matmul(H,Ep))
    do i=1,nc; do j=1,nc
      Hbc(i,j)=Hbar(cas(i),cas(j))
    end do; end do
    call tc_nonsym_tda_eig(Hbc,nc,w,vr,vl,mi,ncx,ierr)
    e_ptc=minval(w)
    deallocate(dets,cas,H,T2op,Hbar,Em,Ep,S2,Hc,Hbc,S2c,w,vr,vl)
  end subroutine s0_at

end program tc_gamma_scan
