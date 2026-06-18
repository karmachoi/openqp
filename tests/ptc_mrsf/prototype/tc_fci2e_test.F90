!> Validate the fast 2-electron FCI engine (ptc_fci2e) against the general
!> determinant engine (tc_geminal_engine) on H2/6-311G+p: same MO integrals must
!> give the same FCI spectrum, and the geminal pTC the same downfolded states.
program tc_fci2e_test
  use precision, only: dp
  use ptc_ao
  use ptc_fci2e
  use tc_geminal_engine, only: ao2mo_1e, ao2mo_2e, build_dets, build_fci_H, &
       build_s2, cas22_compact, build_geminal_T2, expm_nilpotent, sym_eig_vec
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none
  integer, parameter :: MS=16, MP=3
  integer  :: nsh, shl_l(MS), shl_np(MS), nao, nat, i, j, nfail
  real(dp) :: shl_e(MP,MS), shl_c(MP,MS), shl_r(3,MS), zat(2), rat(3,2), R, gamma, enuc
  real(dp), allocatable :: S(:,:), Hc(:,:), eri(:,:,:,:), Gao(:,:,:,:)
  real(dp), allocatable :: Cmo(:,:), eps(:), h1mo(:,:), eri_mo(:,:,:,:), Gmo(:,:,:,:)
  real(dp) :: e_gen, e_fast, s0_gen, s0_fast

  nfail=0; R=1.4_dp; gamma=1.0_dp; enuc=1.0_dp/R
  call h2_basis(nsh,shl_l,shl_np,shl_e,shl_c,shl_r,R)
  nat=2; zat=1.0_dp; rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
  nao=0
  do i=1,nsh
    nao=nao+ao_ncart(shl_l(i))
  end do
  allocate(S(nao,nao),Hc(nao,nao),eri(nao,nao,nao,nao),Gao(nao,nao,nao,nao))
  call build_ints(nsh,shl_l(1:nsh),shl_np(1:nsh),shl_e,shl_c,shl_r(:,1:nsh),nat,zat,rat,nao,S,Hc,eri)
  call build_geminal_ao(nsh,shl_l(1:nsh),shl_np(1:nsh),shl_e,shl_c,shl_r(:,1:nsh),gamma,nao,S,Gao)
  allocate(Cmo(nao,nao),eps(nao),h1mo(nao,nao),eri_mo(nao,nao,nao,nao),Gmo(nao,nao,nao,nao))
  call rohf_triplet(nao,S,Hc,eri,Cmo,eps)
  call ao2mo_1e(Hc,Cmo,nao,h1mo)
  call ao2mo_2e(eri,Cmo,nao,eri_mo)
  call ao2mo_2e(Gao,Cmo,nao,Gmo)

  write(*,'(a,i0,a)') '=== fast vs general 2e-FCI (H2/6-311G+p, ',nao,' MOs) ==='

  ! ---- FCI ground state: general vs fast ----
  block
    integer, allocatable :: dets(:)
    integer :: dim, hfidx
    real(dp), allocatable :: Hg(:,:), wg(:)
    real(dp) :: H2(nao*nao,nao*nao), w2(nao*nao)
    call build_dets(nao,1,1,dets,dim,hfidx)
    allocate(Hg(dim,dim),wg(dim)); call build_fci_H(h1mo,eri_mo,enuc,nao,dets,dim,Hg)
    call sym_eig_vec(Hg,dim,wg); e_gen=wg(1)
    call build_H2e(h1mo,eri_mo,enuc,nao,H2); call sym_eig_vec(H2,nao*nao,w2); e_fast=w2(1)
    deallocate(Hg,wg,dets)
  end block
  call chk('FCI ground E (gen vs fast)', e_gen, e_fast, 1e-9_dp, nfail)

  ! ---- pTC S0 (lowest singlet): general vs fast ----
  s0_gen  = ptc_s0_general(nao, h1mo, eri_mo, Gmo, enuc, gamma)
  s0_fast = ptc_s0_fast(nao, h1mo, eri_mo, Gmo, enuc, gamma)
  call chk('pTC S0 (gen vs fast)      ', s0_gen, s0_fast, 1e-8_dp, nfail)

  if (nfail==0) then
    write(*,'(a)') 'ALL PASS: fast 2e-FCI engine matches the general engine.'
  else
    error stop 1
  end if

contains

  real(dp) function ptc_s0_general(n, h1, e2, G, enc, gam) result(s0)
    integer, intent(in) :: n
    real(dp), intent(in) :: h1(n,n), e2(n,n,n,n), G(n,n,n,n), enc, gam
    integer, allocatable :: dets(:), cas(:)
    integer :: dim,hfidx,nc,iact(2),iext(MS),nextn,i2,j2
    real(dp), allocatable :: H(:,:),T2(:,:),Hb(:,:),Em(:,:),Ep(:,:),S2(:,:),Hbc(:,:),S2c(:,:)
    real(dp) :: w(4),vr(4,4),vl(4,4),mi,ss(4),den
    integer :: ncx,ierr,k
    call build_dets(n,1,1,dets,dim,hfidx)
    allocate(H(dim,dim),T2(dim,dim),Hb(dim,dim),Em(dim,dim),Ep(dim,dim),S2(dim,dim))
    call build_fci_H(h1,e2,enc,n,dets,dim,H); call build_s2(n,dets,dim,S2)
    iact=[1,2]; nextn=n-2
    do i2=1,nextn
      iext(i2)=i2+2
    end do
    call cas22_compact(dets,dim,n,iact,2,cas,nc)
    allocate(Hbc(nc,nc),S2c(nc,nc))
    call build_geminal_T2(G,n,iact,2,iext(1:nextn),nextn,1.0_dp,gam,dets,dim,T2)
    call expm_nilpotent(-1.0_dp,T2,dim,Em); call expm_nilpotent(1.0_dp,T2,dim,Ep)
    Hb=matmul(Em,matmul(H,Ep))
    do i2=1,nc; do j2=1,nc
      Hbc(i2,j2)=Hb(cas(i2),cas(j2)); S2c(i2,j2)=S2(cas(i2),cas(j2))
    end do; end do
    call tc_nonsym_tda_eig(Hbc,nc,w,vr,vl,mi,ncx,ierr)
    do k=1,nc
      den=dot_product(vl(:,k),vr(:,k)); ss(k)=dot_product(vl(:,k),matmul(S2c,vr(:,k)))/den
    end do
    s0=huge(1.0_dp)
    do k=1,nc
      if (ss(k)<1.0_dp .and. w(k)<s0) s0=w(k)
    end do
  end function ptc_s0_general

  real(dp) function ptc_s0_fast(n, h1, e2, G, enc, gam) result(s0)
    integer, intent(in) :: n
    real(dp), intent(in) :: h1(n,n), e2(n,n,n,n), G(n,n,n,n), enc, gam
    integer, allocatable :: cas(:)
    integer :: dim,nc,iact(2),iext(MS),nextn,i2,j2,k,ncx,ierr
    real(dp), allocatable :: H(:,:),T2(:,:),Hb(:,:),Em(:,:),Ep(:,:),S2(:,:),Hbc(:,:),S2c(:,:)
    real(dp) :: w(4),vr(4,4),vl(4,4),mi,ss(4),den
    dim=n*n
    allocate(H(dim,dim),T2(dim,dim),Hb(dim,dim),Em(dim,dim),Ep(dim,dim),S2(dim,dim))
    call build_H2e(h1,e2,enc,n,H); call build_S2_2e(n,S2)
    iact=[1,2]; nextn=n-2
    do i2=1,nextn
      iext(i2)=i2+2
    end do
    call cas22_2e(n,cas,nc)
    allocate(Hbc(nc,nc),S2c(nc,nc))
    call build_geminal_T2_2e(G,n,iact,2,iext(1:nextn),nextn,gam,T2)
    call expm_nilpotent(-1.0_dp,T2,dim,Em); call expm_nilpotent(1.0_dp,T2,dim,Ep)
    Hb=matmul(Em,matmul(H,Ep))
    do i2=1,nc; do j2=1,nc
      Hbc(i2,j2)=Hb(cas(i2),cas(j2)); S2c(i2,j2)=S2(cas(i2),cas(j2))
    end do; end do
    call tc_nonsym_tda_eig(Hbc,nc,w,vr,vl,mi,ncx,ierr)
    do k=1,nc
      den=dot_product(vl(:,k),vr(:,k)); ss(k)=dot_product(vl(:,k),matmul(S2c,vr(:,k)))/den
    end do
    s0=huge(1.0_dp)
    do k=1,nc
      if (ss(k)<1.0_dp .and. w(k)<s0) s0=w(k)
    end do
  end function ptc_s0_fast

  subroutine chk(name,x,y,tol,nf)
    character(*), intent(in) :: name
    real(dp), intent(in) :: x,y,tol
    integer, intent(inout) :: nf
    write(*,'(a,2f16.8,es10.1,a)') name,x,y,abs(x-y),merge('  PASS','  FAIL',abs(x-y)<=tol)
    if (abs(x-y)>tol) nf=nf+1
  end subroutine chk

  subroutine h2_basis(nsh,shl_l,shl_np,shl_e,shl_c,shl_r,R)
    integer, intent(out) :: nsh,shl_l(:),shl_np(:)
    real(dp), intent(out) :: shl_e(:,:),shl_c(:,:),shl_r(:,:)
    real(dp), intent(in) :: R
    real(dp) :: cen(3,2)
    integer :: at
    cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.0_dp,0.0_dp,R]; nsh=0
    do at=1,2
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
        shl_e(1:3,nsh)=[33.8650_dp,5.094790_dp,1.158790_dp]
        shl_c(1:3,nsh)=[0.0254938_dp,0.190373_dp,0.852161_dp]; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.325840_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.102741_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.75_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
    end do
  end subroutine h2_basis

  subroutine rohf_triplet(n,S,Hc,eri,Cmo,eps)
    integer, intent(in) :: n
    real(dp), intent(in) :: S(n,n),Hc(n,n),eri(n,n,n,n)
    real(dp), intent(out) :: Cmo(n,n),eps(n)
    real(dp) :: Fa(n,n),Ca(n,n),Da(n,n),Ja(n,n),Ka(n,n),Scp(n,n),eold,e
    integer :: it,mu,nu,la,si,i
    Fa=Hc; e=0.0_dp
    do it=1,400
      Scp=S; call geig(Fa,Scp,n,eps,Ca)
      Da=0.0_dp
      do mu=1,n;do nu=1,n;do i=1,2
        Da(mu,nu)=Da(mu,nu)+Ca(mu,i)*Ca(nu,i)
      end do;end do;end do
      Ja=0.0_dp; Ka=0.0_dp
      do mu=1,n;do nu=1,n;do la=1,n;do si=1,n
        Ja(mu,nu)=Ja(mu,nu)+Da(la,si)*eri(mu,nu,la,si); Ka(mu,nu)=Ka(mu,nu)+Da(la,si)*eri(mu,la,nu,si)
      end do;end do;end do;end do
      Fa=Hc+Ja-Ka; eold=e; e=0.0_dp
      do mu=1,n;do nu=1,n
        e=e+Da(mu,nu)*(Hc(mu,nu)+Fa(mu,nu))
      end do;end do
      if (abs(e-eold)<1e-11_dp .and. it>1) exit
    end do
    Cmo=Ca
  end subroutine rohf_triplet

  subroutine geig(F,S,n,w,C)
    integer, intent(in) :: n
    real(dp), intent(inout) :: F(n,n),S(n,n)
    real(dp), intent(out) :: w(n),C(n,n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info,lw
    call dsygv(1,'V','U',n,F,n,S,n,w,wq,-1,info); lw=int(wq(1)); allocate(wk(lw))
    call dsygv(1,'V','U',n,F,n,S,n,w,wk,lw,info); C=F; deallocate(wk)
  end subroutine geig

end program tc_fci2e_test
