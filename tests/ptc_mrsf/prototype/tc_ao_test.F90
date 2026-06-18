!> Validate the general-L AO integral builder (ptc_ao): for an s-only H2 basis the
!> RHF energy must equal the s-only engine's; adding p functions must lower it.
program tc_ao_test
  use precision, only: dp
  use ptc_ao
  use ptc_geminal, only: ptc_s_ao_1e, ptc_s_ao_tensor, PTC_1E_OVERLAP, &
                         PTC_1E_KINETIC, PTC_1E_NUCLEAR, PTC_OP_ERI
  implicit none
  integer, parameter :: MP = 3
  real(dp) :: R, enuc, e_ao, e_sref, e_sp
  integer  :: nfail
  nfail = 0
  R = 1.4_dp
  enuc = 1.0_dp/R

  write(*,'(a)') '=== general-L AO builder: validation (H2) ==='

  ! (1) s-only 6-311G via the general-L AO builder
  call h2_rhf_ao(.false., R, e_ao)
  ! (2) same via the s-only engine
  call h2_rhf_sengine(R, e_sref)
  write(*,'(a,f14.8)') 'H2/6-311G RHF (general-L AO) = ', e_ao
  write(*,'(a,f14.8)') 'H2/6-311G RHF (sh-only engine)= ', e_sref
  if (abs(e_ao - e_sref) < 1e-8_dp) then
    write(*,'(a)') 'PASS: general-L AO reproduces the sh-only engine'
  else
    write(*,'(a,es10.2)') 'FAIL: dE = ', abs(e_ao-e_sref); nfail=nfail+1
  end if

  ! (3) add a p shell per atom -> energy must drop (polarization)
  call h2_rhf_ao(.true., R, e_sp)
  write(*,'(a,f14.8)') 'H2/6-311G+p RHF (general-L)  = ', e_sp
  if (e_sp < e_ao - 1e-4_dp) then
    write(*,'(a,f8.2,a)') 'PASS: p functions lower RHF by ', (e_ao-e_sp)*1000, ' mHa'
  else
    write(*,'(a)') 'FAIL: p functions did not lower the energy'; nfail=nfail+1
  end if

  if (nfail == 0) then
    write(*,'(a)') 'ALL PASS: general-L AO integral builder validated.'
  else
    error stop 1
  end if

contains

  subroutine h2_basis(addp, nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    logical, intent(in) :: addp
    integer, intent(out) :: nsh, shl_l(:), shl_np(:)
    real(dp), intent(out) :: shl_e(:,:), shl_c(:,:), shl_r(:,:)
    real(dp), intent(in) :: R
    real(dp) :: cen(3,2)
    integer :: at, base
    cen(:,1) = [0.0_dp,0.0_dp,0.0_dp]; cen(:,2) = [0.0_dp,0.0_dp,R]
    nsh = 0
    do at = 1, 2
      ! 6-311G H: contracted 3s + 2 single s
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
        shl_e(1:3,nsh)=[33.8650_dp,5.094790_dp,1.158790_dp]
        shl_c(1:3,nsh)=[0.0254938_dp,0.190373_dp,0.852161_dp]; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.325840_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.102741_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      if (addp) then
        nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.75_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      end if
    end do
    base = 0
  end subroutine h2_basis

  subroutine h2_rhf_ao(addp, R, e)
    logical, intent(in) :: addp
    real(dp), intent(in) :: R
    real(dp), intent(out) :: e
    integer, parameter :: MS = 16
    integer :: nsh, shl_l(MS), shl_np(MS), nao, nat
    real(dp) :: shl_e(MP,MS), shl_c(MP,MS), shl_r(3,MS), zat(2), rat(3,2)
    real(dp), allocatable :: S(:,:), Hc(:,:), eri(:,:,:,:)
    call h2_basis(addp, nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    nat = 2; zat = 1.0_dp; rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
    nao = 0
    block
      integer :: sh
      do sh=1,nsh
        nao = nao + ao_ncart(shl_l(sh))
      end do
    end block
    allocate(S(nao,nao), Hc(nao,nao), eri(nao,nao,nao,nao))
    call build_ints(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                    nat, zat, rat, nao, S, Hc, eri)
    call rhf(nao, 1, S, Hc, eri, 1.0_dp/R, e)
    deallocate(S,Hc,eri)
  end subroutine h2_rhf_ao

  subroutine h2_rhf_sengine(R, e)
    real(dp), intent(in) :: R
    real(dp), intent(out) :: e
    integer, parameter :: NS=6
    integer :: npr(NS), sh
    real(dp) :: exs(MP,NS), cos_(MP,NS), cns(3,NS), zat(2), rat(3,2)
    real(dp) :: S(NS,NS), T(NS,NS), V(NS,NS), Hc(NS,NS), M(NS,NS,NS,NS), eri(NS,NS,NS,NS)
    integer :: i,j,k,l
    do sh=1,NS
      npr(sh)=1
    end do
    do sh=1,NS,3
      npr(sh)=3; exs(:,sh)=[33.8650_dp,5.094790_dp,1.158790_dp]; cos_(:,sh)=[0.0254938_dp,0.190373_dp,0.852161_dp]
    end do
    do sh=2,NS,3
      exs(1,sh)=0.325840_dp; cos_(1,sh)=1.0_dp
    end do
    do sh=3,NS,3
      exs(1,sh)=0.102741_dp; cos_(1,sh)=1.0_dp
    end do
    cns(:,1:3)=spread([0.0_dp,0.0_dp,0.0_dp],2,3); cns(:,4:6)=spread([0.0_dp,0.0_dp,R],2,3)
    zat=1.0_dp; rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
    call ptc_s_ao_1e(NS,npr,exs,cos_,cns,PTC_1E_OVERLAP,2,zat,rat,S)
    call ptc_s_ao_1e(NS,npr,exs,cos_,cns,PTC_1E_KINETIC,2,zat,rat,T)
    call ptc_s_ao_1e(NS,npr,exs,cos_,cns,PTC_1E_NUCLEAR,2,zat,rat,V)
    Hc=T+V
    call ptc_s_ao_tensor(NS,npr,exs,cos_,cns,PTC_OP_ERI,0.0_dp,M)
    do i=1,NS;do j=1,NS;do k=1,NS;do l=1,NS
      eri(i,j,k,l)=M(i,k,j,l)
    end do;end do;end do;end do
    call rhf(NS, 1, S, Hc, eri, 1.0_dp/R, e)
  end subroutine h2_rhf_sengine

  subroutine rhf(n, nocc, S, Hc, eri, enuc, e_rhf)
    integer, intent(in) :: n, nocc
    real(dp), intent(in) :: S(n,n), Hc(n,n), eri(n,n,n,n), enuc
    real(dp), intent(out) :: e_rhf
    real(dp) :: F(n,n), C(n,n), P(n,n), G(n,n), Scp(n,n), eps(n), eold
    integer :: it, mu, nu, la, si, i
    F = Hc; e_rhf = 0.0_dp
    do it=1,300
      Scp=S; call geig(F,Scp,n,eps,C)
      P=0.0_dp
      do mu=1,n;do nu=1,n;do i=1,nocc
        P(mu,nu)=P(mu,nu)+2.0_dp*C(mu,i)*C(nu,i)
      end do;end do;end do
      G=0.0_dp
      do mu=1,n;do nu=1,n;do la=1,n;do si=1,n
        G(mu,nu)=G(mu,nu)+P(la,si)*(eri(mu,nu,la,si)-0.5_dp*eri(mu,la,nu,si))
      end do;end do;end do;end do
      F=Hc+G
      eold=e_rhf; e_rhf=0.0_dp
      do mu=1,n;do nu=1,n
        e_rhf=e_rhf+0.5_dp*P(mu,nu)*(Hc(mu,nu)+F(mu,nu))
      end do;end do
      if (abs(e_rhf-eold)<1e-11_dp .and. it>1) exit
    end do
    e_rhf=e_rhf+enuc
  end subroutine rhf

  subroutine geig(F,S,n,w,C)
    integer, intent(in) :: n
    real(dp), intent(inout) :: F(n,n), S(n,n)
    real(dp), intent(out) :: w(n), C(n,n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsygv(1,'V','U',n,F,n,S,n,w,wq,-1,info)
    lw=int(wq(1)); allocate(wk(lw))
    call dsygv(1,'V','U',n,F,n,S,n,w,wk,lw,info)
    C=F; deallocate(wk)
  end subroutine geig

end program tc_ao_test
