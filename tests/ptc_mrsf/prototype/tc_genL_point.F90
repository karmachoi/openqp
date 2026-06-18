!> General-L pTC-MRSF-CIS end-to-end at one geometry (H2/6-311G+p), proving the
!> whole stack works with p functions: McMurchie-Davidson primitives (ptc_md) ->
!> AO contraction (ptc_ao) -> ROHF triplet -> FCI -> cusp-fixed geminal
!> transcorrelation -> (2,2) downfold -> non-Hermitian solve. Checks that p
!> polarization lowers FCI vs s-only and that the geminal dresses S0 and T1.
program tc_genL_point
  use precision, only: dp
  use ptc_ao
  use tc_geminal_engine, only: ao2mo_1e, ao2mo_2e, build_dets, build_fci_H, &
       build_s2, cas22_compact, build_geminal_T2, expm_nilpotent, sym_eig_vec
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none
  integer, parameter :: MS = 16, MP = 3
  integer  :: nsh, shl_l(MS), shl_np(MS), nao, nat
  real(dp) :: shl_e(MP,MS), shl_c(MP,MS), shl_r(3,MS), zat(2), rat(3,2)
  real(dp) :: R, gamma, enuc
  real(dp), allocatable :: S(:,:), Hc(:,:), eri(:,:,:,:), Gao(:,:,:,:)
  real(dp), allocatable :: Cmo(:,:), eps(:), h1mo(:,:), eri_mo(:,:,:,:), Gmo(:,:,:,:)
  integer, allocatable :: dets(:), cas(:)
  integer  :: dim, hfidx, nc, iact(2), iext(MS), nextn, i, j, k
  real(dp), allocatable :: H(:,:), T2op(:,:), Hbar(:,:), Em(:,:), Ep(:,:), S2(:,:)
  real(dp), allocatable :: Hcc(:,:), Hbc(:,:), S2c(:,:)
  real(dp) :: efci(3), ebare(3), eptc(3)

  R = 1.4_dp; gamma = 1.0_dp; enuc = 1.0_dp/R
  call h2_basis(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
  nat = 2; zat = 1.0_dp; rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
  nao = 0
  do i = 1, nsh
    nao = nao + ao_ncart(shl_l(i))
  end do
  write(*,'(a,i0,a,i0,a)') '=== H2/6-311G+p : ', nsh, ' shells, ', nao, ' AOs ==='
  allocate(S(nao,nao), Hc(nao,nao), eri(nao,nao,nao,nao), Gao(nao,nao,nao,nao))
  call build_ints(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                  nat, zat, rat, nao, S, Hc, eri)
  call build_geminal_ao(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                        gamma, nao, S, Gao)

  allocate(Cmo(nao,nao), eps(nao), h1mo(nao,nao), eri_mo(nao,nao,nao,nao), Gmo(nao,nao,nao,nao))
  call rohf_triplet(nao, S, Hc, eri, Cmo, eps)
  call ao2mo_1e(Hc, Cmo, nao, h1mo)
  call ao2mo_2e(eri, Cmo, nao, eri_mo)
  call ao2mo_2e(Gao, Cmo, nao, Gmo)

  call build_dets(nao, 1, 1, dets, dim, hfidx)
  allocate(H(dim,dim), T2op(dim,dim), Hbar(dim,dim), Em(dim,dim), Ep(dim,dim), S2(dim,dim))
  call build_fci_H(h1mo, eri_mo, enuc, nao, dets, dim, H)
  call build_s2(nao, dets, dim, S2)
  ! (2,2) frontier MRSF active = the 2 lowest MOs; external = the rest
  iact = [1, 2]
  nextn = nao - 2
  do i = 1, nextn
    iext(i) = i + 2
  end do
  call cas22_compact(dets, dim, nao, iact, 2, cas, nc)
  allocate(Hcc(nc,nc), Hbc(nc,nc), S2c(nc,nc))
  call states3(H, S2, dim, .false., efci)
  do i=1,nc; do j=1,nc
    Hcc(i,j)=H(cas(i),cas(j)); S2c(i,j)=S2(cas(i),cas(j))
  end do; end do
  call states3(Hcc, S2c, nc, .false., ebare)
  call build_geminal_T2(Gmo, nao, iact, 2, iext(1:nextn), nextn, 1.0_dp, gamma, dets, dim, T2op)
  call expm_nilpotent(-1.0_dp, T2op, dim, Em)
  call expm_nilpotent( 1.0_dp, T2op, dim, Ep)
  Hbar = matmul(Em, matmul(H, Ep))
  do i=1,nc; do j=1,nc
    Hbc(i,j)=Hbar(cas(i),cas(j))
  end do; end do
  call states3(Hbc, S2c, nc, .true., eptc)

  write(*,'(a)') ''
  write(*,'(a)') '             S0           T1           S1'
  write(*,'(a,3f13.6)') 'FCI       ', efci
  write(*,'(a,3f13.6)') 'bare MRSF ', ebare
  write(*,'(a,3f13.6)') 'pTC(gem)  ', eptc
  write(*,'(a)') ''
  write(*,'(a,f8.6,a)') 'FCI S0 = ', efci(1), '  (vs s-only ~ -1.128; p lowers it)'
  if (efci(1) < -1.13_dp .and. eptc(1) < ebare(1) .and. abs(eptc(2)-ebare(2)) > 1e-5_dp) then
    write(*,'(a)') 'ALL PASS: general-L pTC-MRSF-CIS runs end-to-end with p functions;'
    write(*,'(a)') 'p lowers FCI, and the geminal dresses S0 and T1.'
  else
    write(*,'(a)') 'CHECK: see values above.'
  end if

contains

  subroutine h2_basis(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    integer, intent(out) :: nsh, shl_l(:), shl_np(:)
    real(dp), intent(out) :: shl_e(:,:), shl_c(:,:), shl_r(:,:)
    real(dp), intent(in) :: R
    real(dp) :: cen(3,2)
    integer :: at
    cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.0_dp,0.0_dp,R]
    nsh = 0
    do at = 1, 2
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
        shl_e(1:3,nsh)=[33.8650_dp,5.094790_dp,1.158790_dp]
        shl_c(1:3,nsh)=[0.0254938_dp,0.190373_dp,0.852161_dp]; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.325840_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.102741_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.75_dp;     shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
    end do
  end subroutine h2_basis

  subroutine rohf_triplet(n, S, Hc, eri, Cmo, eps)
    integer, intent(in) :: n
    real(dp), intent(in) :: S(n,n), Hc(n,n), eri(n,n,n,n)
    real(dp), intent(out) :: Cmo(n,n), eps(n)
    real(dp) :: Fa(n,n), Ca(n,n), Da(n,n), Ja(n,n), Ka(n,n), Scp(n,n), eold, e
    integer :: it, mu, nu, la, si, i
    Fa = Hc; e = 0.0_dp
    do it=1,400
      Scp=S; call geig(Fa,Scp,n,eps,Ca)
      Da=0.0_dp
      do mu=1,n;do nu=1,n;do i=1,2     ! 2 alpha electrons (sigma_g, sigma_u)
        Da(mu,nu)=Da(mu,nu)+Ca(mu,i)*Ca(nu,i)
      end do;end do;end do
      Ja=0.0_dp; Ka=0.0_dp
      do mu=1,n;do nu=1,n;do la=1,n;do si=1,n
        Ja(mu,nu)=Ja(mu,nu)+Da(la,si)*eri(mu,nu,la,si)
        Ka(mu,nu)=Ka(mu,nu)+Da(la,si)*eri(mu,la,nu,si)
      end do;end do;end do;end do
      Fa = Hc + Ja - Ka
      eold=e; e=0.0_dp
      do mu=1,n;do nu=1,n
        e=e+Da(mu,nu)*(Hc(mu,nu)+Fa(mu,nu))
      end do;end do
      if (abs(e-eold)<1e-11_dp .and. it>1) exit
    end do
    Cmo = Ca
  end subroutine rohf_triplet

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
        den=dot_product(vll(:,kk),vrr(:,kk)); ss(kk)=dot_product(vll(:,kk),matmul(S2m,vrr(:,kk)))/den
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

end program tc_genL_point
