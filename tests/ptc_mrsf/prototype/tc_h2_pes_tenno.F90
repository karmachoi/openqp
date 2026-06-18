!> H2 dissociation PES with the GENUINE Ten-no transcorrelation through the MRSF
!> response. ROHF high-spin (triplet) reference (dissociates correctly), the (2,2)
!> frontier active space, and the genuine non-Hermitian H_bar built from Ten-no's
!> cusp-fixed geminal (c=1/2; the H2 pair is antiparallel in the Ms=0 sector). For
!> each bond length we report S0/T1/S1 for (a) the in-basis FCI, (b) bare
!> MRSF-CIS(2,2), and (c) TC-MRSF-CIS(2,2). Output -> pes_tenno.dat for plotting.
!>
!> This replaces the earlier geminal-DOUBLES PES with the genuine first-quantized
!> H_bar = e^{-tau} H e^tau (validated to reach the exact CBS limit on H2).
program tc_h2_pes_tenno
  use precision, only: dp
  use ptc_ao,        only: ao_ncart, build_ints
  use tc_boyshandy,  only: build_tc2e_ao
  use tc_geminal_engine, only: ao2mo_1e, ao2mo_2e
  implicit none
  integer, parameter :: MS = 16, MP = 6
  integer  :: nsh, shl_l(MS), shl_np(MS), nao, nat, i, ir
  real(dp) :: shl_e(MP,MS), shl_c(MP,MS), shl_r(3,MS), zat(2), rat(3,2)
  real(dp) :: R, gamma, enuc
  real(dp), allocatable :: S(:,:), Hc(:,:), eri(:,:,:,:), Lin(:,:,:,:), Quad(:,:,:,:)
  real(dp), allocatable :: Cmo(:,:), eps(:), h1mo(:,:), eri_mo(:,:,:,:), Lin_mo(:,:,:,:), Quad_mo(:,:,:,:)
  real(dp) :: efci(3), ebare(3), etc(3)
  integer  :: u

  gamma = 0.7_dp            ! Slater-geminal exponent (cusp recovery near-optimal for cc-pVDZ H2)
  nat = 2; zat = 1.0_dp
  open(newunit=u, file='pes_tenno.dat', status='replace', action='write')
  write(u,'(a)') '# R(bohr)  FCI:S0 T1 S1   bareMRSF:S0 T1 S1   TC-MRSF:S0 T1 S1   (Hartree)'
  write(*,'(a)') '=== Genuine Ten-no TC-MRSF-CIS : H2 dissociation (cc-pVDZ, gamma=0.7) ==='
  write(*,'(a)') '   R       FCI_S0    bare_S0    TC_S0   |  FCI(T1-S0) bare(T1-S0) TC(T1-S0)  | FCI(S1-S0) TC(S1-S0)'

  do ir = 0, 24
    R = 0.8_dp + 0.18_dp*real(ir,dp)     ! 0.8 .. 5.12 bohr
    enuc = 1.0_dp/R
    rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
    call h2_ccpvdz(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    nao = 0
    do i = 1, nsh
      nao = nao + ao_ncart(shl_l(i))
    end do
    if (allocated(S)) deallocate(S,Hc,eri,Lin,Quad,Cmo,eps,h1mo,eri_mo,Lin_mo,Quad_mo)
    allocate(S(nao,nao),Hc(nao,nao),eri(nao,nao,nao,nao),Lin(nao,nao,nao,nao),Quad(nao,nao,nao,nao))
    allocate(Cmo(nao,nao),eps(nao),h1mo(nao,nao),eri_mo(nao,nao,nao,nao),Lin_mo(nao,nao,nao,nao),Quad_mo(nao,nao,nao,nao))

    call build_ints(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), nat, zat, rat, nao, S, Hc, eri)
    call build_tc2e_ao(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), gamma, nao, Lin, Quad)
    call rohf_triplet(nao, S, Hc, eri, Cmo, eps)
    call ao2mo_1e(Hc, Cmo, nao, h1mo)
    call ao2mo_2e(eri, Cmo, nao, eri_mo)
    call ao2mo_2e(Lin, Cmo, nao, Lin_mo)
    call ao2mo_2e(Quad, Cmo, nao, Quad_mo)

    call solve3(nao, h1mo, eri_mo, Lin_mo, Quad_mo, enuc, 0.0_dp, .false., efci)   ! full FCI, c=0
    call solve3(nao, h1mo, eri_mo, Lin_mo, Quad_mo, enuc, 0.0_dp, .true.,  ebare)  ! (2,2), c=0
    call solve3(nao, h1mo, eri_mo, Lin_mo, Quad_mo, enuc, 0.5_dp, .true.,  etc)    ! (2,2), c=1/2 (genuine TC)

    write(u,'(f7.3,9(1x,f12.7))') R, efci, ebare, etc
    write(*,'(f7.3,3(1x,f10.6),3x,3(1x,f9.4),3x,2(1x,f9.4))') R, efci(1), ebare(1), etc(1), &
         efci(2)-efci(1), ebare(2)-ebare(1), etc(2)-etc(1), efci(3)-efci(1), etc(3)-etc(1)
  end do
  close(u)
  write(*,'(a)') 'wrote pes_tenno.dat'

contains

  subroutine h2_ccpvdz(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    integer, intent(out) :: nsh, shl_l(:), shl_np(:)
    real(dp), intent(out) :: shl_e(:,:), shl_c(:,:), shl_r(:,:)
    real(dp), intent(in) :: R
    real(dp) :: cen(3,2)
    integer :: at
    cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.0_dp,0.0_dp,R]
    nsh = 0
    do at = 1, 2
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
        shl_e(1:3,nsh)=[13.01_dp,1.962_dp,0.4446_dp]
        shl_c(1:3,nsh)=[0.019685_dp,0.137977_dp,0.478148_dp]; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.122_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.727_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
    end do
  end subroutine h2_ccpvdz

  !> ROHF high-spin (triplet Ms=+1): two alpha electrons in the two lowest MOs.
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
      do mu=1,n;do nu=1,n;do i=1,2
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

  !> Build the 2-electron (Ms=0, alpha-beta) TC Hamiltonian and S^2, optionally
  !> downfolded to the (2,2) frontier, and return [S0, T1, S1] classified by <S^2>.
  !> c=0 gives the bare operator; c=1/2 the genuine Ten-no transcorrelation (the H2
  !> pair is antiparallel). Non-Hermitian -> DGEEV with biorthogonal <S^2>.
  subroutine solve3(n, h1, eri_c, Linm, Quadm, enuc, c, cas22, out)
    integer, intent(in) :: n
    real(dp), intent(in) :: h1(n,n), eri_c(n,n,n,n), Linm(n,n,n,n), Quadm(n,n,n,n), enuc, c
    logical, intent(in) :: cas22
    real(dp), intent(out) :: out(3)
    integer :: dfull, d, p, q, pp, qq, r, col, map(4), idx
    real(dp), allocatable :: H(:,:), S2(:,:)
    dfull = n*n
    ! full H and S2 in the (alpha,beta) basis
    block
      real(dp) :: Hf(dfull,dfull), S2f(dfull,dfull)
      Hf = 0.0_dp; S2f = 0.0_dp
      do p=1,n; do q=1,n
        col = (p-1)*n+q
        do pp=1,n; do qq=1,n
          r = (pp-1)*n+qq
          Hf(r,col) = eri_c(pp,p,qq,q) + c*Linm(pp,qq,p,q) + c*c*Quadm(pp,qq,p,q)
          if (qq==q) Hf(r,col) = Hf(r,col) + h1(pp,p)
          if (pp==p) Hf(r,col) = Hf(r,col) + h1(qq,q)
        end do; end do
        Hf(col,col) = Hf(col,col) + enuc
      end do; end do
      do p=1,n; do q=1,n
        if (p/=q) then
          S2f((p-1)*n+q,(p-1)*n+q) = 1.0_dp
          S2f((q-1)*n+p,(p-1)*n+q) = -1.0_dp
        end if
      end do; end do
      if (cas22) then
        d = 4
        map = [ (1-1)*n+1, (1-1)*n+2, (2-1)*n+1, (2-1)*n+2 ]  ! (p,q) in {1,2}
        allocate(H(d,d), S2(d,d))
        do p=1,4; do q=1,4
          H(p,q) = Hf(map(p),map(q)); S2(p,q) = S2f(map(p),map(q))
        end do; end do
      else
        d = dfull
        allocate(H(d,d), S2(d,d)); H = Hf; S2 = S2f
      end if
    end block
    call states_from(H, S2, d, out)
    deallocate(H, S2)
  end subroutine solve3

  !> DGEEV solve; sort by real eigenvalue; classify by <S^2> (biorthogonal for the
  !> non-Hermitian case); return lowest singlet (S0), lowest triplet (T1), 2nd singlet (S1).
  subroutine states_from(Hm, S2m, m, out)
    integer, intent(in) :: m
    real(dp), intent(inout) :: Hm(m,m)
    real(dp), intent(in) :: S2m(m,m)
    real(dp), intent(out) :: out(3)
    real(dp) :: wr(m), wi(m), VL(m,m), VR(m,m), s2val(m), den, ord(m)
    real(dp), allocatable :: work(:)
    integer :: info, lwork, k, i2, j2, o(m), tmp, ns0, ns1, nt0
    lwork = 8*m; allocate(work(lwork))
    call dgeev('V','V', m, Hm, m, wr, wi, VL, m, VR, m, work, lwork, info)
    do k=1,m
      den = dot_product(VL(:,k), VR(:,k))
      s2val(k) = dot_product(VL(:,k), matmul(S2m, VR(:,k))) / den
    end do
    do i2=1,m; o(i2)=i2; end do
    do i2=1,m-1; do j2=i2+1,m
      if (wr(o(j2)) < wr(o(i2))) then
        tmp=o(i2); o(i2)=o(j2); o(j2)=tmp
      end if
    end do; end do
    out = huge(1.0_dp); ns0=0; nt0=0; ns1=0
    do i2=1,m
      k=o(i2)
      if (s2val(k) < 1.0_dp) then          ! singlet
        if (ns0==0) then; out(1)=wr(k); ns0=1
        else if (ns1==0) then; out(3)=wr(k); ns1=1; end if
      else                                  ! triplet
        if (nt0==0) then; out(2)=wr(k); nt0=1; end if
      end if
    end do
    deallocate(work)
  end subroutine states_from

end program tc_h2_pes_tenno
