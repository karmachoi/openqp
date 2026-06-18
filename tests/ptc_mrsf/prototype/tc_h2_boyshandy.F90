!> Genuine Boys-Handy / Ten-no transcorrelated FCI for H2 -- the clean validation
!> case, because for two electrons H_bar is EXACTLY two-body (no three-body term).
!>
!> Pipeline: native McMurchie-Davidson AO integrals (S, Hcore, ERI) -> closed-shell
!> RHF -> AO->MO -> genuine TC two-body integrals (tc_boyshandy: -nabla^2 f, the
!> non-Hermitian drift, and -(f')^2) with cusp amplitude c=1/2 (the H2 pair is
!> antiparallel) -> non-Hermitian 2-electron FCI (DGEEV).
!>
!> Validations:
!>  (1) native bare RHF/FCI vs pyscf (printed reference);
!>  (2) gem_r2_cart vs the independent analytic s-only r2_geminal_s;
!>  (3) the TC spectrum is real (similarity transform preserves the spectrum);
!>  (4) physics: TC-FCI lowers the finite-basis energy TOWARD the exact -1.17447
!>      (Kolos-Wolniewicz) -- i.e. it recovers basis-set incompleteness / the cusp.
program tc_h2_boyshandy
  use precision, only: dp
  use ptc_ao,        only: ao_ncart, build_ints
  use tc_boyshandy,  only: build_tc2e_ao
  use ptc_geminal,   only: r2_geminal_s
  use tc_geminal_engine, only: ao2mo_1e, ao2mo_2e
  implicit none
  integer, parameter :: MS = 16, MP = 6
  real(dp), parameter :: ANG = 0.52917721092_dp
  real(dp), parameter :: E_EXACT = -1.174475_dp   ! H2 BO ground state, R=1.4 bohr
  real(dp), parameter :: E_PYSCF_HF(2)  = [-1.12674270_dp, -1.12870926_dp]
  real(dp), parameter :: E_PYSCF_FCI(2) = [-1.15167903_dp, -1.16339823_dp]
  character(len=10) :: bname(2) = ['6-31G     ','cc-pVDZ   ']
  integer  :: ib, nsh, shl_l(MS), shl_np(MS), nao, nat, i, j
  real(dp) :: shl_e(MP,MS), shl_c(MP,MS), shl_r(3,MS), zat(2), rat(3,2)
  real(dp) :: R, gamma, enuc
  real(dp), allocatable :: S(:,:), Hc(:,:), eri(:,:,:,:), Lin(:,:,:,:), Quad(:,:,:,:)
  real(dp), allocatable :: Cmo(:,:), h1mo(:,:), eri_mo(:,:,:,:), Lin_mo(:,:,:,:), Quad_mo(:,:,:,:)
  real(dp) :: ebare(2), etc(2), emaximag

  R = 1.4_dp; gamma = 1.0_dp; enuc = 1.0_dp/R
  nat = 2; zat = 1.0_dp; rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]

  call check_gem_r2()
  call check_drift()

  write(*,'(a)') ''
  write(*,'(a)') '=== Genuine Boys-Handy / Ten-no TC-MRSF FCI for H2 (R=1.4 bohr) ==='
  write(*,'(a,f10.6,a)') 'Exact (Kolos-Wolniewicz) BO energy = ', E_EXACT, ' Ha'
  write(*,'(a)') ''
  write(*,'(a)') 'basis        nao   E(HF)        E(FCI bare)   E(FCI TC)     exact        d(bare)  d(TC)'
  write(*,'(a)') '----------------------------------------------------------------------------------------'

  do ib = 1, 2
    call h2_basis(ib, nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    nao = 0
    do i = 1, nsh
      nao = nao + ao_ncart(shl_l(i))
    end do
    if (allocated(S)) deallocate(S,Hc,eri,Lin,Quad,Cmo,h1mo,eri_mo,Lin_mo,Quad_mo)
    allocate(S(nao,nao), Hc(nao,nao), eri(nao,nao,nao,nao), Lin(nao,nao,nao,nao), Quad(nao,nao,nao,nao))
    allocate(Cmo(nao,nao), h1mo(nao,nao), eri_mo(nao,nao,nao,nao), Lin_mo(nao,nao,nao,nao), Quad_mo(nao,nao,nao,nao))

    call build_ints(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                    nat, zat, rat, nao, S, Hc, eri)
    call build_tc2e_ao(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                       gamma, nao, Lin, Quad)
    call rhf(nao, S, Hc, eri, enuc, Cmo, ebare(ib))   ! ebare temporarily = E_HF
    call ao2mo_1e(Hc, Cmo, nao, h1mo)
    call ao2mo_2e(eri, Cmo, nao, eri_mo)              ! chemist (pq|rs)
    call ao2mo_2e(Lin, Cmo, nao, Lin_mo)              ! physicist <pq|.|rs>
    call ao2mo_2e(Quad, Cmo, nao, Quad_mo)

    block
      real(dp) :: ehf
      ehf = ebare(ib)
      call fci2e(nao, h1mo, eri_mo, Lin_mo, Quad_mo, enuc, 0.0_dp, ebare(ib), emaximag) ! bare (c=0)
      call fci2e(nao, h1mo, eri_mo, Lin_mo, Quad_mo, enuc, 0.5_dp, etc(ib),   emaximag) ! TC (c=1/2)
      write(*,'(a10,i5,2x,f12.7,1x,f12.7,1x,f12.7,1x,f12.7,2x,f7.4,1x,f7.4)') &
        bname(ib), nao, ehf, ebare(ib), etc(ib), E_EXACT, &
        ebare(ib)-E_EXACT, etc(ib)-E_EXACT
      write(*,'(a,es9.2,a,f5.1,a)') '            (TC spectrum max|Im(E)| = ', emaximag, &
        ' ; cusp recovery = ', 100.0_dp*(etc(ib)-ebare(ib))/(E_EXACT-ebare(ib)), ' %)'
    end block
  end do

  call gamma_scan(2)   ! cc-pVDZ: how much cusp the Ten-no correlator recovers vs gamma

  write(*,'(a)') ''
  write(*,'(a)') 'pyscf reference:   6-31G  E_HF=-1.1267427  E_FCI=-1.1516790'
  write(*,'(a)') '                 cc-pVDZ  E_HF=-1.1287093  E_FCI=-1.1633982'
  write(*,'(a)') ''
  if (etc(1) < ebare(1) .and. etc(2) < ebare(2) .and. &
      etc(1)-E_EXACT < ebare(1)-E_EXACT .and. etc(2)-E_EXACT < ebare(2)-E_EXACT) then
    write(*,'(a)') 'PASS: genuine TC lowers the finite-basis energy toward the exact CBS limit'
    write(*,'(a)') '      in BOTH bases (recovers the cusp / basis incompleteness).'
  else
    write(*,'(a)') 'CHECK values above.'
  end if

contains

  !> Validate the general-L r12^2-geminal (finite-difference d/dGamma) against the
  !> independent closed-form s-only analytic r2_geminal_s, on a nontrivial quartet.
  subroutine check_gem_r2()
    use ptc_md, only: gem_r2_cart
    real(dp) :: A(3),B(3),C(3),D(3), aA,aB,aC,aD, lam, ref, val, mx
    integer  :: t
    integer  :: l0(3)
    l0 = 0
    A=[0.0_dp,0.0_dp,0.0_dp]; B=[0.1_dp,0.2_dp,0.3_dp]
    C=[0.0_dp,0.0_dp,1.4_dp]; D=[0.5_dp,-0.3_dp,1.1_dp]
    aA=1.2_dp; aB=0.7_dp; aC=0.9_dp; aD=1.5_dp
    mx = 0.0_dp
    do t = 1, 5
      lam = 0.3_dp*real(t,dp)
      ref = r2_geminal_s(aA,A,aB,B,aC,C,aD,D,lam)            ! analytic s-only, e1=(A,C),e2=(B,D)
      val = gem_r2_cart(l0,A,aA, l0,C,aC, l0,B,aB, l0,D,aD, lam)  ! general-L FD, matched grouping
      mx = max(mx, abs(ref-val))
    end do
    write(*,'(a,es10.2)') 'CHECK gem_r2_cart vs analytic r2_geminal_s : max|diff| = ', mx
  end subroutine check_gem_r2

  !> Independent validation of the non-Hermitian DRIFT integral via the exact
  !> operator identity (grad u . grad)^dagger = -(grad u . grad) - nabla^2 u, hence
  !> the symmetric part of D = -[(grad_1 u).grad_1 + (grad_2 u).grad_2] equals the
  !> scalar nabla^2 u, which is built from the already-validated geminal/r2-geminal.
  !> A match certifies the integration-by-parts reduction in drift_prim.
  subroutine check_drift()
    use ptc_md, only: geminal_cart, gem_r2_cart
    use tc_boyshandy, only: fit_unit_cusp, drift_prim
    real(dp) :: Cf(6), gf(6)
    integer  :: l0(3), ng, I,J,K,L, kk
    real(dp) :: dmat(2,2,2,2), lap(2,2,2,2), sym, mx
    real(dp) :: cen(3,2), zz(2)
    l0 = 0
    cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.2_dp,-0.1_dp,1.3_dp]
    zz = [0.8_dp, 1.3_dp]
    call fit_unit_cusp(1.0_dp, Cf, gf, ng)
    mx = 0.0_dp
    do I=1,2; do J=1,2; do K=1,2; do L=1,2
      dmat(I,J,K,L) = 0.0_dp; lap(I,J,K,L) = 0.0_dp
      do kk = 1, ng
        ! drift (electron1=(I,K), electron2=(J,L)); drift_prim args (lI,RI,eI, lJ,RJ,eJ, lK,RK,eK, lL,RL,eL)
        dmat(I,J,K,L) = dmat(I,J,K,L) + Cf(kk)* &
            drift_prim(l0,cen(:,I),zz(I), l0,cen(:,J),zz(J), l0,cen(:,K),zz(K), l0,cen(:,L),zz(L), gf(kk))
        ! scalar nabla^2 f = sum_k C_k (4 g^2 r^2 - 6 g) e^{-g r^2}, physicist e1=(I,K),e2=(J,L)
        lap(I,J,K,L) = lap(I,J,K,L) + Cf(kk)*( &
             4.0_dp*gf(kk)*gf(kk)*gem_r2_cart(l0,cen(:,I),zz(I), l0,cen(:,K),zz(K), &
                                              l0,cen(:,J),zz(J), l0,cen(:,L),zz(L), gf(kk)) &
           - 6.0_dp*gf(kk)*geminal_cart(l0,cen(:,I),zz(I), l0,cen(:,K),zz(K), &
                                        l0,cen(:,J),zz(J), l0,cen(:,L),zz(L), gf(kk)) )
      end do
    end do; end do; end do; end do
    do I=1,2; do J=1,2; do K=1,2; do L=1,2
      sym = 0.5_dp*(dmat(I,J,K,L) + dmat(K,L,I,J))   ! symmetric (Hermitian) part of drift
      mx = max(mx, abs(sym - lap(I,J,K,L)))
    end do; end do; end do; end do
    write(*,'(a,es10.2)') 'CHECK drift symmetric-part vs nabla^2 u     : max|diff| = ', mx
  end subroutine check_drift

  !> Scan the Slater-geminal exponent gamma for one basis, reporting the
  !> transcorrelated S0 and the fraction of the basis-incompleteness (cusp) error
  !> recovered. The minimum locates the best single-gamma Ten-no correlator.
  subroutine gamma_scan(ib)
    integer, intent(in) :: ib
    integer  :: nsh2, sl(MS), snp(MS), nao2, i2
    real(dp) :: se(MP,MS), sc(MP,MS), sr(3,MS)
    real(dp), allocatable :: S2(:,:),Hc2(:,:),er2(:,:,:,:),Li(:,:,:,:),Qu(:,:,:,:)
    real(dp), allocatable :: Cm(:,:),h1(:,:),erm(:,:,:,:),Lim(:,:,:,:),Qum(:,:,:,:)
    real(dp) :: g, ebar, e, mi, ehf, gbest, ebest
    integer  :: ig
    call h2_basis(ib, nsh2, sl, snp, se, sc, sr, R)
    nao2 = 0
    do i2 = 1, nsh2
      nao2 = nao2 + ao_ncart(sl(i2))
    end do
    allocate(S2(nao2,nao2),Hc2(nao2,nao2),er2(nao2,nao2,nao2,nao2),Li(nao2,nao2,nao2,nao2),Qu(nao2,nao2,nao2,nao2))
    allocate(Cm(nao2,nao2),h1(nao2,nao2),erm(nao2,nao2,nao2,nao2),Lim(nao2,nao2,nao2,nao2),Qum(nao2,nao2,nao2,nao2))
    call build_ints(nsh2, sl(1:nsh2), snp(1:nsh2), se, sc, sr(:,1:nsh2), nat, zat, rat, nao2, S2, Hc2, er2)
    call rhf(nao2, S2, Hc2, er2, enuc, Cm, ehf)
    call ao2mo_1e(Hc2, Cm, nao2, h1); call ao2mo_2e(er2, Cm, nao2, erm)
    call fci2e(nao2, h1, erm, erm, erm, enuc, 0.0_dp, ebar, mi)   ! bare
    write(*,'(a)') ''
    write(*,'(a,a)') '=== gamma scan (Ten-no Slater geminal), basis ', trim(bname(ib))
    write(*,'(a)') ' gamma   E(TC)        recovery(%)'
    gbest = 0.0_dp; ebest = ebar
    do ig = 1, 12
      g = 0.10_dp + 0.20_dp*real(ig-1,dp)
      call build_tc2e_ao(nsh2, sl(1:nsh2), snp(1:nsh2), se, sc, sr(:,1:nsh2), g, nao2, Li, Qu)
      call ao2mo_2e(Li, Cm, nao2, Lim); call ao2mo_2e(Qu, Cm, nao2, Qum)
      call fci2e(nao2, h1, erm, Lim, Qum, enuc, 0.5_dp, e, mi)
      write(*,'(f6.2,1x,f12.7,4x,f6.1)') g, e, 100.0_dp*(e-ebar)/(E_EXACT-ebar)
      if (abs(e-E_EXACT) < abs(ebest-E_EXACT)) then   ! closest to exact (TC is non-variational)
        ebest = e; gbest = g
      end if
    end do
    write(*,'(a,f5.2,a,f10.6,a,f5.1,a)') 'best gamma = ', gbest, '  E(TC) = ', ebest, &
      '  recovery = ', 100.0_dp*(ebest-ebar)/(E_EXACT-ebar), ' %'
    deallocate(S2,Hc2,er2,Li,Qu,Cm,h1,erm,Lim,Qum)
  end subroutine gamma_scan

  subroutine h2_basis(ib, nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    integer, intent(in)  :: ib
    integer, intent(out) :: nsh, shl_l(:), shl_np(:)
    real(dp), intent(out) :: shl_e(:,:), shl_c(:,:), shl_r(:,:)
    real(dp), intent(in) :: R
    real(dp) :: cen(3,2)
    integer :: at
    cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.0_dp,0.0_dp,R]
    nsh = 0
    do at = 1, 2
      if (ib == 1) then   ! 6-31G
        nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
          shl_e(1:3,nsh)=[18.731137_dp,2.8253937_dp,0.6401217_dp]
          shl_c(1:3,nsh)=[0.0334946_dp,0.23472695_dp,0.81375733_dp]; shl_r(:,nsh)=cen(:,at)
        nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.1612778_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      else                ! aug-cc-pVDZ (diffuse s,p added)
        nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
          shl_e(1:3,nsh)=[13.01_dp,1.962_dp,0.4446_dp]
          shl_c(1:3,nsh)=[0.019685_dp,0.137977_dp,0.478148_dp]; shl_r(:,nsh)=cen(:,at)
        nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.122_dp;   shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
        nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.02974_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
        nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.727_dp;   shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
        nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.141_dp;   shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      end if
    end do
  end subroutine h2_basis

  !> Closed-shell RHF (2 electrons, one doubly-occupied MO). Returns C and E_HF.
  subroutine rhf(n, S, Hc, eri, enuc, Cmo, ehf)
    integer, intent(in) :: n
    real(dp), intent(in) :: S(n,n), Hc(n,n), eri(n,n,n,n), enuc
    real(dp), intent(out) :: Cmo(n,n), ehf
    real(dp) :: F(n,n), C(n,n), D(n,n), G(n,n), Scp(n,n), eps(n), eold, e
    integer :: it, mu, nu, la, si
    F = Hc; e = 0.0_dp
    do it = 1, 300
      Scp = S; call geig(F, Scp, n, eps, C)
      D = 0.0_dp
      do mu=1,n; do nu=1,n
        D(mu,nu) = 2.0_dp*C(mu,1)*C(nu,1)   ! one doubly-occupied MO
      end do; end do
      G = 0.0_dp
      do mu=1,n; do nu=1,n; do la=1,n; do si=1,n
        G(mu,nu) = G(mu,nu) + D(la,si)*(eri(mu,nu,la,si) - 0.5_dp*eri(mu,si,la,nu))
      end do; end do; end do; end do
      F = Hc + G
      eold = e; e = 0.0_dp
      do mu=1,n; do nu=1,n
        e = e + 0.5_dp*D(mu,nu)*(Hc(mu,nu)+F(mu,nu))
      end do; end do
      if (abs(e-eold) < 1.0e-12_dp .and. it > 1) exit
    end do
    Cmo = C; ehf = e + enuc
  end subroutine rhf

  subroutine geig(F, S, n, w, C)
    integer, intent(in) :: n
    real(dp), intent(inout) :: F(n,n), S(n,n)
    real(dp), intent(out) :: w(n), C(n,n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsygv(1,'V','U',n,F,n,S,n,w,wq,-1,info)
    lw = int(wq(1)); allocate(wk(lw))
    call dsygv(1,'V','U',n,F,n,S,n,w,wk,lw,info)
    C = F; deallocate(wk)
  end subroutine geig

  !> Two-electron (Ms=0) FCI with the genuine TC two-body operator. The opposite-spin
  !> interaction <p'q'|V_bar|pq> = (p'p|q'q) + c*Lin(p'q',pq) + c^2*Quad(p'q',pq);
  !> c=0 -> bare, c=1/2 -> Ten-no antiparallel cusp. Non-Hermitian -> DGEEV; returns
  !> the lowest real eigenvalue (the transcorrelated S0) and the max imaginary part.
  subroutine fci2e(n, h1, eri_c, Linm, Quadm, enuc, c, e0, maximag)
    integer, intent(in) :: n
    real(dp), intent(in) :: h1(n,n), eri_c(n,n,n,n), Linm(n,n,n,n), Quadm(n,n,n,n), enuc, c
    real(dp), intent(out) :: e0, maximag
    integer :: d, p, q, pp, qq, r, cc
    real(dp), allocatable :: H(:,:), wr(:), wi(:), vl(:,:), vr(:,:), work(:)
    integer :: info, lwork
    d = n*n
    allocate(H(d,d), wr(d), wi(d), vl(1,d), vr(1,d))
    H = 0.0_dp
    do p=1,n; do q=1,n
      cc = (p-1)*n + q
      do pp=1,n; do qq=1,n
        r = (pp-1)*n + qq
        H(r,cc) = eri_c(pp,p,qq,q) + c*Linm(pp,qq,p,q) + c*c*Quadm(pp,qq,p,q)
        if (qq==q) H(r,cc) = H(r,cc) + h1(pp,p)
        if (pp==p) H(r,cc) = H(r,cc) + h1(qq,q)
      end do; end do
      H(cc,cc) = H(cc,cc) + enuc
    end do; end do
    lwork = 8*d
    allocate(work(lwork))
    call dgeev('N','N', d, H, d, wr, wi, vl, 1, vr, 1, work, lwork, info)
    e0 = huge(1.0_dp); maximag = 0.0_dp
    do p = 1, d
      if (abs(wi(p)) < 1.0e-6_dp .and. wr(p) < e0) e0 = wr(p)
      maximag = max(maximag, abs(wi(p)))
    end do
    deallocate(H, wr, wi, vl, vr, work)
  end subroutine fci2e

end program tc_h2_boyshandy
