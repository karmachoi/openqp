!> Ten-no projective transcorrelation pTC(2) for H2 -- F12 route.
!>   E_pTC(2) = E_MP2 + V_gg^gg ,
!>   V_gg^gg = <gg|f12/r12|gg> - sum_{pq in GBS} F_gg^pq g_pq^gg - 2 sum_{x in CABS} F_gg^gx g_gx^gg.
!> f12 = -(1/gamma) e^{-gamma r12} (STG-NG). This stage implements the GBS-RI part
!> (seed minus the GBS pair-sum via the AO RI sum_p|p><p| = S^-1); the CABS term is
!> the next refinement. Validates the V-term assembly and sign (V<0). The seed and F
!> use the validated geminal-Coulomb (vgem_cart) and geminal integrals.
program tc_h2_ptc2
  use precision, only: dp
  use ptc_md,    only: vgem_cart, geminal_cart, cart_norm
  use ptc_ao,    only: ao_ncart, build_ints, build_geminal_ao, stg6
  use tc_geminal_engine, only: ao2mo_2e
  implicit none
  real(dp), parameter :: E_EXACT = -1.174475_dp
  integer, parameter  :: MS = 16, MP = 3
  integer  :: nsh, shl_l(MS), shl_np(MS), nao, nat, i, j, la, sg
  real(dp) :: shl_e(MP,MS), shl_c(MP,MS), shl_r(3,MS), zat(2), rat(3,2)
  real(dp) :: R, gamma, enuc, ehf, emp2c, emp2, seed, gbs_term, V_gbsri, eptc2
  real(dp), allocatable :: S(:,:), Hc(:,:), eri(:,:,:,:), Gao(:,:,:,:), FRao(:,:,:,:)
  real(dp), allocatable :: Cmo(:,:), eps(:), Sinv(:,:), Cg(:)
  real(dp), allocatable :: Mf(:,:), Meri(:,:), T1(:,:)
  real(dp) :: cc(6), omg(6)
  integer  :: ng_stg

  R = 1.4_dp; gamma = 1.0_dp; enuc = 1.0_dp/R
  nat = 2; zat = 1.0_dp; rat(:,1)=[0.0_dp,0.0_dp,0.0_dp]; rat(:,2)=[0.0_dp,0.0_dp,R]
  call h2_aug(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
  nao = 0
  do i = 1, nsh; nao = nao + ao_ncart(shl_l(i)); end do
  allocate(S(nao,nao), Hc(nao,nao), eri(nao,nao,nao,nao), Gao(nao,nao,nao,nao), FRao(nao,nao,nao,nao))
  allocate(Cmo(nao,nao), eps(nao), Sinv(nao,nao), Cg(nao), Mf(nao,nao), Meri(nao,nao), T1(nao,nao))

  call build_ints(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                  nat, zat, rat, nao, S, Hc, eri)
  ! f12 geminal and f12/r12 (geminal-Coulomb) AO tensors, physicist <IJ|.|KL>,
  ! with f12 = -(1/gamma) e^{-gamma r12} (STG-NG; c_k = -a_k/gamma, gem exponent w_k).
  call build_geminal_ao(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), gamma, nao, S, Gao)
  Gao = -(1.0_dp/gamma) * Gao          ! build_geminal_ao fits e^{-gamma r}; scale to f12
  call stg6(gamma, cc, omg, ng_stg)
  call build_fr_ao(nsh, shl_l(1:nsh), shl_np(1:nsh), shl_e, shl_c, shl_r(:,1:nsh), &
                   cc, omg, ng_stg, gamma, nao, S, FRao)

  call rhf(nao, S, Hc, eri, enuc, Cmo, eps, ehf)
  Cg = Cmo(:,1)                          ! sigma_g occupied MO
  call mp2_2e(nao, eri, Cmo, eps, emp2c)
  emp2 = ehf + emp2c
  call inv_sym(S, nao, Sinv)             ! GBS RI = S^-1

  sg = 1
  ! seed = <gg|f12/r12|gg>  (FRao is physicist <IJ|f12/r12|KL>, e1=(I,K),e2=(J,L))
  seed = mo4(FRao, Cg, Cg, Cg, Cg, nao)
  ! Mf(la,sig) = <gg|f12|la sig> = sum_{mu nu} Cg_mu Cg_nu Gao(mu,nu,la,sig)
  ! Meri(la,sig) = <gg|1/r12|la sig> = (g la | g sig) chemist = sum Cg_mu Cg_nu eri(mu,la,nu,sig)
  do la=1,nao; do j=1,nao
    Mf(la,j)   = gg_bra(Gao, Cg, la, j, nao, .true.)    ! physicist tensor
    Meri(la,j) = gg_bra_eri(eri, Cg, la, j, nao)        ! chemist ERI
  end do; end do
  ! GBS term = sum Mf(la,sig) Sinv(la,la') Sinv(sig,sig') Meri(la',sig')
  T1 = matmul(Sinv, matmul(Mf, Sinv))
  gbs_term = 0.0_dp
  do i=1,nao; do j=1,nao; gbs_term = gbs_term + T1(i,j)*Meri(i,j); end do; end do

  V_gbsri = seed - gbs_term
  eptc2   = emp2 + V_gbsri

  write(*,'(a)') '=== Ten-no pTC(2) for H2 (aug-cc-pVDZ GBS, GBS-RI stage) ==='
  write(*,'(a,i0)')   'GBS AOs = ', nao
  write(*,'(a,f12.7)') 'E_HF     = ', ehf
  write(*,'(a,f12.7,a,f8.4,a)') 'E_MP2    = ', emp2, '   (corr ', emp2c*1000, ' mEh)'
  write(*,'(a)') ''
  write(*,'(a,es13.5)') '  seed <gg|f12/r12|gg>   = ', seed
  write(*,'(a,es13.5)') '  GBS-RI pair sum        = ', gbs_term
  write(*,'(a,f12.7,a)') '  V_gg^gg (GBS-RI)       = ', V_gbsri, '   (must be < 0)'
  write(*,'(a)') ''
  write(*,'(a,f12.7)') 'E_pTC(2) [GBS-RI] ~     = ', eptc2
  write(*,'(a,f12.7)') 'E_MP2                   = ', emp2
  write(*,'(a,f12.7)') 'exact (K-W)             = ', E_EXACT
  write(*,'(a)') ''
  if (V_gbsri < 0.0_dp) then
    write(*,'(a)') 'OK: V_gg^gg < 0 -- the explicitly-correlated term lowers MP2 (correct sign).'
    write(*,'(a)') '    (CABS complement term is the next refinement toward the CBS limit.)'
  else
    write(*,'(a)') 'CHECK: V_gg^gg should be negative.'
  end if

contains

  subroutine h2_aug(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, R)
    integer, intent(out) :: nsh, shl_l(:), shl_np(:)
    real(dp), intent(out) :: shl_e(:,:), shl_c(:,:), shl_r(:,:)
    real(dp), intent(in) :: R
    real(dp) :: cen(3,2); integer :: at
    cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.0_dp,0.0_dp,R]
    nsh = 0
    do at = 1, 2
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=3
        shl_e(1:3,nsh)=[13.01_dp,1.962_dp,0.4446_dp]; shl_c(1:3,nsh)=[0.019685_dp,0.137977_dp,0.478148_dp]; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.122_dp;   shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=0; shl_np(nsh)=1; shl_e(1,nsh)=0.02974_dp; shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.727_dp;   shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
      nsh=nsh+1; shl_l(nsh)=1; shl_np(nsh)=1; shl_e(1,nsh)=0.141_dp;   shl_c(1,nsh)=1.0_dp; shl_r(:,nsh)=cen(:,at)
    end do
  end subroutine h2_aug

  !> f12/r12 geminal-Coulomb AO tensor, physicist <IJ|f12/r12|KL> (e1=(I,K),e2=(J,L)),
  !> f12 = -(1/gamma) sum_k a_k e^{-w_k r^2} -> sum_k (-a_k/gamma) vgem(w_k). Normalized.
  subroutine build_fr_ao(nsh, sl, snp, se, sc, sr, cc, omg, ng, gamma, nao, S, FR)
    integer,  intent(in) :: nsh, sl(nsh), snp(nsh), ng, nao
    real(dp), intent(in) :: se(:,:), sc(:,:), sr(3,nsh), cc(:), omg(:), gamma, S(:,:)
    real(dp), intent(out) :: FR(:,:,:,:)
    integer  :: aos(nao), aol(3,nao)
    real(dp) :: nrm(nao)
    integer  :: ss,ll,lx,ly,lz,n, I,J,K,L
    integer  :: si,sj,sk,slh, pi,pj,pk,pl, g
    real(dp) :: ci,cj,ck,cl, ei,ej,ek,el, vacc
    n=0                                          ! enumerate Cartesian AOs
    do ss=1,nsh; ll=sl(ss)
      do lx=ll,0,-1; do ly=ll-lx,0,-1; lz=ll-lx-ly
        n=n+1; aos(n)=ss; aol(1,n)=lx; aol(2,n)=ly; aol(3,n)=lz
      end do; end do
    end do
    do I=1,nao; nrm(I)=1.0_dp/sqrt(S(I,I)); end do
    !$omp parallel do collapse(2) default(shared) schedule(dynamic) &
    !$omp   private(I,J,K,L,si,sj,sk,slh,pi,pj,pk,pl,g,ci,cj,ck,cl,ei,ej,ek,el,vacc)
    do I=1,nao; do J=1,nao; do K=1,nao; do L=1,nao
      si=aos(I); sj=aos(J); sk=aos(K); slh=aos(L)
      vacc=0.0_dp
      do pi=1,snp(si); ei=se(pi,si); ci=sc(pi,si)*cart_norm(aol(:,I),ei)
       do pj=1,snp(sj); ej=se(pj,sj); cj=sc(pj,sj)*cart_norm(aol(:,J),ej)
        do pk=1,snp(sk); ek=se(pk,sk); ck=sc(pk,sk)*cart_norm(aol(:,K),ek)
         do pl=1,snp(slh); el=se(pl,slh); cl=sc(pl,slh)*cart_norm(aol(:,L),el)
          do g=1,ng
            vacc = vacc + (-cc(g)/gamma)*ci*cj*ck*cl* &
                vgem_cart(aol(:,I),sr(:,si),ei, aol(:,K),sr(:,sk),ek, &
                          aol(:,J),sr(:,sj),ej, aol(:,L),sr(:,slh),el, omg(g))
          end do
         end do
        end do
       end do
      end do
      FR(I,J,K,L) = vacc*nrm(I)*nrm(J)*nrm(K)*nrm(L)
    end do; end do; end do; end do
  end subroutine build_fr_ao

  !> contract a physicist tensor T(I,J,K,L)=<IJ|.|KL> with four MO vectors.
  real(dp) function mo4(T, a, b, c, d, n) result(v)
    integer, intent(in) :: n
    real(dp), intent(in) :: T(n,n,n,n), a(n), b(n), c(n), d(n)
    integer :: I,J,K,L
    v=0.0_dp
    do I=1,n; do J=1,n; do K=1,n; do L=1,n
      v = v + a(I)*b(J)*c(K)*d(L)*T(I,J,K,L)
    end do; end do; end do; end do
  end function mo4

  !> <gg|.|la sig> for a physicist tensor: sum_{mu nu} Cg_mu Cg_nu T(mu,nu,la,sig).
  real(dp) function gg_bra(T, Cg, la, sig, n, phys) result(v)
    integer, intent(in) :: la, sig, n
    real(dp), intent(in) :: T(n,n,n,n), Cg(n)
    logical, intent(in) :: phys
    integer :: mu, nu
    v=0.0_dp
    do mu=1,n; do nu=1,n; v = v + Cg(mu)*Cg(nu)*T(mu,nu,la,sig); end do; end do
  end function gg_bra

  !> <gg|1/r12|la sig> = (g la | g sig) chemist = sum Cg_mu Cg_nu eri(mu,la,nu,sig).
  real(dp) function gg_bra_eri(eri, Cg, la, sig, n) result(v)
    integer, intent(in) :: la, sig, n
    real(dp), intent(in) :: eri(n,n,n,n), Cg(n)
    integer :: mu, nu
    v=0.0_dp
    do mu=1,n; do nu=1,n; v = v + Cg(mu)*Cg(nu)*eri(mu,la,nu,sig); end do; end do
  end function gg_bra_eri

  subroutine rhf(n, S, Hc, eri, enuc, Cmo, eps, ehf)
    integer, intent(in) :: n
    real(dp), intent(in) :: S(n,n), Hc(n,n), eri(n,n,n,n), enuc
    real(dp), intent(out) :: Cmo(n,n), eps(n), ehf
    real(dp) :: F(n,n), C(n,n), D(n,n), G(n,n), Scp(n,n), eold, e
    integer :: it,mu,nu,la,si
    F=Hc; e=0.0_dp
    do it=1,300
      Scp=S; call geig(F,Scp,n,eps,C)
      D=0.0_dp; do mu=1,n;do nu=1,n; D(mu,nu)=2.0_dp*C(mu,1)*C(nu,1); end do;end do
      G=0.0_dp
      do mu=1,n;do nu=1,n;do la=1,n;do si=1,n
        G(mu,nu)=G(mu,nu)+D(la,si)*(eri(mu,nu,la,si)-0.5_dp*eri(mu,si,la,nu))
      end do;end do;end do;end do
      F=Hc+G; eold=e; e=0.0_dp
      do mu=1,n;do nu=1,n; e=e+0.5_dp*D(mu,nu)*(Hc(mu,nu)+F(mu,nu)); end do;end do
      if (abs(e-eold)<1e-12_dp .and. it>1) exit
    end do
    Cmo=C; ehf=e+enuc
  end subroutine rhf

  !> closed-shell MP2 correlation for 1 occupied MO (H2): sum_{ab} (ia|jb)[2(ia|jb)-(ib|ja)]/D, i=j=1.
  subroutine mp2_2e(n, eri, Cmo, eps, ec)
    integer, intent(in) :: n
    real(dp), intent(in) :: eri(n,n,n,n), Cmo(n,n), eps(n)
    real(dp), intent(out) :: ec
    real(dp) :: emo(n,n,n,n)
    integer :: a,b
    call ao2mo_2e(eri, Cmo, n, emo)           ! chemist (pq|rs) in MO
    ec=0.0_dp
    do a=2,n; do b=2,n
      ec = ec + emo(1,a,1,b)*(2.0_dp*emo(1,a,1,b)-emo(1,b,1,a))/(2.0_dp*eps(1)-eps(a)-eps(b))
    end do; end do
  end subroutine mp2_2e

  subroutine geig(F,S,n,w,C)
    integer, intent(in) :: n
    real(dp), intent(inout) :: F(n,n), S(n,n)
    real(dp), intent(out) :: w(n), C(n,n)
    real(dp) :: wq(1); real(dp), allocatable :: wk(:); integer :: info, lw
    call dsygv(1,'V','U',n,F,n,S,n,w,wq,-1,info); lw=int(wq(1)); allocate(wk(lw))
    call dsygv(1,'V','U',n,F,n,S,n,w,wk,lw,info); C=F; deallocate(wk)
  end subroutine geig

  subroutine inv_sym(A, n, Ainv)
    integer, intent(in) :: n
    real(dp), intent(in) :: A(n,n)
    real(dp), intent(out) :: Ainv(n,n)
    real(dp) :: M(n,n), wq(1); real(dp), allocatable :: wk(:); integer :: ipiv(n), info, lw
    M=A; call dgetrf(n,n,M,n,ipiv,info)
    call dgetri(n,M,n,ipiv,wq,-1,info); lw=int(wq(1)); allocate(wk(lw))
    call dgetri(n,M,n,ipiv,wk,lw,info); Ainv=M; deallocate(wk)
  end subroutine inv_sym

end program tc_h2_ptc2
