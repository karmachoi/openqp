!> Shared engine for native genuine-geminal pTC-MRSF-CIS on H_n (s-only):
!>   ROHF high-spin reference, native integrals (module ptc_geminal), determinant
!>   FCI, the cusp-fixed GEMINAL F12 transcorrelation (spin-resolved 1/2 opposite-
!>   spin, 1/4 same-spin amplitudes from the Slater geminal, NOT an MP2 proxy),
!>   downfold into the (2,2) frontier MRSF space, and the production non-Hermitian
!>   solver tc_nonsym_tda_eig.
module tc_geminal_engine
  use precision, only: dp
  use ptc_geminal
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none
  private
  public :: rohf_highspin, ao2mo_1e, ao2mo_2e, build_dets, build_fci_H
  public :: build_geminal_T2, build_conv_T2, build_msz_T2, expm_nilpotent, sym_eig_vec
  public :: cas22_compact, build_s2, geminal_mo, eri_chem
  public :: ann, cre, findidx

contains

  !> High-spin ROHF (na alpha, nb beta) via UHF; returns alpha MO coeffs/energies
  !> (the reference orbitals) plus AO core/overlap and the chemist ERI tensor.
  subroutine rohf_highspin(n, na, nb, npr, exs, cos_, cns, nat, rat, &
                           Cmo, eps, e_scf, enuc, h1ao, eri_c)
    integer,  intent(in)  :: n, na, nb, npr(n), nat
    real(dp), intent(in)  :: exs(:,:), cos_(:,:), cns(3,n), rat(:,:)
    real(dp), intent(out) :: Cmo(n,n), eps(n), e_scf, enuc, h1ao(n,n), eri_c(n,n,n,n)
    real(dp) :: S(n,n), T(n,n), V(n,n), Hc(n,n), Fa(n,n), Fb(n,n)
    real(dp) :: Da(n,n), Db(n,n), Ca(n,n), Cb(n,n), epb(n), M(n,n,n,n)
    real(dp) :: Scp(n,n), zat(nat), eold, Ja(n,n), Ka(n,n), Jb(n,n), Kb(n,n)
    integer  :: i, j, k, l, it, mu, nu, la, si
    zat = 1.0_dp
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_OVERLAP, nat, zat, rat, S)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_KINETIC, nat, zat, rat, T)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_NUCLEAR, nat, zat, rat, V)
    Hc = T + V; h1ao = Hc
    call ptc_s_ao_tensor(n, npr, exs, cos_, cns, PTC_OP_ERI, 0.0_dp, M)
    do i=1,n; do j=1,n; do k=1,n; do l=1,n
      eri_c(i,j,k,l) = M(i,k,j,l)
    end do; end do; end do; end do
    enuc = 0.0_dp
    do i=1,nat; do j=i+1,nat
      enuc = enuc + 1.0_dp/sqrt(sum((rat(:,i)-rat(:,j))**2))
    end do; end do
    Fa = Hc; Fb = Hc; e_scf = 0.0_dp
    do it = 1, 400
      Scp = Fa     ! placeholder, reset below
      Scp = S; call gen_eig(Fa, Scp, n, eps, Ca)
      Scp = S; call gen_eig(Fb, Scp, n, epb, Cb)
      Da = 0.0_dp; Db = 0.0_dp
      do mu=1,n; do nu=1,n
        do i=1,na
          Da(mu,nu) = Da(mu,nu) + Ca(mu,i)*Ca(nu,i)
        end do
        do i=1,nb
          Db(mu,nu) = Db(mu,nu) + Cb(mu,i)*Cb(nu,i)
        end do
      end do; end do
      Ja=0; Ka=0; Jb=0; Kb=0
      do mu=1,n; do nu=1,n; do la=1,n; do si=1,n
        Ja(mu,nu) = Ja(mu,nu) + (Da(la,si)+Db(la,si))*eri_c(mu,nu,la,si)
        Ka(mu,nu) = Ka(mu,nu) + Da(la,si)*eri_c(mu,la,nu,si)
        Kb(mu,nu) = Kb(mu,nu) + Db(la,si)*eri_c(mu,la,nu,si)
      end do; end do; end do; end do
      Fa = Hc + Ja - Ka
      Fb = Hc + Ja - Kb
      eold = e_scf
      e_scf = 0.0_dp
      do mu=1,n; do nu=1,n
        e_scf = e_scf + 0.5_dp*((Da(mu,nu)+Db(mu,nu))*Hc(mu,nu) &
                              + Da(mu,nu)*Fa(mu,nu) + Db(mu,nu)*Fb(mu,nu))
      end do; end do
      if (abs(e_scf-eold) < 1.0e-11_dp .and. it > 1) exit
    end do
    e_scf = e_scf + enuc
    Cmo = Ca                ! alpha orbitals are the reference set
  end subroutine rohf_highspin

  subroutine gen_eig(F, S, n, w, C)
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
  end subroutine gen_eig

  subroutine ao2mo_1e(h, C, n, hmo)
    integer, intent(in) :: n
    real(dp), intent(in) :: h(n,n), C(n,n)
    real(dp), intent(out) :: hmo(n,n)
    hmo = matmul(transpose(C), matmul(h, C))
  end subroutine ao2mo_1e

  subroutine ao2mo_2e(eri, C, n, emo)
    integer, intent(in) :: n
    real(dp), intent(in) :: eri(n,n,n,n), C(n,n)
    real(dp), intent(out) :: emo(n,n,n,n)
    real(dp) :: a1(n,n,n,n), a2(n,n,n,n)
    integer :: p,q,r,s,mu
    a1 = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      a1(p,q,r,s) = a1(p,q,r,s) + C(mu,p)*eri(mu,q,r,s)
    end do; end do; end do; end do; end do
    a2 = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      a2(p,q,r,s) = a2(p,q,r,s) + C(mu,q)*a1(p,mu,r,s)
    end do; end do; end do; end do; end do
    a1 = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      a1(p,q,r,s) = a1(p,q,r,s) + C(mu,r)*a2(p,q,mu,s)
    end do; end do; end do; end do; end do
    emo = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      emo(p,q,r,s) = emo(p,q,r,s) + C(mu,s)*a1(p,q,r,mu)
    end do; end do; end do; end do; end do
  end subroutine ao2mo_2e

  !> chemist ERI (ij|kl) from the physicist geminal-tensor layout M(i,k,j,l).
  subroutine eri_chem(M, n, eri_c)
    integer, intent(in) :: n
    real(dp), intent(in) :: M(n,n,n,n)
    real(dp), intent(out) :: eri_c(n,n,n,n)
    integer :: i,j,k,l
    do i=1,n; do j=1,n; do k=1,n; do l=1,n
      eri_c(i,j,k,l) = M(i,k,j,l)
    end do; end do; end do; end do
  end subroutine eri_chem

  !> MO geminal tensor Gmo(i,j,a,b) = <ij|f12|ab> for the Slater factor
  !> f12 = exp(-gamma r12) (STG-6G), in physicist layout (e1=(i,a), e2=(j,b)).
  subroutine geminal_mo(n, npr, exs, cos_, cns, gamma, C, Gmo)
    integer,  intent(in)  :: n, npr(n)
    real(dp), intent(in)  :: exs(:,:), cos_(:,:), cns(3,n), gamma, C(n,n)
    real(dp), intent(out) :: Gmo(n,n,n,n)
    real(dp) :: Gao(n,n,n,n), Mk(n,n,n,n), cc(6), omg(6)
    integer  :: ng, k
    call stg_ng(gamma, cc, omg, ng)
    Gao = 0.0_dp
    do k = 1, ng
      call ptc_s_ao_tensor(n, npr, exs, cos_, cns, PTC_OP_GEMINAL, omg(k), Mk)
      Gao = Gao + cc(k)*Mk
    end do
    call ao2mo_2e(Gao, C, n, Gmo)
  end subroutine geminal_mo

  ! ---- determinant engine (spin-orbital bitmasks, P = 2*orb + spin) ----
  subroutine build_dets(norb, na, nb, dets, dim, hfidx)
    integer, intent(in) :: norb, na, nb
    integer, allocatable, intent(out) :: dets(:)
    integer, intent(out) :: dim, hfidx
    integer :: da, db, cnt, hf, naa, nbb
    integer, allocatable :: al(:), bl(:)
    call combos(norb, na, al, naa)
    call combos(norb, nb, bl, nbb)
    dim = naa*nbb; allocate(dets(dim)); cnt = 0
    do da=1,naa; do db=1,nbb
      cnt = cnt+1
      dets(cnt) = ior(spread_spin(al(da),0), spread_spin(bl(db),1))
    end do; end do
    hf = ior(spread_spin(2**na-1,0), spread_spin(2**nb-1,1))
    hfidx = 0
    do cnt=1,dim
      if (dets(cnt)==hf) hfidx = cnt
    end do
  end subroutine build_dets

  integer function spread_spin(omask, spin) result(p)
    integer, intent(in) :: omask, spin
    integer :: o
    p = 0
    do o=0,30
      if (btest(omask,o)) p = ibset(p, 2*o+spin)
    end do
  end function spread_spin

  subroutine combos(norb, k, list, ncomb)
    integer, intent(in) :: norb, k
    integer, allocatable, intent(out) :: list(:)
    integer, intent(out) :: ncomb
    integer :: m, cnt
    integer, allocatable :: tmp(:)
    allocate(tmp(2**norb)); cnt=0
    do m=0,2**norb-1
      if (popcnt_(m)==k) then
        cnt=cnt+1; tmp(cnt)=m
      end if
    end do
    ncomb=cnt; allocate(list(cnt)); list=tmp(1:cnt)
  end subroutine combos

  integer function popcnt_(m) result(c)
    integer, intent(in) :: m
    integer :: i
    c=0
    do i=0,30
      if (btest(m,i)) c=c+1
    end do
  end function popcnt_

  integer function parity_below(d, x) result(s)
    integer, intent(in) :: d, x
    integer :: i, c
    c=0
    do i=0,x-1
      if (btest(d,i)) c=c+1
    end do
    s = 1 - 2*mod(c,2)
  end function parity_below

  integer function ann(det, p, nd) result(sgn)
    integer, intent(in) :: det, p
    integer, intent(out) :: nd
    sgn=0; nd=det
    if (.not. btest(det,p)) return
    sgn=parity_below(det,p); nd=ibclr(det,p)
  end function ann

  integer function cre(det, p, nd) result(sgn)
    integer, intent(in) :: det, p
    integer, intent(out) :: nd
    sgn=0; nd=det
    if (btest(det,p)) return
    sgn=parity_below(det,p); nd=ibset(det,p)
  end function cre

  integer function findidx(dets, dim, d) result(idx)
    integer, intent(in) :: dets(:), dim, d
    integer :: i
    idx=0
    do i=1,dim
      if (dets(i)==d) then
        idx=i; return
      end if
    end do
  end function findidx

  subroutine build_fci_H(h1, eri, ecore, norb, dets, dim, H)
    integer, intent(in) :: norb, dets(:), dim
    real(dp), intent(in) :: h1(norb,norb), eri(norb,norb,norb,norb), ecore
    real(dp), intent(out) :: H(dim,dim)
    integer :: col,p,q,r,u,s1,s2,d1,d2,dd,g,g2,jj,det
    real(dp) :: v
    H=0.0_dp
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(det,p,q,r,u,s1,s2,d1,d2,dd,g,g2,jj,v)
    do col=1,dim
      det=dets(col)
      do p=1,norb; do q=1,norb
        if (abs(h1(p,q))<1e-14_dp) cycle
        do s1=0,1
          g=ann(det,2*(q-1)+s1,d1); if(g==0) cycle
          g2=cre(d1,2*(p-1)+s1,d2); if(g2==0) cycle
          jj=findidx(dets,dim,d2); if(jj>0) H(jj,col)=H(jj,col)+h1(p,q)*g*g2
        end do
      end do; end do
      do p=1,norb; do q=1,norb; do r=1,norb; do u=1,norb
        v=eri(p,q,r,u); if(abs(v)<1e-14_dp) cycle
        do s1=0,1; do s2=0,1
          g=ann(det,2*(q-1)+s1,d1); if(g==0) cycle
          g2=ann(d1,2*(u-1)+s2,dd); if(g2==0) cycle
          g=g*g2
          g2=cre(dd,2*(r-1)+s2,d1); if(g2==0) cycle
          g=g*g2
          g2=cre(d1,2*(p-1)+s1,d2); if(g2==0) cycle
          g=g*g2
          jj=findidx(dets,dim,d2); if(jj>0) H(jj,col)=H(jj,col)+0.5_dp*v*g
        end do; end do
      end do; end do; end do; end do
      H(col,col)=H(col,col)+ecore
    end do
  end subroutine build_fci_H

  !> Spin-resolved geminal F12 doubles operator T = sum t_ij^ab E_ai E_bj over
  !> active (i,j) -> external (a,b), amplitude = cusp(spin) * Gmo(i,j,a,b), with
  !> cusp = 1/2 for opposite-spin, 1/4 for same-spin (the fixed-amplitude F12
  !> cusp conditions). nact active orbitals (1-based list iact), next external.
  !> gamma: Slater exponent. The F12 correlation factor is f12 = -(1/gamma)
  !> exp(-gamma r12) (negative -> correlation lowers the energy); its cusp slope
  !> is 1 at coalescence, scaled to 1/2 (opposite-spin) or 1/4 (same-spin).
  subroutine build_geminal_T2(Gmo, norb, iact, nact, iext, next, scale, gamma, dets, dim, T2op)
    integer,  intent(in)  :: norb, iact(:), nact, iext(:), next, dets(:), dim
    real(dp), intent(in)  :: Gmo(norb,norb,norb,norb), scale, gamma
    real(dp), intent(out) :: T2op(dim,dim)
    integer :: col,ii,jj2,aa,bb,i,j,a,b,det,d1,d1b,d2,d2f,jx,sa,sb,g1,g2,g3,g4
    real(dp) :: amp, cusp
    T2op=0.0_dp
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp private(det,ii,jj2,aa,bb,i,j,a,b,amp,cusp,sa,sb,g1,g2,g3,g4,d1,d1b,d2,d2f,jx)
    do col=1,dim
      det=dets(col)
      do ii=1,nact; do jj2=1,nact; do aa=1,next; do bb=1,next
        i=iact(ii); j=iact(jj2); a=iext(aa); b=iext(bb)
        amp = scale*Gmo(i,j,a,b)
        if (abs(amp)<1e-14_dp) cycle
        do sb=0,1
          g1=ann(det,2*(j-1)+sb,d1); if(g1==0) cycle
          g2=cre(d1,2*(b-1)+sb,d1b); if(g2==0) cycle
          do sa=0,1
            ! F12 factor -(1/gamma) with cusp slope 1/2 (opp-spin) or 1/4 (same)
            cusp = -merge(0.25_dp, 0.5_dp, sa==sb)/gamma
            g3=ann(d1b,2*(i-1)+sa,d2); if(g3==0) cycle
            g4=cre(d2,2*(a-1)+sa,d2f); if(g4==0) cycle
            jx=findidx(dets,dim,d2f)
            if(jx>0) T2op(jx,col)=T2op(jx,col)+cusp*amp*g1*g2*g3*g4
          end do
        end do
      end do; end do; end do; end do
    end do
  end subroutine build_geminal_T2

  !> Conventional doubles operator: T = sum t_ij^ab E_ai E_bj over active (i,j) ->
  !> external (a,b), with MP2 amplitudes t_ij^ab = (ia|jb)/(e_i+e_j-e_a-e_b). This
  !> is the bulk (in-basis) dynamic correlation; for a 2-electron system the
  !> doubles ARE the full correlation, so bare MRSF-CIS + this -> FCI. The geminal
  !> F12 piece (build_geminal_T2) is added on top as the explicit-correlation cusp.
  subroutine build_conv_T2(eri_mo, eps, norb, iact, nact, iext, next, scale, dets, dim, T2op)
    integer,  intent(in)  :: norb, iact(:), nact, iext(:), next, dets(:), dim
    real(dp), intent(in)  :: eri_mo(norb,norb,norb,norb), eps(norb), scale
    real(dp), intent(out) :: T2op(dim,dim)
    integer :: col,ii,jj2,aa,bb,i,j,a,b,det,d1,d1b,d2,d2f,jx,sa,sb,g1,g2,g3,g4
    real(dp) :: amp, den
    T2op=0.0_dp
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp private(det,ii,jj2,aa,bb,i,j,a,b,amp,den,sa,sb,g1,g2,g3,g4,d1,d1b,d2,d2f,jx)
    do col=1,dim
      det=dets(col)
      do ii=1,nact; do jj2=1,nact; do aa=1,next; do bb=1,next
        i=iact(ii); j=iact(jj2); a=iext(aa); b=iext(bb)
        den = eps(i)+eps(j)-eps(a)-eps(b)
        if (abs(den) < 1e-10_dp) cycle
        amp = scale*eri_mo(i,a,j,b)/den          ! MP2 amplitude (ia|jb)/Delta
        if (abs(amp)<1e-14_dp) cycle
        do sb=0,1
          g1=ann(det,2*(j-1)+sb,d1); if(g1==0) cycle
          g2=cre(d1,2*(b-1)+sb,d1b); if(g2==0) cycle
          do sa=0,1
            g3=ann(d1b,2*(i-1)+sa,d2); if(g3==0) cycle
            g4=cre(d2,2*(a-1)+sa,d2f); if(g4==0) cycle
            jx=findidx(dets,dim,d2f)
            if(jx>0) T2op(jx,col)=T2op(jx,col)+amp*g1*g2*g3*g4
          end do
        end do
      end do; end do; end do; end do
    end do
  end subroutine build_conv_T2

  ! kept for reference: MP2-style spatial T2 (not used in the geminal path)
  subroutine build_msz_T2(t2, nocc, nvir, norb, dets, dim, T2op)
    integer, intent(in) :: nocc, nvir, norb, dets(:), dim
    real(dp), intent(in) :: t2(nocc,nocc,nvir,nvir)
    real(dp), intent(out) :: T2op(dim,dim)
    integer :: col,i,j,a,b,det,d1,d1b,d2,d2f,jj,sa,sb,g1,g2,g3,g4
    real(dp) :: amp
    T2op=0.0_dp
    do col=1,dim
      det=dets(col)
      do i=1,nocc; do j=1,nocc; do a=1,nvir; do b=1,nvir
        amp=0.5_dp*t2(i,j,a,b); if(abs(amp)<1e-14_dp) cycle
        do sb=0,1
          g1=ann(det,2*(j-1)+sb,d1); if(g1==0) cycle
          g2=cre(d1,2*(nocc+b-1)+sb,d1b); if(g2==0) cycle
          do sa=0,1
            g3=ann(d1b,2*(i-1)+sa,d2); if(g3==0) cycle
            g4=cre(d2,2*(nocc+a-1)+sa,d2f); if(g4==0) cycle
            jj=findidx(dets,dim,d2f)
            if(jj>0) T2op(jj,col)=T2op(jj,col)+amp*g1*g2*g3*g4
          end do
        end do
      end do; end do; end do; end do
    end do
  end subroutine build_msz_T2

  subroutine build_s2(norb, dets, dim, S2)
    integer, intent(in) :: norb, dets(:), dim
    real(dp), intent(out) :: S2(dim,dim)
    integer, allocatable :: pd(:)
    integer :: np2, p, g, g2, d1, d2, col, row, a, b, cnt
    integer, allocatable :: tmp(:)
    real(dp), allocatable :: Sp(:,:)
    allocate(tmp(norb*norb)); cnt=0
    do a=0,norb-1; do b=a+1,norb-1
      cnt=cnt+1; tmp(cnt)=ior(ibset(0,2*a),ibset(0,2*b))
    end do; end do
    np2=cnt; allocate(pd(np2)); pd=tmp(1:np2)
    allocate(Sp(np2,dim)); Sp=0.0_dp
    do col=1,dim
      do p=1,norb
        g=ann(dets(col),2*(p-1)+1,d1); if(g==0) cycle
        g2=cre(d1,2*(p-1)+0,d2); if(g2==0) cycle
        row=findidx(pd,np2,d2)
        if(row>0) Sp(row,col)=Sp(row,col)+g*g2
      end do
    end do
    S2 = matmul(transpose(Sp), Sp)
    deallocate(Sp,pd,tmp)
  end subroutine build_s2

  subroutine cas22_compact(dets, dim, norb, iact, nact, cas, nc)
    integer, intent(in) :: dets(:), dim, norb, iact(:), nact
    integer, allocatable, intent(out) :: cas(:)
    integer, intent(out) :: nc
    integer :: i, p, det, cnt, ok, q, inact
    integer, allocatable :: tmp(:)
    allocate(tmp(dim)); cnt=0
    do i=1,dim
      det=dets(i); ok=1
      do p=0,norb-1
        if (btest(det,2*p) .or. btest(det,2*p+1)) then
          inact=0
          do q=1,nact
            if (iact(q)==p+1) inact=1
          end do
          if (inact==0) ok=0
        end if
      end do
      if (ok==1) then
        cnt=cnt+1; tmp(cnt)=i
      end if
    end do
    nc=cnt; allocate(cas(nc)); cas=tmp(1:nc)
  end subroutine cas22_compact

  subroutine sym_eig_vec(A, n, w)
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n)
    real(dp), intent(out) :: w(n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsyev('V','U',n,A,n,w,wq,-1,info)
    lw=int(wq(1)); allocate(wk(lw))
    call dsyev('V','U',n,A,n,w,wk,lw,info)
    deallocate(wk)
  end subroutine sym_eig_vec

  subroutine expm_nilpotent(scale, A, n, E)
    integer, intent(in) :: n
    real(dp), intent(in) :: scale, A(n,n)
    real(dp), intent(out) :: E(n,n)
    real(dp) :: term(n,n)
    integer :: k
    E=0.0_dp; term=0.0_dp
    do k=1,n
      E(k,k)=1.0_dp; term(k,k)=1.0_dp
    end do
    do k=1,18
      term=matmul(term,scale*A)/real(k,dp)
      E=E+term
      if (maxval(abs(term))<1e-15_dp) exit
    end do
  end subroutine expm_nilpotent

end module tc_geminal_engine
