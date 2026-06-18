!> Native (pyscf-free) H2 dissociation curves: bare MRSF-CIS vs pTC-MRSF-CIS vs
!> FCI, in an all-s 6-311G basis (so the native s-only integral engine applies).
!> For each bond length the (2,2) frontier spin-flip space gives bare MRSF-CIS
!> (CASCI(2,2)); the MP2-T2 transcorrelated H_bar downfolded into the same space
!> gives pTC-MRSF-CIS; full FCI (6 orbitals) is the benchmark. The point is that
!> bare MRSF-CIS misses the dynamic correlation from the 4 external orbitals while
!> pTC recovers it -- three DISTINCT curves per state, unlike a FCI-only plot.
!>
!> Writes pes_h2_native.dat: R  S0_fci T1_fci S1_fci  S0_bare T1_bare S1_bare
!>                              S0_ptc T1_ptc S1_ptc.
!>
!> Build (OpenMP):
!>   gfortran -O2 -fopenmp source/precision.F90 source/modules/ptc_geminal.F90 \
!>     source/modules/tdhf_mrsf_ptc.F90 \
!>     tests/ptc_mrsf/prototype/tc_h2_pes_native.F90 -llapack -lblas -o /tmp/ptc_pes && /tmp/ptc_pes
program tc_h2_pes_native
  use precision, only: dp
  use ptc_geminal
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none

  integer,  parameter :: NSH = 6, NOCC = 1, MP = 3   ! H2/6-311G: 3 shells/atom
  integer  :: npr(NSH), i, nR, ir
  real(dp) :: exs(MP,NSH), cos_(MP,NSH), cns(3,NSH), zat(2), rat(3,2)
  real(dp) :: R, e_fci(3), e_bare(3), e_ptc(3)
  real(dp), allocatable :: Rs(:)
  integer  :: u

  ! 6-311G hydrogen: one 3-primitive contracted s + two single s per atom.
  call set_h_basis()

  ! R grid (bohr): dense near equilibrium, out to dissociation
  nR = 28
  allocate(Rs(nR))
  do i = 1, nR
    Rs(i) = 1.0_dp + (8.0_dp - 1.0_dp)*real(i-1,dp)/real(nR-1,dp)
  end do

  open(newunit=u, file='pes_h2_native.dat', status='replace')
  write(u,'(a)') '# R   S0_fci T1_fci S1_fci   S0_bare T1_bare S1_bare   S0_ptc T1_ptc S1_ptc'
  do ir = 1, nR
    R = Rs(ir)
    rat(:,1) = [0.0_dp,0.0_dp,0.0_dp]; rat(:,2) = [0.0_dp,0.0_dp,R]
    cns(:,1:3) = spread(rat(:,1), 2, 3)   ! shells 1-3 on atom 1
    cns(:,4:6) = spread(rat(:,2), 2, 3)   ! shells 4-6 on atom 2
    call one_point(R, e_fci, e_bare, e_ptc)
    write(u,'(f8.3, 3f12.6, 3x, 3f12.6, 3x, 3f12.6)') R, e_fci, e_bare, e_ptc
    write(*,'(a,f6.3,a,3f11.6,a,3f11.6,a,3f11.6)') 'R=',R, &
      '  FCI(S0,T1,S1)=', e_fci, '  bare=', e_bare, '  pTC=', e_ptc
  end do
  close(u)
  write(*,'(a)') 'wrote pes_h2_native.dat'

contains

  subroutine set_h_basis()
    integer :: s
    do s = 1, NSH
      npr(s) = 1
    end do
    ! shell 1 & 4: contracted 3-prim
    do s = 1, NSH, 3
      npr(s) = 3
      exs(:,s)  = [33.8650_dp, 5.094790_dp, 1.158790_dp]
      cos_(:,s) = [0.0254938_dp, 0.190373_dp, 0.852161_dp]
    end do
    ! shells 2,5: single s 0.325840 ; shells 3,6: single s 0.102741
    do s = 2, NSH, 3
      exs(1,s) = 0.325840_dp; cos_(1,s) = 1.0_dp
    end do
    do s = 3, NSH, 3
      exs(1,s) = 0.102741_dp; cos_(1,s) = 1.0_dp
    end do
  end subroutine set_h_basis

  subroutine one_point(R, efci, ebare, eptc)
    real(dp), intent(in)  :: R
    real(dp), intent(out) :: efci(3), ebare(3), eptc(3)
    integer  :: n, nvir
    real(dp) :: Cmo(NSH,NSH), eps(NSH), e_rhf, enuc, h1ao(NSH,NSH), eri_c(NSH,NSH,NSH,NSH)
    real(dp) :: h1mo(NSH,NSH), eri_mo(NSH,NSH,NSH,NSH)
    real(dp), allocatable :: t2(:,:,:,:)
    integer, allocatable :: dets(:), cas(:)
    integer  :: dim, hfidx, nc, i, j, a, b, k
    real(dp), allocatable :: H(:,:), T2op(:,:), Hbar(:,:), Em(:,:), Ep(:,:)
    real(dp), allocatable :: S2(:,:), Hc(:,:), Hbc(:,:), S2c(:,:)
    real(dp), allocatable :: ee(:), vr(:,:), vl(:,:)
    real(dp) :: denom, gia, maxim
    integer  :: ncx, ierr
    n = NSH; nvir = n - NOCC
    allocate(t2(NOCC,NOCC,nvir,nvir))
    call rhf(n, NOCC, npr, exs, cos_, cns, 2, zat, rat, Cmo, eps, e_rhf, enuc, h1ao, eri_c)
    call ao2mo_1e(h1ao, Cmo, n, h1mo)
    call ao2mo_2e(eri_c, Cmo, n, eri_mo)
    ! MP2 amplitudes
    do i = 1, NOCC; do j = 1, NOCC; do a = 1, nvir; do b = 1, nvir
      denom = eps(i)+eps(j)-eps(NOCC+a)-eps(NOCC+b)
      gia = eri_mo(i, NOCC+a, j, NOCC+b)
      t2(i,j,a,b) = gia/denom
    end do; end do; end do; end do
    call build_dets(n, NOCC, dets, dim, hfidx)
    allocate(H(dim,dim), T2op(dim,dim), Hbar(dim,dim), Em(dim,dim), Ep(dim,dim), S2(dim,dim))
    call build_fci_H(h1mo, eri_mo, enuc, n, dets, dim, H)
    call build_s2(n, dets, dim, S2)
    call build_T2op(t2, NOCC, nvir, n, dets, dim, T2op)
    call expm_nilpotent(-1.0_dp, T2op, dim, Em)
    call expm_nilpotent( 1.0_dp, T2op, dim, Ep)
    Hbar = matmul(Em, matmul(H, Ep))
    ! FCI states (full space)
    call states_from(H, S2, dim, .false., efci)
    ! (2,2) frontier CAS: orbitals {NOCC, NOCC+1}
    call cas22_compact(dets, dim, n, NOCC, cas, nc)
    allocate(Hc(nc,nc), Hbc(nc,nc), S2c(nc,nc))
    do i = 1, nc; do j = 1, nc
      Hc(i,j)  = H(cas(i), cas(j))
      Hbc(i,j) = Hbar(cas(i), cas(j))
      S2c(i,j) = S2(cas(i), cas(j))
    end do; end do
    call states_from(Hc, S2c, nc, .false., ebare)
    call states_from(Hbc, S2c, nc, .true.,  eptc)
    deallocate(t2, dets, H, T2op, Hbar, Em, Ep, S2, cas, Hc, Hbc, S2c)
  end subroutine one_point

  !> Extract [S0, T1, S1] = lowest singlet, lowest triplet, 2nd-lowest singlet,
  !> labeling by <S^2> (singlet < 1, triplet >= 1). nonherm=.true. uses the
  !> non-Hermitian solver and biorthonormal <S^2>; else symmetric.
  subroutine states_from(Hm, S2m, m, nonherm, out)
    integer,  intent(in)  :: m
    real(dp), intent(in)  :: Hm(m,m), S2m(m,m)
    logical,  intent(in)  :: nonherm
    real(dp), intent(out) :: out(3)
    real(dp) :: w(m), V(m,m), vrr(m,m), vll(m,m), ss(m), mi, den
    integer  :: k, ncx, ierr, order(m), i, j, tmp
    integer  :: ns0, nt0, ns1
    if (nonherm) then
      call tc_nonsym_tda_eig(Hm, m, w, vrr, vll, mi, ncx, ierr)
      do k = 1, m
        den = dot_product(vll(:,k), vrr(:,k))
        ss(k) = dot_product(vll(:,k), matmul(S2m, vrr(:,k)))/den
      end do
    else
      V = Hm
      call sym_eig_vec(V, m, w)
      do k = 1, m
        ss(k) = dot_product(V(:,k), matmul(S2m, V(:,k)))
      end do
    end if
    ! sort by energy
    do i = 1, m
      order(i) = i
    end do
    do i = 1, m-1
      do j = i+1, m
        if (w(order(j)) < w(order(i))) then
          tmp = order(i); order(i) = order(j); order(j) = tmp
        end if
      end do
    end do
    out = 0.0_dp
    ns0 = 0; nt0 = 0; ns1 = 0
    out(1) = huge(1.0_dp); out(2) = huge(1.0_dp); out(3) = huge(1.0_dp)
    do i = 1, m
      k = order(i)
      if (ss(k) < 1.0_dp) then       ! singlet
        if (ns0 == 0) then
          out(1) = w(k); ns0 = 1
        else if (ns1 == 0) then
          out(3) = w(k); ns1 = 1
        end if
      else                            ! triplet
        if (nt0 == 0) then
          out(2) = w(k); nt0 = 1
        end if
      end if
    end do
  end subroutine states_from

  ! ===== engine (ported from tc_mrsf_native_test.F90, validated) =====

  subroutine rhf(n, nocc, npr, exs, cos_, cns, nat, zat, rat, Cmo, eps, e_rhf, enuc, h1ao, eri_c)
    integer,  intent(in)  :: n, nocc, npr(n), nat
    real(dp), intent(in)  :: exs(:,:), cos_(:,:), cns(3,n), zat(:), rat(:,:)
    real(dp), intent(out) :: Cmo(n,n), eps(n), e_rhf, enuc, h1ao(n,n), eri_c(n,n,n,n)
    real(dp) :: S(n,n), T(n,n), V(n,n), Hcore(n,n), F(n,n), P(n,n), G(n,n)
    real(dp) :: M(n,n,n,n), Scopy(n,n), eold
    integer  :: i, j, k, l, it, mu, nu, la, si
    real(dp) :: zatt(nat), ratt(3,nat)
    zatt = 1.0_dp
    ratt = rat(:,1:nat)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_OVERLAP, nat, zatt, ratt, S)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_KINETIC, nat, zatt, ratt, T)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_NUCLEAR, nat, zatt, ratt, V)
    Hcore = T + V; h1ao = Hcore
    call ptc_s_ao_tensor(n, npr, exs, cos_, cns, PTC_OP_ERI, 0.0_dp, M)
    do i=1,n; do j=1,n; do k=1,n; do l=1,n
      eri_c(i,j,k,l) = M(i,k,j,l)
    end do; end do; end do; end do
    enuc = 0.0_dp
    do i = 1, nat; do j = i+1, nat
      enuc = enuc + zatt(i)*zatt(j)/sqrt(sum((ratt(:,i)-ratt(:,j))**2))
    end do; end do
    F = Hcore; e_rhf = 0.0_dp
    do it = 1, 300
      Scopy = S
      call gen_eig(F, Scopy, n, eps, Cmo)
      P = 0.0_dp
      do mu=1,n; do nu=1,n; do i=1,nocc
        P(mu,nu) = P(mu,nu) + 2.0_dp*Cmo(mu,i)*Cmo(nu,i)
      end do; end do; end do
      G = 0.0_dp
      do mu=1,n; do nu=1,n; do la=1,n; do si=1,n
        G(mu,nu) = G(mu,nu) + P(la,si)*(eri_c(mu,nu,la,si) - 0.5_dp*eri_c(mu,la,nu,si))
      end do; end do; end do; end do
      F = Hcore + G
      eold = e_rhf; e_rhf = 0.0_dp
      do mu=1,n; do nu=1,n
        e_rhf = e_rhf + 0.5_dp*P(mu,nu)*(Hcore(mu,nu)+F(mu,nu))
      end do; end do
      if (abs(e_rhf-eold) < 1.0e-11_dp .and. it > 1) exit
    end do
    e_rhf = e_rhf + enuc
  end subroutine rhf

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
    real(dp) :: t1(n,n,n,n), t2t(n,n,n,n)
    integer :: p,q,r,s,mu
    t1 = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      t1(p,q,r,s) = t1(p,q,r,s) + C(mu,p)*eri(mu,q,r,s)
    end do; end do; end do; end do; end do
    t2t = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      t2t(p,q,r,s) = t2t(p,q,r,s) + C(mu,q)*t1(p,mu,r,s)
    end do; end do; end do; end do; end do
    t1 = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      t1(p,q,r,s) = t1(p,q,r,s) + C(mu,r)*t2t(p,q,mu,s)
    end do; end do; end do; end do; end do
    emo = 0.0_dp
    do p=1,n; do q=1,n; do r=1,n; do s=1,n; do mu=1,n
      emo(p,q,r,s) = emo(p,q,r,s) + C(mu,s)*t1(p,q,r,mu)
    end do; end do; end do; end do; end do
  end subroutine ao2mo_2e

  subroutine build_dets(norb, nocc, dets, dim, hfidx)
    integer, intent(in) :: norb, nocc
    integer, allocatable, intent(out) :: dets(:)
    integer, intent(out) :: dim, hfidx
    integer :: amask, bmask, da, db, cnt, hf, na
    integer, allocatable :: alist(:)
    call combos(norb, nocc, alist, na)
    dim = na*na; allocate(dets(dim)); cnt = 0
    do da = 1, na; do db = 1, na
      cnt = cnt + 1
      amask = spread_spin(alist(da), 0); bmask = spread_spin(alist(db), 1)
      dets(cnt) = ior(amask, bmask)
    end do; end do
    hf = ior(spread_spin(2**nocc-1, 0), spread_spin(2**nocc-1, 1))
    hfidx = 0
    do cnt = 1, dim
      if (dets(cnt) == hf) hfidx = cnt
    end do
  end subroutine build_dets

  integer function spread_spin(omask, spin) result(p)
    integer, intent(in) :: omask, spin
    integer :: o
    p = 0
    do o = 0, 30
      if (btest(omask,o)) p = ibset(p, 2*o+spin)
    end do
  end function spread_spin

  subroutine combos(norb, k, list, ncomb)
    integer, intent(in) :: norb, k
    integer, allocatable, intent(out) :: list(:)
    integer, intent(out) :: ncomb
    integer :: m, cnt
    integer, allocatable :: tmp(:)
    allocate(tmp(2**norb)); cnt = 0
    do m = 0, 2**norb-1
      if (popcnt_(m) == k) then
        cnt = cnt + 1; tmp(cnt) = m
      end if
    end do
    ncomb = cnt; allocate(list(cnt)); list = tmp(1:cnt)
  end subroutine combos

  integer function popcnt_(m) result(c)
    integer, intent(in) :: m
    integer :: i
    c = 0
    do i = 0, 30
      if (btest(m,i)) c = c + 1
    end do
  end function popcnt_

  integer function parity_below(d, x) result(s)
    integer, intent(in) :: d, x
    integer :: i, c
    c = 0
    do i = 0, x-1
      if (btest(d,i)) c = c + 1
    end do
    s = 1 - 2*mod(c,2)
  end function parity_below

  integer function ann(det, p, newdet) result(sgn)
    integer, intent(in) :: det, p
    integer, intent(out) :: newdet
    sgn = 0; newdet = det
    if (.not. btest(det,p)) return
    sgn = parity_below(det,p); newdet = ibclr(det,p)
  end function ann

  integer function cre(det, p, newdet) result(sgn)
    integer, intent(in) :: det, p
    integer, intent(out) :: newdet
    sgn = 0; newdet = det
    if (btest(det,p)) return
    sgn = parity_below(det,p); newdet = ibset(det,p)
  end function cre

  integer function findidx(dets, dim, d) result(idx)
    integer, intent(in) :: dets(:), dim, d
    integer :: i
    idx = 0
    do i = 1, dim
      if (dets(i) == d) then
        idx = i; return
      end if
    end do
  end function findidx

  subroutine build_fci_H(h1, eri, ecore, norb, dets, dim, H)
    integer, intent(in) :: norb, dets(:), dim
    real(dp), intent(in) :: h1(norb,norb), eri(norb,norb,norb,norb), ecore
    real(dp), intent(out) :: H(dim,dim)
    integer :: col, p, q, r, u, s1, s2, d1, d2, dd, g, g2, jj, det
    real(dp) :: v
    H = 0.0_dp
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(det,p,q,r,u,s1,s2,d1,d2,dd,g,g2,jj,v)
    do col = 1, dim
      det = dets(col)
      do p = 1, norb; do q = 1, norb
        if (abs(h1(p,q)) < 1.0e-14_dp) cycle
        do s1 = 0, 1
          g = ann(det, 2*(q-1)+s1, d1); if (g==0) cycle
          g2 = cre(d1, 2*(p-1)+s1, d2); if (g2==0) cycle
          jj = findidx(dets, dim, d2)
          if (jj>0) H(jj,col) = H(jj,col) + h1(p,q)*g*g2
        end do
      end do; end do
      do p = 1, norb; do q = 1, norb; do r = 1, norb; do u = 1, norb
        v = eri(p,q,r,u); if (abs(v) < 1.0e-14_dp) cycle
        do s1 = 0, 1; do s2 = 0, 1
          g = ann(det, 2*(q-1)+s1, d1); if (g==0) cycle
          g2 = ann(d1, 2*(u-1)+s2, dd); if (g2==0) cycle
          g = g*g2
          g2 = cre(dd, 2*(r-1)+s2, d1); if (g2==0) cycle
          g = g*g2
          g2 = cre(d1, 2*(p-1)+s1, d2); if (g2==0) cycle
          g = g*g2
          jj = findidx(dets, dim, d2)
          if (jj>0) H(jj,col) = H(jj,col) + 0.5_dp*v*g
        end do; end do
      end do; end do; end do; end do
      H(col,col) = H(col,col) + ecore
    end do
  end subroutine build_fci_H

  subroutine build_T2op(t2, nocc, nvir, norb, dets, dim, T2op)
    integer, intent(in) :: nocc, nvir, norb, dets(:), dim
    real(dp), intent(in) :: t2(nocc,nocc,nvir,nvir)
    real(dp), intent(out) :: T2op(dim,dim)
    integer :: col, i, j, a, b, det, d1, d1b, d2, d2f, jj, sa, sb, g1, g2, g3, g4
    real(dp) :: amp
    T2op = 0.0_dp
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(det,i,j,a,b,amp,sa,sb,g1,g2,g3,g4,d1,d1b,d2,d2f,jj)
    do col = 1, dim
      det = dets(col)
      do i=1,nocc; do j=1,nocc; do a=1,nvir; do b=1,nvir
        amp = 0.5_dp*t2(i,j,a,b); if (abs(amp) < 1.0e-14_dp) cycle
        do sb = 0, 1
          g1 = ann(det, 2*(j-1)+sb, d1); if (g1==0) cycle
          g2 = cre(d1, 2*(nocc+b-1)+sb, d1b); if (g2==0) cycle
          do sa = 0, 1
            g3 = ann(d1b, 2*(i-1)+sa, d2); if (g3==0) cycle
            g4 = cre(d2, 2*(nocc+a-1)+sa, d2f); if (g4==0) cycle
            jj = findidx(dets, dim, d2f)
            if (jj>0) T2op(jj,col) = T2op(jj,col) + amp*g1*g2*g3*g4
          end do
        end do
      end do; end do; end do; end do
    end do
  end subroutine build_T2op

  !> <S^2> matrix in the Ms=0 determinant basis via S+ (maps Ms=0 -> Ms=+1):
  !> S2 = (S+)^dag S+  (= <S^2> at Ms=0).
  subroutine build_s2(norb, dets, dim, S2)
    integer,  intent(in)  :: norb, dets(:), dim
    real(dp), intent(out) :: S2(dim,dim)
    integer, allocatable :: pdets(:)
    integer :: np2, i, j, p, g, g2, d1, d2, col, row
    real(dp), allocatable :: Sp(:,:)
    ! Ms=+1 basis: (nocc+1) alpha, (nocc-1) beta ; here nocc=1 -> 2 alpha, 0 beta
    call build_msplus(norb, pdets, np2)
    allocate(Sp(np2, dim)); Sp = 0.0_dp
    do col = 1, dim
      do p = 1, norb
        g = ann(dets(col), 2*(p-1)+1, d1); if (g==0) cycle   ! annihilate beta
        g2 = cre(d1, 2*(p-1)+0, d2); if (g2==0) cycle         ! create alpha
        row = findidx(pdets, np2, d2)
        if (row > 0) Sp(row,col) = Sp(row,col) + g*g2
      end do
    end do
    S2 = matmul(transpose(Sp), Sp)
    deallocate(Sp, pdets)
  end subroutine build_s2

  !> Ms=+1 determinants with 2 alpha electrons (nocc=1 case): all alpha pairs.
  subroutine build_msplus(norb, pdets, np2)
    integer, intent(in) :: norb
    integer, allocatable, intent(out) :: pdets(:)
    integer, intent(out) :: np2
    integer :: a, b, cnt
    integer, allocatable :: tmp(:)
    allocate(tmp(norb*norb)); cnt = 0
    do a = 0, norb-1
      do b = a+1, norb-1
        cnt = cnt + 1
        tmp(cnt) = ior(ibset(0, 2*a), ibset(0, 2*b))   ! two alpha
      end do
    end do
    np2 = cnt; allocate(pdets(cnt)); pdets = tmp(1:cnt)
  end subroutine build_msplus

  !> (2,2) frontier CAS: determinants with all occupied spatial orbitals in
  !> {nocc, nocc+1} (1-based).
  subroutine cas22_compact(dets, dim, norb, nocc, cas, nc)
    integer, intent(in) :: dets(:), dim, norb, nocc
    integer, allocatable, intent(out) :: cas(:)
    integer, intent(out) :: nc
    integer :: i, p, det, cnt, ok
    integer, allocatable :: tmp(:)
    allocate(tmp(dim)); cnt = 0
    do i = 1, dim
      det = dets(i); ok = 1
      do p = 0, norb-1
        if (btest(det, 2*p) .or. btest(det, 2*p+1)) then
          if (.not. (p == nocc-1 .or. p == nocc)) ok = 0   ! orbital p (0-based)
        end if
      end do
      if (ok == 1) then
        cnt = cnt + 1; tmp(cnt) = i
      end if
    end do
    nc = cnt; allocate(cas(nc)); cas = tmp(1:nc)
  end subroutine cas22_compact

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

  subroutine sym_eig_vec(A, n, w)
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n)
    real(dp), intent(out) :: w(n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsyev('V','U',n,A,n,w,wq,-1,info)
    lw = int(wq(1)); allocate(wk(lw))
    call dsyev('V','U',n,A,n,w,wk,lw,info)
    deallocate(wk)
  end subroutine sym_eig_vec

  subroutine expm_nilpotent(scale, A, n, E)
    integer, intent(in) :: n
    real(dp), intent(in) :: scale, A(n,n)
    real(dp), intent(out) :: E(n,n)
    real(dp) :: term(n,n)
    integer :: k
    E = 0.0_dp; term = 0.0_dp
    do k = 1, n
      E(k,k) = 1.0_dp; term(k,k) = 1.0_dp
    end do
    do k = 1, 16
      term = matmul(term, scale*A)/real(k,dp)
      E = E + term
      if (maxval(abs(term)) < 1.0e-15_dp) exit
    end do
  end subroutine expm_nilpotent

end program tc_h2_pes_native
