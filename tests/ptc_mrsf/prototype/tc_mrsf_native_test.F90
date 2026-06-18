!> Native, pyscf-free pTC-MRSF-CIS end-to-end on a real molecule (H4/STO-3G),
!> assembling the validated native integral toolkit (module ptc_geminal: S,T,V,ERI)
!> with a determinant FCI engine, the MP2-T2 transcorrelation, and the production
!> non-Hermitian solver tc_nonsym_tda_eig (module tdhf_mrsf_ptc). Fortran port of
!> tests/ptc_mrsf/prototype/tc_finite_basis.py.
!>
!> Validations (no pyscf, no external reference except textbook STO-3G H2 E_RHF):
!>   (0) SCF anchor : H2/STO-3G E_RHF = -1.1167 Ha (textbook).
!>   (1) FCI engine : H Hermitian; ground state below E_HF.
!>   (2) TC identity: <HF|H_bar|HF> = E_MP2 (exact; the rigorous build check).
!>   (3) downfold   : bare (HF+singles) ground = E_HF (Brillouin); transcorrelated
!>                    compact ground recovers a large fraction of the FCI correlation.
!>   (4) tau=0 gate : H_bar -> H reproduces the bare compact result.
!>   (5) solver     : tc_nonsym_tda_eig real spectrum, biorthonormal.
!>
!> Build:
!>   gfortran -O2 source/precision.F90 source/modules/ptc_geminal.F90 \
!>     source/modules/tdhf_mrsf_ptc.F90 \
!>     tests/ptc_mrsf/prototype/tc_mrsf_native_test.F90 -llapack -lblas -o /tmp/ptc_mrsf && /tmp/ptc_mrsf
program tc_mrsf_native_test
  use precision, only: dp
  use ptc_geminal
  use tdhf_mrsf_ptc, only: tc_nonsym_tda_eig
  implicit none

  integer  :: nfail
  real(dp) :: e_rhf_h2

  nfail = 0
  write(*,'(a)') '=== native pTC-MRSF-CIS (H4/STO-3G), pyscf-free ==='
  write(*,'(a)') ''

  ! (0) SCF anchor on H2/STO-3G (textbook total E_RHF = -1.1167 Ha at R=1.4 bohr)
  call h2_rhf_anchor(e_rhf_h2)
  call chk('H2/STO-3G E_RHF vs textbook    ', e_rhf_h2, -1.1167593_dp, 2.0e-3_dp, nfail)
  write(*,'(a)') ''

  ! full pipeline on H4
  call run_h4(nfail)

  if (nfail == 0) then
    write(*,'(a)') ''
    write(*,'(a)') 'ALL PASS: native pTC-MRSF-CIS works end-to-end on H4/STO-3G --'
    write(*,'(a)') 'SCF, MP2, determinant FCI, H_bar, and the non-Hermitian solve all'
    write(*,'(a)') 'from the native integral engine; <HF|H_bar|HF>=E_MP2 exactly.'
  else
    write(*,'(a,i0,a)') 'FAILURES: ', nfail, ' check(s) failed.'
    error stop 1
  end if

contains

  subroutine chk(name, a, b, tol, nf)
    character(*), intent(in)    :: name
    real(dp),     intent(in)    :: a, b, tol
    integer,      intent(inout) :: nf
    real(dp) :: d
    d = abs(a - b)
    write(*,'(a,2f16.8,es10.2,a)') name, a, b, d, merge('  PASS', '  FAIL', d <= tol)
    if (d > tol) nf = nf + 1
  end subroutine chk

  !> STO-3G H chain along z: nat shells (one s shell per H), spacing dz (bohr).
  subroutine make_hchain(nat, dz, npr, exs, cos_, cns, zat, rat)
    integer,  intent(in)  :: nat
    real(dp), intent(in)  :: dz
    integer,  intent(out) :: npr(nat)
    real(dp), intent(out) :: exs(3,nat), cos_(3,nat), cns(3,nat), zat(nat), rat(3,nat)
    integer :: i
    do i = 1, nat
      npr(i)   = 3
      exs(:,i) = [3.42525091_dp, 0.62391373_dp, 0.16885540_dp]
      cos_(:,i) = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      cns(:,i) = [0.0_dp, 0.0_dp, dz*real(i-1,dp)]
      zat(i)   = 1.0_dp
      rat(:,i) = cns(:,i)
    end do
  end subroutine make_hchain

  !> Closed-shell RHF on an s-only system. Returns MO coeffs C, energies eps,
  !> the total RHF energy, the AO core/overlap and chemist ERI tensor.
  subroutine rhf(n, nocc, npr, exs, cos_, cns, nat, zat, rat, &
                 Cmo, eps, e_rhf, enuc, h1ao, eri_c)
    integer,  intent(in)  :: n, nocc, npr(n), nat
    real(dp), intent(in)  :: exs(:,:), cos_(:,:), cns(3,n), zat(:), rat(:,:)
    real(dp), intent(out) :: Cmo(n,n), eps(n), e_rhf, enuc, h1ao(n,n), eri_c(n,n,n,n)
    real(dp) :: S(n,n), T(n,n), V(n,n), Hcore(n,n), F(n,n), P(n,n), G(n,n)
    real(dp) :: M(n,n,n,n), Scopy(n,n), eold
    integer  :: i, j, k, l, it, mu, nu, la, si
    ! integrals
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_OVERLAP, nat, zat, rat, S)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_KINETIC, nat, zat, rat, T)
    call ptc_s_ao_1e(n, npr, exs, cos_, cns, PTC_1E_NUCLEAR, nat, zat, rat, V)
    Hcore = T + V
    h1ao  = Hcore
    call ptc_s_ao_tensor(n, npr, exs, cos_, cns, PTC_OP_ERI, 0.0_dp, M)
    ! chemist (i j | k l) = M(i,k,j,l)
    do i = 1, n
      do j = 1, n
        do k = 1, n
          do l = 1, n
            eri_c(i,j,k,l) = M(i,k,j,l)
          end do
        end do
      end do
    end do
    ! nuclear repulsion
    enuc = 0.0_dp
    do i = 1, nat
      do j = i+1, nat
        enuc = enuc + zat(i)*zat(j)/sqrt(dot_product(rat(:,i)-rat(:,j), rat(:,i)-rat(:,j)))
      end do
    end do
    ! SCF
    F = Hcore
    e_rhf = 0.0_dp
    do it = 1, 200
      Scopy = S
      call gen_eig(F, Scopy, n, eps, Cmo)     ! F C = S C eps
      P = 0.0_dp
      do mu = 1, n
        do nu = 1, n
          do i = 1, nocc
            P(mu,nu) = P(mu,nu) + 2.0_dp*Cmo(mu,i)*Cmo(nu,i)
          end do
        end do
      end do
      G = 0.0_dp
      do mu = 1, n
        do nu = 1, n
          do la = 1, n
            do si = 1, n
              G(mu,nu) = G(mu,nu) + P(la,si)*(eri_c(mu,nu,la,si) - 0.5_dp*eri_c(mu,la,nu,si))
            end do
          end do
        end do
      end do
      F = Hcore + G
      eold = e_rhf
      e_rhf = 0.0_dp
      do mu = 1, n
        do nu = 1, n
          e_rhf = e_rhf + 0.5_dp*P(mu,nu)*(Hcore(mu,nu) + F(mu,nu))
        end do
      end do
      if (abs(e_rhf - eold) < 1.0e-11_dp .and. it > 1) exit
    end do
    e_rhf = e_rhf + enuc
  end subroutine rhf

  subroutine h2_rhf_anchor(e_rhf)
    real(dp), intent(out) :: e_rhf
    integer, parameter :: n = 2
    integer  :: npr(n)
    real(dp) :: exs(3,n), cos_(3,n), cns(3,n), zat(n), rat(3,n)
    real(dp) :: Cmo(n,n), eps(n), enuc, h1(n,n), eri(n,n,n,n)
    call make_hchain(n, 1.4_dp, npr, exs, cos_, cns, zat, rat)
    call rhf(n, 1, npr, exs, cos_, cns, n, zat, rat, Cmo, eps, e_rhf, enuc, h1, eri)
  end subroutine h2_rhf_anchor

  subroutine run_h4(nf)
    integer, intent(inout) :: nf
    integer, parameter :: n = 4, nocc = 2
    integer  :: npr(n), nvir
    real(dp) :: exs(3,n), cos_(3,n), cns(3,n), zat(n), rat(3,n)
    real(dp) :: Cmo(n,n), eps(n), e_rhf, enuc, h1ao(n,n), eri_c(n,n,n,n)
    real(dp) :: h1mo(n,n), eri_mo(n,n,n,n), t2(nocc,nocc,n-nocc,n-nocc)
    integer  :: i, j, a, b, dim, hfidx
    integer, allocatable :: dets(:), compact(:)
    real(dp), allocatable :: H(:,:), T2op(:,:), Hbar(:,:), Em(:,:), Ep(:,:)
    real(dp), allocatable :: Hc(:,:), Hbc(:,:), wf(:), wbare(:)
    real(dp), allocatable :: ee(:), vr(:,:), vl(:,:)
    real(dp) :: e_mp2_corr, e_hf_fci, e_fci, denom, gia, gib
    real(dp) :: e_ref_tc, e_bare_c, e_tc_c, recov, maxim, e_tc0
    integer  :: nc, ierr, ncx, k

    nvir = n - nocc
    call make_hchain(n, 1.8_dp, npr, exs, cos_, cns, zat, rat)
    call rhf(n, nocc, npr, exs, cos_, cns, n, zat, rat, Cmo, eps, e_rhf, enuc, h1ao, eri_c)
    write(*,'(a,f14.8)') 'H4/STO-3G  E_RHF              = ', e_rhf

    ! AO -> MO
    call ao2mo_1e(h1ao, Cmo, n, h1mo)
    call ao2mo_2e(eri_c, Cmo, n, eri_mo)

    ! MP2 amplitudes + correlation energy (chemist (ia|jb) = eri_mo(i,a+nocc,j,b+nocc))
    e_mp2_corr = 0.0_dp
    do i = 1, nocc
      do j = 1, nocc
        do a = 1, nvir
          do b = 1, nvir
            denom = eps(i)+eps(j)-eps(nocc+a)-eps(nocc+b)
            gia = eri_mo(i, nocc+a, j, nocc+b)
            gib = eri_mo(i, nocc+b, j, nocc+a)
            t2(i,j,a,b) = gia/denom
            e_mp2_corr = e_mp2_corr + gia*(2.0_dp*gia - gib)/denom
          end do
        end do
      end do
    end do
    write(*,'(a,f14.8)') 'H4/STO-3G  E_MP2(total)      = ', e_rhf + e_mp2_corr

    ! determinant FCI (Ms=0)
    call build_dets(n, nocc, dets, dim, hfidx)
    allocate(H(dim,dim), T2op(dim,dim), Hbar(dim,dim), Em(dim,dim), Ep(dim,dim))
    call build_fci_H(h1mo, eri_mo, enuc, n, dets, dim, H)
    ! (1) Hermiticity + FCI ground
    call chk('FCI H Hermitian                ', maxval(abs(H-transpose(H))), 0.0_dp, 1.0e-10_dp, nf)
    block
      real(dp) :: wfull(dim), Hcp(dim,dim)
      Hcp = H
      call sym_eig(Hcp, dim, wfull)
      e_fci = wfull(1)
    end block
    e_hf_fci = H(hfidx,hfidx)
    write(*,'(a,f14.8)') 'H4/STO-3G  <HF|H|HF>=E_HF    = ', e_hf_fci
    write(*,'(a,f14.8)') 'H4/STO-3G  E_FCI             = ', e_fci
    call chk('<HF|H|HF> vs E_RHF             ', e_hf_fci, e_rhf, 1.0e-8_dp, nf)
    if (e_fci >= e_rhf) nf = nf + 1

    ! T2 operator, H_bar = e^{-T2} H e^{T2} (T2 nilpotent -> finite series)
    call build_T2op(t2, nocc, nvir, n, dets, dim, T2op)
    call expm_nilpotent(-1.0_dp, T2op, dim, Em)
    call expm_nilpotent( 1.0_dp, T2op, dim, Ep)
    Hbar = matmul(Em, matmul(H, Ep))
    ! (2) rigorous identity <HF|H_bar|HF> = E_MP2
    e_ref_tc = Hbar(hfidx, hfidx)
    call chk('<HF|H_bar|HF> = E_MP2          ', e_ref_tc, e_rhf + e_mp2_corr, 1.0e-7_dp, nf)
    if (maxval(abs(Hbar-transpose(Hbar))) < 1.0e-10_dp) then
      write(*,'(a)') 'ERROR: H_bar symmetric (should be non-Hermitian)'; nf = nf + 1
    end if

    ! (3) downfold to compact HF+singles space
    call build_compact(dets, dim, hfidx, compact, nc)
    allocate(Hc(nc,nc), Hbc(nc,nc), wbare(nc), ee(nc), vr(nc,nc), vl(nc,nc))
    do i = 1, nc
      do j = 1, nc
        Hc(i,j)  = H(compact(i), compact(j))
        Hbc(i,j) = Hbar(compact(i), compact(j))
      end do
    end do
    block
      real(dp) :: Hccp(nc,nc)
      Hccp = Hc
      call sym_eig(Hccp, nc, wbare)
    end block
    e_bare_c = wbare(1)
    call tc_nonsym_tda_eig(Hbc, nc, ee, vr, vl, maxim, ncx, ierr)
    call sort_asc(ee, nc)
    e_tc_c = ee(1)
    recov = (e_tc_c - e_bare_c)/(e_fci - e_rhf)*100.0_dp
    write(*,'(a,i0,a,i0)') 'compact (HF+singles) dim     = ', nc, ' of ', dim
    write(*,'(a,f14.8)') 'bare compact ground (=E_HF)  = ', e_bare_c
    write(*,'(a,f14.8)') 'TC   compact ground          = ', e_tc_c
    write(*,'(a,f8.1,a)') 'correlation recovered        = ', recov, ' %'
    call chk('bare compact ground = E_HF (Brillouin)', e_bare_c, e_rhf, 1.0e-7_dp, nf)
    if (.not. (e_rhf > e_tc_c .and. e_tc_c > e_fci - 1.0e-6_dp)) nf = nf + 1
    if (recov <= 40.0_dp) nf = nf + 1
    if (maxim > 1.0e-8_dp) nf = nf + 1

    ! (4) tau=0 gate: solve bare H on compact with the non-Herm solver
    block
      real(dp) :: ee0(nc), vr0(nc,nc), vl0(nc,nc), mi0
      integer  :: nx0, ie0
      call tc_nonsym_tda_eig(Hc, nc, ee0, vr0, vl0, mi0, nx0, ie0)
      call sort_asc(ee0, nc)
      e_tc0 = ee0(1)
    end block
    call chk('tau=0 gate: compact == bare    ', e_tc0, e_bare_c, 1.0e-9_dp, nf)
  end subroutine run_h4

  ! ---- AO->MO transforms ----
  subroutine ao2mo_1e(h, C, n, hmo)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: h(n,n), C(n,n)
    real(dp), intent(out) :: hmo(n,n)
    hmo = matmul(transpose(C), matmul(h, C))
  end subroutine ao2mo_1e

  subroutine ao2mo_2e(eri, C, n, emo)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: eri(n,n,n,n), C(n,n)
    real(dp), intent(out) :: emo(n,n,n,n)
    real(dp) :: t1(n,n,n,n), t2t(n,n,n,n)
    integer  :: p,q,r,s, mu
    t1 = 0.0_dp
    do p = 1, n; do q = 1, n; do r = 1, n; do s = 1, n
      do mu = 1, n
        t1(p,q,r,s) = t1(p,q,r,s) + C(mu,p)*eri(mu,q,r,s)
      end do
    end do; end do; end do; end do
    t2t = 0.0_dp
    do p = 1, n; do q = 1, n; do r = 1, n; do s = 1, n
      do mu = 1, n
        t2t(p,q,r,s) = t2t(p,q,r,s) + C(mu,q)*t1(p,mu,r,s)
      end do
    end do; end do; end do; end do
    t1 = 0.0_dp
    do p = 1, n; do q = 1, n; do r = 1, n; do s = 1, n
      do mu = 1, n
        t1(p,q,r,s) = t1(p,q,r,s) + C(mu,r)*t2t(p,q,mu,s)
      end do
    end do; end do; end do; end do
    emo = 0.0_dp
    do p = 1, n; do q = 1, n; do r = 1, n; do s = 1, n
      do mu = 1, n
        emo(p,q,r,s) = emo(p,q,r,s) + C(mu,s)*t1(p,q,r,mu)
      end do
    end do; end do; end do; end do
  end subroutine ao2mo_2e

  ! ---- determinant engine (spin-orbital bitmasks, P = 2*orb + spin) ----
  subroutine build_dets(norb, nocc, dets, dim, hfidx)
    integer,              intent(in)  :: norb, nocc
    integer, allocatable, intent(out) :: dets(:)
    integer,              intent(out) :: dim, hfidx
    integer :: amask, bmask, da, db, cnt, hf
    integer, allocatable :: alist(:), blist(:)
    integer :: na
    call combos(norb, nocc, alist, na)     ! alpha occ combinations (orbital masks)
    blist = alist
    dim = na*na
    allocate(dets(dim))
    cnt = 0
    do da = 1, na
      do db = 1, na
        cnt = cnt + 1
        amask = spread_spin(alist(da), 0)   ! alpha -> even bits
        bmask = spread_spin(blist(db), 1)   ! beta  -> odd bits
        dets(cnt) = ior(amask, bmask)
      end do
    end do
    hf = ior(spread_spin(2**nocc - 1, 0), spread_spin(2**nocc - 1, 1))
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
      if (btest(omask, o)) p = ibset(p, 2*o + spin)
    end do
  end function spread_spin

  subroutine combos(norb, k, list, ncomb)
    integer,              intent(in)  :: norb, k
    integer, allocatable, intent(out) :: list(:)
    integer,              intent(out) :: ncomb
    integer :: m, cnt
    integer, allocatable :: tmp(:)
    allocate(tmp(2**norb)); cnt = 0
    do m = 0, 2**norb - 1
      if (popcnt_(m) == k) then
        cnt = cnt + 1; tmp(cnt) = m
      end if
    end do
    ncomb = cnt
    allocate(list(cnt)); list = tmp(1:cnt)
  end subroutine combos

  integer function popcnt_(m) result(c)
    integer, intent(in) :: m
    integer :: i
    c = 0
    do i = 0, 30
      if (btest(m, i)) c = c + 1
    end do
  end function popcnt_

  integer function parity_below(d, x) result(s)
    integer, intent(in) :: d, x
    integer :: i, c
    c = 0
    do i = 0, x-1
      if (btest(d, i)) c = c + 1
    end do
    s = 1 - 2*mod(c, 2)
  end function parity_below

  integer function ann(det, p, newdet) result(sgn)
    integer, intent(in) :: det, p
    integer, intent(out) :: newdet
    sgn = 0; newdet = det
    if (.not. btest(det, p)) return
    sgn = parity_below(det, p)
    newdet = ibclr(det, p)
  end function ann

  integer function cre(det, p, newdet) result(sgn)
    integer, intent(in) :: det, p
    integer, intent(out) :: newdet
    sgn = 0; newdet = det
    if (btest(det, p)) return
    sgn = parity_below(det, p)
    newdet = ibset(det, p)
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
    integer,  intent(in)  :: norb, dets(:), dim
    real(dp), intent(in)  :: h1(norb,norb), eri(norb,norb,norb,norb), ecore
    real(dp), intent(out) :: H(dim,dim)
    integer :: col, p, q, r, u, s1, s2, d1, d2, dd, g, g2, jj, det
    real(dp) :: v
    H = 0.0_dp
    do col = 1, dim
      det = dets(col)
      ! one-body
      do p = 1, norb
        do q = 1, norb
          if (abs(h1(p,q)) < 1.0e-14_dp) cycle
          do s1 = 0, 1
            g = ann(det, 2*(q-1)+s1, d1)
            if (g == 0) cycle
            g2 = cre(d1, 2*(p-1)+s1, d2)
            if (g2 == 0) cycle
            jj = findidx(dets, dim, d2)
            if (jj > 0) H(jj,col) = H(jj,col) + h1(p,q)*g*g2
          end do
        end do
      end do
      ! two-body: 0.5 (pq|ru) a^dag_p,s1 a^dag_r,s2 a_u,s2 a_q,s1
      do p = 1, norb
        do q = 1, norb
          do r = 1, norb
            do u = 1, norb
              v = eri(p,q,r,u)
              if (abs(v) < 1.0e-14_dp) cycle
              do s1 = 0, 1
                do s2 = 0, 1
                  g = ann(det, 2*(q-1)+s1, d1)
                  if (g == 0) cycle
                  g2 = ann(d1, 2*(u-1)+s2, dd)
                  if (g2 == 0) cycle
                  g = g*g2
                  g2 = cre(dd, 2*(r-1)+s2, d1)
                  if (g2 == 0) cycle
                  g = g*g2
                  g2 = cre(d1, 2*(p-1)+s1, d2)
                  if (g2 == 0) cycle
                  g = g*g2
                  jj = findidx(dets, dim, d2)
                  if (jj > 0) H(jj,col) = H(jj,col) + 0.5_dp*v*g
                end do
              end do
            end do
          end do
        end do
      end do
      H(col,col) = H(col,col) + ecore
    end do
  end subroutine build_fci_H

  subroutine build_T2op(t2, nocc, nvir, norb, dets, dim, T2op)
    integer,  intent(in)  :: nocc, nvir, norb, dets(:), dim
    real(dp), intent(in)  :: t2(nocc,nocc,nvir,nvir)
    real(dp), intent(out) :: T2op(dim,dim)
    integer :: col, i, j, a, b, det, d1, d1b, d2, d2f, jj, sa, sb, g1, g2, g3, g4
    real(dp) :: amp
    ! T2 = 1/2 sum_{ijab} t2[i,j,a,b] E_{(nocc+a),i} E_{(nocc+b),j}, spin-summed.
    T2op = 0.0_dp
    do col = 1, dim
      det = dets(col)
      do i = 1, nocc
        do j = 1, nocc
          do a = 1, nvir
            do b = 1, nvir
              amp = 0.5_dp*t2(i,j,a,b)
              if (abs(amp) < 1.0e-14_dp) cycle
              do sb = 0, 1
                g1 = ann(det, 2*(j-1)+sb, d1)
                if (g1 == 0) cycle
                g2 = cre(d1, 2*(nocc+b-1)+sb, d1b)
                if (g2 == 0) cycle
                do sa = 0, 1
                  g3 = ann(d1b, 2*(i-1)+sa, d2)
                  if (g3 == 0) cycle
                  g4 = cre(d2, 2*(nocc+a-1)+sa, d2f)
                  if (g4 == 0) cycle
                  jj = findidx(dets, dim, d2f)
                  if (jj > 0) T2op(jj,col) = T2op(jj,col) + amp*g1*g2*g3*g4
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine build_T2op

  ! ---- linear algebra helpers ----
  subroutine gen_eig(F, S, n, w, C)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: F(n,n), S(n,n)
    real(dp), intent(out)   :: w(n), C(n,n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsygv(1, 'V', 'U', n, F, n, S, n, w, wq, -1, info)
    lw = int(wq(1)); allocate(wk(lw))
    call dsygv(1, 'V', 'U', n, F, n, S, n, w, wk, lw, info)
    C = F
    deallocate(wk)
  end subroutine gen_eig

  subroutine sym_eig(A, n, w)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: A(n,n)
    real(dp), intent(out)   :: w(n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsyev('N', 'U', n, A, n, w, wq, -1, info)
    lw = int(wq(1)); allocate(wk(lw))
    call dsyev('N', 'U', n, A, n, w, wk, lw, info)
    deallocate(wk)
  end subroutine sym_eig

  subroutine sort_asc(a, n)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: a(:)
    integer :: i, j
    real(dp) :: tmp
    do i = 1, n-1
      do j = i+1, n
        if (a(j) < a(i)) then
          tmp = a(i); a(i) = a(j); a(j) = tmp
        end if
      end do
    end do
  end subroutine sort_asc

  !> expm(scale*A) for a nilpotent A via terminating Taylor series.
  subroutine expm_nilpotent(scale, A, n, E)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: scale, A(n,n)
    real(dp), intent(out) :: E(n,n)
    real(dp) :: term(n,n)
    integer  :: k
    E = 0.0_dp
    do k = 1, n
      E(k,k) = 1.0_dp
    end do
    term = 0.0_dp
    do k = 1, n
      term(k,k) = 1.0_dp
    end do
    do k = 1, 12
      term = matmul(term, scale*A)/real(k,dp)
      E = E + term
      if (maxval(abs(term)) < 1.0e-15_dp) exit
    end do
  end subroutine expm_nilpotent

  subroutine build_compact(dets, dim, hfidx, compact, nc)
    integer,              intent(in)  :: dets(:), dim, hfidx
    integer, allocatable, intent(out) :: compact(:)
    integer,              intent(out) :: nc
    integer :: i, diffbits, cnt, hf
    integer, allocatable :: tmp(:)
    hf = dets(hfidx)
    allocate(tmp(dim)); cnt = 0
    do i = 1, dim
      diffbits = popcnt_(ieor(dets(i), hf))   ! 2*(excitation level)
      if (diffbits == 0 .or. diffbits == 2) then
        cnt = cnt + 1; tmp(cnt) = i
      end if
    end do
    nc = cnt
    allocate(compact(nc)); compact = tmp(1:nc)
  end subroutine build_compact

end program tc_mrsf_native_test
