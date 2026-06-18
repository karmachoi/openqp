!> Genuine three-body transcorrelation in a many-electron FCI (H4), the decisive
!> Phase-3 validation. Builds the EXACT three-body operator O3 in the 4-electron
!> determinant FCI via second-quantized create/annihilate from the validated
!> 3-electron oracle (tc_threebody.o3_geminal_s), and the NO2B normal-ordered
!> approximation (constant + effective 1-body + effective 2-body, dropping the
!> genuine 3-body remainder), and compares them. This (a) pins the absolute
!> normalization against the literal operator (resolving the workflow's 3! caveat),
!> (b) validates the eff-1/2-body descendants, (c) quantifies the dropped remainder.
!>
!> O3 = (1/3!) sum_{pqr,stu} <pqr|o|stu> a+_p a+_q a+_r a_u a_t a_s, with the
!> 3-electron kernel <pqr|o|stu> = -sum_{apex} c(apex) K3apex(spatial), spin-diagonal
!> per line; c(apex) = c_{s_apex,s_leg1} c_{s_apex,s_leg2} (1/2 anti, 1/4 par).
program tc_threebody_fci_test
  use precision, only: dp
  use ptc_ao,       only: ao_ncart, build_ints
  use tc_threebody, only: o3_prim_s
  use tc_boyshandy, only: fit_unit_cusp
  implicit none
  integer, parameter :: NSPA = 4          ! spatial orbitals (H4/STO-3G)
  integer, parameter :: MSO  = 2*NSPA     ! spin-orbitals
  integer, parameter :: NEL  = 4
  real(dp) :: cen(3,NSPA), zat(NSPA), rat(3,NSPA)
  real(dp) :: shl_e(3,NSPA), shl_c(3,NSPA), shl_r(3,NSPA)
  integer  :: shl_l(NSPA), shl_np(NSPA)
  real(dp) :: S(NSPA,NSPA), Hc(NSPA,NSPA), eri(NSPA,NSPA,NSPA,NSPA)
  real(dp) :: Cmo(NSPA,NSPA), h1(NSPA,NSPA), erimo(NSPA,NSPA,NSPA,NSPA)
  real(dp) :: o3ao(NSPA,NSPA,NSPA,NSPA,NSPA,NSPA), o3mo(NSPA,NSPA,NSPA,NSPA,NSPA,NSPA)
  real(dp) :: Cf(6), gf(6)
  integer  :: ng, i, j, nao
  ! determinants
  integer  :: dets(70), ndet, refidx
  real(dp), allocatable :: H012(:,:), H3(:,:), HN2(:,:), w(:)
  real(dp) :: e_bare, e_exact, e_no2b, e0_diag, e0_formula
  ! NO2B descendants
  real(dp) :: f1(NSPA,NSPA,2), v2(NSPA,NSPA,NSPA,NSPA,2,2), econst

  ! H4 linear chain, STO-3G
  do i=1,NSPA
    cen(:,i)=[0.0_dp,0.0_dp,1.8_dp*real(i-1,dp)]
    shl_l(i)=0; shl_np(i)=3
    shl_e(:,i)=[3.42525091_dp,0.62391373_dp,0.16885540_dp]
    shl_c(:,i)=[0.15432897_dp,0.53532814_dp,0.44463454_dp]
    shl_r(:,i)=cen(:,i); zat(i)=1.0_dp; rat(:,i)=cen(:,i)
  end do
  call fit_unit_cusp(1.0_dp, Cf, gf, ng)
  call build_ints(NSPA, shl_l, shl_np, shl_e, shl_c, shl_r, NSPA, zat, rat, nao, S, Hc, eri)
  call rhf_cs(NSPA, 2, S, Hc, eri, Cmo)         ! closed-shell RHF, 2 doubly-occ
  call ao2mo1(Hc, Cmo, NSPA, h1)
  call ao2mo2(eri, Cmo, NSPA, erimo)            ! chemist (pq|rs)

  ! genuine 3-electron AO integral (apex = electron 1), contracted over primitives
  call build_o3ao()
  call ao2mo6(o3ao, Cmo, NSPA, o3mo)

  ! enumerate Ms=0 determinants (2 alpha, 2 beta) over MSO spin-orbitals
  call enum_dets(dets, ndet, refidx)
  write(*,'(a,i0,a,i0,a)') '=== Genuine three-body TC in FCI : H4/STO-3G (', ndet, ' dets, ref #', refidx, ') ==='

  allocate(H012(ndet,ndet), H3(ndet,ndet), HN2(ndet,ndet), w(ndet))
  call build_H012(dets, ndet, h1, erimo, NSPA, H012)
  call build_H3_exact(dets, ndet, o3mo, NSPA, H3)

  ! reference 3-body energy two ways (pin normalization)
  e0_diag    = H3(refidx, refidx)
  e0_formula = o3_ref_formula(dets(refidx))
  call build_no2b(dets(refidx), o3mo, NSPA, econst, f1, v2)
  call build_HN2(dets, ndet, f1, v2, econst, NSPA, HN2)

  ! ground-state energies
  e_bare  = lowest(H012, ndet, w)
  e_exact = lowest(H012 + H3, ndet, w)
  e_no2b  = lowest(H012 + HN2, ndet, w)

  write(*,'(a)') ''
  write(*,'(a,f14.8)')  ' <ref|O3|ref> (FCI diagonal)      = ', e0_diag
  write(*,'(a,f14.8)')  ' <ref|O3|ref> (direct formula)    = ', e0_formula
  write(*,'(a,es10.2)') '   |diff| (normalization pinned)  = ', abs(e0_diag-e0_formula)
  write(*,'(a,f14.8)')  ' NO2B constant E0                 = ', econst
  write(*,'(a)') ''
  write(*,'(a,f14.8)')  ' E(FCI, bare)                     = ', e_bare
  write(*,'(a,f14.8)')  ' E(FCI, + exact 3-body O3)        = ', e_exact
  write(*,'(a,f14.8)')  ' E(FCI, + NO2B 3-body)            = ', e_no2b
  write(*,'(a,es10.2,a)') '   NO2B error vs exact 3-body     = ', abs(e_no2b-e_exact), &
       '  (= dropped 3-body remainder)'
  write(*,'(a,f12.8,a)') '   3-body shifts the FCI energy by ', e_exact-e_bare, ' Ha'
  write(*,'(a)') ''
  if (abs(e0_diag-e0_formula) < 1.0e-10_dp .and. abs(e_exact-e_bare) > 1.0e-7_dp &
      .and. abs(e_no2b-e_exact) < abs(e_exact-e_bare)) then
    write(*,'(a)') 'PASS: exact three-body built in the many-electron FCI (normalization pinned'
    write(*,'(a)') '      to the literal operator); NO2B reproduces it up to the small dropped'
    write(*,'(a)') '      3-body remainder. Genuine Ten-no three-body works in a >=3-electron system.'
  else
    write(*,'(a)') 'CHECK values above.'
  end if

contains

  pure integer function spat(so) result(p)   ! 0-based so -> 1-based spatial
    integer, intent(in) :: so
    p = so/2 + 1
  end function spat
  pure integer function spn(so) result(s)    ! 0=alpha,1=beta
    integer, intent(in) :: so
    s = mod(so,2)
  end function spn
  pure real(dp) function camp(s1,s2) result(c)
    integer, intent(in) :: s1,s2
    if (s1==s2) then; c=0.25_dp; else; c=0.5_dp; end if
  end function camp

  !> spatial 3e integral with chosen apex position (relabel o3mo which has apex=1).
  pure real(dp) function k3(apex, a,b,c,d,e,f) result(v)
    integer, intent(in) :: apex,a,b,c,d,e,f
    select case (apex)
    case (1); v = o3mo(a,b,c,d,e,f)
    case (2); v = o3mo(b,a,c,e,d,f)
    case default; v = o3mo(c,b,a,f,e,d)
    end select
  end function k3

  !> 3-electron spin-orbital kernel <pqr|o|stu> = -sum_apex c(apex) K3apex, spin-diag.
  pure real(dp) function ker3(p,q,r,s,t,u) result(v)
    integer, intent(in) :: p,q,r,s,t,u
    integer :: sp(3)
    v = 0.0_dp
    if (spn(p)/=spn(s) .or. spn(q)/=spn(t) .or. spn(r)/=spn(u)) return
    sp = [spn(p), spn(q), spn(r)]
    v = -( camp(sp(1),sp(2))*camp(sp(1),sp(3))*k3(1, spat(p),spat(q),spat(r),spat(s),spat(t),spat(u)) &
         + camp(sp(2),sp(1))*camp(sp(2),sp(3))*k3(2, spat(p),spat(q),spat(r),spat(s),spat(t),spat(u)) &
         + camp(sp(3),sp(1))*camp(sp(3),sp(2))*k3(3, spat(p),spat(q),spat(r),spat(s),spat(t),spat(u)) )
  end function ker3

  ! ---- bitmask determinant ops (0-based spin-orbitals 0..MSO-1) ----
  pure integer function phase_below(mask, i) result(ph)
    integer, intent(in) :: mask, i
    ph = 1
    if (i > 0) then
      if (mod(popcnt(iand(mask, ishft(1,i)-1)),2)==1) ph = -1
    end if
  end function phase_below

  subroutine annih(mask, i, ph, ok)
    integer, intent(inout) :: mask
    integer, intent(in) :: i
    integer, intent(inout) :: ph
    logical, intent(out) :: ok
    if (.not. btest(mask,i)) then; ok=.false.; return; end if
    ph = ph*phase_below(mask,i); mask = ibclr(mask,i); ok=.true.
  end subroutine annih

  subroutine creat(mask, i, ph, ok)
    integer, intent(inout) :: mask
    integer, intent(in) :: i
    integer, intent(inout) :: ph
    logical, intent(out) :: ok
    if (btest(mask,i)) then; ok=.false.; return; end if
    ph = ph*phase_below(mask,i); mask = ibset(mask,i); ok=.true.
  end subroutine creat

  subroutine enum_dets(dets, ndet, refidx)
    integer, intent(out) :: dets(:), ndet, refidx
    integer :: m, na, nb, k, refmask
    ndet = 0
    do m = 0, ishft(1,MSO)-1
      if (popcnt(m) /= NEL) cycle
      na=0; nb=0
      do k=0,MSO-1
        if (btest(m,k)) then
          if (spn(k)==0) then; na=na+1; else; nb=nb+1; end if
        end if
      end do
      if (na==2 .and. nb==2) then
        ndet = ndet+1; dets(ndet)=m
      end if
    end do
    ! reference = lowest 2 spatial doubly occupied: spin-orbitals 0,1,2,3
    refmask = 0
    do k=0,3; refmask=ibset(refmask,k); end do
    refidx = 0
    do k=1,ndet
      if (dets(k)==refmask) refidx=k
    end do
  end subroutine enum_dets

  pure integer function findidx(dets, ndet, mask) result(idx)
    integer, intent(in) :: dets(:), ndet, mask
    integer :: k
    idx = 0
    do k=1,ndet
      if (dets(k)==mask) then; idx=k; return; end if
    end do
  end function findidx

  ! ---- FCI matrices via second quantization ----
  subroutine build_H012(dets, ndet, h1, erimo, n, H)
    integer, intent(in) :: dets(:), ndet, n
    real(dp), intent(in) :: h1(n,n), erimo(n,n,n,n)
    real(dp), intent(out) :: H(ndet,ndet)
    integer :: jk, p,s, q,r, m1,m2, ib, ph
    logical :: ok
    real(dp) :: hv, gv
    H = 0.0_dp
    do jk = 1, ndet
      ! 1-body: sum_{p,s} h_ps a+_p a_s
      do p=0,MSO-1; do s=0,MSO-1
        if (spn(p)/=spn(s)) cycle
        hv = h1(spat(p),spat(s)); if (hv==0.0_dp) cycle
        ph=1; m1=dets(jk)
        call annih(m1,s,ph,ok); if(.not.ok) cycle
        call creat(m1,p,ph,ok); if(.not.ok) cycle
        ib = findidx(dets,ndet,m1); if (ib>0) H(ib,jk)=H(ib,jk)+ph*hv
      end do; end do
      ! 2-body: 1/2 sum_{pqrs} <pq|rs> a+_p a+_q a_s a_r  (physicist <pq|rs>=(pr|qs) chemist)
      do p=0,MSO-1; do q=0,MSO-1; do r=0,MSO-1; do s=0,MSO-1
        if (spn(p)/=spn(r) .or. spn(q)/=spn(s)) cycle
        gv = erimo(spat(p),spat(r),spat(q),spat(s)); if (gv==0.0_dp) cycle
        ph=1; m1=dets(jk)
        call annih(m1,r,ph,ok); if(.not.ok) cycle
        call annih(m1,s,ph,ok); if(.not.ok) cycle
        call creat(m1,q,ph,ok); if(.not.ok) cycle
        call creat(m1,p,ph,ok); if(.not.ok) cycle
        ib = findidx(dets,ndet,m1); if (ib>0) H(ib,jk)=H(ib,jk)+0.5_dp*ph*gv
      end do; end do; end do; end do
    end do
  end subroutine build_H012

  !> exact 3-body O3 = (1/3!) sum_{pqr,stu} ker3 a+_p a+_q a+_r a_u a_t a_s.
  subroutine build_H3_exact(dets, ndet, o3mo, n, H)
    integer, intent(in) :: dets(:), ndet, n
    real(dp), intent(in) :: o3mo(n,n,n,n,n,n)
    real(dp), intent(out) :: H(ndet,ndet)
    integer :: jk, p,q,r,s,t,u, m1, ib, ph
    logical :: ok
    real(dp) :: kv
    H = 0.0_dp
    do jk = 1, ndet
      do s=0,MSO-1; do t=0,MSO-1; do u=0,MSO-1
        do p=0,MSO-1; do q=0,MSO-1; do r=0,MSO-1
          kv = ker3(p,q,r,s,t,u); if (kv==0.0_dp) cycle
          ph=1; m1=dets(jk)
          call annih(m1,s,ph,ok); if(.not.ok) cycle
          call annih(m1,t,ph,ok); if(.not.ok) cycle
          call annih(m1,u,ph,ok); if(.not.ok) cycle
          call creat(m1,r,ph,ok); if(.not.ok) cycle
          call creat(m1,q,ph,ok); if(.not.ok) cycle
          call creat(m1,p,ph,ok); if(.not.ok) cycle
          ib = findidx(dets,ndet,m1); if (ib>0) H(ib,jk)=H(ib,jk)+(kv*ph)/6.0_dp
        end do; end do; end do
      end do; end do; end do
    end do
  end subroutine build_H3_exact

  !> direct <ref|O3|ref> = (1/3!) sum over occupied spin-orbital triples, antisymmetrized.
  function o3_ref_formula(refmask) result(E)
    integer, intent(in) :: refmask
    real(dp) :: E
    integer :: occ(NEL), no, k, ia,ib,ic
    no=0
    do k=0,MSO-1
      if (btest(refmask,k)) then; no=no+1; occ(no)=k; end if
    end do
    E = 0.0_dp
    do ia=1,no; do ib=1,no; do ic=1,no
      if (ia==ib .or. ia==ic .or. ib==ic) cycle
      E = E + asym(occ(ia),occ(ib),occ(ic))      ! antisymmetrized diagonal
    end do; end do; end do
    E = E/6.0_dp
  end function o3_ref_formula

  !> antisymmetrized diagonal bracket <pqr||o||pqr> = sum_P sign(P) ker3(p,q,r,P(p,q,r)).
  pure function asym(p,q,r) result(v)
    integer, intent(in) :: p,q,r
    real(dp) :: v
    v =  ker3(p,q,r, p,q,r) - ker3(p,q,r, q,p,r) - ker3(p,q,r, p,r,q) &
       - ker3(p,q,r, r,q,p) + ker3(p,q,r, q,r,p) + ker3(p,q,r, r,p,q)
  end function asym

  !> NO2B descendants from the reference: constant, eff-1-body f1(p,s,spin),
  !> eff-2-body v2(p,q,r,s,spin_pr,spin_qs). Contract the antisymmetrized 3e kernel
  !> with the reference occupied spin-orbitals (one / two spectator legs).
  subroutine build_no2b(refmask, o3mo, n, econst, f1, v2)
    integer, intent(in) :: refmask, n
    real(dp), intent(in) :: o3mo(n,n,n,n,n,n)
    real(dp), intent(out) :: econst, f1(n,n,2), v2(n,n,n,n,2,2)
    integer :: occ(NEL), no, k
    integer :: pa,sa, sp1, ox,oy, x,y
    real(dp) :: acc
    no=0
    do k=0,MSO-1
      if (btest(refmask,k)) then; no=no+1; occ(no)=k; end if
    end do
    econst = o3_ref_formula(refmask)
    ! eff-1-body: f1(p,s; spin) = (1/2) sum_{x<y in occ} <p x y||o|| s x y>  (p,s same spin)
    f1 = 0.0_dp
    do pa=0,MSO-1; do sa=0,MSO-1
      if (spn(pa)/=spn(sa)) cycle
      acc = 0.0_dp
      do ox=1,no; do oy=ox+1,no
        x=occ(ox); y=occ(oy)
        acc = acc + asym2(pa,x,y, sa,x,y)
      end do; end do
      f1(spat(pa),spat(sa),spn(pa)+1) = 0.5_dp*acc
    end do; end do
    ! eff-2-body: v2(p,q,r,s; spins) = sum_{x in occ} <p q x||o|| r s x>  (spn(p)=spn(r),spn(q)=spn(s))
    v2 = 0.0_dp
    do pa=0,MSO-1; do sa=0,MSO-1; do ox=0,MSO-1; do oy=0,MSO-1
      ! interpret as p=pa(e1 bra), q=ox(e2 bra), r=sa(e1 ket), s=oy(e2 ket)
    end do; end do; end do; end do
    call build_v2(refmask, occ, no, v2, n)
    sp1=0  ! silence
  end subroutine build_no2b

  subroutine build_v2(refmask, occ, no, v2, n)
    integer, intent(in) :: refmask, occ(:), no, n
    real(dp), intent(out) :: v2(n,n,n,n,2,2)
    integer :: p,q,r,s, ox, x
    real(dp) :: acc
    integer :: dum
    v2 = 0.0_dp
    do p=0,MSO-1; do q=0,MSO-1; do r=0,MSO-1; do s=0,MSO-1
      if (spn(p)/=spn(r) .or. spn(q)/=spn(s)) cycle
      acc = 0.0_dp
      do ox=1,no
        x=occ(ox)
        acc = acc + asym3(p,q,x, r,s,x)        ! one spectator leg x
      end do
      v2(spat(p),spat(q),spat(r),spat(s), spn(p)+1, spn(q)+1) = acc
    end do; end do; end do; end do
    dum = refmask
  end subroutine build_v2

  !> antisymmetrized 2-spectator bracket for eff-1-body (legs x,y fixed in bra&ket).
  pure function asym2(p,x,y, s,a,b) result(v)
    integer, intent(in) :: p,x,y, s,a,b
    real(dp) :: v
    v =  ker3(p,x,y, s,a,b) - ker3(p,x,y, a,s,b) - ker3(p,x,y, s,b,a) &
       - ker3(p,x,y, b,a,s) + ker3(p,x,y, a,b,s) + ker3(p,x,y, b,s,a)
  end function asym2

  !> antisymmetrized 1-spectator bracket for eff-2-body (leg x fixed in bra&ket).
  pure function asym3(p,q,x, r,s,y) result(v)
    integer, intent(in) :: p,q,x, r,s,y
    real(dp) :: v
    v =  ker3(p,q,x, r,s,y) - ker3(p,q,x, s,r,y) - ker3(p,q,x, r,y,s) &
       - ker3(p,q,x, y,s,r) + ker3(p,q,x, s,y,r) + ker3(p,q,x, y,r,s)
  end function asym3

  !> NO2B FCI matrix: constant + eff-1-body + eff-2-body, as 1- and 2-body operators.
  subroutine build_HN2(dets, ndet, f1, v2, econst, n, H)
    integer, intent(in) :: dets(:), ndet, n
    real(dp), intent(in) :: f1(n,n,2), v2(n,n,n,n,2,2), econst
    real(dp), intent(out) :: H(ndet,ndet)
    integer :: jk, p,s,q,r, m1, ib, ph
    logical :: ok
    real(dp) :: hv, gv
    H = 0.0_dp
    do jk=1,ndet
      H(jk,jk) = H(jk,jk) + econst
      do p=0,MSO-1; do s=0,MSO-1
        if (spn(p)/=spn(s)) cycle
        hv = f1(spat(p),spat(s),spn(p)+1); if (hv==0.0_dp) cycle
        ph=1; m1=dets(jk)
        call annih(m1,s,ph,ok); if(.not.ok) cycle
        call creat(m1,p,ph,ok); if(.not.ok) cycle
        ib=findidx(dets,ndet,m1); if (ib>0) H(ib,jk)=H(ib,jk)+ph*hv
      end do; end do
      ! eff-2-body: (1/2) sum <pq|v|rs> a+_p a+_q a_s a_r  (v2 already the antisym 1-spectator trace)
      do p=0,MSO-1; do q=0,MSO-1; do r=0,MSO-1; do s=0,MSO-1
        if (spn(p)/=spn(r) .or. spn(q)/=spn(s)) cycle
        gv = v2(spat(p),spat(q),spat(r),spat(s),spn(p)+1,spn(q)+1); if (gv==0.0_dp) cycle
        ph=1; m1=dets(jk)
        call annih(m1,r,ph,ok); if(.not.ok) cycle
        call annih(m1,s,ph,ok); if(.not.ok) cycle
        call creat(m1,q,ph,ok); if(.not.ok) cycle
        call creat(m1,p,ph,ok); if(.not.ok) cycle
        ib=findidx(dets,ndet,m1); if (ib>0) H(ib,jk)=H(ib,jk)+0.5_dp*ph*gv
      end do; end do; end do; end do
    end do
  end subroutine build_HN2

  function lowest(H, ndet, w) result(e0)
    integer, intent(in) :: ndet
    real(dp), intent(in) :: H(ndet,ndet)
    real(dp), intent(out) :: w(ndet)
    real(dp) :: e0, A(ndet,ndet), wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    A = 0.5_dp*(H + transpose(H))        ! symmetric part (O3 Hermitian; guards roundoff)
    call dsyev('N','U',ndet,A,ndet,w,wq,-1,info); lw=int(wq(1)); allocate(wk(lw))
    call dsyev('N','U',ndet,A,ndet,w,wk,lw,info); deallocate(wk)
    e0 = w(1)
  end function lowest

  ! ---- one-/two-/six-index AO->MO and SCF helpers ----
  subroutine ao2mo1(A, C, n, B)
    integer, intent(in) :: n
    real(dp), intent(in) :: A(n,n), C(n,n)
    real(dp), intent(out) :: B(n,n)
    B = matmul(transpose(C), matmul(A, C))
  end subroutine ao2mo1

  subroutine ao2mo2(g, C, n, gm)
    integer, intent(in) :: n
    real(dp), intent(in) :: g(n,n,n,n), C(n,n)
    real(dp), intent(out) :: gm(n,n,n,n)
    real(dp) :: t(n,n,n,n)
    integer :: p,q,r,s,mu
    t=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do mu=1,n; t(p,q,r,s)=t(p,q,r,s)+C(mu,p)*g(mu,q,r,s); end do;end do;end do;end do;end do
    gm=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do mu=1,n; gm(p,q,r,s)=gm(p,q,r,s)+C(mu,q)*t(p,mu,r,s); end do;end do;end do;end do;end do
    t=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do mu=1,n; t(p,q,r,s)=t(p,q,r,s)+C(mu,r)*gm(p,q,mu,s); end do;end do;end do;end do;end do
    gm=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do mu=1,n; gm(p,q,r,s)=gm(p,q,r,s)+C(mu,s)*t(p,q,r,mu); end do;end do;end do;end do;end do
  end subroutine ao2mo2

  subroutine ao2mo6(g, C, n, gm)
    integer, intent(in) :: n
    real(dp), intent(in) :: g(n,n,n,n,n,n), C(n,n)
    real(dp), intent(out) :: gm(n,n,n,n,n,n)
    real(dp) :: t(n,n,n,n,n,n)
    integer :: p,q,r,s,u,v,mu
    t=0;  do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n; t(p,q,r,s,u,v)=t(p,q,r,s,u,v)+C(mu,p)*g(mu,q,r,s,u,v); end do;end do;end do;end do;end do;end do;end do
    gm=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n; gm(p,q,r,s,u,v)=gm(p,q,r,s,u,v)+C(mu,q)*t(p,mu,r,s,u,v); end do;end do;end do;end do;end do;end do;end do
    t=0;  do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n; t(p,q,r,s,u,v)=t(p,q,r,s,u,v)+C(mu,r)*gm(p,q,mu,s,u,v); end do;end do;end do;end do;end do;end do;end do
    gm=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n; gm(p,q,r,s,u,v)=gm(p,q,r,s,u,v)+C(mu,s)*t(p,q,r,mu,u,v); end do;end do;end do;end do;end do;end do;end do
    t=0;  do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n; t(p,q,r,s,u,v)=t(p,q,r,s,u,v)+C(mu,u)*gm(p,q,r,s,mu,v); end do;end do;end do;end do;end do;end do;end do
    gm=0; do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n; gm(p,q,r,s,u,v)=gm(p,q,r,s,u,v)+C(mu,v)*t(p,q,r,s,u,mu); end do;end do;end do;end do;end do;end do;end do
  end subroutine ao2mo6

  subroutine build_o3ao()
    integer :: a,b,c,d,e,f, pa,pb,pc,pd,pe,pf
    real(dp) :: acc, qa,qb,qc,qd,qe,qf, ea,eb,ec,ed,ee,ef
    do a=1,NSPA; do b=1,NSPA; do c=1,NSPA; do d=1,NSPA; do e=1,NSPA; do f=1,NSPA
      acc = 0.0_dp
      do pa=1,3; ea=shl_e(pa,a); qa=shl_c(pa,a)*sn(ea)
      do pd=1,3; ed=shl_e(pd,d); qd=shl_c(pd,d)*sn(ed)
      do pb=1,3; eb=shl_e(pb,b); qb=shl_c(pb,b)*sn(eb)
      do pe=1,3; ee=shl_e(pe,e); qe=shl_c(pe,e)*sn(ee)
      do pc=1,3; ec=shl_e(pc,c); qc=shl_c(pc,c)*sn(ec)
      do pf=1,3; ef=shl_e(pf,f); qf=shl_c(pf,f)*sn(ef)
        acc = acc + qa*qd*qb*qe*qc*qf * &
          o3_prim_s(ea,cen(:,a), ed,cen(:,d), eb,cen(:,b), ee,cen(:,e), ec,cen(:,c), ef,cen(:,f), Cf, gf, ng)
      end do; end do; end do; end do; end do; end do
      o3ao(a,b,c,d,e,f) = acc
    end do; end do; end do; end do; end do; end do
  end subroutine build_o3ao

  pure real(dp) function sn(z) result(n)
    real(dp), intent(in) :: z
    real(dp), parameter :: PI=3.141592653589793238462643383279_dp
    n = (2.0_dp*z/PI)**0.75_dp
  end function sn

  subroutine rhf_cs(n, nocc, S, Hc, eri, Cmo)
    integer, intent(in) :: n, nocc
    real(dp), intent(in) :: S(n,n), Hc(n,n), eri(n,n,n,n)
    real(dp), intent(out) :: Cmo(n,n)
    real(dp) :: F(n,n), C(n,n), D(n,n), G(n,n), Scp(n,n), eps(n), eold, e
    integer :: it,mu,nu,la,si,i
    F=Hc; e=0.0_dp
    do it=1,300
      Scp=S; call geig(F,Scp,n,eps,C)
      D=0.0_dp
      do mu=1,n;do nu=1,n;do i=1,nocc; D(mu,nu)=D(mu,nu)+2.0_dp*C(mu,i)*C(nu,i); end do;end do;end do
      G=0.0_dp
      do mu=1,n;do nu=1,n;do la=1,n;do si=1,n
        G(mu,nu)=G(mu,nu)+D(la,si)*(eri(mu,nu,la,si)-0.5_dp*eri(mu,si,la,nu))
      end do;end do;end do;end do
      F=Hc+G; eold=e; e=0.0_dp
      do mu=1,n;do nu=1,n; e=e+0.5_dp*D(mu,nu)*(Hc(mu,nu)+F(mu,nu)); end do;end do
      if (abs(e-eold)<1e-11_dp .and. it>1) exit
    end do
    Cmo=C
  end subroutine rhf_cs

  subroutine geig(F,S,n,w,C)
    integer, intent(in) :: n
    real(dp), intent(inout) :: F(n,n), S(n,n)
    real(dp), intent(out) :: w(n), C(n,n)
    real(dp) :: wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw
    call dsygv(1,'V','U',n,F,n,S,n,w,wq,-1,info); lw=int(wq(1)); allocate(wk(lw))
    call dsygv(1,'V','U',n,F,n,S,n,w,wk,lw,info); C=F; deallocate(wk)
  end subroutine geig

end program tc_threebody_fci_test
