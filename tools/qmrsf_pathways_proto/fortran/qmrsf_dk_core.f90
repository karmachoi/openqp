! QMRSF-DK dressed-kernel pole search -- LIVE-shaped core (standalone, build-testable).
!
! Faithful Fortran port of the validated NumPy LIVE prototype
! (qmrsf_dk_live_proto.py). Reads the live-shaped model from dk_ref.dat:
!     A0       (Ns x Ns)  adiabatic single-spin-flip response block (symmetric),
!     omega_d  (Nd,)      bare 0OS closed-shell double-spin-flip energies,
!     V        (Ns x Nd)  single<->0OS couplings,
! then INDEPENDENTLY:
!   (1) assembles the frequency-dependent dressed kernel
!         g_xc(omega)_{c c'} = sum_d V(c,d) V(c',d) / (omega - omega_d)
!   (2) forms the pole-cancelled secular function (= augmented characteristic poly)
!         fsec(omega) = det[ omega I - A0 - g_xc(omega) ] * prod_d (omega - omega_d)
!       (a smooth, pole-free function whose Ns+Nd roots are the exact spectrum),
!   (3) root-finds ALL Ns+Nd roots by a fine sign scan + bisection over the
!       Gershgorin window of the augmented matrix,
!   (4) cross-checks against (a) the NumPy reference spectra in dk_ref.dat AND
!       (b) its OWN exact diagonalization (dsyev) of the explicit augmented
!       matrix [[A0,V],[V^T,diag(omega_d)]] -- so it is self-validating,
!   (5) demonstrates that the ADIABATIC kernel (g_xc=0) yields only Ns roots
!       (eigvals of A0) and MISSES the Nd 0OS doubles.
!
! This validates the dressed-kernel ALGEBRA in isolation; the production module
! source/modules/tdhf_qmrsf_dk.F90 will wire the same assembly (steps 1-3, the
! reusable kernel) to live OpenQP tagarray quantities (A0 from the backbone
! singles block, omega_d / V from the active-space integrals).
!
! Build & run (standalone, does NOT touch liboqp):
!   gfortran-15 -O2 qmrsf_dk_core.f90 -o /tmp/dk_core -framework Accelerate
!   cd <this dir> && /tmp/dk_core
!
! Fortran case-insensitivity note: local vars never collide (case-only) with
! dummy args; the dressed-matrix and gxc use distinct names (amat, gx).
module qmrsf_dk_core
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains

  !> Dressed kernel g_xc(omega)_{c c'} = sum_d V(c,d) V(c',d) / (omega - omega_d).
  subroutine build_gxc(ns, nd, omega, V, omega_d, gx)
    integer,  intent(in)  :: ns, nd
    real(dp), intent(in)  :: omega, V(ns,nd), omega_d(nd)
    real(dp), intent(out) :: gx(ns,ns)
    integer  :: c, cp, d
    real(dp) :: denom
    gx = 0.0_dp
    do d = 1, nd
      denom = omega - omega_d(d)
      do cp = 1, ns
        do c = 1, ns
          gx(c,cp) = gx(c,cp) + V(c,d)*V(cp,d)/denom
        end do
      end do
    end do
  end subroutine build_gxc

  !> det of a general ns x ns matrix via LU (LAPACK dgetrf). Matrix destroyed.
  real(dp) function det_lu(ns, amat) result(det)
    integer,  intent(in)    :: ns
    real(dp), intent(inout) :: amat(ns,ns)
    integer :: ipiv(ns), info, i
    call dgetrf(ns, ns, amat, ns, ipiv, info)
    if (info < 0) stop 'dgetrf illegal arg'
    det = 1.0_dp
    do i = 1, ns
      det = det * amat(i,i)
      if (ipiv(i) /= i) det = -det     ! row swap flips the sign
    end do
  end function det_lu

  !> Pole-cancelled secular function:
  !>   fsec(omega) = det[ omega I - A0 - g_xc(omega) ] * prod_d (omega - omega_d).
  real(dp) function fsec(ns, nd, omega, A0, V, omega_d) result(fval)
    integer,  intent(in) :: ns, nd
    real(dp), intent(in) :: omega, A0(ns,ns), V(ns,nd), omega_d(nd)
    real(dp) :: gx(ns,ns), amat(ns,ns), poleprod
    integer  :: i, d
    call build_gxc(ns, nd, omega, V, omega_d, gx)
    amat = -A0 - gx
    do i = 1, ns
      amat(i,i) = amat(i,i) + omega
    end do
    poleprod = 1.0_dp
    do d = 1, nd
      poleprod = poleprod * (omega - omega_d(d))
    end do
    fval = det_lu(ns, amat) * poleprod
  end function fsec

  !> Bisect a root of fsec in [a,b] given a sign change (fa = fsec(a)).
  real(dp) function root_bisect(ns, nd, a0in, b0in, faIn, A0, V, omega_d) result(rt)
    integer,  intent(in) :: ns, nd
    real(dp), intent(in) :: a0in, b0in, faIn, A0(ns,ns), V(ns,nd), omega_d(nd)
    real(dp) :: a, b, fa, m, fm
    integer  :: it
    a = a0in; b = b0in; fa = faIn
    do it = 1, 300
      m = 0.5_dp*(a+b)
      fm = fsec(ns, nd, m, A0, V, omega_d)
      if (b - a < 1.0d-14 .or. fm == 0.0_dp) exit
      if (fm*fa > 0.0_dp) then
        a = m; fa = fm
      else
        b = m
      end if
    end do
    rt = 0.5_dp*(a+b)
  end function root_bisect

  subroutine dsortvec(a, m)
    real(dp), intent(inout) :: a(*)
    integer,  intent(in)    :: m
    integer  :: i, j
    real(dp) :: t
    do i = 1, m-1
      do j = 1, m-i
        if (a(j) > a(j+1)) then; t=a(j); a(j)=a(j+1); a(j+1)=t; end if
      end do
    end do
  end subroutine dsortvec

end module qmrsf_dk_core

program test_dk_core
  use qmrsf_dk_core
  implicit none
  integer  :: ns, nd, ntot, i, j, d, info, lwork, nr, nscan, ngrid
  real(dp), allocatable :: A0(:,:), V(:,:), omega_d(:)
  real(dp), allocatable :: refExact(:), refDressed(:), dressed(:), exactF(:)
  real(dp), allocatable :: Haug(:,:), evx(:), work(:), adiab(:), Acopy(:,:), awork(:)
  real(dp), allocatable :: xs(:), fv(:), dblw(:)
  real(dp) :: lo, hi, fa, dmaxRef, dmaxExF, admaxRef
  real(dp) :: worst_adiab_gap, worst_dressed_gap, EE, dblgap, dedup
  logical  :: pass1, pass2

  ! ---- read the live-shaped model + reference spectra ----
  open(10, file="dk_ref.dat", status="old", action="read")
  read(10,*) ns, nd
  ntot = ns + nd
  allocate(A0(ns,ns), V(ns,nd), omega_d(nd))
  allocate(refExact(ntot), refDressed(ntot), dressed(ntot), exactF(ntot))
  do i = 1, ns; read(10,*) (A0(i,j), j=1,ns); end do
  read(10,*) (omega_d(d), d=1,nd)
  do i = 1, ns; read(10,*) (V(i,d), d=1,nd); end do
  read(10,*) (refExact(i), i=1,ntot)
  read(10,*) (refDressed(i), i=1,ntot)
  close(10)

  ! ===================================================================
  ! (4b) self-check: OWN exact diagonalization of the augmented matrix
  !      Haug = [[A0, V],[V^T, diag(omega_d)]]
  ! ===================================================================
  allocate(Haug(ntot,ntot), evx(ntot))
  Haug = 0.0_dp
  Haug(1:ns,1:ns) = A0
  do d = 1, nd
    Haug(1:ns, ns+d) = V(:,d)
    Haug(ns+d, 1:ns) = V(:,d)
    Haug(ns+d, ns+d) = omega_d(d)
  end do
  lwork = max(64*ntot, 256)
  allocate(work(lwork), dblw(ntot))
  call dsyev('V','U', ntot, Haug, ntot, evx, work, lwork, info)   ! 'V' -> Haug holds eigenvecs
  if (info /= 0) stop 'dsyev Haug failed'
  exactF = evx                      ! Fortran's own exact spectrum (ascending)
  ! doubles-sector weight of each eigenvector (rows ns+1..ntot), columns = states
  do i = 1, ntot
    dblw(i) = 0.0_dp
    do d = 1, nd
      dblw(i) = dblw(i) + Haug(ns+d, i)**2
    end do
  end do

  ! ===================================================================
  ! (5) adiabatic kernel: g_xc = 0  ->  spectrum = eigvals(A0), only ns roots
  ! ===================================================================
  allocate(adiab(ns), Acopy(ns,ns), awork(max(64*ns,256)))
  Acopy = A0
  call dsyev('N','U', ns, Acopy, ns, adiab, awork, size(awork), info)
  if (info /= 0) stop 'dsyev A0 failed'

  ! ===================================================================
  ! (1-3) DRESSED roots via the pole-cancelled secular function
  ! ===================================================================
  ! Gershgorin window of Haug brackets the whole spectrum.
  lo =  1.0d30; hi = -1.0d30
  ! rebuild Haug values (dsyev destroyed it) just to take Gershgorin bounds
  block
    real(dp) :: H2(ntot,ntot), centre, rad
    integer  :: ii, jj
    H2 = 0.0_dp
    H2(1:ns,1:ns) = A0
    do d = 1, nd
      H2(1:ns, ns+d) = V(:,d); H2(ns+d, 1:ns) = V(:,d); H2(ns+d, ns+d) = omega_d(d)
    end do
    do ii = 1, ntot
      rad = 0.0_dp
      do jj = 1, ntot
        if (jj /= ii) rad = rad + abs(H2(ii,jj))
      end do
      centre = H2(ii,ii)
      lo = min(lo, centre - rad)
      hi = max(hi, centre + rad)
    end do
  end block
  lo = lo - 1.0_dp; hi = hi + 1.0_dp

  ! fine uniform scan + bisection
  nscan = max(20000, 2000*ntot)
  ngrid = nscan
  allocate(xs(ngrid), fv(ngrid))
  do i = 1, ngrid
    xs(i) = lo + (hi - lo)*real(i-1,dp)/real(ngrid-1,dp)
    fv(i) = fsec(ns, nd, xs(i), A0, V, omega_d)
  end do

  nr = 0
  dedup = 1.0d-10
  do i = 1, ngrid-1
    if (fv(i) == 0.0_dp) then
      call push_root(xs(i))
    else if (fv(i)*fv(i+1) < 0.0_dp) then
      fa = fv(i)
      call push_root(root_bisect(ns, nd, xs(i), xs(i+1), fa, A0, V, omega_d))
    end if
  end do

  call dsortvec(dressed, nr)

  ! ---- compare dressed vs reference dressed, reference exact, own exact ----
  dmaxRef = 0.0_dp; dmaxExF = 0.0_dp; admaxRef = 0.0_dp
  if (nr == ntot) then
    do i = 1, ntot
      dmaxRef = max(dmaxRef, abs(dressed(i) - refDressed(i)))
      dmaxExF = max(dmaxExF, abs(dressed(i) - exactF(i)))
    end do
  end if
  do i = 1, ntot
    admaxRef = max(admaxRef, abs(refExact(i) - exactF(i)))   ! NumPy exact vs Fortran exact
  end do

  ! ===================================================================
  ! GATE 2 metric: adiabatic misses the doubles (mirrors NumPy gate 2)
  ! ===================================================================
  ! Pick the Nd EXACT states with the LARGEST doubles-sector weight (the
  ! predominantly-double 0OS states). For THOSE, the nearest ADIABATIC root
  ! (eigval of A0) is far (>1e-3 -> missed), while the nearest DRESSED root is
  ! at machine precision (injected). A0 has only ns roots, so ntot-ns states are
  ! structurally absent from the adiabatic spectrum regardless.
  worst_adiab_gap = 0.0_dp; worst_dressed_gap = 0.0_dp
  block
    integer  :: k, jj, nmiss, pick
    real(dp) :: wtmp(ntot), ga, gd, biggest
    logical  :: used(ntot)
    ! sort states by descending doubles-weight (simple selection)
    wtmp = dblw
    used = .false.
    nmiss = 0
    do k = 1, nd                       ! the Nd most-doubly-excited states
      biggest = -1.0_dp; pick = 0
      do jj = 1, ntot
        if (.not. used(jj) .and. wtmp(jj) > biggest) then
          biggest = wtmp(jj); pick = jj
        end if
      end do
      used(pick) = .true.
      EE = exactF(pick)
      ga = 1.0d30
      do j = 1, ns; ga = min(ga, abs(adiab(j) - EE)); end do
      gd = 1.0d30
      do j = 1, nr; gd = min(gd, abs(dressed(j) - EE)); end do
      worst_adiab_gap   = max(worst_adiab_gap, ga)
      worst_dressed_gap = max(worst_dressed_gap, gd)
      if (ga > 1.0d-3) nmiss = nmiss + 1
    end do
    dblgap = real(nmiss, dp)           ! how many of the Nd doubles the adiabatic kernel misses
  end block

  ! ===================================================================
  ! report
  ! ===================================================================
  print '(a)',        "==== QMRSF-DK dressed-kernel pole search (Fortran, LIVE-shaped) ===="
  print '(a,i0,a,i0,a,i0)', "  Ns(single-spin-flip)=", ns, "  Nd(0OS doubles)=", nd, &
       "  total roots=", ntot
  print '(a)',        "  bare 0OS doubles omega_d:"
  print '(6f12.6)',   (omega_d(d), d=1,nd)
  print '(a)',        "  EXACT spectrum (Fortran dsyev of augmented matrix):"
  print '(6f12.6)',   (exactF(i), i=1,ntot)
  print '(a)',        "  DRESSED roots (pole-cancelled secular search):"
  print '(6f12.6)',   (dressed(i), i=1,nr)
  print '(a)',        "  ADIABATIC roots (g_xc=0 -> eigvals of A0 only):"
  print '(6f12.6)',   (adiab(i), i=1,ns)
  print '(a,i0,a,i0)', "  root count: dressed=", nr, "  expected=", ntot
  print '(a,es10.2)', "  max|dressed - refDressed(NumPy)| = ", dmaxRef
  print '(a,es10.2)', "  max|dressed - exact(Fortran)|    = ", dmaxExF
  print '(a,es10.2)', "  max|refExact(NumPy) - exact(Fort)|= ", admaxRef
  print '(a,i0,a,i0,a)', "  of the ", nd, " most-doubly-excited 0OS states, adiabatic MISSES ", &
       int(dblgap), " (>1e-3 from any A0 eigval; A0 has only ns roots)"
  print '(a,es10.2)', "  worst 0OS-double gap to nearest ADIABATIC root = ", worst_adiab_gap
  print '(a,es10.2)', "  worst 0OS-double gap to nearest DRESSED   root = ", worst_dressed_gap

  ! ---- gates ----
  pass1 = (nr == ntot) .and. (dmaxRef < 1.0d-9) .and. (dmaxExF < 1.0d-9) &
          .and. (admaxRef < 1.0d-9)
  pass2 = (int(dblgap) == nd) .and. (worst_adiab_gap > 1.0d-3) &
          .and. (worst_dressed_gap < 1.0d-9)

  print '(a)', "  ----------------------------------------------------------------"
  if (pass1) then
    print '(a,es10.2,a)', "  GATE 1: PASS  (dressed == exact == NumPy ref to ", &
         max(dmaxRef, dmaxExF, admaxRef), ", all roots)"
  else
    print '(a)', "  GATE 1: FAIL"
  end if
  if (pass2) then
    print '(a,i0,a)', "  GATE 2: PASS  (adiabatic misses all ", nd, &
         " 0OS doubles; dressed injects them <1e-9)"
  else
    print '(a)', "  GATE 2: FAIL"
  end if
  if (pass1 .and. pass2) then
    print '(a)', "  RESULT: PASS  (Fortran DK matches NumPy + own exact to <1e-9;"
    print '(a)', "                 adiabatic kernel structurally misses the 0OS doubles)"
  else
    print '(a)', "  RESULT: FAIL"
  end if

contains
  !> append a root to `dressed`, skipping near-duplicates (within dedup).
  subroutine push_root(r)
    real(dp), intent(in) :: r
    if (nr > 0) then
      if (abs(r - dressed(nr)) < dedup) return
    end if
    if (nr < ntot) then
      nr = nr + 1
      dressed(nr) = r
    end if
  end subroutine push_root
end program test_dk_core
