!> @file tdhf_qmrsf_dk.F90
!> @brief QMRSF-DK dressed-kernel dynamic-correlation pathway (Pathway II) -- LIVE.
!>
!> @details
!>   QMRSF = Quintet Mixed-Reference Spin-Flip.  This module implements the
!>   *density-functional-picture* dynamic-correlation layer (the dressed /
!>   frequency-dependent quadratic exchange-correlation kernel g_xc(omega)) on
!>   top of the CAS(4,4) M_s=0 backbone of the quintet (S=2) ROHF reference.
!>
!>   PHYSICS (validated as prototypes -- see
!>   tools/qmrsf_pathways_proto/qmrsf_dk_live_proto.py,
!>   tools/qmrsf_pathways_proto/fortran/qmrsf_dk_core.f90):
!>
!>     The CAS(4,4) M_s=0 space has 36 = C(4,2)^2 determinants.  SIX of them are
!>     CLOSED-SHELL (0OS, "zero open shells": the alpha and beta spatial
!>     occupations coincide) -- these are the double-spin-flip configurations that
!>     an ADIABATIC (frequency-independent) kernel cannot reach.  The other THIRTY
!>     are open-shell (single / double-spin-flip with unpaired electrons) and form
!>     the adiabatic single-spin-flip response block A0.
!>
!>     We partition the CAS Hamiltonian H (Slater-Condon, spin-orbital basis):
!>         A0   = H restricted to the 30 open-shell determinants        (Ns=30)
!>         Wdd  = H restricted to the  6 closed-shell (0OS) determinants (Nd=6)
!>         Vc   = H coupling  (open-shell <-> 0OS)                       (30 x 6)
!>     The 0OS block Wdd is NOT diagonal, so we diagonalize it,
!>         Wdd U = U diag(omega_d),
!>     giving the SIX bare 0OS double-excitation energies omega_d and rotating the
!>     coupling  V = Vc U  into the 0OS eigenbasis.  Then the explicit augmented
!>     matrix  [[A0, V],[V^T, diag(omega_d)]]  is an EXACT (orthogonal) similarity
!>     transform of H, so it has the same 36 eigenvalues as the CAS Hamiltonian.
!>
!>     The dressed-kernel route downfolds the 0OS sector onto the singles sector
!>     analytically through the frequency-dependent quadratic kernel
!>         g_xc(omega)_{c c'} = sum_d V_{c,d} V_{c',d} / (omega - omega_d),
!>     and solves the pole-cancelled secular function
!>         fsec(omega) = det[ omega I - A0 - g_xc(omega) ] * prod_d (omega-omega_d)
!>     for all Ns+Nd roots.  Each pole of g_xc at omega_d INJECTS one 0OS
!>     double-excitation root the adiabatic kernel (g_xc=0) structurally misses.
!>
!>   CONSISTENCY (the gate): on BARE HF integrals the dressed-kernel DK spectrum
!>   equals the full CAS(4,4) spectrum to machine precision, because g_xc with the
!>   exact omega_d, V is the exact Feshbach downfold of the augmented matrix.  This
!>   establishes the LIVE pathway mechanism.
!>
!>   PHYSICS CAVEAT (scope): on HF integrals DK == CAS (a pure consistency check).
!>   The genuine DK *value* -- a KS/ROKS reference, an adiabatic-dressed singles
!>   block A0, and a DFT-DERIVED g_xc(omega) (the quadratic xc kernel, NOT the bare
!>   Coulomb/exchange coupling used here) -- is the NEXT layer and is NOT in scope.
!>   This module implements the HF-integral mechanism that establishes and
!>   validates the live pathway; swapping the integral source (KS A0 + DFT g_xc)
!>   for the next layer leaves the secular machinery below unchanged.
!>
!>   CONVENTIONS MIRRORED FROM:
!>     - source/modules/tdhf_qmrsf_icpt2.F90  (active-space determination from the
!>       quintet, qmrsf_active_integrals call, log/dump discipline).
!>     - source/modules/qmrsf_cas.F90         (the validated CAS det machinery,
!>       replicated here as the DK partition needs per-block access).
!>     - source/modules/tdhf_mrsf_energy.F90  (C-binding pattern, information
!>       handle via c_interop, print_module_info, log-file discipline).
module tdhf_qmrsf_dk_mod
  use precision, only: dp

  implicit none

  private
  public :: tdhf_qmrsf_dk_C
  public :: tdhf_qmrsf_dk

  character(len=*), parameter :: module_name = "tdhf_qmrsf_dk_mod"

  integer, parameter :: QMRSF_NACT = 4          !< CAS(4,4) active orbitals (SOMOs)
  integer, parameter :: NSO        = 8          !< spin orbitals (4 spatial x 2 spin)
  integer, parameter :: QMRSF_NDET = 36         !< C(4,2)^2 M_s=0 determinants
  integer, parameter :: NCLOSED    = 6          !< closed-shell (0OS) determinants
  integer, parameter :: NOPEN      = 30         !< open-shell determinants (Ns)

contains

  !> @brief C-bound entry point (matches the `void f(struct oqp_handle_t*)` ABI
  !>        declared in include/oqp.h and parsed by pyoqp via cffi).
  subroutine tdhf_qmrsf_dk_C(c_handle) bind(C, name="tdhf_qmrsf_dk")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_qmrsf_dk(inf)
  end subroutine tdhf_qmrsf_dk_C

  !> @brief Inner driver for the QMRSF dressed-kernel pathway (LIVE).
  !> @param[inout] infos  OpenQP run container.
  subroutine tdhf_qmrsf_dk(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use printing, only: print_module_info
    use qmrsf_ao2mo_mod, only: qmrsf_active_integrals
    use qmrsf_cas_mod, only: qmrsf_cas_solve

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_qmrsf_dk"

    type(information), target, intent(inout) :: infos

    integer :: nact, ncore, nbf, i, j, p, q, r, s
    integer, allocatable :: act(:)
    real(dp), allocatable :: h_act(:,:), eri_act(:,:,:,:)
    real(dp) :: ho1(QMRSF_NACT,QMRSF_NACT)
    real(dp) :: eri4(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp) :: ecore

    ! CAS partition
    real(dp) :: Hcas(QMRSF_NDET,QMRSF_NDET)
    integer  :: dets(4,QMRSF_NDET)
    integer  :: idx_open(NOPEN), idx_closed(NCLOSED)
    real(dp) :: A0(NOPEN,NOPEN), Vc(NOPEN,NCLOSED), Wdd(NCLOSED,NCLOSED)
    real(dp) :: Uw(NCLOSED,NCLOSED), omega_d(NCLOSED), V(NOPEN,NCLOSED)

    ! spectra
    real(dp) :: dressed(QMRSF_NDET), exactF(QMRSF_NDET)
    real(dp) :: adiab(NOPEN), cas_ref(QMRSF_NDET)
    real(dp) :: dblw(QMRSF_NDET)

    ! gate metrics
    real(dp) :: gate1_dk_cas, gate1_dk_exact, herm
    real(dp) :: worst_adiab_gap, worst_dressed_gap
    integer  :: nmiss, nr
    logical  :: pass1, pass2

    ! --- Open the main log file (append), matching the backbone discipline. ---
    open(unit=iw, file=infos%log_filename, position="append")

    call print_module_info('QMRSF_DK', &
         'Frequency-dependent dressed xc kernel g_xc(omega): 0OS doubles injection')

    ! ---- Active space from the quintet (S=2) reference ------------------------
    !   nelec_A = ncore + 4 , nelec_B = ncore. The four SOMOs are MOs
    !   ncore+1..ncore+4 (CAS(4,4) = the M_s=0 active space).
    nbf   = int(infos%basis%nbf)
    nact  = QMRSF_NACT
    ncore = int(infos%mol_prop%nelec_B)
    if (int(infos%mol_prop%nelec_A) - int(infos%mol_prop%nelec_B) /= nact) then
      write(iw,'(/,5x,a)') 'QMRSF-DK ERROR: reference is not a quintet (S=2) '// &
           'high-spin ROHF (need nelec_A - nelec_B = 4). Aborting pathway.'
      call flush(iw); close(iw); return
    end if

    allocate(act(nact))
    do i = 1, nact
      act(i) = ncore + i
    end do

    write(iw,'(/,5x,a,i0)') 'QMRSF-DK: basis functions          = ', nbf
    write(iw,'(5x,a,i0)')   'QMRSF-DK: frozen-core MOs          = ', ncore
    write(iw,'(5x,a,i0,a)') 'QMRSF-DK: CAS active (SOMOs)       = ', nact, ' (MOs ncore+1..ncore+4)'

    ! ---- Active-space integrals (same int2-reuse path as icPT2) --------------
    allocate(h_act(nact,nact), eri_act(nact,nact,nact,nact))
    call qmrsf_active_integrals(infos, nact, act, ncore, h_act, eri_act, ecore)
    ho1  = h_act
    eri4 = eri_act

    ! ---- CAS reference spectrum (the VALIDATION REFERENCE) -------------------
    call qmrsf_cas_solve(ho1, eri4, cas_ref, herm=herm)

    ! ---- Build the CAS Hamiltonian, classify 0OS, partition into A0/Vc/Wdd ---
    call dk_build_cas_partition(ho1, eri4, Hcas, dets, idx_open, idx_closed, &
                                A0, Vc, Wdd)

    ! ---- Diagonalize the 0OS block -> bare double energies omega_d + rotation -
    call dk_diag_sym(NCLOSED, Wdd, omega_d, Uw)
    !  V = Vc * Uw  (rotate the coupling into the 0OS eigenbasis)
    V = matmul(Vc, Uw)

    write(iw,'(/,5x,a,i0,a,i0,a,i0)') &
         'QMRSF-DK: Ns(open-shell singles block) = ', NOPEN, &
         ' ; Nd(0OS doubles) = ', NCLOSED, ' ; total roots = ', QMRSF_NDET
    write(iw,'(5x,a)') 'QMRSF-DK: bare 0OS double energies omega_d (electronic):'
    write(iw,'(7x,6f14.8)') (omega_d(i), i=1,NCLOSED)

    ! ---- (internal arbiter) explicit augmented diagonalization --------------
    !   [[A0, V],[V^T, diag(omega_d)]] -- exact, with doubles-sector weights.
    call dk_augmented_spectrum(A0, V, omega_d, exactF, dblw)

    ! ---- (deliverable) DRESSED secular / pole-cancelled root search ----------
    call dk_dressed_roots(A0, V, omega_d, dressed, nr)

    ! ---- ADIABATIC spectrum (g_xc = 0 -> eigvals of A0 only) -----------------
    call dk_adiabatic_spectrum(A0, adiab)

    ! ---- GATE 1: DK dressed roots == CAS == augmented-exact ------------------
    gate1_dk_cas   = 0.0_dp
    gate1_dk_exact = 0.0_dp
    if (nr == QMRSF_NDET) then
      do i = 1, QMRSF_NDET
        gate1_dk_cas   = max(gate1_dk_cas,   abs(dressed(i) - cas_ref(i)))
        gate1_dk_exact = max(gate1_dk_exact, abs(dressed(i) - exactF(i)))
      end do
    end if

    ! ---- GATE 2: adiabatic MISSES the 0OS doubles ---------------------------
    !   Among the Nd EXACT states with the largest doubles-sector weight, the
    !   nearest adiabatic root is far (>1e-3 -> missed) while the nearest dressed
    !   root is at machine precision (injected). A0 has only Ns roots, so Nd
    !   states are structurally absent from the adiabatic spectrum regardless.
    call dk_gate2_metrics(exactF, dblw, adiab, dressed, nr, &
                          nmiss, worst_adiab_gap, worst_dressed_gap)

    ! ---- report -------------------------------------------------------------
    write(iw,'(/,5x,a,f18.10)') 'QMRSF-DK: E_core (nuc + frozen core)   = ', ecore
    write(iw,'(5x,a,es10.2)')   'QMRSF-DK: CAS Hamiltonian |H-H^T|       = ', herm
    write(iw,'(5x,a,i0)')       'QMRSF-DK: dressed-kernel root count     = ', nr

    write(iw,'(/,5x,a)') 'QMRSF-DK: state    E_CAS(total)      E_DK(total)        E_adiab(total)'
    do i = 1, QMRSF_NDET
      if (i <= NOPEN) then
        write(iw,'(7x,i3,3f18.10)') i-1, cas_ref(i)+ecore, dressed(i)+ecore, adiab(i)+ecore
      else
        write(iw,'(7x,i3,2f18.10,a)') i-1, cas_ref(i)+ecore, dressed(i)+ecore, &
             '        (no adiabatic root: 0OS double)'
      end if
    end do

    write(iw,'(/,5x,a)') '================  QMRSF-DK validation gates  ================'
    write(iw,'(5x,a,es10.2)') 'GATE 1: max|E_DK - E_CAS|             = ', gate1_dk_cas
    write(iw,'(5x,a,es10.2)') 'GATE 1: max|E_DK - E_augmented-exact| = ', gate1_dk_exact
    pass1 = (nr == QMRSF_NDET) .and. (gate1_dk_cas < 1.0d-9) .and. (gate1_dk_exact < 1.0d-9)
    if (pass1) then
      write(iw,'(5x,a)') 'GATE 1: PASS  (dressed-kernel DK spectrum == CAS(4,4) to <1e-9)'
    else
      write(iw,'(5x,a)') 'GATE 1: FAIL'
    end if

    write(iw,'(5x,a,i0,a,i0,a)') 'GATE 2: adiabatic MISSES ', nmiss, ' of ', NCLOSED, &
         ' most-doubly-excited 0OS states (>1e-3 from any A0 eigval)'
    write(iw,'(5x,a,es10.2)') 'GATE 2: worst 0OS-double gap to nearest ADIABATIC root = ', worst_adiab_gap
    write(iw,'(5x,a,es10.2)') 'GATE 2: worst 0OS-double gap to nearest DRESSED   root = ', worst_dressed_gap
    pass2 = (nmiss == NCLOSED) .and. (worst_adiab_gap > 1.0d-3) .and. (worst_dressed_gap < 1.0d-9)
    if (pass2) then
      write(iw,'(5x,a,i0,a)') 'GATE 2: PASS  (adiabatic structurally misses all ', NCLOSED, &
           ' 0OS doubles; dressed kernel injects them <1e-9)'
    else
      write(iw,'(5x,a)') 'GATE 2: FAIL'
    end if
    write(iw,'(5x,a)') '------------------------------------------------------------'
    if (pass1 .and. pass2) then
      write(iw,'(5x,a)') 'QMRSF-DK RESULT: PASS (live HF-integral dressed-kernel pathway established)'
    else
      write(iw,'(5x,a)') 'QMRSF-DK RESULT: FAIL'
    end if
    write(iw,'(5x,a)') 'NOTE: on HF integrals DK==CAS by construction (consistency check).'
    write(iw,'(5x,a)') '      The DFT-dressed value (KS A0 + DFT-derived g_xc) is the next layer.'

    ! ---- validation dump (parsed by pyoqp -> JSON + log table) ---------------
    call dk_write_dump(ho1, eri4, ecore, omega_d, cas_ref, dressed, adiab, &
                       gate1_dk_cas, gate1_dk_exact, worst_adiab_gap, &
                       worst_dressed_gap, nmiss)
    write(iw,'(/,5x,a)') 'QMRSF-DK: wrote validation dump (qmrsf_dk_full_live.dat).'

    deallocate(act, h_act, eri_act)

    call flush(iw)
    close(iw)

  end subroutine tdhf_qmrsf_dk

!===============================================================================
! CAS(4,4) determinant machinery (replicated from qmrsf_cas.F90, validated;
! the DK partition needs explicit per-block access the public solver does not give)
!===============================================================================

  !> Build the full CAS Hamiltonian H, classify 0OS vs open-shell, and slice the
  !> A0 (open block), Vc (open<->0OS coupling), Wdd (0OS block) partition.
  subroutine dk_build_cas_partition(h_act, eri_act, Hmat, dets, &
                                    idx_open, idx_closed, A0, Vc, Wdd)
    real(dp), intent(in)  :: h_act(QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(in)  :: eri_act(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(out) :: Hmat(QMRSF_NDET,QMRSF_NDET)
    integer,  intent(out) :: dets(4,QMRSF_NDET)
    integer,  intent(out) :: idx_open(NOPEN), idx_closed(NCLOSED)
    real(dp), intent(out) :: A0(NOPEN,NOPEN), Vc(NOPEN,NCLOSED), Wdd(NCLOSED,NCLOSED)

    real(dp) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer  :: i, j, no, nc

    call dk_build_spinorb(h_act, eri_act, H1, g)
    call dk_gen_dets(dets)
    call dk_build_H(dets, H1, g, Hmat)

    ! classify: closed-shell (0OS) = alpha spatial set == beta spatial set.
    no = 0; nc = 0
    do i = 1, QMRSF_NDET
      if (dk_is_closed(dets(:,i))) then
        nc = nc + 1
        idx_closed(nc) = i
      else
        no = no + 1
        idx_open(no) = i
      end if
    end do

    ! slice the partition
    do i = 1, NOPEN
      do j = 1, NOPEN
        A0(i,j) = Hmat(idx_open(i), idx_open(j))
      end do
      do j = 1, NCLOSED
        Vc(i,j) = Hmat(idx_open(i), idx_closed(j))
      end do
    end do
    do i = 1, NCLOSED
      do j = 1, NCLOSED
        Wdd(i,j) = Hmat(idx_closed(i), idx_closed(j))
      end do
    end do
  end subroutine dk_build_cas_partition

  !> Is determinant D closed-shell (0OS): each occupied alpha spatial orbital
  !> also has its beta partner occupied (alpha spatial set == beta spatial set)?
  pure logical function dk_is_closed(D) result(closed)
    integer, intent(in) :: D(4)
    integer :: i, j, na, nb, sa(4), sb(4)
    na = 0; nb = 0
    do i = 1, 4
      if (D(i) <= QMRSF_NACT) then
        na = na + 1; sa(na) = D(i)
      else
        nb = nb + 1; sb(nb) = D(i) - QMRSF_NACT
      end if
    end do
    closed = .false.
    if (na == 2 .and. nb == 2) then
      ! M_s=0 here always has 2 alpha + 2 beta; closed iff the two spatial sets match
      closed = .true.
      do i = 1, 2
        if (.not. any(sb(1:2) == sa(i))) closed = .false.
      end do
    end if
  end function dk_is_closed

  subroutine dk_build_spinorb(h_act, eri_act, H1, g)
    real(dp), intent(in)  :: h_act(QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(in)  :: eri_act(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(out) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: P,Q,R,S, spat(NSO), spin(NSO), i
    real(dp) :: a, b
    do i = 1, NSO
      if (i <= QMRSF_NACT) then; spat(i) = i;             spin(i) = 0
      else;                      spat(i) = i - QMRSF_NACT; spin(i) = 1; end if
    end do
    H1 = 0.0_dp
    do P = 1, NSO; do Q = 1, NSO
      if (spin(P) == spin(Q)) H1(P,Q) = h_act(spat(P), spat(Q))
    end do; end do
    g = 0.0_dp
    do P = 1, NSO; do Q = 1, NSO; do R = 1, NSO; do S = 1, NSO
      a = 0.0_dp; b = 0.0_dp
      if (spin(P)==spin(R) .and. spin(Q)==spin(S)) a = eri_act(spat(P),spat(R),spat(Q),spat(S))
      if (spin(P)==spin(S) .and. spin(Q)==spin(R)) b = eri_act(spat(P),spat(S),spat(Q),spat(R))
      g(P,Q,R,S) = a - b
    end do; end do; end do; end do
  end subroutine dk_build_spinorb

  subroutine dk_gen_dets(dets)
    integer, intent(out) :: dets(4,QMRSF_NDET)
    integer :: a1,a2,b1,b2,k, t(4), i,j,tmp
    k = 0
    do a1 = 1, QMRSF_NACT-1; do a2 = a1+1, QMRSF_NACT
      do b1 = 1, QMRSF_NACT-1; do b2 = b1+1, QMRSF_NACT
        k = k + 1
        t = (/ a1, a2, b1+QMRSF_NACT, b2+QMRSF_NACT /)
        do i = 1,3; do j = i+1,4
          if (t(j) < t(i)) then; tmp=t(i); t(i)=t(j); t(j)=tmp; end if
        end do; end do
        dets(:,k) = t
      end do; end do
    end do; end do
  end subroutine dk_gen_dets

  pure logical function dk_inset(x, D)
    integer, intent(in) :: x, D(4)
    dk_inset = any(D == x)
  end function dk_inset

  real(dp) function dk_melem(D1, D2, H1, g)
    integer,  intent(in) :: D1(4), D2(4)
    real(dp), intent(in) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: holes(4), parts(4), common(4), nh, np, nc
    integer :: occ(4), nocc, i, idx, k
    integer :: p1,p2,ho1,ho2, Pp, Hh, Qc
    real(dp) :: sgn, val, e
    nh=0; np=0; nc=0
    do i=1,4
      if (.not. dk_inset(D2(i), D1)) then; nh=nh+1; holes(nh)=D2(i); end if
    end do
    do i=1,4
      if (.not. dk_inset(D1(i), D2)) then; np=np+1; parts(np)=D1(i); end if
      if (      dk_inset(D1(i), D2)) then; nc=nc+1; common(nc)=D1(i); end if
    end do
    if (nh > 2) then; dk_melem = 0.0_dp; return; end if
    occ = D2; nocc = 4; sgn = 1.0_dp
    do k = 1, nh
      idx = 0
      do i = 1, nocc
        if (occ(i) == holes(k)) then; idx = i; exit; end if
      end do
      if (mod(idx-1,2) == 1) sgn = -sgn
      do i = idx, nocc-1; occ(i) = occ(i+1); end do
      nocc = nocc - 1
    end do
    do k = np, 1, -1
      idx = 1
      do i = 1, nocc
        if (occ(i) < parts(k)) idx = idx + 1
      end do
      if (mod(idx-1,2) == 1) sgn = -sgn
      do i = nocc, idx, -1; occ(i+1) = occ(i); end do
      occ(idx) = parts(k); nocc = nocc + 1
    end do
    if (nh == 0) then
      e = 0.0_dp
      do i = 1,4; e = e + H1(D1(i),D1(i)); end do
      do i = 1,3
        do k = i+1,4; e = e + g(D1(i),D1(k),D1(i),D1(k)); end do
      end do
      dk_melem = e
    else if (nh == 1) then
      Pp = parts(1); Hh = holes(1)
      val = H1(Pp,Hh)
      do i = 1, nc; Qc = common(i); val = val + g(Pp,Qc,Hh,Qc); end do
      dk_melem = sgn * val
    else
      p1=parts(1); p2=parts(2); ho1=holes(1); ho2=holes(2)
      dk_melem = sgn * g(p1,p2,ho1,ho2)
    end if
  end function dk_melem

  subroutine dk_build_H(dets, H1, g, Hmat)
    integer,  intent(in)  :: dets(4,QMRSF_NDET)
    real(dp), intent(in)  :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp), intent(out) :: Hmat(QMRSF_NDET,QMRSF_NDET)
    integer :: i, j
    do i = 1, QMRSF_NDET
      do j = i, QMRSF_NDET
        Hmat(i,j) = dk_melem(dets(:,i), dets(:,j), H1, g)
        Hmat(j,i) = Hmat(i,j)
      end do
    end do
  end subroutine dk_build_H

!===============================================================================
! Dressed-kernel secular machinery (ported from qmrsf_dk_core.f90, validated)
!===============================================================================

  !> Symmetric eigensolve via the OpenQP eigen module (ILP64-safe), ascending.
  subroutine dk_diag_sym(n, Ain, eval, evec)
    use eigen, only: diag_symm_full
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: Ain(n,n)
    real(dp), intent(out) :: eval(n), evec(n,n)
    integer :: ierr
    evec = Ain
    call diag_symm_full(1, n, evec, n, eval, ierr)   ! mode 1 -> eigenvectors in evec
  end subroutine dk_diag_sym

  !> Explicit augmented matrix [[A0,V],[V^T,diag(omega_d)]]: exact spectrum (the
  !> internal arbiter) + doubles-sector weight of each eigenvector.
  subroutine dk_augmented_spectrum(A0, V, omega_d, eval, dblw)
    real(dp), intent(in)  :: A0(NOPEN,NOPEN), V(NOPEN,NCLOSED), omega_d(NCLOSED)
    real(dp), intent(out) :: eval(QMRSF_NDET), dblw(QMRSF_NDET)
    real(dp) :: Haug(QMRSF_NDET,QMRSF_NDET), evec(QMRSF_NDET,QMRSF_NDET)
    integer  :: d, i
    Haug = 0.0_dp
    Haug(1:NOPEN,1:NOPEN) = A0
    do d = 1, NCLOSED
      Haug(1:NOPEN, NOPEN+d) = V(:,d)
      Haug(NOPEN+d, 1:NOPEN) = V(:,d)
      Haug(NOPEN+d, NOPEN+d) = omega_d(d)
    end do
    call dk_diag_sym(QMRSF_NDET, Haug, eval, evec)
    do i = 1, QMRSF_NDET
      dblw(i) = 0.0_dp
      do d = 1, NCLOSED
        dblw(i) = dblw(i) + evec(NOPEN+d, i)**2
      end do
    end do
  end subroutine dk_augmented_spectrum

  !> Adiabatic spectrum: g_xc = 0 -> eigvals of A0 (only Ns roots).
  subroutine dk_adiabatic_spectrum(A0, adiab)
    use eigen, only: diag_symm_full
    real(dp), intent(in)  :: A0(NOPEN,NOPEN)
    real(dp), intent(out) :: adiab(NOPEN)
    real(dp) :: tmp(NOPEN,NOPEN)
    integer  :: ierr
    tmp = A0
    call diag_symm_full(0, NOPEN, tmp, NOPEN, adiab, ierr)
  end subroutine dk_adiabatic_spectrum

  !> Dressed kernel g_xc(omega)_{c c'} = sum_d V(c,d) V(c',d) / (omega - omega_d).
  subroutine dk_build_gxc(omega, V, omega_d, gx)
    real(dp), intent(in)  :: omega, V(NOPEN,NCLOSED), omega_d(NCLOSED)
    real(dp), intent(out) :: gx(NOPEN,NOPEN)
    integer  :: c, cp, d
    real(dp) :: denom
    gx = 0.0_dp
    do d = 1, NCLOSED
      denom = omega - omega_d(d)
      do cp = 1, NOPEN
        do c = 1, NOPEN
          gx(c,cp) = gx(c,cp) + V(c,d)*V(cp,d)/denom
        end do
      end do
    end do
  end subroutine dk_build_gxc

  !> det of a general n x n matrix via LU (LAPACK dgetrf). Matrix destroyed.
  real(dp) function dk_det_lu(n, amat) result(det)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: amat(n,n)
    integer :: ipiv(n), info, i
    call dgetrf(n, n, amat, n, ipiv, info)
    if (info < 0) then; det = 0.0_dp; return; end if
    det = 1.0_dp
    do i = 1, n
      det = det * amat(i,i)
      if (ipiv(i) /= i) det = -det     ! row swap flips the sign
    end do
  end function dk_det_lu

  !> Pole-cancelled secular function:
  !>   fsec(omega) = det[ omega I - A0 - g_xc(omega) ] * prod_d (omega - omega_d).
  real(dp) function dk_fsec(omega, A0, V, omega_d) result(fval)
    real(dp), intent(in) :: omega, A0(NOPEN,NOPEN), V(NOPEN,NCLOSED), omega_d(NCLOSED)
    real(dp) :: gx(NOPEN,NOPEN), amat(NOPEN,NOPEN), poleprod
    integer  :: i, d
    call dk_build_gxc(omega, V, omega_d, gx)
    amat = -A0 - gx
    do i = 1, NOPEN
      amat(i,i) = amat(i,i) + omega
    end do
    poleprod = 1.0_dp
    do d = 1, NCLOSED
      poleprod = poleprod * (omega - omega_d(d))
    end do
    fval = dk_det_lu(NOPEN, amat) * poleprod
  end function dk_fsec

  !> Bisect a root of fsec in [a,b] given a sign change (fa = fsec(a)).
  real(dp) function dk_root_bisect(a0in, b0in, faIn, A0, V, omega_d) result(rt)
    real(dp), intent(in) :: a0in, b0in, faIn, A0(NOPEN,NOPEN), V(NOPEN,NCLOSED), omega_d(NCLOSED)
    real(dp) :: a, b, fa, m, fm
    integer  :: it
    a = a0in; b = b0in; fa = faIn
    do it = 1, 300
      m = 0.5_dp*(a+b)
      fm = dk_fsec(m, A0, V, omega_d)
      if (b - a < 1.0d-14 .or. fm == 0.0_dp) exit
      if (fm*fa > 0.0_dp) then
        a = m; fa = fm
      else
        b = m
      end if
    end do
    rt = 0.5_dp*(a+b)
  end function dk_root_bisect

  !> Solve the dressed eigenvalue condition for all Ns+Nd roots via the
  !> pole-cancelled secular function (fine sign scan + bisection over the
  !> Gershgorin window of the augmented matrix). Mirrors qmrsf_dk_core.f90.
  subroutine dk_dressed_roots(A0, V, omega_d, dressed, nr)
    real(dp), intent(in)  :: A0(NOPEN,NOPEN), V(NOPEN,NCLOSED), omega_d(NCLOSED)
    real(dp), intent(out) :: dressed(QMRSF_NDET)
    integer,  intent(out) :: nr

    real(dp) :: lo, hi, centre, rad, fa
    real(dp) :: x_prev, f_prev, x_cur, f_cur
    real(dp), parameter :: dedup = 1.0d-10
    integer  :: i, j, d, ngrid
    real(dp) :: Haug(QMRSF_NDET,QMRSF_NDET)

    ! Gershgorin window of the augmented matrix brackets the whole spectrum.
    Haug = 0.0_dp
    Haug(1:NOPEN,1:NOPEN) = A0
    do d = 1, NCLOSED
      Haug(1:NOPEN, NOPEN+d) = V(:,d)
      Haug(NOPEN+d, 1:NOPEN) = V(:,d)
      Haug(NOPEN+d, NOPEN+d) = omega_d(d)
    end do
    lo =  1.0d30; hi = -1.0d30
    do i = 1, QMRSF_NDET
      rad = 0.0_dp
      do j = 1, QMRSF_NDET
        if (j /= i) rad = rad + abs(Haug(i,j))
      end do
      centre = Haug(i,i)
      lo = min(lo, centre - rad)
      hi = max(hi, centre + rad)
    end do
    lo = lo - 1.0_dp; hi = hi + 1.0_dp

    ! fine uniform scan + bisection (resolution finer than the closest spacing)
    ngrid = max(20000, 2000*QMRSF_NDET)
    nr = 0
    x_prev = lo
    f_prev = dk_fsec(x_prev, A0, V, omega_d)
    do i = 2, ngrid
      x_cur = lo + (hi - lo)*real(i-1,dp)/real(ngrid-1,dp)
      f_cur = dk_fsec(x_cur, A0, V, omega_d)
      if (f_prev == 0.0_dp) then
        call dk_push_root(dressed, nr, x_prev, dedup)
      else if (f_prev*f_cur < 0.0_dp) then
        fa = f_prev
        call dk_push_root(dressed, nr, &
             dk_root_bisect(x_prev, x_cur, fa, A0, V, omega_d), dedup)
      end if
      x_prev = x_cur; f_prev = f_cur
    end do

    call dk_sortvec(dressed, nr)
  end subroutine dk_dressed_roots

  subroutine dk_push_root(dressed, nr, r, dedup)
    real(dp), intent(inout) :: dressed(QMRSF_NDET)
    integer,  intent(inout) :: nr
    real(dp), intent(in)    :: r, dedup
    if (nr > 0) then
      if (abs(r - dressed(nr)) < dedup) return
    end if
    if (nr < QMRSF_NDET) then
      nr = nr + 1
      dressed(nr) = r
    end if
  end subroutine dk_push_root

  subroutine dk_sortvec(a, m)
    real(dp), intent(inout) :: a(*)
    integer,  intent(in)    :: m
    integer  :: i, j
    real(dp) :: t
    do i = 1, m-1
      do j = 1, m-i
        if (a(j) > a(j+1)) then; t=a(j); a(j)=a(j+1); a(j+1)=t; end if
      end do
    end do
  end subroutine dk_sortvec

  !> GATE 2 metric: among the Nd EXACT states with the largest doubles-sector
  !> weight, count how many the adiabatic kernel misses (>1e-3 from any A0
  !> eigval) and report the worst adiabatic / dressed gaps.
  subroutine dk_gate2_metrics(exactF, dblw, adiab, dressed, nr, &
                              nmiss, worst_adiab_gap, worst_dressed_gap)
    real(dp), intent(in)  :: exactF(QMRSF_NDET), dblw(QMRSF_NDET)
    real(dp), intent(in)  :: adiab(NOPEN), dressed(QMRSF_NDET)
    integer,  intent(in)  :: nr
    integer,  intent(out) :: nmiss
    real(dp), intent(out) :: worst_adiab_gap, worst_dressed_gap
    logical  :: used(QMRSF_NDET)
    real(dp) :: biggest, EE, ga, gd
    integer  :: k, jj, j, pick
    worst_adiab_gap = 0.0_dp; worst_dressed_gap = 0.0_dp; nmiss = 0
    used = .false.
    do k = 1, NCLOSED                       ! the Nd most-doubly-excited states
      biggest = -1.0_dp; pick = 0
      do jj = 1, QMRSF_NDET
        if (.not. used(jj) .and. dblw(jj) > biggest) then
          biggest = dblw(jj); pick = jj
        end if
      end do
      used(pick) = .true.
      EE = exactF(pick)
      ga = 1.0d30
      do j = 1, NOPEN; ga = min(ga, abs(adiab(j) - EE)); end do
      gd = 1.0d30
      do j = 1, nr;    gd = min(gd, abs(dressed(j) - EE)); end do
      worst_adiab_gap   = max(worst_adiab_gap,   ga)
      worst_dressed_gap = max(worst_dressed_gap, gd)
      if (ga > 1.0d-3) nmiss = nmiss + 1
    end do
  end subroutine dk_gate2_metrics

  !> Validation dump consumed by pyoqp (oqp/library/qmrsf_results.py).
  !> Format mirrors qmrsf_icpt2_full_live.dat but is DK-specific.
  subroutine dk_write_dump(h_act, eri4, ecore, omega_d, cas_ref, dressed, adiab, &
                           g1_cas, g1_exact, gap_adiab, gap_dressed, nmiss)
    real(dp), intent(in) :: h_act(QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(in) :: eri4(QMRSF_NACT,QMRSF_NACT,QMRSF_NACT,QMRSF_NACT)
    real(dp), intent(in) :: ecore, omega_d(NCLOSED)
    real(dp), intent(in) :: cas_ref(QMRSF_NDET), dressed(QMRSF_NDET), adiab(NOPEN)
    real(dp), intent(in) :: g1_cas, g1_exact, gap_adiab, gap_dressed
    integer,  intent(in) :: nmiss
    integer :: u, p, q, r, s, i
    u = 94
    open(unit=u, file='qmrsf_dk_full_live.dat', status='replace', action='write')
    !  header: nact, ndet, nopen, nclosed
    write(u,'(i0,1x,i0,1x,i0,1x,i0)') QMRSF_NACT, QMRSF_NDET, NOPEN, NCLOSED
    !  active integrals (for the AO oracle gate, same layout as icPT2)
    do p = 1, QMRSF_NACT; write(u,'(*(es24.16))') (h_act(p,q), q=1,QMRSF_NACT); end do
    do p = 1, QMRSF_NACT; do q = 1, QMRSF_NACT; do r = 1, QMRSF_NACT
      write(u,'(*(es24.16))') (eri4(p,q,r,s), s=1,QMRSF_NACT)
    end do; end do; end do
    write(u,'(es24.16)') ecore
    write(u,'(*(es24.16))') (omega_d(i),  i=1,NCLOSED)     ! bare 0OS double energies
    write(u,'(*(es24.16))') (cas_ref(i),  i=1,QMRSF_NDET)  ! CAS reference (electronic)
    write(u,'(*(es24.16))') (dressed(i),  i=1,QMRSF_NDET)  ! DK dressed (electronic)
    write(u,'(*(es24.16))') (adiab(i),    i=1,NOPEN)       ! adiabatic (electronic)
    !  gate metrics
    write(u,'(*(es24.16))') g1_cas, g1_exact, gap_adiab, gap_dressed
    write(u,'(i0)') nmiss
    close(u)
  end subroutine dk_write_dump

end module tdhf_qmrsf_dk_mod
