! QMRSF backbone core (standalone, build-testable) -- the DSF-generated CAS(4,4) Ms=0
! determinant Hamiltonian by the Slater-Condon rules, a faithful Fortran port of the
! validated NumPy reference (qmrsf_icpt2_ppp_proto.py: gen_dets/spinorb/melem/build_H).
! Validated against that reference on identical integrals (qmrsf_cas_ref.dat).
! This is the in-space backbone that QMRSF-icPT2 (external-Q downfold) and the CSF
! spin-adaptation build on. NO OpenQP dependency: integrals are read from file so the
! Slater-Condon machinery can be unit-tested in isolation before OpenQP integral wiring.
module qmrsf_backbone
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: NACT = 4, NSO = 8, NDET = 36
contains

  subroutine build_spinorb(h_act, eri_act, H1, g)
    real(dp), intent(in)  :: h_act(NACT,NACT), eri_act(NACT,NACT,NACT,NACT)
    real(dp), intent(out) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: P,Q,R,S, spat(NSO), spin(NSO), i
    real(dp) :: a, b
    do i = 1, NSO
      if (i <= NACT) then; spat(i) = i;        spin(i) = 0
      else;                spat(i) = i - NACT; spin(i) = 1; end if
    end do
    H1 = 0.0_dp
    do P = 1, NSO; do Q = 1, NSO
      if (spin(P) == spin(Q)) H1(P,Q) = h_act(spat(P), spat(Q))
    end do; end do
    ! antisymmetrized g(P,Q,R,S) = <PQ||RS> = (PR|QS) - (PS|QR), chemist eri, spin deltas
    g = 0.0_dp
    do P = 1, NSO; do Q = 1, NSO; do R = 1, NSO; do S = 1, NSO
      a = 0.0_dp; b = 0.0_dp
      if (spin(P)==spin(R) .and. spin(Q)==spin(S)) a = eri_act(spat(P),spat(R),spat(Q),spat(S))
      if (spin(P)==spin(S) .and. spin(Q)==spin(R)) b = eri_act(spat(P),spat(S),spat(Q),spat(R))
      g(P,Q,R,S) = a - b
    end do; end do; end do; end do
  end subroutine build_spinorb

  subroutine gen_dets(dets)
    integer, intent(out) :: dets(4,NDET)
    integer :: a1,a2,b1,b2,k, t(4), i,j,tmp
    k = 0
    do a1 = 1, NACT-1; do a2 = a1+1, NACT
      do b1 = 1, NACT-1; do b2 = b1+1, NACT
        k = k + 1
        t = (/ a1, a2, b1+NACT, b2+NACT /)
        do i = 1,3; do j = i+1,4
          if (t(j) < t(i)) then; tmp=t(i); t(i)=t(j); t(j)=tmp; end if
        end do; end do
        dets(:,k) = t
      end do; end do
    end do; end do
  end subroutine gen_dets

  pure logical function inset(x, D)
    integer, intent(in) :: x, D(4)
    inset = any(D == x)
  end function inset

  real(dp) function melem(D1, D2, H1, g)
    integer,  intent(in) :: D1(4), D2(4)
    real(dp), intent(in) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: holes(4), parts(4), common(4), nh, np, nc
    integer :: occ(4), nocc, i, idx, k
    integer :: p1,p2,ho1,ho2, Pp, Hh, Qc
    real(dp) :: sgn, val, e
    nh=0; np=0; nc=0
    do i=1,4
      if (.not. inset(D2(i), D1)) then; nh=nh+1; holes(nh)=D2(i); end if   ! in D2 not D1
    end do
    do i=1,4
      if (.not. inset(D1(i), D2)) then; np=np+1; parts(np)=D1(i); end if   ! in D1 not D2
      if (      inset(D1(i), D2)) then; nc=nc+1; common(nc)=D1(i); end if
    end do
    ! holes,parts,common are ascending (D1,D2 sorted)
    if (nh > 2) then; melem = 0.0_dp; return; end if
    ! phase: annihilate holes (ascending) then create parts (descending) on D2
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
    do k = np, 1, -1                                   ! create parts descending
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
      melem = e
    else if (nh == 1) then
      Pp = parts(1); Hh = holes(1)
      val = H1(Pp,Hh)
      do i = 1, nc; Qc = common(i); val = val + g(Pp,Qc,Hh,Qc); end do
      melem = sgn * val
    else
      p1=parts(1); p2=parts(2); ho1=holes(1); ho2=holes(2)
      melem = sgn * g(p1,p2,ho1,ho2)
    end if
  end function melem

  subroutine build_H(dets, H1, g, Hmat)
    integer,  intent(in)  :: dets(4,NDET)
    real(dp), intent(in)  :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    real(dp), intent(out) :: Hmat(NDET,NDET)
    integer :: i, j
    do i = 1, NDET
      do j = i, NDET
        Hmat(i,j) = melem(dets(:,i), dets(:,j), H1, g)
        Hmat(j,i) = Hmat(i,j)
      end do
    end do
  end subroutine build_H

end module qmrsf_backbone


program test_backbone
  use qmrsf_backbone
  implicit none
  real(dp) :: h_act(NACT,NACT), eri_act(NACT,NACT,NACT,NACT)
  real(dp) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO), Hmat(NDET,NDET)
  real(dp) :: ev(NDET), ref(NDET), work(512)
  integer  :: dets(4,NDET), p,q,r,s, i, nev, nr, info, lwork
  real(dp) :: herm, dmax

  open(10, file="qmrsf_cas_ref.dat", status="old", action="read")
  read(10,*) nr
  do p=1,NACT; read(10,*) (h_act(p,q), q=1,NACT); end do
  do p=1,NACT; do q=1,NACT; do r=1,NACT
     read(10,*) (eri_act(p,q,r,s), s=1,NACT)
  end do; end do; end do
  read(10,*) nev
  read(10,*) (ref(i), i=1,nev)
  close(10)

  call build_spinorb(h_act, eri_act, H1, g)
  call gen_dets(dets)
  call build_H(dets, H1, g, Hmat)

  herm = 0.0_dp
  do i=1,NDET; do p=1,NDET; herm = max(herm, abs(Hmat(i,p)-Hmat(p,i))); end do; end do

  lwork = 512
  call dsyev('N','U', NDET, Hmat, NDET, ev, work, lwork, info)

  dmax = 0.0_dp
  do i=1,NDET; dmax = max(dmax, abs(ev(i)-ref(i))); end do

  print '(a)', "==== QMRSF backbone core (Fortran) vs NumPy reference ===="
  print '(a,i0)',      "  ndet                 = ", NDET
  print '(a,es12.3)',  "  Hermiticity |H-H^T|  = ", herm
  print '(a,i0)',      "  dsyev info           = ", info
  print '(a,f16.10,a,f16.10)', "  ground state  Fortran= ", ev(1), "   ref= ", ref(1)
  print '(a,es12.3)',  "  max|ev - ref| (36)   = ", dmax
  if (dmax < 1.0d-9 .and. herm < 1.0d-10 .and. info == 0) then
     print '(a)', "  RESULT: PASS  (Fortran determinant-CI matches NumPy to <1e-9)"
  else
     print '(a)', "  RESULT: FAIL"
  end if
end program test_backbone
