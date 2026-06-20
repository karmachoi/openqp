! QMRSF-icPT2 external-Q downfold (standalone, build-testable) -- Stage A (WFT dressing).
! Reads the P/Q blocks of the CAS(4,4)+external model (dump_icpt2_ref.py) and forms the
! state-specific second-order self-energy downfold on the spin-pure CAS backbone ground state:
!   H_eff = H_PP + Sigma,   Sigma = sum_{q in Q} |<q|H|Psi_P>|^2 / (E_P - H_qq)   (Epstein-Nesbet)
! Validated against the NumPy icPT2 reference. This is the WFT-pathway dressing of the backbone;
! decoupled from OpenQP (blocks read from file) so the downfold algebra is unit-tested in isolation.
program qmrsf_icpt2
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: nP, nQ, i, j, info, lwork
  real(dp) :: ECASref, Eicref, sigref, ECAS, sigma, denom
  real(dp), allocatable :: HPP(:,:), HPQ(:,:), Hqq(:), ev(:), c(:), coup(:), work(:)

  open(10, file="qmrsf_icpt2_ref.dat", status="old", action="read")
  read(10,*) nP, nQ
  allocate(HPP(nP,nP), HPQ(nP,nQ), Hqq(nQ), ev(nP), c(nP), coup(nQ))
  do i=1,nP; read(10,*) (HPP(i,j), j=1,nP); end do
  do i=1,nP; read(10,*) (HPQ(i,j), j=1,nQ); end do
  read(10,*) (Hqq(j), j=1,nQ)
  read(10,*) ECASref, Eicref, sigref
  close(10)

  lwork = max(1, 3*nP); allocate(work(lwork))
  call dsyev('V','U', nP, HPP, nP, ev, work, lwork, info)   ! HPP -> eigenvectors
  ECAS = ev(1)                                              ! ground CAS energy
  c = HPP(:,1)                                              ! ground CAS eigenvector (over P)
  coup = matmul(transpose(HPQ), c)                          ! <q|H|Psi_P> = (H_QP c)_q
  sigma = 0.0_dp
  do j=1,nQ
    denom = ECAS - Hqq(j)
    if (abs(denom) < 1.0d-6) denom = sign(1.0d-6, denom)    ! intruder guard (none expected for ground)
    sigma = sigma + coup(j)**2 / denom
  end do

  print '(a)',        "==== QMRSF-icPT2 external-Q downfold (Fortran) vs NumPy ===="
  print '(a,i0,a,i0)',"  nP = ", nP, "   nQ = ", nQ
  print '(a,i0)',     "  dsyev info = ", info
  print '(a,f16.10,a,f16.10)', "  E_CAS  Fortran= ", ECAS,  "   ref= ", ECASref
  print '(a,f16.10,a,f16.10)', "  sigma  Fortran= ", sigma, "   ref= ", sigref
  print '(a,f16.10,a,f16.10)', "  icPT2  Fortran= ", ECAS+sigma, "   ref= ", Eicref
  print '(a,es12.3)', "  |icPT2 - ref| = ", abs(ECAS+sigma - Eicref)
  if (abs(ECAS-ECASref) < 1.0d-9 .and. abs(ECAS+sigma - Eicref) < 1.0d-9 .and. info==0) then
     print '(a)', "  RESULT: PASS  (Fortran external-Q downfold matches NumPy icPT2 to <1e-9)"
  else
     print '(a)', "  RESULT: FAIL"
  end if
end program
