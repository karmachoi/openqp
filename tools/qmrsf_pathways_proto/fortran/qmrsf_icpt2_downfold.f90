! QMRSF-icPT2 external-Q self-energy downfold (standalone, build-testable).
! Faithful Fortran port of the validated NumPy multistate prototype
! (qmrsf_icpt2_multistate.py: icpt2_multistate / dyall_denoms). Reads the raw
! matrices from icpt2_downfold_ref.dat, independently
!   (1) diagonalizes H_PP                 -> CAS roots eP, vectors cP
!   (2) contracts couplings  coup = H_QP cP[:, :nP]      (internal contraction)
!   (3) builds per-root Q denominators    d_qk = eP[k] - H0_qk   (EN or Dyall)
!   (4) assembles the des Cloizeaux symmetric effective Hamiltonian
!         H_eff[k,l] = delta_kl eP[k]
!                      + 1/2 sum_q coup_qk coup_ql (1/d_qk + 1/d_ql)
!   (5) diagonalizes H_eff -> dressed energies
! and compares the EN and Dyall dressed spectra to the NumPy reference.
!
! This validates the downfold ALGEBRA in isolation; the production module wires
! the same assembly (steps 4-5, the reusable kernel) to live OpenQP perturbers.
module qmrsf_icpt2_downfold
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains

  !> Near-singular-safe reciprocal matching the prototype's guard.
  elemental real(dp) function safe_inv(d) result(r)
    real(dp), intent(in) :: d
    real(dp) :: dd
    dd = d
    if (abs(dd) < 1.0e-6_dp) dd = sign(1.0e-6_dp, dd) + 1.0e-30_dp
    r = 1.0_dp / dd
  end function safe_inv

  !> des Cloizeaux symmetric multistate effective Hamiltonian + spectrum.
  !>   eP(nP)        CAS roots (lowest nP)
  !>   coup(nQ,nP)   internally-contracted couplings <q|H|Psi_P^k>
  !>   invd(nQ,nP)   1/(eP[k]-H0_qk)
  subroutine icpt2_eff_hamiltonian(nP, nQ, eP, coup, invd, Edressed, herm)
    integer,  intent(in)  :: nP, nQ
    real(dp), intent(in)  :: eP(nP), coup(nQ,nP), invd(nQ,nP)
    real(dp), intent(out) :: Edressed(nP)
    real(dp), intent(out) :: herm
    real(dp) :: Heff(nP,nP), s, work(64*64)
    integer  :: k, l, q, lwork, info
    Heff = 0.0_dp
    do k = 1, nP
      Heff(k,k) = eP(k)
    end do
    do k = 1, nP
      do l = 1, nP
        s = 0.0_dp
        do q = 1, nQ
          s = s + coup(q,k)*coup(q,l)*0.5_dp*(invd(q,k)+invd(q,l))
        end do
        Heff(k,l) = Heff(k,l) + s
      end do
    end do
    herm = 0.0_dp
    do k = 1, nP
      do l = 1, nP
        herm = max(herm, abs(Heff(k,l)-Heff(l,k)))
      end do
    end do
    Heff = 0.5_dp*(Heff + transpose(Heff))
    lwork = size(work)
    call dsyev('N','U', nP, Heff, nP, Edressed, work, lwork, info)
    if (info /= 0) stop 'dsyev Heff failed'
  end subroutine icpt2_eff_hamiltonian

end module qmrsf_icpt2_downfold

program test_icpt2_downfold
  use qmrsf_icpt2_downfold
  implicit none
  integer :: dimP, dimQ, nP, i, j, k, info, lwork
  real(dp), allocatable :: HPP(:,:), HQP(:,:), Hqq(:), H0dy(:,:)
  real(dp), allocatable :: eP_all(:), cP(:,:), eP(:), coup(:,:), invd(:,:)
  real(dp), allocatable :: refEN(:), refDy(:), edEN(:), edDy(:), work(:)
  real(dp) :: hermEN, hermDy, dEN, dDy

  open(10, file="icpt2_downfold_ref.dat", status="old", action="read")
  read(10,*) dimP, dimQ, nP
  allocate(HPP(dimP,dimP), HQP(dimQ,dimP), Hqq(dimQ), H0dy(dimQ,nP))
  allocate(eP_all(dimP), cP(dimP,dimP), eP(nP), coup(dimQ,nP), invd(dimQ,nP))
  allocate(refEN(nP), refDy(nP), edEN(nP), edDy(nP))
  do i = 1, dimP; read(10,*) (HPP(i,j), j=1,dimP); end do
  do i = 1, dimQ; read(10,*) (HQP(i,j), j=1,dimP); end do
  read(10,*) (Hqq(i), i=1,dimQ)
  do i = 1, dimQ; read(10,*) (H0dy(i,k), k=1,nP); end do
  read(10,*) (refEN(i), i=1,nP)
  read(10,*) (refDy(i), i=1,nP)
  close(10)

  ! (1) diagonalize H_PP -> eP_all, cP (eigenvectors in columns)
  cP = HPP
  lwork = 64*dimP
  allocate(work(lwork))
  call dsyev('V','U', dimP, cP, dimP, eP_all, work, lwork, info)
  if (info /= 0) stop 'dsyev HPP failed'
  eP = eP_all(1:nP)

  ! (2) contract couplings coup(q,k) = sum_j HQP(q,j) cP(j,k)
  coup = matmul(HQP, cP(:,1:nP))

  ! (3) EN denominators, (4-5) downfold
  do k = 1, nP
    invd(:,k) = safe_inv(eP(k) - Hqq(:))
  end do
  call icpt2_eff_hamiltonian(nP, dimQ, eP, coup, invd, edEN, hermEN)

  ! (3') Dyall denominators, (4-5) downfold
  do k = 1, nP
    invd(:,k) = safe_inv(eP(k) - H0dy(:,k))
  end do
  call icpt2_eff_hamiltonian(nP, dimQ, eP, coup, invd, edDy, hermDy)

  dEN = maxval(abs(edEN - refEN))
  dDy = maxval(abs(edDy - refDy))

  print '(a)', "==== QMRSF-icPT2 downfold (Fortran) vs NumPy multistate reference ===="
  print '(a,i0,a,i0,a,i0)', "  dimP=", dimP, "  dimQ=", dimQ, "  nP(dressed)=", nP
  print '(a,4f14.8)', "  CAS roots eP      = ", eP
  print '(a,4f14.8)', "  dressed EN (Fort) = ", edEN
  print '(a,4f14.8)', "  dressed EN (ref)  = ", refEN
  print '(a,es10.2)', "  max|EN - ref|     = ", dEN
  print '(a,4f14.8)', "  dressed Dy (Fort) = ", edDy
  print '(a,4f14.8)', "  dressed Dy (ref)  = ", refDy
  print '(a,es10.2)', "  max|Dy - ref|     = ", dDy
  print '(a,es10.2,es10.2)', "  Heff herm (EN,Dy) = ", hermEN, hermDy
  if (dEN < 1.0d-9 .and. dDy < 1.0d-9) then
     print '(a)', "  RESULT: PASS  (Fortran downfold matches NumPy to <1e-9, EN and Dyall)"
  else
     print '(a)', "  RESULT: FAIL"
  end if
end program test_icpt2_downfold
