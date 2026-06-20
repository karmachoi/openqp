!> @file qmrsf_icpt2_downfold.F90
!> @brief QMRSF-icPT2 external-Q self-energy downfold (des Cloizeaux multistate).
!>
!> @details
!>   Given the CAS roots eP, the internally-contracted external-Q couplings
!>   coup_qk = <q|H|Psi_P^k>, and per-root zeroth-order denominators
!>   1/(eP_k - H0_qk) (Epstein-Nesbet or Dyall), assemble the state-independent
!>   Hermitian effective Hamiltonian (des Cloizeaux / QDPT symmetric average)
!>
!>     H_eff[k,l] = delta_kl eP_k
!>                  + 1/2 sum_q coup_qk coup_ql ( 1/(eP_k-H0_qk) + 1/(eP_l-H0_ql) )
!>
!>   and diagonalize it to give the dressed (icPT2-corrected) state energies.
!>   This is the reusable downfold KERNEL: it is agnostic to how the couplings
!>   and denominators were produced (full-CI brute force on a small system, or
!>   the contracted NEVPT2-style perturber generation for production).
!>
!>   Validated to <1e-13 (EN and Dyall) against the NumPy multistate prototype
!>   (tools/qmrsf_pathways_proto/qmrsf_icpt2_multistate.py) by the standalone
!>   tools/qmrsf_pathways_proto/fortran/qmrsf_icpt2_downfold.f90.
module qmrsf_icpt2_downfold_mod
  use precision, only: dp
  implicit none
  private
  public :: icpt2_eff_hamiltonian, icpt2_safe_inv

contains

  !> @brief Near-singular-safe reciprocal (intruder guard), matching the proto.
  elemental real(dp) function icpt2_safe_inv(d) result(r)
    real(dp), intent(in) :: d
    real(dp) :: dd
    dd = d
    if (abs(dd) < 1.0e-6_dp) dd = sign(1.0e-6_dp, dd) + 1.0e-30_dp
    r = 1.0_dp / dd
  end function icpt2_safe_inv

  !> @brief des Cloizeaux symmetric multistate effective Hamiltonian + spectrum.
  !> @param[in]  nP        number of dressed CAS roots
  !> @param[in]  nQ        number of external-Q perturbers
  !> @param[in]  eP(nP)    CAS roots (the lowest nP)
  !> @param[in]  coup(nQ,nP)  contracted couplings <q|H|Psi_P^k>
  !> @param[in]  invd(nQ,nP)  1/(eP_k - H0_qk)  (EN or Dyall; intruder-guarded)
  !> @param[out] Edressed(nP) dressed energies (ascending)
  !> @param[out] herm      max |H_eff - H_eff^T| before the defensive symmetrize
  subroutine icpt2_eff_hamiltonian(nP, nQ, eP, coup, invd, Edressed, herm, evec)
    use eigen, only: diag_symm_full
    integer,  intent(in)  :: nP, nQ
    real(dp), intent(in)  :: eP(nP), coup(nQ,nP), invd(nQ,nP)
    real(dp), intent(out) :: Edressed(nP)
    real(dp), intent(out) :: herm
    real(dp), intent(out), optional :: evec(nP,nP)   !< H_eff eigenvectors (CAS-root basis)
    real(dp) :: Heff(nP,nP), s
    integer  :: k, l, q, ierr
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
    call diag_symm_full(0, nP, Heff, nP, Edressed, ierr)   ! Heff overwritten by eigenvectors
    if (present(evec)) evec = Heff
  end subroutine icpt2_eff_hamiltonian

end module qmrsf_icpt2_downfold_mod
