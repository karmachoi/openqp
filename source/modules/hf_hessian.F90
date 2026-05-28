module hf_hessian_mod

  use precision, only: dp

  implicit none

  character(len=*), parameter :: module_name = "hf_hessian_mod"

contains

!###############################################################################

  subroutine hf_hessian_C(c_handle) bind(C, name="hf_hessian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call hf_hessian(inf)
  end subroutine hf_hessian_C

!###############################################################################

  subroutine hf_hessian(infos)
    use types, only: information
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos

    real(kind=dp), allocatable :: hessian(:,:)
    integer :: natom, ndim, ok

    natom = int(infos%mol_prop%natom)
    ndim = 3 * natom
    allocate(hessian(ndim, ndim), source=0.0_dp, stat=ok)
    if (ok /= 0) call show_message('Cannot allocate HF/DFT Hessian work matrix.', WITH_ABORT)

    call build_nuclear_repulsion_hessian(infos, hessian)

    ! native_openqp_hf_nuclear_repulsion_only partial_kernel
    ! The exact nuclear-repulsion block is now built natively, but the full
    ! HF/DFT analytic Hessian still requires one-/two-electron derivative and
    ! CPHF/CPKS response terms.  Keep runtime guarded until those terms are
    ! implemented and finite-difference validated.
    call show_message(&
      'Native HF/DFT Hessian partial_kernel reached: nuclear repulsion only; implementation is not complete yet.', &
      WITH_ABORT)
  end subroutine hf_hessian

!###############################################################################

  subroutine build_nuclear_repulsion_hessian(infos, hessian)
    use types, only: information

    implicit none

    type(information), target, intent(in) :: infos
    real(kind=dp), intent(inout) :: hessian(:,:)

    integer :: natom
    integer :: iatom, jatom, idir, jdir
    integer :: ia, ib, ja, jb
    real(kind=dp) :: qi, qj, qij
    real(kind=dp) :: rij(3), r2, r, r3, r5
    real(kind=dp) :: block(3,3)

    natom = int(infos%mol_prop%natom)
    hessian = 0.0_dp

    do iatom = 1, natom - 1
      qi = infos%atoms%zn(iatom) - infos%basis%ecp_zn_num(iatom)
      do jatom = iatom + 1, natom
        qj = infos%atoms%zn(jatom) - infos%basis%ecp_zn_num(jatom)
        qij = qi * qj

        rij = infos%atoms%xyz(:,iatom) - infos%atoms%xyz(:,jatom)
        r2 = dot_product(rij, rij)
        r = sqrt(r2)
        r3 = r2 * r
        r5 = r3 * r2

        block = 0.0_dp
        do idir = 1, 3
          do jdir = 1, 3
            block(idir,jdir) = 3.0_dp * rij(idir) * rij(jdir) / r5
          end do
          block(idir,idir) = block(idir,idir) - 1.0_dp / r3
        end do
        block = qij * block

        do idir = 1, 3
          ia = 3 * (iatom - 1) + idir
          ib = 3 * (jatom - 1) + idir
          do jdir = 1, 3
            ja = 3 * (iatom - 1) + jdir
            jb = 3 * (jatom - 1) + jdir

            hessian(ia,ja) = hessian(ia,ja) + block(idir,jdir)
            hessian(ib,jb) = hessian(ib,jb) + block(idir,jdir)
            hessian(ia,jb) = hessian(ia,jb) - block(idir,jdir)
            hessian(ib,ja) = hessian(ib,ja) - block(idir,jdir)
          end do
        end do
      end do
    end do
  end subroutine build_nuclear_repulsion_hessian

end module hf_hessian_mod
