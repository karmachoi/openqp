program test_1e_der2_fd
  use precision, only: dp
  use basis_tools, only: basis_set
  use mod_shell_tools, only: shell_t, shpair_t
  use mod_1e_primitives, only: comp_overlap_der1, comp_overlap_der2, &
      comp_kinetic_der1, comp_kinetic_der2, comp_coulomb_der1, comp_coulomb_der2
  implicit none

  call run_fixture(1, 1, [0.73_dp])
  call run_fixture(2, 1, [0.31_dp, -0.27_dp, 0.19_dp])

contains

  subroutine run_fixture(iang, jang, dij_values)
    integer, intent(in) :: iang, jang
    real(dp), intent(in) :: dij_values(:)
    integer :: inao, jnao, n
    real(dp), allocatable :: dij(:,:)
    real(dp) :: ri(3), rj(3), c(3), znuc

    ri = [0.20_dp, -0.30_dp, 0.40_dp]
    rj = [-0.45_dp, 0.25_dp, -0.15_dp]
    c = [0.31_dp, -0.22_dp, 0.73_dp]
    znuc = 1.7_dp

    inao = nao_for_ang(iang)
    jnao = nao_for_ang(jang)
    if (size(dij_values) /= inao*jnao) error stop 'Bad dij fixture size'

    allocate(dij(inao,jnao))
    n = 0
    do jnao = 1, size(dij, 2)
      do inao = 1, size(dij, 1)
        n = n + 1
        dij(inao,jnao) = dij_values(n)
      end do
    end do

    call check_operator('overlap', iang, jang, ri, rj, c, znuc, dij, 1.0e-7_dp, 1.0e-7_dp)
    call check_operator('kinetic', iang, jang, ri, rj, c, znuc, dij, 1.0e-7_dp, 1.0e-7_dp)
    call check_operator('coulomb', iang, jang, ri, rj, c, znuc, dij, 1.0e-6_dp, 1.0e-6_dp)
  end subroutine

  subroutine check_operator(op, iang, jang, ri, rj, c, znuc, dij, atol, rtol)
    character(len=*), intent(in) :: op
    integer, intent(in) :: iang, jang
    real(dp), intent(in) :: ri(3), rj(3), c(3), znuc, dij(:,:), atol, rtol
    type(shpair_t) :: cp
    real(dp) :: analytic(3,3), fd(3,3), diff(3,3), scale(3,3)
    real(dp) :: err

    call make_pair(cp, iang, jang, ri, rj)
    analytic = 0.0_dp
    select case (op)
    case ('overlap')
      call comp_overlap_der2(cp, dij, analytic)
      call finite_difference_der1(op, iang, jang, ri, rj, c, znuc, dij, fd)
    case ('kinetic')
      call comp_kinetic_der2(cp, dij, analytic)
      call finite_difference_der1(op, iang, jang, ri, rj, c, znuc, dij, fd)
    case ('coulomb')
      call comp_coulomb_der2(cp, c, znuc, dij, analytic)
      call finite_difference_der1(op, iang, jang, ri, rj, c, znuc, dij, fd)
    case default
      error stop 'Unknown operator'
    end select

    diff = abs(analytic - fd)
    scale = atol + rtol * max(abs(analytic), abs(fd))
    err = maxval(diff / scale)
    if (err > 1.0_dp) then
      write(*,'(a)') 'Finite-difference derivative mismatch: '//trim(op)
      write(*,'(a,es24.16)') 'scaled max error = ', err
      write(*,'(a)') 'analytic:'
      call print_matrix(analytic)
      write(*,'(a)') 'finite difference:'
      call print_matrix(fd)
      write(*,'(a)') 'absolute diff:'
      call print_matrix(diff)
      error stop 'one-electron der2 finite-difference check failed'
    end if
  end subroutine

  subroutine finite_difference_der1(op, iang, jang, ri, rj, c, znuc, dij, fd)
    character(len=*), intent(in) :: op
    integer, intent(in) :: iang, jang
    real(dp), intent(in) :: ri(3), rj(3), c(3), znuc, dij(:,:)
    real(dp), intent(out) :: fd(3,3)
    real(dp), parameter :: h = 1.0e-6_dp
    real(dp) :: ri_plus(3), ri_minus(3), g_plus(3), g_minus(3)
    integer :: beta

    do beta = 1, 3
      ri_plus = ri
      ri_minus = ri
      ri_plus(beta) = ri_plus(beta) + h
      ri_minus(beta) = ri_minus(beta) - h
      call der1_at(op, iang, jang, ri_plus, rj, c, znuc, dij, g_plus)
      call der1_at(op, iang, jang, ri_minus, rj, c, znuc, dij, g_minus)
      fd(:, beta) = (g_plus - g_minus) / (2.0_dp * h)
    end do
  end subroutine

  subroutine der1_at(op, iang, jang, ri, rj, c, znuc, dij, grad)
    character(len=*), intent(in) :: op
    integer, intent(in) :: iang, jang
    real(dp), intent(in) :: ri(3), rj(3), c(3), znuc, dij(:,:)
    real(dp), intent(out) :: grad(3)
    type(shpair_t) :: cp

    call make_pair(cp, iang, jang, ri, rj)
    grad = 0.0_dp
    select case (op)
    case ('overlap')
      call comp_overlap_der1(cp, dij, grad)
    case ('kinetic')
      call comp_kinetic_der1(cp, dij, grad)
    case ('coulomb')
      call comp_coulomb_der1(cp, c, znuc, dij, grad)
    case default
      error stop 'Unknown operator'
    end select
  end subroutine

  subroutine make_pair(cp, iang, jang, ri, rj)
    type(shpair_t), intent(inout) :: cp
    integer, intent(in) :: iang, jang
    real(dp), intent(in) :: ri(3), rj(3)
    type(basis_set) :: basis
    type(shell_t) :: shi, shj

    basis%mxcontr = 1
    allocate(basis%ex(2), basis%cc(2))
    basis%ex = [0.7_dp, 1.2_dp]
    basis%cc = [1.0_dp, -0.8_dp]

    shi%shid = 1
    shi%atid = 1
    shi%ig1 = 1
    shi%ig2 = 1
    shi%ang = iang
    shi%locao = 1
    shi%nao = nao_for_ang(iang)
    shi%r = ri

    shj%shid = 2
    shj%atid = 2
    shj%ig1 = 2
    shj%ig2 = 2
    shj%ang = jang
    shj%locao = 1 + shi%nao
    shj%nao = nao_for_ang(jang)
    shj%r = rj

    call cp%alloc(basis)
    call cp%shell_pair(basis, shi, shj, huge(1.0_dp), dup=.false.)
    if (cp%numpairs /= 1) error stop 'Unexpected screened shell pair'
  end subroutine

  integer function nao_for_ang(ang) result(nao)
    integer, intent(in) :: ang
    nao = ang * (ang + 1) / 2
  end function

  subroutine print_matrix(mat)
    real(dp), intent(in) :: mat(3,3)
    integer :: i
    do i = 1, 3
      write(*,'(3es24.16)') mat(i,1), mat(i,2), mat(i,3)
    end do
  end subroutine

end program test_1e_der2_fd
