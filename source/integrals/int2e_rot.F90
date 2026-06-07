!> @brief  M5: "rot" two-electron backend -- libintRot spherical-resolution
!>         engine (tabulated radial values + spherical-harmonic rotations) wired
!>         into OpenQP as an alternative to rotspd/libint/rys.
!>
!> Opt-in via environment variables (default OpenQP behaviour is unchanged):
!>   OQP_ERI_ROT=1            enable the rot backend
!>   OQP_ROT_RADIAL=<path>    radial_table.bin
!>   OQP_ROT_GAUNT=<path>     gaunt_table.bin
!>
!> The C side (src/openqp_iface.cpp in libintRot) returns a raw-monomial
!> Cartesian shell-quartet block in OpenQP's component order and ints(nd,nc,nb,na)
!> layout; the caller then applies normalize_ints, exactly as for rys/libint.
module int2e_rot

  use iso_c_binding, only: c_int, c_double, c_char, c_null_char
  use precision, only: dp
  use basis_tools, only: basis_set
  implicit none

  private
  public :: rot_active, rot_static_init, rot_compute_eri

  logical, save :: rot_active = .false.
  logical, save :: rot_inited = .false.

  interface
    function librot_init(radial, gaunt) bind(C, name="librot_init") result(ierr)
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: radial(*)
      character(kind=c_char), intent(in) :: gaunt(*)
      integer(c_int) :: ierr
    end function

    function librot_ready() bind(C, name="librot_ready") result(r)
      import :: c_int
      integer(c_int) :: r
    end function

    function librot_eri_cart(la,lb,lc,ld, A,B,C,D, &
                             ea,ca,Ka, eb,cb,Kb, ec,cc,Kc, ed,cd,Kd, out) &
        bind(C, name="librot_eri_cart") result(ierr)
      import :: c_int, c_double
      integer(c_int), value :: la,lb,lc,ld, Ka,Kb,Kc,Kd
      real(c_double), intent(in)  :: A(*),B(*),C(*),D(*)
      real(c_double), intent(in)  :: ea(*),ca(*),eb(*),cb(*),ec(*),cc(*),ed(*),cd(*)
      real(c_double), intent(out) :: out(*)
      integer(c_int) :: ierr
    end function
  end interface

contains

  !> Read environment, load tables once.  Sets rot_active.
  subroutine rot_static_init()
    character(len=512) :: sval, rad, gau
    integer :: ln, ierr
    if (rot_inited) return
    rot_inited = .true.

    call get_environment_variable("OQP_ERI_ROT", sval, ln)
    if (ln <= 0) return
    if (.not.(sval(1:1) == '1' .or. sval(1:1) == 'y' .or. sval(1:1) == 'Y' &
              .or. sval(1:1) == 't' .or. sval(1:1) == 'T')) return

    call get_environment_variable("OQP_ROT_RADIAL", rad, ln)
    if (ln <= 0) rad = "radial_table.bin"
    call get_environment_variable("OQP_ROT_GAUNT", gau, ln)
    if (ln <= 0) gau = "gaunt_table.bin"

    ierr = librot_init(trim(rad)//c_null_char, trim(gau)//c_null_char)
    if (ierr == 0 .and. librot_ready() == 1) then
      rot_active = .true.
    end if
  end subroutine

  !> Compute the Cartesian shell-quartet block for shells ids(1:4) into ints,
  !> in OpenQP's ints(nd,nc,nb,na) layout (un-normalized; caller normalizes).
  subroutine rot_compute_eri(basis, ids, ints, nbf, ok)
    type(basis_set), intent(in) :: basis
    integer, intent(in)  :: ids(4)
    real(dp), intent(out) :: ints(*)
    integer, intent(out) :: nbf(4)
    logical, intent(out) :: ok

    integer :: am(4), s, k(4), nc(4)
    integer(c_int) :: cam(4), cnc(4), ierr
    real(c_double) :: A(3), B(3), C(3), D(3)

    am = basis%am(ids)
    nbf = (am+1)*(am+2)/2
    do s = 1, 4
      k(s)  = basis%g_offset(ids(s))
      nc(s) = basis%ncontr(ids(s))
    end do
    cam = int(am, c_int)         ! OpenQP default integer is 8-byte; C entry uses c_int
    cnc = int(nc, c_int)
    A = basis%shell_centers(ids(1),1:3)
    B = basis%shell_centers(ids(2),1:3)
    C = basis%shell_centers(ids(3),1:3)
    D = basis%shell_centers(ids(4),1:3)

    ierr = librot_eri_cart(cam(1),cam(2),cam(3),cam(4), A,B,C,D, &
       basis%ex(k(1)), basis%cc(k(1)), cnc(1), &
       basis%ex(k(2)), basis%cc(k(2)), cnc(2), &
       basis%ex(k(3)), basis%cc(k(3)), cnc(3), &
       basis%ex(k(4)), basis%cc(k(4)), cnc(4), ints)
    ok = (ierr == 0)
  end subroutine

end module int2e_rot
