!> Genuine three-body term of the Ten-no / Boys-Handy transcorrelated Hamiltonian.
!>
!> The (1/2)(grad_i tau)^2 piece of H_bar, expanded with grad_i tau = sum_{j!=i}
!> grad_i u(r_ij), has a j=k diagonal (two-body, already in tc_boyshandy) and a
!> j!=k cross term that is a genuine THREE-electron operator:
!>
!>   O3 = -(1/2) sum_i sum_{j!=i} sum_{k!=i,k!=j} (grad_i u_ij).(grad_i u_ik).
!>
!> O3 is multiplicative (no derivative acts on the wavefunction) and Hermitian.
!> Electron i is the apex (differentiated against both partners j,k). For the
!> Gaussian geminal u = sum_m C_m exp(-g_m r^2),
!>   (grad_i u_ij).(grad_i u_ik) = sum_{m,n} C_m C_n 4 g_m g_n (r_i-r_j).(r_i-r_k)
!>        exp(-g_m r_ij^2 - g_n r_ik^2).
!>
!> This module provides o3_geminal_s / o3_prim_s: the exact s-primitive 3-electron
!> integral <a b c| (grad_1 u_12).(grad_1 u_13) |a' b' c'>, obtained by integrating
!> the two leg electrons (2,3) analytically against their geminal-gradient (each
!> leaves a Gaussian x linear-in-r1 factor at the apex), then doing the apex
!> (electron-1) Gaussian-polynomial integral in closed form. Validation oracle for
!> the normal-ordered effective operators. (Names: exponent e, centre Q, const K --
!> kept letter-distinct to avoid Fortran's case-insensitive p/P clash.)
module tc_threebody
  use precision, only: dp
  implicit none
  private
  public :: o3_prim_s, o3_geminal_s

  real(dp), parameter :: PI = 3.141592653589793238462643383279_dp

contains

  !> Charge distribution of two s primitives a(r)=e^{-za|r-A|^2}, b(r)=e^{-zb|r-B|^2}:
  !> a*b = K exp(-e |r-Q|^2). Returns exponent e, centre Q, prefactor K.
  pure subroutine sprod(za, A, zb, B, e, Q, K)
    real(dp), intent(in)  :: za, zb, A(3), B(3)
    real(dp), intent(out) :: e, Q(3), K
    e = za + zb
    Q = (za*A + zb*B)/e
    K = exp(-za*zb/e * dot_product(A-B, A-B))
  end subroutine sprod

  !> Genuine 3-electron integral for ONE Gaussian-geminal pair (gm on the 1-2 leg,
  !> gn on the 1-3 leg), s primitives only:
  !>   I(m,n) = 4 gm gn * <a(1)a'(1) b(2)b'(2) c(3)c'(3)| (r1-r2).(r1-r3)
  !>            e^{-gm r12^2 - gn r13^2} >.
  !> electron 1 = (a,a'), electron 2 = (b,b'), electron 3 = (c,c').
  function o3_geminal_s(zA,A, zAp,Ap, zB,B, zBp,Bp, zC,C, zCp,Cp, gm, gn) result(val)
    real(dp), intent(in) :: zA,A(3), zAp,Ap(3), zB,B(3), zBp,Bp(3), zC,C(3), zCp,Cp(3), gm, gn
    real(dp) :: val
    real(dp) :: e1,Q1(3),K1, e2,Q2(3),K2, e3,Q3(3),K3
    real(dp) :: mu2, mu3, pref2, pref3, etot, Rt(3), gconst, poly
    integer  :: d
    call sprod(zA,A, zAp,Ap, e1,Q1,K1)   ! apex (electron 1)
    call sprod(zB,B, zBp,Bp, e2,Q2,K2)   ! leg 2
    call sprod(zC,C, zCp,Cp, e3,Q3,K3)   ! leg 3
    ! integrate leg 2: int Omega2(r2)(x1-x2)e^{-gm|r1-r2|^2} dr2
    !   = K2 (pi/(e2+gm))^{3/2} [e2/(e2+gm)] (x1-Q2) e^{-mu2 |r1-Q2|^2}
    mu2 = e2*gm/(e2+gm); pref2 = K2*(PI/(e2+gm))**1.5_dp * (e2/(e2+gm))
    mu3 = e3*gn/(e3+gn); pref3 = K3*(PI/(e3+gn))**1.5_dp * (e3/(e3+gn))
    ! apex integral: int Omega1(r1) [pref2 (x1-Q2)] [pref3 (x1-Q3)] summed over dir
    etot = e1 + mu2 + mu3
    Rt   = (e1*Q1 + mu2*Q2 + mu3*Q3)/etot
    gconst = exp( -( e1*mu2*dot_product(Q1-Q2,Q1-Q2) &
                   + e1*mu3*dot_product(Q1-Q3,Q1-Q3) &
                   + mu2*mu3*dot_product(Q2-Q3,Q2-Q3) )/etot )
    ! int e^{-etot|r1-Rt|^2} (x1-Q2)_d (x1-Q3)_d dr1
    !   = (pi/etot)^{3/2} [ (Rt_d-Q2_d)(Rt_d-Q3_d) + 1/(2 etot) ]   (per dir)
    poly = 0.0_dp
    do d = 1, 3
      poly = poly + ( (Rt(d)-Q2(d))*(Rt(d)-Q3(d)) + 0.5_dp/etot )
    end do
    val = 4.0_dp*gm*gn * K1*pref2*pref3 * gconst * (PI/etot)**1.5_dp * poly
  end function o3_geminal_s

  !> Full s-primitive 3-electron O3 integral summed over the Gaussian expansion of u
  !> (coefficients C_m, exponents g_m):
  !>   <a b c| (grad_1 u_12).(grad_1 u_13) |a' b' c'> = sum_{m,n} C_m C_n o3_geminal_s.
  !> The leading -(1/2) of O3 and the spin amplitudes are applied by the caller.
  function o3_prim_s(zA,A, zAp,Ap, zB,B, zBp,Bp, zC,C, zCp,Cp, Cf, gf, ng) result(val)
    real(dp), intent(in) :: zA,A(3), zAp,Ap(3), zB,B(3), zBp,Bp(3), zC,C(3), zCp,Cp(3)
    real(dp), intent(in) :: Cf(:), gf(:)
    integer,  intent(in) :: ng
    real(dp) :: val
    integer  :: m, n
    val = 0.0_dp
    do m = 1, ng
      do n = 1, ng
        val = val + Cf(m)*Cf(n)* &
          o3_geminal_s(zA,A, zAp,Ap, zB,B, zBp,Bp, zC,C, zCp,Cp, gf(m), gf(n))
      end do
    end do
  end function o3_prim_s

end module tc_threebody
