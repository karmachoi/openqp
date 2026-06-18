!> @brief  pTC-MRSF-CIS Phase 2: genuine F12 Gaussian-geminal two-electron
!>         integrals and the F12 intermediates (X, B, V) built from them.
!>
!> @details
!>   Ten-no's projective transcorrelation uses the Slater correlation factor
!>   f12 = exp(-gamma r12), represented as a fixed linear combination of Gaussian
!>   geminals exp(-omega r12^2) (Tew-Klopper STG-NG fit). Every F12 matrix element
!>   the transcorrelated Hamiltonian needs reduces to the closed-form four-centre
!>   Gaussian-geminal integral over s primitives and its omega-derivative:
!>
!>     geminal  (ab| e^{-omega r12^2} |cd)
!>                = Kac Kbd (pi^2/Dl)^{3/2} exp(-(pq omega/Dl) |P-Q|^2),
!>                  Dl = pq + (p+q) omega,   p=aA+aC (->P,Kac), q=aB+aD (->Q,Kbd)
!>     X (metric)   <ab| f12^2     |cd>  = geminal(2 omega)
!>     B (kinetic)  <ab| (grad f12)^2|cd> = 8 omega^2 <ab| r12^2 e^{-2omega r12^2}|cd>
!>                  with <ab| r12^2 e^{-lam r12^2}|cd> = -d/dlam geminal(lam)
!>     V (Coulomb)  <ab| f12 / r12 |cd>  = (2/sqrt(pi)) int_0^inf geminal(omega+t^2) dt
!>                  (-> the standard 1/r12 ERI, Boys F0, as omega -> 0)
!>
!>   This module is the native (pyscf-free) port of the validated prototypes
!>   tests/ptc_mrsf/prototype/{r12_geminal,f12_intermediates}.py. It is the
!>   primitive layer that tc_build_eff_integrals (tdhf_mrsf_ptc.F90, Phase 2
!>   assembly) will contract over shell pairs and the STG-NG expansion.
!>
!>   Self-test: tests/ptc_mrsf/prototype/ptc_geminal_test.F90 validates every
!>   routine against a pyscf-free oracle (overlap-product omega->0 limit,
!>   6D numerical quadrature, Boys-F0 ERI).
!
module ptc_geminal

  use precision, only: dp
  implicit none
  private

  public :: gem_overlap_s      ! one-electron s overlap (a|c)
  public :: gaussian_geminal_s ! (ab| e^{-omega r12^2} |cd)
  public :: r2_geminal_s       ! <ab| r12^2 e^{-lam r12^2} |cd> = -d/dlam geminal
  public :: f12_X_s            ! X metric   = geminal(2 omega)
  public :: f12_B_s            ! B kinetic  = 8 omega^2 r2_geminal(2 omega)
  public :: eri_s              ! (ab| 1/r12 |cd) closed form (Boys F0)
  public :: v_geminal_s        ! V Coulomb  = <ab| e^{-omega r12^2}/r12 |cd>
  public :: boys0
  public :: stg_ng             ! STG-NG fit of exp(-gamma r12) -> {c_k, omega_k}
  public :: s_norm             ! primitive s normalization (2a/pi)^{3/4}
  public :: ptc_s_ao_tensor    ! contracted-shell AO 4-index tensor for an operator

  ! Operator selectors for ptc_s_ao_tensor
  integer, parameter, public :: PTC_OP_GEMINAL = 1  ! (ab| e^{-w r12^2} |cd)
  integer, parameter, public :: PTC_OP_X       = 2  ! geminal(2w)
  integer, parameter, public :: PTC_OP_B       = 3  ! 8 w^2 r2_geminal(2w)
  integer, parameter, public :: PTC_OP_V       = 4  ! <ab| e^{-w r12^2}/r12 |cd>
  integer, parameter, public :: PTC_OP_ERI     = 5  ! (ab| 1/r12 |cd)

  real(dp), parameter :: PI = 3.141592653589793238462643383279_dp

  ! Tew-Klopper 6-term STG-6G fit of exp(-r12): coefficients c_k and Gaussian
  ! exponents e_k, so exp(-gamma r12) ~ sum_k c_k exp(-(e_k gamma^2) r12^2).
  integer, parameter :: NSTG = 6
  real(dp), parameter :: STG6_C(NSTG) = &
    [0.3144_dp, 0.3037_dp, 0.1681_dp, 0.09811_dp, 0.06024_dp, 0.03726_dp]
  real(dp), parameter :: STG6_E(NSTG) = &
    [0.2209_dp, 1.004_dp, 3.622_dp, 12.16_dp, 45.87_dp, 254.4_dp]

contains

  !> Gaussian product: exponent p, centre Pc, prefactor K for a*c (un-normalized).
  pure subroutine gem_product(aA, A, aC, C, p, Pc, K)
    real(dp), intent(in)  :: aA, aC, A(3), C(3)
    real(dp), intent(out) :: p, Pc(3), K
    p  = aA + aC
    Pc = (aA*A + aC*C)/p
    K  = exp(-aA*aC/p * dot_product(A - C, A - C))
  end subroutine gem_product

  !> One-electron overlap of two s primitives (a|c), un-normalized.
  pure function gem_overlap_s(aA, A, aC, C) result(s)
    real(dp), intent(in) :: aA, aC, A(3), C(3)
    real(dp) :: s, p
    p = aA + aC
    s = (PI/p)**1.5_dp * exp(-aA*aC/p * dot_product(A - C, A - C))
  end function gem_overlap_s

  !> Four-centre Gaussian-geminal integral (ab| exp(-omega r12^2) |cd) over s
  !> primitives. omega -> 0 reduces to the overlap product (a|c)(b|d).
  pure function gaussian_geminal_s(aA,A,aB,B,aC,C,aD,D,omega) result(g)
    real(dp), intent(in) :: aA,aB,aC,aD, A(3),B(3),C(3),D(3), omega
    real(dp) :: g, p, q, Pc(3), Qc(3), Kac, Kbd, denom, pref, ex
    call gem_product(aA,A,aC,C, p,Pc,Kac)
    call gem_product(aB,B,aD,D, q,Qc,Kbd)
    denom = p*q + p*omega + q*omega
    pref  = (PI**2/denom)**1.5_dp
    ex    = exp(-(p*q*omega/denom) * dot_product(Pc - Qc, Pc - Qc))
    g     = Kac*Kbd*pref*ex
  end function gaussian_geminal_s

  !> <ab| r12^2 exp(-lam r12^2) |cd> = -d/dlam geminal(lam)  (the B kernel core).
  pure function r2_geminal_s(aA,A,aB,B,aC,C,aD,D,lam) result(r)
    real(dp), intent(in) :: aA,aB,aC,aD, A(3),B(3),C(3),D(3), lam
    real(dp) :: r, p, q, Pc(3), Qc(3), Kac, Kbd, Dl, pq, PQ2, gem
    call gem_product(aA,A,aC,C, p,Pc,Kac)
    call gem_product(aB,B,aD,D, q,Qc,Kbd)
    Dl  = p*q + (p+q)*lam
    pq  = p*q
    PQ2 = dot_product(Pc - Qc, Pc - Qc)
    gem = Kac*Kbd*(PI**2/Dl)**1.5_dp * exp(-(pq*lam/Dl)*PQ2)
    r   = gem * (1.5_dp*(p+q)/Dl + PQ2*(pq/Dl)**2)
  end function r2_geminal_s

  !> F12 metric intermediate X = <ab| f12^2 |cd> = geminal(2 omega).
  pure function f12_X_s(aA,A,aB,B,aC,C,aD,D,omega) result(x)
    real(dp), intent(in) :: aA,aB,aC,aD, A(3),B(3),C(3),D(3), omega
    real(dp) :: x
    x = gaussian_geminal_s(aA,A,aB,B,aC,C,aD,D, 2.0_dp*omega)
  end function f12_X_s

  !> F12 kinetic intermediate B = <ab| (grad f12)^2 |cd>
  !>   = 8 omega^2 <ab| r12^2 e^{-2 omega r12^2} |cd>.
  pure function f12_B_s(aA,A,aB,B,aC,C,aD,D,omega) result(bval)
    real(dp), intent(in) :: aA,aB,aC,aD, A(3),B(3),C(3),D(3), omega
    real(dp) :: bval
    bval = 8.0_dp*omega*omega * r2_geminal_s(aA,A,aB,B,aC,C,aD,D, 2.0_dp*omega)
  end function f12_B_s

  !> Boys function F0(t) = (1/2) sqrt(pi/t) erf(sqrt(t)), F0(0) = 1.
  pure function boys0(t) result(f)
    real(dp), intent(in) :: t
    real(dp) :: f
    if (t < 1.0e-12_dp) then
      f = 1.0_dp
    else
      f = 0.5_dp*sqrt(PI/t)*erf(sqrt(t))
    end if
  end function boys0

  !> Closed-form electron-repulsion integral (ab| 1/r12 |cd) for s primitives.
  pure function eri_s(aA,A,aB,B,aC,C,aD,D) result(v)
    real(dp), intent(in) :: aA,aB,aC,aD, A(3),B(3),C(3),D(3)
    real(dp) :: v, p, q, Pc(3), Qc(3), Kac, Kbd, alpha, PQ2
    call gem_product(aA,A,aC,C, p,Pc,Kac)
    call gem_product(aB,B,aD,D, q,Qc,Kbd)
    alpha = p*q/(p+q)
    PQ2   = dot_product(Pc - Qc, Pc - Qc)
    v = (2.0_dp*PI**2.5_dp/(p*q*sqrt(p+q))) * Kac*Kbd * boys0(alpha*PQ2)
  end function eri_s

  !> V Coulomb intermediate <ab| exp(-omega r12^2)/r12 |cd>
  !>   = (2/sqrt(pi)) int_0^inf geminal(omega + t^2) dt,  mapped Gauss-Legendre
  !> on t = x/(1-x), x in (0,1). omega = 0 reproduces eri_s to ~1e-10.
  function v_geminal_s(aA,A,aB,B,aC,C,aD,D,omega, nq) result(v)
    real(dp), intent(in) :: aA,aB,aC,aD, A(3),B(3),C(3),D(3), omega
    integer, intent(in), optional :: nq
    real(dp) :: v
    integer :: n, i
    real(dp), allocatable :: xg(:), wg(:)
    real(dp) :: x, t, jac, acc
    n = 128
    if (present(nq)) n = nq
    allocate(xg(n), wg(n))
    call gauss_legendre01(n, xg, wg)
    acc = 0.0_dp
    do i = 1, n
      x   = xg(i)
      t   = x/(1.0_dp - x)
      jac = 1.0_dp/(1.0_dp - x)**2
      acc = acc + wg(i)*jac*gaussian_geminal_s(aA,A,aB,B,aC,C,aD,D, omega + t*t)
    end do
    v = (2.0_dp/sqrt(PI)) * acc
    deallocate(xg, wg)
  end function v_geminal_s

  !> Gauss-Legendre nodes/weights on [0,1] (Newton iteration on P_n).
  subroutine gauss_legendre01(n, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: x(n), w(n)
    integer  :: i, j, it
    real(dp) :: z, z1, p1, p2, p3, pp
    real(dp), parameter :: EPS = 1.0e-15_dp
    do i = 1, (n+1)/2
      z = cos(PI*(real(i,dp) - 0.25_dp)/(real(n,dp) + 0.5_dp))
      do it = 1, 100
        p1 = 1.0_dp
        p2 = 0.0_dp
        do j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp*j - 1.0_dp)*z*p2 - (j - 1.0_dp)*p3)/real(j,dp)
        end do
        pp = real(n,dp)*(z*p1 - p2)/(z*z - 1.0_dp)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= EPS) exit
      end do
      ! [-1,1] node z, weight 2/((1-z^2) pp^2); map to [0,1] (Jacobian 1/2).
      x(i)       = 0.5_dp*(1.0_dp - z)
      x(n+1-i)   = 0.5_dp*(1.0_dp + z)
      w(i)       = 1.0_dp/((1.0_dp - z*z)*pp*pp)
      w(n+1-i)   = w(i)
    end do
  end subroutine gauss_legendre01

  !> STG-NG fit of the Slater factor exp(-gamma r12): on return c(k),omg(k) give
  !> exp(-gamma r12) ~ sum_k c(k) exp(-omg(k) r12^2), k=1..NSTG.
  pure subroutine stg_ng(gamma, c, omg, ng)
    real(dp), intent(in)  :: gamma
    real(dp), intent(out) :: c(:), omg(:)
    integer,  intent(out) :: ng
    integer :: k
    ng = NSTG
    do k = 1, NSTG
      c(k)   = STG6_C(k)
      omg(k) = STG6_E(k)*gamma*gamma
    end do
  end subroutine stg_ng

  !> Primitive s-Gaussian normalization N = (2a/pi)^{3/4}.
  pure function s_norm(a) result(n)
    real(dp), intent(in) :: a
    real(dp) :: n
    n = (2.0_dp*a/PI)**0.75_dp
  end function s_norm

  !> Contracted-shell AO 4-index tensor M(i,j,k,l) = <ij|op|kl> (physicist
  !> convention: electron 1 = shells i,k; electron 2 = shells j,l) for one of the
  !> geminal operators, over nsh contracted s shells. Shell s has np(s) primitives
  !> with exponents ex(:,s), contraction coefficients co(:,s), centre cn(:,s).
  !> The primitive s normalization (2a/pi)^{3/4} is applied internally, so co are
  !> the bare contraction coefficients. This is the contraction machinery the
  !> Phase-2 assembly (tc_build_eff_integrals) uses; generalizing op past s shells
  !> or routing 1/r12 through the Rys engine is the remaining production work.
  subroutine ptc_s_ao_tensor(nsh, np, ex, co, cn, op, omega, M)
    integer,  intent(in)  :: nsh
    integer,  intent(in)  :: np(nsh)
    real(dp), intent(in)  :: ex(:,:), co(:,:), cn(3,nsh)
    integer,  intent(in)  :: op
    real(dp), intent(in)  :: omega
    real(dp), intent(out) :: M(nsh,nsh,nsh,nsh)
    integer  :: i, j, k, l, pa, pb, pc, pd
    real(dp) :: acc, ci, cj, ck, cl, ea, eb, ec, ed, val
    real(dp) :: Ai(3), Aj(3), Ak(3), Al(3)
    do i = 1, nsh
      Ai = cn(:,i)
      do j = 1, nsh
        Aj = cn(:,j)
        do k = 1, nsh
          Ak = cn(:,k)
          do l = 1, nsh
            Al = cn(:,l)
            acc = 0.0_dp
            do pa = 1, np(i)
              ea = ex(pa,i); ci = co(pa,i)*s_norm(ea)
              do pb = 1, np(j)
                eb = ex(pb,j); cj = co(pb,j)*s_norm(eb)
                do pc = 1, np(k)
                  ec = ex(pc,k); ck = co(pc,k)*s_norm(ec)
                  do pd = 1, np(l)
                    ed = ex(pd,l); cl = co(pd,l)*s_norm(ed)
                    select case (op)
                    case (PTC_OP_GEMINAL)
                      val = gaussian_geminal_s(ea,Ai,eb,Aj,ec,Ak,ed,Al, omega)
                    case (PTC_OP_X)
                      val = f12_X_s(ea,Ai,eb,Aj,ec,Ak,ed,Al, omega)
                    case (PTC_OP_B)
                      val = f12_B_s(ea,Ai,eb,Aj,ec,Ak,ed,Al, omega)
                    case (PTC_OP_V)
                      val = v_geminal_s(ea,Ai,eb,Aj,ec,Ak,ed,Al, omega)
                    case (PTC_OP_ERI)
                      val = eri_s(ea,Ai,eb,Aj,ec,Ak,ed,Al)
                    case default
                      val = 0.0_dp
                    end select
                    acc = acc + ci*cj*ck*cl*val
                  end do
                end do
              end do
            end do
            M(i,j,k,l) = acc
          end do
        end do
      end do
    end do
  end subroutine ptc_s_ao_tensor

end module ptc_geminal
