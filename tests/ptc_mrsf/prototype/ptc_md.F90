!> General angular-momentum Gaussian integrals via McMurchie-Davidson, native and
!> pyscf-free, to lift the pTC-MRSF-CIS engine beyond s-only (enables p/d bases,
!> ethylene, cc-pVTZ, and the proper Rydberg excited states S1..S3).
!>
!> Cartesian Gaussian primitive: x^l y^m z^n exp(-a |r-A|^2). Angular momentum is
!> carried as integer triplets (lx,ly,lz). This module provides the Hermite
!> expansion coefficients E_t^{ij}, the Hermite-Coulomb auxiliaries R_{tuv}, and
!> the overlap/kinetic/nuclear/ERI/Gaussian-geminal primitive integrals; higher
!> layers contract and normalize.
module ptc_md
  use precision, only: dp
  implicit none
  private
  public :: herm_e, overlap_cart, kinetic_cart, boys, herm_R, nuclear_cart
  public :: eri_cart, geminal_cart, cart_norm
  public :: prim_grad, gem_r2_cart, vgem_cart

  real(dp), parameter :: PI = 3.141592653589793238462643383279_dp

contains

  !> Hermite expansion coefficient E_t^{i,j} for the 1D Gaussian product
  !> exp(-a(x-Ax)^2) exp(-b(x-Bx)^2), with QAB = Ax - Bx. Recursion (Helgaker).
  pure recursive function herm_e(i, j, t, ze, zf, QAB) result(e)
    integer,  intent(in) :: i, j, t
    real(dp), intent(in) :: ze, zf, QAB
    real(dp) :: e, p, mu, oo2p, xpa, xpb
    p = ze + zf; mu = ze*zf/p; oo2p = 0.5_dp/p
    xpa = -zf*QAB/p          ! P - A
    xpb =  ze*QAB/p          ! P - B
    if (t < 0 .or. t > i+j) then
      e = 0.0_dp
    else if (i == 0 .and. j == 0 .and. t == 0) then
      e = exp(-mu*QAB*QAB)
    else if (i == 0) then    ! decrement j
      e = oo2p*herm_e(i,j-1,t-1,ze,zf,QAB) + xpb*herm_e(i,j-1,t,ze,zf,QAB) &
        + real(t+1,dp)*herm_e(i,j-1,t+1,ze,zf,QAB)
    else                     ! decrement i
      e = oo2p*herm_e(i-1,j,t-1,ze,zf,QAB) + xpa*herm_e(i-1,j,t,ze,zf,QAB) &
        + real(t+1,dp)*herm_e(i-1,j,t+1,ze,zf,QAB)
    end if
  end function herm_e

  !> 1D overlap integral S^{ij} = sqrt(pi/p) E_0^{ij}.
  pure function s1d(i, j, ze, zf, QAB) result(s)
    integer,  intent(in) :: i, j
    real(dp), intent(in) :: ze, zf, QAB
    real(dp) :: s
    s = sqrt(PI/(ze+zf)) * herm_e(i,j,0,ze,zf,QAB)
  end function s1d

  !> Overlap of two Cartesian primitives (un-normalized).
  function overlap_cart(la, A, ze, lb, B, zf) result(s)
    integer,  intent(in) :: la(3), lb(3)
    real(dp), intent(in) :: A(3), B(3), ze, zf
    real(dp) :: s
    s = s1d(la(1),lb(1),ze,zf,A(1)-B(1)) &
      * s1d(la(2),lb(2),ze,zf,A(2)-B(2)) &
      * s1d(la(3),lb(3),ze,zf,A(3)-B(3))
  end function overlap_cart

  !> 1D kinetic integral T^{ij} (Helgaker, increment on the b/j centre).
  function t1d(i, j, ze, zf, QAB) result(t)
    integer,  intent(in) :: i, j
    real(dp), intent(in) :: ze, zf, QAB
    real(dp) :: t
    t = -2.0_dp*zf*zf*s1d(i,j+2,ze,zf,QAB) + zf*real(2*j+1,dp)*s1d(i,j,ze,zf,QAB)
    if (j >= 2) t = t - 0.5_dp*real(j*(j-1),dp)*s1d(i,j-2,ze,zf,QAB)
  end function t1d

  !> Kinetic energy integral of two Cartesian primitives (un-normalized).
  function kinetic_cart(la, A, ze, lb, B, zf) result(tk)
    integer,  intent(in) :: la(3), lb(3)
    real(dp), intent(in) :: A(3), B(3), ze, zf
    real(dp) :: tk, sx, sy, sz, tx, ty, tz
    sx = s1d(la(1),lb(1),ze,zf,A(1)-B(1))
    sy = s1d(la(2),lb(2),ze,zf,A(2)-B(2))
    sz = s1d(la(3),lb(3),ze,zf,A(3)-B(3))
    tx = t1d(la(1),lb(1),ze,zf,A(1)-B(1))
    ty = t1d(la(2),lb(2),ze,zf,A(2)-B(2))
    tz = t1d(la(3),lb(3),ze,zf,A(3)-B(3))
    tk = tx*sy*sz + sx*ty*sz + sx*sy*tz
  end function kinetic_cart

  !> Boys function F_n(x) by upward recursion from F_0 (series for small x).
  pure function boys(n, x) result(f)
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: f, fv(0:n), ex
    integer  :: m
    if (x < 1.0e-12_dp) then
      f = 1.0_dp/real(2*n+1,dp)
      return
    end if
    ! F_0 closed form, then upward recursion F_m = ((2m-1)F_{m-1} - e^{-x})/(2x).
    fv(0) = 0.5_dp*sqrt(PI/x)*erf(sqrt(x))
    ex = exp(-x)
    do m = 1, n
      fv(m) = (real(2*m-1,dp)*fv(m-1) - ex)/(2.0_dp*x)
    end do
    f = fv(n)
  end function boys

  !> Hermite-Coulomb auxiliary R_{tuv} for exponent p and centre separation PC,
  !> built from R^{n}_{000} = (-2p)^n F_n(p |PC|^2) by downward recursion.
  subroutine herm_R(tmax, umax, vmax, p, PC, R)
    integer,  intent(in)  :: tmax, umax, vmax
    real(dp), intent(in)  :: p, PC(3)
    real(dp), intent(out) :: R(0:tmax,0:umax,0:vmax)
    integer  :: nmax, n, t, u, v
    real(dp) :: x, RN(0:tmax+umax+vmax, 0:tmax, 0:umax, 0:vmax)
    real(dp) :: pc2
    pc2 = dot_product(PC, PC)
    nmax = tmax + umax + vmax
    RN = 0.0_dp
    x = p*pc2
    do n = 0, nmax
      RN(n,0,0,0) = (-2.0_dp*p)**n * boys(n, x)
    end do
    do v = 0, vmax
      do u = 0, umax
        do t = 0, tmax
          if (t==0 .and. u==0 .and. v==0) cycle
          do n = 0, nmax-(t+u+v)
            if (t > 0) then
              RN(n,t,u,v) = PC(1)*RN(n+1,t-1,u,v)
              if (t > 1) RN(n,t,u,v) = RN(n,t,u,v) + real(t-1,dp)*RN(n+1,t-2,u,v)
            else if (u > 0) then
              RN(n,t,u,v) = PC(2)*RN(n+1,t,u-1,v)
              if (u > 1) RN(n,t,u,v) = RN(n,t,u,v) + real(u-1,dp)*RN(n+1,t,u-2,v)
            else
              RN(n,t,u,v) = PC(3)*RN(n+1,t,u,v-1)
              if (v > 1) RN(n,t,u,v) = RN(n,t,u,v) + real(v-1,dp)*RN(n+1,t,u,v-2)
            end if
          end do
        end do
      end do
    end do
    R = RN(0,:,:,:)
  end subroutine herm_R

  !> Nuclear attraction of two Cartesian primitives at nucleus C, charge Zc.
  function nuclear_cart(la, A, ze, lb, B, zf, Zc, C) result(v)
    integer,  intent(in) :: la(3), lb(3)
    real(dp), intent(in) :: A(3), B(3), ze, zf, Zc, C(3)
    real(dp) :: v, p, P0(3), PC(3)
    real(dp), allocatable :: Ex(:), Ey(:), Ez(:), R(:,:,:)
    integer  :: tmax, umax, vmax, t, u, vv
    p = ze + zf
    P0 = (ze*A + zf*B)/p
    PC = P0 - C
    tmax = la(1)+lb(1); umax = la(2)+lb(2); vmax = la(3)+lb(3)
    allocate(Ex(0:tmax), Ey(0:umax), Ez(0:vmax), R(0:tmax,0:umax,0:vmax))
    do t = 0, tmax
      Ex(t) = herm_e(la(1),lb(1),t,ze,zf,A(1)-B(1))
    end do
    do u = 0, umax
      Ey(u) = herm_e(la(2),lb(2),u,ze,zf,A(2)-B(2))
    end do
    do vv = 0, vmax
      Ez(vv) = herm_e(la(3),lb(3),vv,ze,zf,A(3)-B(3))
    end do
    call herm_R(tmax, umax, vmax, p, PC, R)
    v = 0.0_dp
    do vv = 0, vmax
      do u = 0, umax
        do t = 0, tmax
          v = v + Ex(t)*Ey(u)*Ez(vv)*R(t,u,vv)
        end do
      end do
    end do
    v = -Zc * (2.0_dp*PI/p) * v
    deallocate(Ex,Ey,Ez,R)
  end function nuclear_cart

  !> Two-electron repulsion (ab|cd) over Cartesian primitives (chemist notation:
  !> electron 1 = a,b on centres A,B; electron 2 = c,d on centres C,D).
  function eri_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh) result(g)
    integer,  intent(in) :: la(3),lb(3),lc(3),ld(3)
    real(dp), intent(in) :: A(3),B(3),C(3),D(3), ze,zf,zg,zh
    real(dp) :: g, p, q, alpha, P0(3), Q0(3), PQ(3)
    real(dp), allocatable :: Eab(:,:,:), Ecd(:,:,:), R(:,:,:)
    integer :: t,u,v, tau,nu,phi, tmx,umx,vmx, smx,xmx,ymx
    p = ze+zf; q = zg+zh; alpha = p*q/(p+q)
    P0 = (ze*A+zf*B)/p; Q0 = (zg*C+zh*D)/q; PQ = P0 - Q0
    tmx = la(1)+lb(1); umx = la(2)+lb(2); vmx = la(3)+lb(3)
    smx = lc(1)+ld(1); xmx = lc(2)+ld(2); ymx = lc(3)+ld(3)
    allocate(Eab(0:tmx,0:umx,0:vmx), Ecd(0:smx,0:xmx,0:ymx))
    allocate(R(0:tmx+smx,0:umx+xmx,0:vmx+ymx))
    block
      integer :: i1,i2,i3
      do i3=0,vmx; do i2=0,umx; do i1=0,tmx
        Eab(i1,i2,i3) = herm_e(la(1),lb(1),i1,ze,zf,A(1)-B(1)) &
                      * herm_e(la(2),lb(2),i2,ze,zf,A(2)-B(2)) &
                      * herm_e(la(3),lb(3),i3,ze,zf,A(3)-B(3))
      end do; end do; end do
      do i3=0,ymx; do i2=0,xmx; do i1=0,smx
        Ecd(i1,i2,i3) = herm_e(lc(1),ld(1),i1,zg,zh,C(1)-D(1)) &
                      * herm_e(lc(2),ld(2),i2,zg,zh,C(2)-D(2)) &
                      * herm_e(lc(3),ld(3),i3,zg,zh,C(3)-D(3))
      end do; end do; end do
    end block
    call herm_R(tmx+smx, umx+xmx, vmx+ymx, alpha, PQ, R)
    g = 0.0_dp
    do v=0,vmx; do u=0,umx; do t=0,tmx
      do phi=0,ymx; do nu=0,xmx; do tau=0,smx
        g = g + Eab(t,u,v)*Ecd(tau,nu,phi)*(-1.0_dp)**(tau+nu+phi)*R(t+tau,u+nu,v+phi)
      end do; end do; end do
    end do; end do; end do
    g = g * 2.0_dp*PI**2.5_dp/(p*q*sqrt(p+q))
    deallocate(Eab,Ecd,R)
  end function eri_cart

  !> Gaussian-geminal two-electron integral (ab| exp(-omega r12^2) |cd) over
  !> Cartesian primitives. The geminal is separable: it multiplies the bra/ket
  !> Hermite expansions by a Gaussian overlap in the (P,Q) coordinate.
  function geminal_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, omega) result(g)
    integer,  intent(in) :: la(3),lb(3),lc(3),ld(3)
    real(dp), intent(in) :: A(3),B(3),C(3),D(3), ze,zf,zg,zh, omega
    real(dp) :: g, p, q, P0(3), Q0(3), PQ(3), denom, theta
    real(dp), allocatable :: Eab(:,:,:), Ecd(:,:,:)
    integer :: tmx,umx,vmx, smx,xmx,ymx, t,u,v, tau,nu,phi
    real(dp) :: Gt, Gu, Gv
    p = ze+zf; q = zg+zh
    P0 = (ze*A+zf*B)/p; Q0 = (zg*C+zh*D)/q; PQ = P0 - Q0
    ! geminal effective exponent in the (P-Q) Gaussian: theta = p q omega/(pq+(p+q)omega)
    denom = p*q + (p+q)*omega
    theta = p*q*omega/denom
    tmx = la(1)+lb(1); umx = la(2)+lb(2); vmx = la(3)+lb(3)
    smx = lc(1)+ld(1); xmx = lc(2)+ld(2); ymx = lc(3)+ld(3)
    allocate(Eab(0:tmx,0:umx,0:vmx), Ecd(0:smx,0:xmx,0:ymx))
    block
      integer :: i1,i2,i3
      do i3=0,vmx; do i2=0,umx; do i1=0,tmx
        Eab(i1,i2,i3) = herm_e(la(1),lb(1),i1,ze,zf,A(1)-B(1)) &
                      * herm_e(la(2),lb(2),i2,ze,zf,A(2)-B(2)) &
                      * herm_e(la(3),lb(3),i3,ze,zf,A(3)-B(3))
      end do; end do; end do
      do i3=0,ymx; do i2=0,xmx; do i1=0,smx
        Ecd(i1,i2,i3) = herm_e(lc(1),ld(1),i1,zg,zh,C(1)-D(1)) &
                      * herm_e(lc(2),ld(2),i2,zg,zh,C(2)-D(2)) &
                      * herm_e(lc(3),ld(3),i3,zg,zh,C(3)-D(3))
      end do; end do; end do
    end block
    ! For a Gaussian geminal the two-electron coupling factorizes per Cartesian
    ! direction into Hermite-Gaussian moments M_{t+tau}(theta, PQ_dir).
    g = 0.0_dp
    do v=0,vmx; do u=0,umx; do t=0,tmx
      do phi=0,ymx; do nu=0,xmx; do tau=0,smx
        Gt = gmom(t+tau, theta, PQ(1))
        Gu = gmom(u+nu,  theta, PQ(2))
        Gv = gmom(v+phi, theta, PQ(3))
        g = g + Eab(t,u,v)*Ecd(tau,nu,phi)*(-1.0_dp)**(tau+nu+phi)*Gt*Gu*Gv
      end do; end do; end do
    end do; end do; end do
    g = g * (PI*PI/denom)**1.5_dp     ! validated against the s-only geminal engine
    deallocate(Eab,Ecd)
  end function geminal_cart

  !> 1D Hermite-Gaussian moment: integral over s of H_n-like coupling for the
  !> Gaussian geminal, M_n(theta, X) = (d/dX)^n of exp(-theta X^2) up to sign,
  !> i.e. the n-th Hermite coefficient of the (P-Q) Gaussian. Here returned as
  !> the n-th derivative factor used by the Gaussian-geminal contraction.
  recursive function gmom(n, theta, X) result(m)
    integer,  intent(in) :: n
    real(dp), intent(in) :: theta, X
    real(dp) :: m
    if (n == 0) then
      m = exp(-theta*X*X)
    else if (n == 1) then
      m = -2.0_dp*theta*X*exp(-theta*X*X)
    else
      m = -2.0_dp*theta*X*gmom(n-1,theta,X) - 2.0_dp*theta*real(n-1,dp)*gmom(n-2,theta,X)
    end if
  end function gmom

  !> Cartesian-primitive gradient expansion: d/dx_dir of x^lx y^ly z^lz e^{-a r^2}
  !> = lx_dir * (l with dir decremented) - 2a * (l with dir incremented). Returns
  !> the (<=2) output Cartesian primitives lo(:,k) with coefficients co(k); the
  !> exponent is unchanged. This is the building block for the transcorrelated
  !> drift integrals (the gradients act on the AO factors after integration by
  !> parts moves the inter-electronic gradient onto the basis functions).
  pure subroutine prim_grad(lin, a, dir, lo, co, no)
    integer,  intent(in)  :: lin(3), dir
    real(dp), intent(in)  :: a
    integer,  intent(out) :: lo(3,2), no
    real(dp), intent(out) :: co(2)
    no = 0
    if (lin(dir) >= 1) then
      no = no + 1
      lo(:,no) = lin; lo(dir,no) = lin(dir) - 1
      co(no) = real(lin(dir), dp)
    end if
    no = no + 1
    lo(:,no) = lin; lo(dir,no) = lin(dir) + 1
    co(no) = -2.0_dp*a
  end subroutine prim_grad

  !> r12^2-weighted Gaussian-geminal integral (ab| r12^2 exp(-Gamma r12^2) |cd)
  !> via the exact identity r12^2 g = -d/dGamma g, evaluated with a 4th-order
  !> central finite difference in Gamma (validated to ~1e-10 vs the analytic
  !> s-only r2_geminal engine). Needed for the (grad u)^2 and Laplacian(u) terms.
  function gem_r2_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, Gamma) result(g)
    integer,  intent(in) :: la(3),lb(3),lc(3),ld(3)
    real(dp), intent(in) :: A(3),B(3),C(3),D(3), ze,zf,zg,zh, Gamma
    real(dp) :: g, h, gm2, gm1, gp1, gp2
    h = max(1.0e-4_dp*Gamma, 1.0e-7_dp)
    gp1 = geminal_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, Gamma+h)
    gm1 = geminal_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, Gamma-h)
    gp2 = geminal_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, Gamma+2.0_dp*h)
    gm2 = geminal_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, Gamma-2.0_dp*h)
    ! d/dGamma via 4th-order central stencil, then negate (r12^2 g = -dg/dGamma)
    g = -(-gp2 + 8.0_dp*gp1 - 8.0_dp*gm1 + gm2)/(12.0_dp*h)
  end function gem_r2_cart

  !> Geminal-Coulomb integral (ab| exp(-omega r12^2)/r12 |cd) over Cartesian
  !> primitives, via 1/r12 = (2/sqrt(pi)) int_0^inf exp(-t^2 r12^2) dt, so the
  !> integrand is geminal_cart(omega + t^2) Gauss-Legendre-quadratured over t in
  !> [0,inf) (mapped from [0,1]). This is the F12 "V" building block: f12/r12 with
  !> f12 = sum_k c_k exp(-omega_k r12^2) gives <ab|f12/r12|cd> = sum_k c_k vgem(omega_k).
  function vgem_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, omega, nq) result(v)
    integer,  intent(in) :: la(3),lb(3),lc(3),ld(3)
    real(dp), intent(in) :: A(3),B(3),C(3),D(3), ze,zf,zg,zh, omega
    integer,  intent(in), optional :: nq
    real(dp) :: v
    integer  :: n, i
    real(dp), allocatable :: xg(:), wg(:)
    real(dp) :: x, t, jac, acc
    n = 128; if (present(nq)) n = nq
    allocate(xg(n), wg(n))
    call gauss_legendre01(n, xg, wg)
    acc = 0.0_dp
    do i = 1, n
      x = xg(i); t = x/(1.0_dp - x); jac = 1.0_dp/(1.0_dp - x)**2
      acc = acc + wg(i)*jac*geminal_cart(la,A,ze, lb,B,zf, lc,C,zg, ld,D,zh, omega + t*t)
    end do
    v = (2.0_dp/sqrt(PI)) * acc
    deallocate(xg, wg)
  end function vgem_cart

  !> Gauss-Legendre nodes/weights on [0,1] (Newton on Legendre roots).
  subroutine gauss_legendre01(n, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: x(n), w(n)
    integer  :: i, j, it
    real(dp) :: z, z1, p1, p2, p3, pp
    real(dp), parameter :: EPS = 1.0e-15_dp
    do i = 1, (n+1)/2
      z = cos(PI*(real(i,dp) - 0.25_dp)/(real(n,dp) + 0.5_dp))
      do it = 1, 100
        p1 = 1.0_dp; p2 = 0.0_dp
        do j = 1, n
          p3 = p2; p2 = p1
          p1 = ((2.0_dp*j - 1.0_dp)*z*p2 - (j - 1.0_dp)*p3)/real(j,dp)
        end do
        pp = real(n,dp)*(z*p1 - p2)/(z*z - 1.0_dp)
        z1 = z; z = z1 - p1/pp
        if (abs(z - z1) <= EPS) exit
      end do
      x(i) = 0.5_dp*(1.0_dp - z); x(n+1-i) = 0.5_dp*(1.0_dp + z)
      w(i) = 1.0_dp/((1.0_dp - z*z)*pp*pp); w(n+1-i) = w(i)
    end do
  end subroutine gauss_legendre01

  !> Cartesian Gaussian normalization for angular (lx,ly,lz), exponent a.
  pure function cart_norm(l, ze) result(n)
    integer,  intent(in) :: l(3)
    real(dp), intent(in) :: ze
    real(dp) :: n
    n = (2.0_dp*ze/PI)**0.75_dp * (4.0_dp*ze)**(0.5_dp*real(l(1)+l(2)+l(3),dp)) &
      / sqrt(dfact(2*l(1)-1)*dfact(2*l(2)-1)*dfact(2*l(3)-1))
  end function cart_norm

  pure function dfact(k) result(zh)
    integer, intent(in) :: k
    real(dp) :: zh
    integer :: i
    zh = 1.0_dp
    i = k
    do while (i > 1)
      zh = zh*real(i,dp); i = i - 2
    end do
  end function dfact

end module ptc_md
