!> Independent validation of the genuine 3-electron integral o3_geminal_s by
!> real-space quadrature, in two layers that together certify the full chain:
!>  (1) LEG: the analytic leg integral J2_d(r1) = int Omega2(r2)(x1-x2)_d
!>      e^{-gm|r1-r2|^2} dr2 vs a 3D box quadrature over r2, at several r1.
!>  (2) APEX: o3_geminal_s vs 4 gm gn * (3D box quadrature over r1 of
!>      Omega1(r1) sum_d J2_d(r1) J3_d(r1)), using the analytic legs.
program tc_threebody_test
  use precision, only: dp
  use tc_threebody, only: o3_geminal_s
  implicit none
  real(dp), parameter :: PI = 3.141592653589793238462643383279_dp
  real(dp) :: zA,A(3), zAp,Ap(3), zB,B(3), zBp,Bp(3), zC,C(3), zCp,Cp(3), gm, gn
  real(dp) :: e1,Q1(3),K1, e2,Q2(3),K2, e3,Q3(3),K3
  real(dp) :: ana, num, r1(3), mx
  integer  :: t, d

  zA=1.1_dp; A=[0.0_dp,0.0_dp,0.0_dp]
  zAp=0.8_dp; Ap=[0.1_dp,-0.2_dp,0.3_dp]
  zB=0.9_dp; B=[0.0_dp,0.0_dp,1.3_dp]
  zBp=1.3_dp; Bp=[0.2_dp,0.4_dp,1.0_dp]
  zC=0.7_dp; C=[-0.3_dp,0.5_dp,0.4_dp]
  zCp=1.0_dp; Cp=[0.1_dp,0.1_dp,-0.6_dp]
  gm=0.6_dp; gn=1.1_dp

  call sprod(zA,A, zAp,Ap, e1,Q1,K1)
  call sprod(zB,B, zBp,Bp, e2,Q2,K2)
  call sprod(zC,C, zCp,Cp, e3,Q3,K3)

  ! ---- Layer (1): leg J2_d at several r1 ----
  mx = 0.0_dp
  do t = 1, 4
    r1 = [0.2_dp*t, -0.1_dp*t, 1.0_dp+0.1_dp*t]
    do d = 1, 3
      ana = leg_analytic(e2,Q2,K2, gm, r1, d)
      num = leg_numeric (e2,Q2,K2, gm, r1, d)
      mx = max(mx, abs(ana-num))
    end do
  end do
  write(*,'(a,es10.2)') 'CHECK leg J2_d analytic vs 3D quadrature  : max|diff| = ', mx

  ! ---- Layer (2): apex assembly ----
  ana = o3_geminal_s(zA,A, zAp,Ap, zB,B, zBp,Bp, zC,C, zCp,Cp, gm, gn)
  num = apex_numeric()
  write(*,'(a,es12.4)') 'o3_geminal_s analytic = ', ana
  write(*,'(a,es12.4)') 'o3_geminal_s quadr    = ', num
  write(*,'(a,es10.2)') 'CHECK 3-electron O3 analytic vs quadrature: rel|diff| = ', abs(ana-num)/abs(ana)

  if (mx < 1.0e-6_dp .and. abs(ana-num)/abs(ana) < 1.0e-4_dp) then
    write(*,'(a)') 'PASS: genuine 3-electron integral o3_geminal_s validated by real-space quadrature.'
  else
    write(*,'(a)') 'CHECK values above.'
  end if

contains

  pure subroutine sprod(za, A, zb, B, e, Q, K)
    real(dp), intent(in)  :: za, zb, A(3), B(3)
    real(dp), intent(out) :: e, Q(3), K
    e = za + zb; Q = (za*A + zb*B)/e
    K = exp(-za*zb/e * dot_product(A-B, A-B))
  end subroutine sprod

  function leg_analytic(e2,Q2,K2, gm, r1, d) result(v)
    real(dp), intent(in) :: e2,Q2(3),K2, gm, r1(3)
    integer,  intent(in) :: d
    real(dp) :: v, mu2
    mu2 = e2*gm/(e2+gm)
    v = K2*(PI/(e2+gm))**1.5_dp*(e2/(e2+gm))*(r1(d)-Q2(d))*exp(-mu2*dot_product(r1-Q2,r1-Q2))
  end function leg_analytic

  function leg_numeric(e2,Q2,K2, gm, r1, d) result(v)
    real(dp), intent(in) :: e2,Q2(3),K2, gm, r1(3)
    integer,  intent(in) :: d
    real(dp) :: v, r2(3), h, lo(3), om2, ge
    integer  :: ix,iy,iz, ng
    ng = 80; h = 12.0_dp/real(ng,dp)
    lo = Q2 - 6.0_dp
    v = 0.0_dp
    do ix=1,ng; r2(1)=lo(1)+(ix-0.5_dp)*h
     do iy=1,ng; r2(2)=lo(2)+(iy-0.5_dp)*h
      do iz=1,ng; r2(3)=lo(3)+(iz-0.5_dp)*h
        om2 = K2*exp(-e2*dot_product(r2-Q2,r2-Q2))
        ge  = exp(-gm*dot_product(r1-r2,r1-r2))
        v = v + om2*(r1(d)-r2(d))*ge
      end do
     end do
    end do
    v = v*h*h*h
  end function leg_numeric

  function apex_numeric() result(v)
    real(dp) :: v, r1g(3), h, lo(3), om1, s, j2, j3
    integer  :: ix,iy,iz, ng, d
    ng = 80; h = 14.0_dp/real(ng,dp)
    lo = (Q1+Q2+Q3)/3.0_dp - 7.0_dp
    v = 0.0_dp
    do ix=1,ng; r1g(1)=lo(1)+(ix-0.5_dp)*h
     do iy=1,ng; r1g(2)=lo(2)+(iy-0.5_dp)*h
      do iz=1,ng; r1g(3)=lo(3)+(iz-0.5_dp)*h
        om1 = K1*exp(-e1*dot_product(r1g-Q1,r1g-Q1))
        s = 0.0_dp
        do d=1,3
          j2 = leg_analytic(e2,Q2,K2, gm, r1g, d)
          j3 = leg_analytic(e3,Q3,K3, gn, r1g, d)
          s = s + j2*j3
        end do
        v = v + om1*s
      end do
     end do
    end do
    v = v*h*h*h*4.0_dp*gm*gn
  end function apex_numeric

end program tc_threebody_test
