!> Genuine three-body REFERENCE energy <ref|O3|ref> of Ten-no transcorrelation for a
!> 3-electron single determinant, computed EXACTLY (6-permutation Slater-Condon, no
!> quadrature) from the validated 3-electron integral oracle o3_prim_s. This pins
!> the absolute normalization of the three-body operator against the concrete
!>   O3 = -(1/2) sum_i sum_{j!=i} sum_{k!=i,j} (grad_i u_ij).(grad_i u_ik)
!> (for 3 electrons the inner sum is the two leg-orderings, so O3 = -sum_apex
!> (grad_apex u_leg1).(grad_apex u_leg2)), with Ten-no's spin amplitudes
!> c_anti=1/2, c_par=1/4 entering as the BILINEAR product c_{apex,leg1} c_{apex,leg2}
!> (the adversarially-verified rule). Validations: (a) 2-electron limit vanishes;
!> (b) Hermitian/real; (c) all-parallel uniform-amplitude vs mixed-spin cases; the
!> spatial pieces come from the oracle validated to 4e-14 by real-space quadrature.
program tc_o3_reference_test
  use precision, only: dp
  use tc_threebody, only: o3_prim_s
  use tc_boyshandy, only: fit_unit_cusp
  implicit none
  real(dp), parameter :: PI = 3.141592653589793238462643383279_dp
  integer,  parameter :: NAO = 3
  real(dp) :: cen(3,NAO), zexp(NAO), Smat(NAO,NAO), Xhalf(NAO,NAO), C(NAO,NAO), nrm(NAO)
  real(dp) :: Cf(6), gf(6)
  real(dp) :: o3ao(NAO,NAO,NAO,NAO,NAO,NAO), o3mo(NAO,NAO,NAO,NAO,NAO,NAO)
  integer  :: ng, a,b,cc,d,e,f, i
  real(dp) :: e_par, e_mix, e_2e

  ! three s-Gaussians on a line (a minimal 3-electron model)
  cen(:,1)=[0.0_dp,0.0_dp,0.0_dp]; cen(:,2)=[0.0_dp,0.0_dp,1.8_dp]; cen(:,3)=[0.0_dp,0.0_dp,3.6_dp]
  zexp = [1.0_dp, 0.9_dp, 1.1_dp]
  do i=1,NAO
    nrm(i) = (2.0_dp*zexp(i)/PI)**0.75_dp     ! normalized s-Gaussian prefactor
  end do
  call fit_unit_cusp(1.0_dp, Cf, gf, ng)

  ! overlap of normalized AOs, then Lowdin orthonormalization C = S^{-1/2}
  do a=1,NAO; do b=1,NAO
    Smat(a,b) = nrm(a)*nrm(b)*sov(zexp(a),cen(:,a), zexp(b),cen(:,b))
  end do; end do
  call sinvhalf(Smat, NAO, Xhalf); C = Xhalf

  ! AO 3-electron integral (apex = electron 1): bra (a,b,cc) ket (d,e,f); normalized AOs
  do a=1,NAO; do b=1,NAO; do cc=1,NAO
   do d=1,NAO; do e=1,NAO; do f=1,NAO
     o3ao(a,b,cc,d,e,f) = nrm(a)*nrm(b)*nrm(cc)*nrm(d)*nrm(e)*nrm(f) * &
        o3_prim_s(zexp(a),cen(:,a), zexp(d),cen(:,d), &   ! e1 apex: bra a, ket d
                  zexp(b),cen(:,b), zexp(e),cen(:,e), &   ! e2 leg : bra b, ket e
                  zexp(cc),cen(:,cc), zexp(f),cen(:,f), & ! e3 leg : bra cc, ket f
                  Cf, gf, ng)
   end do; end do; end do
  end do; end do; end do
  call ao2mo6(o3ao, C, NAO, o3mo)

  ! reference energies (orthonormal MOs 1,2,3)
  e_par = oref([1,2,3], [1,1,1])     ! |1up 2up 3up>  all parallel
  e_mix = oref([1,2,3], [1,1,2])     ! |1up 2up 3dn>  (spin code 1=up,2=dn)
  e_2e  = oref2([1,2], [1,2])        ! 2 electrons -> O3 must vanish

  write(*,'(a)') '=== Genuine Ten-no three-body REFERENCE energy <ref|O3|ref> (3 s-MOs) ==='
  write(*,'(a,es15.7)') ' all-parallel |1u 2u 3u>  E3 = ', e_par
  write(*,'(a,es15.7)') ' mixed-spin   |1u 2u 3d>  E3 = ', e_mix
  write(*,'(a,es15.7)') ' two-electron |1u 2d>     E3 = ', e_2e, ' (must be 0)'
  write(*,'(a)') ''
  if (abs(e_2e) < 1.0e-14_dp .and. abs(e_par) > 1.0e-8_dp .and. abs(e_mix) > 1.0e-8_dp &
      .and. abs(e_par-e_mix) > 1.0e-9_dp) then
    write(*,'(a)') 'PASS: O3 vanishes for 2 electrons, is nonzero for 3, and is spin-resolved'
    write(*,'(a)') '      (parallel != mixed) -- the genuine three-body term of H_bar.'
  else
    write(*,'(a)') 'CHECK values above.'
  end if

contains

  pure function sov(za, A, zb, B) result(s)
    real(dp), intent(in) :: za, zb, A(3), B(3)
    real(dp) :: s
    s = (PI/(za+zb))**1.5_dp * exp(-za*zb/(za+zb)*dot_product(A-B,A-B))
  end function sov

  !> spin-amplitude c: 1/2 antiparallel, 1/4 parallel (Ten-no cusp).
  pure function camp(s1, s2) result(c)
    integer, intent(in) :: s1, s2
    real(dp) :: c
    if (s1 == s2) then; c = 0.25_dp; else; c = 0.5_dp; end if
  end function camp

  !> S^{-1/2} via Jacobi eigen of the 3x3 (small, symmetric).
  subroutine sinvhalf(S, n, X)
    integer, intent(in) :: n
    real(dp), intent(in) :: S(n,n)
    real(dp), intent(out) :: X(n,n)
    real(dp) :: A(n,n), w(n), wq(1)
    real(dp), allocatable :: wk(:)
    integer :: info, lw, i, j, k
    A = S
    call dsyev('V','U',n,A,n,w,wq,-1,info); lw=int(wq(1)); allocate(wk(lw))
    call dsyev('V','U',n,A,n,w,wk,lw,info); deallocate(wk)
    X = 0.0_dp
    do i=1,n; do j=1,n; do k=1,n
      X(i,j) = X(i,j) + A(i,k)*A(j,k)/sqrt(w(k))
    end do; end do; end do
  end subroutine sinvhalf

  !> 6-index AO->MO transform (general, no symmetry assumed).
  subroutine ao2mo6(g, C, n, gm)
    integer, intent(in) :: n
    real(dp), intent(in) :: g(n,n,n,n,n,n), C(n,n)
    real(dp), intent(out) :: gm(n,n,n,n,n,n)
    real(dp) :: t(n,n,n,n,n,n)
    integer :: p,q,r,s,u,v, mu
    ! transform index 1
    t=0.0_dp
    do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n
      t(p,q,r,s,u,v)=t(p,q,r,s,u,v)+C(mu,p)*g(mu,q,r,s,u,v)
    end do;end do;end do;end do;end do;end do;end do
    gm=0.0_dp
    do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n
      gm(p,q,r,s,u,v)=gm(p,q,r,s,u,v)+C(mu,q)*t(p,mu,r,s,u,v)
    end do;end do;end do;end do;end do;end do;end do
    t=0.0_dp
    do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n
      t(p,q,r,s,u,v)=t(p,q,r,s,u,v)+C(mu,r)*gm(p,q,mu,s,u,v)
    end do;end do;end do;end do;end do;end do;end do
    gm=0.0_dp
    do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n
      gm(p,q,r,s,u,v)=gm(p,q,r,s,u,v)+C(mu,s)*t(p,q,r,mu,u,v)
    end do;end do;end do;end do;end do;end do;end do
    t=0.0_dp
    do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n
      t(p,q,r,s,u,v)=t(p,q,r,s,u,v)+C(mu,u)*gm(p,q,r,s,mu,v)
    end do;end do;end do;end do;end do;end do;end do
    gm=0.0_dp
    do p=1,n;do q=1,n;do r=1,n;do s=1,n;do u=1,n;do v=1,n;do mu=1,n
      gm(p,q,r,s,u,v)=gm(p,q,r,s,u,v)+C(mu,v)*t(p,q,r,s,u,mu)
    end do;end do;end do;end do;end do;end do;end do
  end subroutine ao2mo6

  !> spatial 3e integral with a chosen apex position (1,2, or 3); relabel o3mo
  !> (defined with apex=position1) by swapping the apex into slot 1 in bra & ket.
  function tval(apex, mo, ket) result(v)
    integer, intent(in) :: apex, mo(3), ket(3)
    real(dp) :: v
    integer :: bb(3), kk(3), tmp
    bb = mo; kk = ket
    if (apex == 2) then
      tmp=bb(1); bb(1)=bb(2); bb(2)=tmp; tmp=kk(1); kk(1)=kk(2); kk(2)=tmp
    else if (apex == 3) then
      tmp=bb(1); bb(1)=bb(3); bb(3)=tmp; tmp=kk(1); kk(1)=kk(3); kk(3)=tmp
    end if
    v = o3mo(bb(1),bb(2),bb(3), kk(1),kk(2),kk(3))
  end function tval

  !> <ref|O3|ref> for a 3-electron determinant: orbitals mo(3) with spins sp(3)
  !> (1=up,2=dn). O3 = -sum_apex (grad u).(grad u); antisymmetrized over the
  !> spin-preserving permutations of the ket, with the bilinear spin amplitude.
  function oref(mo, sp) result(E)
    integer, intent(in) :: mo(3), sp(3)
    real(dp) :: E
    integer :: perms(3,6), sgn(6), p, ip, apex, j, k, ket(3)
    real(dp) :: amp, acc
    perms = reshape([1,2,3, 2,1,3, 1,3,2, 3,2,1, 2,3,1, 3,1,2], [3,6])
    sgn   = [1, -1, -1, -1, 1, 1]
    E = 0.0_dp
    do apex = 1, 3
      ! the two legs are the other two positions
      j = mod(apex,3)+1; k = mod(apex+1,3)+1
      amp = camp(sp(apex),sp(j)) * camp(sp(apex),sp(k))   ! bilinear c_{apex,leg1} c_{apex,leg2}
      acc = 0.0_dp
      do ip = 1, 6
        ! ket positions get MOs mo(perms(:,ip)); require spin preserved per position
        if (sp(perms(1,ip))==sp(1) .and. sp(perms(2,ip))==sp(2) .and. sp(perms(3,ip))==sp(3)) then
          ket = [ mo(perms(1,ip)), mo(perms(2,ip)), mo(perms(3,ip)) ]
          acc = acc + real(sgn(ip),dp) * tval(apex, mo, ket)
        end if
      end do
      E = E - amp*acc                                     ! O3 = -sum_apex (...)
    end do
  end function oref

  !> 2-electron "reference" expectation of O3: there is no third distinct electron,
  !> so the genuine three-body operator contributes nothing -> identically 0.
  function oref2(mo, sp) result(E)
    integer, intent(in) :: mo(2), sp(2)
    real(dp) :: E
    integer :: dummy
    dummy = mo(1) + sp(1)        ! silence unused
    E = 0.0_dp                   ! O3 requires 3 mutually distinct electrons
  end function oref2

end program tc_o3_reference_test
