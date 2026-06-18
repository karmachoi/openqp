!> Genuine Boys-Handy / Ten-no transcorrelation: the explicit first-quantized
!> similarity transform H_bar = e^{-tau} H e^{tau} with tau = sum_{i<j} u(r_ij).
!>
!> Because the correlation factor u is multiplicative it commutes with every part
!> of H except the kinetic energy, so the Baker-Campbell-Hausdorff series
!> terminates exactly: for a pair (1,2),
!>
!>   H_bar(1,2) = 1/r12  -  [ nabla^2 u + (u')^2 ]        (Hermitian scalar W)
!>                       -  [ (grad_1 u).grad_1 + (grad_2 u).grad_2 ]  (drift D)
!>
!> plus, for >=3 electrons, the genuine three-body term -1/2 sum_i sum_{l/=m/=i}
!> (grad_i u_il).(grad_i u_im) (handled separately by normal ordering). This module
!> builds the TWO-body part exactly -- which IS the complete H_bar for a two-electron
!> system such as H2 -- using Ten-no's Slater geminal u = c * f, f = -(1/gamma)
!> e^{-gamma r12}, with cusp-fixed amplitudes c = 1/2 (antiparallel) and 1/4
!> (parallel). f is represented by STG-6G; only u',u'' enter, so f's sign/offset
!> are immaterial.
!>
!> The drift is reduced to plain Gaussian-geminal integrals over gradient-shifted
!> Cartesian primitives by integrating the inter-electronic gradient by parts:
!>   <ab|(grad_1 u).grad_1|cd> -> sum_k C_k sum_d [ (d_a,d_c|b,d) + (a,dd_c|b,d) ].
!> The scalar W uses the geminal and r12^2-geminal primitives. All integrals carry
!> the unit-cusp factor; the caller scales the linear-in-u part by c and the
!> quadratic (u')^2 part by c^2, so spin amplitudes are applied at assembly.
module tc_boyshandy
  use precision, only: dp
  use ptc_md
  use ptc_ao, only: ao_enum, ao_ncart, stg6
  implicit none
  private
  public :: build_tc2e_ao, fit_unit_cusp, drift_prim

contains

  !> Unit-cusp correlation factor f(r) = -(1/gamma) e^{-gamma r} fit as
  !> sum_k C_k exp(-g_k r^2): from e^{-gamma r} ~ sum_k a_k e^{-w_k r^2} (STG-6G),
  !> C_k = -a_k/gamma, g_k = w_k.
  subroutine fit_unit_cusp(gamma, C, g, ng)
    real(dp), intent(in)  :: gamma
    real(dp), intent(out) :: C(:), g(:)
    integer,  intent(out) :: ng
    real(dp) :: aa(6), ww(6)
    integer  :: k
    call stg6(gamma, aa, ww, ng)
    do k = 1, ng
      C(k) = -aa(k)/gamma
      g(k) = ww(k)
    end do
  end subroutine fit_unit_cusp

  !> Second Cartesian derivative d^2/dx_dir^2 of a primitive (compose prim_grad).
  pure subroutine grad2(lin, a, dir, lo, co, no)
    integer,  intent(in)  :: lin(3), dir
    real(dp), intent(in)  :: a
    integer,  intent(out) :: lo(3,4), no
    real(dp), intent(out) :: co(4)
    integer  :: l1(3,2), n1, k, l2(3,2), n2, j
    real(dp) :: c1(2), c2(2)
    call prim_grad(lin, a, dir, l1, c1, n1)
    no = 0
    do k = 1, n1
      call prim_grad(l1(:,k), a, dir, l2, c2, n2)
      do j = 1, n2
        no = no + 1
        lo(:,no) = l2(:,j); co(no) = c1(k)*c2(j)
      end do
    end do
  end subroutine grad2

  !> Drift integral <IJ|(grad_1 u).grad_1 + (grad_2 u).grad_2|KL> for ONE Gaussian
  !> geminal exp(-Gam r12^2) (the C_k weight and the leading minus sign of D are
  !> applied by the caller). Physicist layout: electron 1 = (I bra, K ket), electron
  !> 2 = (J bra, L ket). Returns term1+term2 from the integration-by-parts reduction.
  function drift_prim(lI,RI,eI, lJ,RJ,eJ, lK,RK,eK, lL,RL,eL, Gam) result(d)
    integer,  intent(in) :: lI(3),lJ(3),lK(3),lL(3)
    real(dp), intent(in) :: RI(3),RJ(3),RK(3),RL(3), eI,eJ,eK,eL, Gam
    real(dp) :: d
    integer  :: dir, i1, i2
    integer  :: gI(3,2), gK(3,2), gJ(3,2), gL(3,2), nI,nK,nJ,nL
    real(dp) :: cI(2), cK(2), cJ(2), cL(2)
    integer  :: hK(3,4), hL(3,4), nKK, nLL
    real(dp) :: dK(4), dL(4)
    d = 0.0_dp
    do dir = 1, 3
      call prim_grad(lI, eI, dir, gI, cI, nI)
      call prim_grad(lK, eK, dir, gK, cK, nK)
      call prim_grad(lJ, eJ, dir, gJ, cJ, nJ)
      call prim_grad(lL, eL, dir, gL, cL, nL)
      call grad2(lK, eK, dir, hK, dK, nKK)
      call grad2(lL, eL, dir, hL, dL, nLL)
      ! term1a: (d_I, d_K | J, L)   (electron-1 gradient lands on both I and K)
      do i1 = 1, nI; do i2 = 1, nK
        d = d + cI(i1)*cK(i2)* &
            geminal_cart(gI(:,i1),RI,eI, gK(:,i2),RK,eK, lJ,RJ,eJ, lL,RL,eL, Gam)
      end do; end do
      ! term1b: (I, dd_K | J, L)
      do i2 = 1, nKK
        d = d + dK(i2)* &
            geminal_cart(lI,RI,eI, hK(:,i2),RK,eK, lJ,RJ,eJ, lL,RL,eL, Gam)
      end do
      ! term2a: (I, K | d_J, d_L)   (electron-2 gradient)
      do i1 = 1, nJ; do i2 = 1, nL
        d = d + cJ(i1)*cL(i2)* &
            geminal_cart(lI,RI,eI, lK,RK,eK, gJ(:,i1),RJ,eJ, gL(:,i2),RL,eL, Gam)
      end do; end do
      ! term2b: (I, K | J, dd_L)
      do i2 = 1, nLL
        d = d + dL(i2)* &
            geminal_cart(lI,RI,eI, lK,RK,eK, lJ,RJ,eJ, hL(:,i2),RL,eL, Gam)
      end do
    end do
  end function drift_prim

  !> Build the genuine two-body transcorrelation tensors in physicist layout
  !> <IJ|.|KL> (electron 1 = I,K ; electron 2 = J,L), normalized identically to
  !> build_ints (AO self-overlap = 1):
  !>   Lin(I,J,K,L)  = <IJ| -nabla^2 f  +  D[f] |KL>     (scales with c)
  !>   Quad(I,J,K,L) = <IJ| -(f')^2 |KL>                 (scales with c^2)
  !> so the spin-resolved correction is c*Lin + c^2*Quad. OpenMP over the bra index.
  subroutine build_tc2e_ao(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, gamma, nao, Lin, Quad)
    integer,  intent(in)  :: nsh, shl_l(nsh), shl_np(nsh), nao
    real(dp), intent(in)  :: shl_e(:,:), shl_c(:,:), shl_r(3,nsh), gamma
    real(dp), intent(out) :: Lin(:,:,:,:), Quad(:,:,:,:)
    integer, allocatable :: ao_sh(:), ao_l(:,:)
    real(dp), allocatable :: nrm(:)
    real(dp) :: Cf(6), gf(6)
    integer  :: maxao, I, J, K, L, ng, nn
    maxao = 0
    do I = 1, nsh
      maxao = maxao + ao_ncart(shl_l(I))
    end do
    allocate(ao_sh(maxao), ao_l(3,maxao), nrm(maxao))
    call ao_enum(nsh, shl_l, nn, ao_sh, ao_l)
    call fit_unit_cusp(gamma, Cf, gf, ng)
    ! AO normalization 1/sqrt(self-overlap), identical to build_ints
    do I = 1, nao
      nrm(I) = 1.0_dp/sqrt(selfov(I))
    end do
    !$omp parallel do collapse(2) default(shared) private(I,J,K,L) schedule(dynamic)
    do I = 1, nao
      do J = 1, nao
        do K = 1, nao
          do L = 1, nao
            Lin(I,J,K,L)  = lin_el(I,J,K,L) *nrm(I)*nrm(J)*nrm(K)*nrm(L)
            Quad(I,J,K,L) = quad_el(I,J,K,L)*nrm(I)*nrm(J)*nrm(K)*nrm(L)
          end do
        end do
      end do
    end do
    deallocate(ao_sh, ao_l, nrm)

  contains

    real(dp) function selfov(I) result(v)
      integer, intent(in) :: I
      integer :: s, p, q
      real(dp) :: e1, e2, c1, c2
      s = ao_sh(I); v = 0.0_dp
      do p = 1, shl_np(s); e1 = shl_e(p,s); c1 = shl_c(p,s)*cart_norm(ao_l(:,I),e1)
        do q = 1, shl_np(s); e2 = shl_e(q,s); c2 = shl_c(q,s)*cart_norm(ao_l(:,I),e2)
          v = v + c1*c2*overlap_cart(ao_l(:,I),shl_r(:,s),e1, ao_l(:,I),shl_r(:,s),e2)
        end do
      end do
    end function selfov

    !> linear-in-u part: -nabla^2 f (scalar) + drift, contracted over primitives.
    real(dp) function lin_el(I,J,K,L) result(v)
      integer, intent(in) :: I,J,K,L
      integer :: si,sj,sk,sl, pi,pj,pk,pl, kk
      real(dp) :: ci,cj,ck,cl, ei,ej,ek,el, RIv(3),RJv(3),RKv(3),RLv(3)
      real(dp) :: g0, gr2, drift, cp
      si=ao_sh(I); sj=ao_sh(J); sk=ao_sh(K); sl=ao_sh(L)
      RIv=shl_r(:,si); RJv=shl_r(:,sj); RKv=shl_r(:,sk); RLv=shl_r(:,sl)
      v = 0.0_dp
      do pi=1,shl_np(si); ei=shl_e(pi,si); ci=shl_c(pi,si)*cart_norm(ao_l(:,I),ei)
       do pj=1,shl_np(sj); ej=shl_e(pj,sj); cj=shl_c(pj,sj)*cart_norm(ao_l(:,J),ej)
        do pk=1,shl_np(sk); ek=shl_e(pk,sk); ck=shl_c(pk,sk)*cart_norm(ao_l(:,K),ek)
         do pl=1,shl_np(sl); el=shl_e(pl,sl); cl=shl_c(pl,sl)*cart_norm(ao_l(:,L),el)
          cp = ci*cj*ck*cl
          do kk = 1, ng
            ! physicist e1=(I,K), e2=(J,L)
            g0  = geminal_cart(ao_l(:,I),RIv,ei, ao_l(:,K),RKv,ek, &
                               ao_l(:,J),RJv,ej, ao_l(:,L),RLv,el, gf(kk))
            gr2 = gem_r2_cart(ao_l(:,I),RIv,ei, ao_l(:,K),RKv,ek, &
                              ao_l(:,J),RJv,ej, ao_l(:,L),RLv,el, gf(kk))
            drift = drift_prim(ao_l(:,I),RIv,ei, ao_l(:,J),RJv,ej, &
                               ao_l(:,K),RKv,ek, ao_l(:,L),RLv,el, gf(kk))
            ! -nabla^2 f = sum_k C_k (6 g_k - 4 g_k^2 r^2) e^{-g_k r^2}; + drift
            v = v + cp*Cf(kk)*( 6.0_dp*gf(kk)*g0 - 4.0_dp*gf(kk)*gf(kk)*gr2 + drift )
          end do
         end do
        end do
       end do
      end do
    end function lin_el

    !> quadratic part: -(f')^2 = -sum_{k,k'} C_k C_k' 4 g_k g_k' r^2 e^{-(g_k+g_k')r^2}.
    real(dp) function quad_el(I,J,K,L) result(v)
      integer, intent(in) :: I,J,K,L
      integer :: si,sj,sk,sl, pi,pj,pk,pl, k1,k2
      real(dp) :: ci,cj,ck,cl, ei,ej,ek,el, RIv(3),RJv(3),RKv(3),RLv(3)
      real(dp) :: gr2, cp
      si=ao_sh(I); sj=ao_sh(J); sk=ao_sh(K); sl=ao_sh(L)
      RIv=shl_r(:,si); RJv=shl_r(:,sj); RKv=shl_r(:,sk); RLv=shl_r(:,sl)
      v = 0.0_dp
      do pi=1,shl_np(si); ei=shl_e(pi,si); ci=shl_c(pi,si)*cart_norm(ao_l(:,I),ei)
       do pj=1,shl_np(sj); ej=shl_e(pj,sj); cj=shl_c(pj,sj)*cart_norm(ao_l(:,J),ej)
        do pk=1,shl_np(sk); ek=shl_e(pk,sk); ck=shl_c(pk,sk)*cart_norm(ao_l(:,K),ek)
         do pl=1,shl_np(sl); el=shl_e(pl,sl); cl=shl_c(pl,sl)*cart_norm(ao_l(:,L),el)
          cp = ci*cj*ck*cl
          do k1 = 1, ng; do k2 = 1, ng
            gr2 = gem_r2_cart(ao_l(:,I),RIv,ei, ao_l(:,K),RKv,ek, &
                              ao_l(:,J),RJv,ej, ao_l(:,L),RLv,el, gf(k1)+gf(k2))
            v = v - cp*Cf(k1)*Cf(k2)*4.0_dp*gf(k1)*gf(k2)*gr2
          end do; end do
         end do
        end do
       end do
      end do
    end function quad_el

  end subroutine build_tc2e_ao

end module tc_boyshandy
