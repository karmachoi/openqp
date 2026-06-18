!> General-L AO integral builder: contracts the McMurchie-Davidson primitive
!> integrals (ptc_md) over shells of arbitrary angular momentum into normalized
!> AO matrices S, Hcore, chemist ERI, and the Gaussian-geminal tensor. This is the
!> contraction layer that lets the pTC-MRSF-CIS pipeline use p/diffuse (and d)
!> bases instead of s-only. Cartesian functions (s,p identical to spherical;
!> d+ would need a spherical transform, added later).
!>
!> Shell arrays: shl_l(nsh) angular momentum, shl_np(nsh) #prims,
!> shl_e(mp,nsh)/shl_c(mp,nsh) primitive exponents/coeffs, shl_r(3,nsh) centre.
module ptc_ao
  use precision, only: dp
  use ptc_md
  implicit none
  private
  public :: ao_enum, ao_ncart, build_ints, build_geminal_ao, stg6

contains

  integer function ao_ncart(l) result(n)   ! # Cartesian functions for shell L
    integer, intent(in) :: l
    n = (l+1)*(l+2)/2
  end function ao_ncart

  !> Enumerate AO functions: ao_sh(nao) parent shell, ao_l(3,nao) = (lx,ly,lz).
  subroutine ao_enum(nsh, shl_l, nao, ao_sh, ao_l)
    integer, intent(in)  :: nsh, shl_l(nsh)
    integer, intent(out) :: nao, ao_sh(*), ao_l(3,*)
    integer :: s, l, lx, ly, lz, k
    nao = 0
    do s = 1, nsh
      l = shl_l(s)
      do lx = l, 0, -1
        do ly = l-lx, 0, -1
          lz = l - lx - ly
          nao = nao + 1
          ao_sh(nao) = s
          ao_l(1,nao) = lx; ao_l(2,nao) = ly; ao_l(3,nao) = lz
        end do
      end do
    end do
    k = nao   ! silence
  end subroutine ao_enum

  !> Build normalized AO S, Hcore (=T+V), and chemist ERI eri_c(I,J,K,L)=(IJ|KL).
  subroutine build_ints(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, nat, zat, rat, &
                        nao, S, Hc, eri_c)
    integer,  intent(in)  :: nsh, shl_l(nsh), shl_np(nsh), nat
    real(dp), intent(in)  :: shl_e(:,:), shl_c(:,:), shl_r(3,nsh), zat(:), rat(:,:)
    integer,  intent(out) :: nao
    real(dp), intent(out) :: S(:,:), Hc(:,:), eri_c(:,:,:,:)
    integer, allocatable :: ao_sh(:), ao_l(:,:)
    real(dp), allocatable :: nrm(:)
    integer :: maxao, I, J, K, L
    maxao = 0
    do I = 1, nsh
      maxao = maxao + ao_ncart(shl_l(I))
    end do
    allocate(ao_sh(maxao), ao_l(3,maxao), nrm(maxao))
    call ao_enum(nsh, shl_l, nao, ao_sh, ao_l)
    ! one-electron
    do I = 1, nao
      do J = 1, nao
        S(I,J)  = c1e(I,J,0)
        Hc(I,J) = c1e(I,J,1)      ! T + V
      end do
    end do
    ! AO normalization so S_II = 1
    do I = 1, nao
      nrm(I) = 1.0_dp/sqrt(S(I,I))
    end do
    do I = 1, nao
      do J = 1, nao
        S(I,J)  = S(I,J)*nrm(I)*nrm(J)
        Hc(I,J) = Hc(I,J)*nrm(I)*nrm(J)
      end do
    end do
    ! two-electron (chemist (IJ|KL)), normalized
    do I = 1, nao; do J = 1, nao; do K = 1, nao; do L = 1, nao
      eri_c(I,J,K,L) = c2e(I,J,K,L)*nrm(I)*nrm(J)*nrm(K)*nrm(L)
    end do; end do; end do; end do
    deallocate(ao_sh, ao_l, nrm)

  contains

    real(dp) function c1e(I, J, mode) result(v)   ! mode 0=S, 1=T+V
      integer, intent(in) :: I, J, mode
      integer :: si, sj, pi, pj, ia
      real(dp) :: ci, cj, ei, ej, RI(3), RJ(3)
      si = ao_sh(I); sj = ao_sh(J); RI = shl_r(:,si); RJ = shl_r(:,sj)
      v = 0.0_dp
      do pi = 1, shl_np(si)
        ei = shl_e(pi,si); ci = shl_c(pi,si)*cart_norm(ao_l(:,I), ei)
        do pj = 1, shl_np(sj)
          ej = shl_e(pj,sj); cj = shl_c(pj,sj)*cart_norm(ao_l(:,J), ej)
          if (mode == 0) then
            v = v + ci*cj*overlap_cart(ao_l(:,I),RI,ei, ao_l(:,J),RJ,ej)
          else
            v = v + ci*cj*kinetic_cart(ao_l(:,I),RI,ei, ao_l(:,J),RJ,ej)
            do ia = 1, nat
              v = v + ci*cj*nuclear_cart(ao_l(:,I),RI,ei, ao_l(:,J),RJ,ej, zat(ia),rat(:,ia))
            end do
          end if
        end do
      end do
    end function c1e

    real(dp) function c2e(I, J, K, L) result(v)   ! chemist (IJ|KL)
      integer, intent(in) :: I, J, K, L
      integer :: si,sj,sk,sl, pi,pj,pk,pl
      real(dp) :: ci,cj,ck,cl, ei,ej,ek,el, RI(3),RJ(3),RK(3),RL(3)
      si=ao_sh(I); sj=ao_sh(J); sk=ao_sh(K); sl=ao_sh(L)
      RI=shl_r(:,si); RJ=shl_r(:,sj); RK=shl_r(:,sk); RL=shl_r(:,sl)
      v = 0.0_dp
      do pi=1,shl_np(si); ei=shl_e(pi,si); ci=shl_c(pi,si)*cart_norm(ao_l(:,I),ei)
        do pj=1,shl_np(sj); ej=shl_e(pj,sj); cj=shl_c(pj,sj)*cart_norm(ao_l(:,J),ej)
          do pk=1,shl_np(sk); ek=shl_e(pk,sk); ck=shl_c(pk,sk)*cart_norm(ao_l(:,K),ek)
            do pl=1,shl_np(sl); el=shl_e(pl,sl); cl=shl_c(pl,sl)*cart_norm(ao_l(:,L),el)
              ! chemist (IJ|KL): e1 = I,J ; e2 = K,L
              v = v + ci*cj*ck*cl*eri_cart(ao_l(:,I),RI,ei, ao_l(:,J),RJ,ej, &
                                           ao_l(:,K),RK,ek, ao_l(:,L),RL,el)
            end do
          end do
        end do
      end do
    end function c2e

  end subroutine build_ints

  !> Gaussian-geminal AO tensor for the Slater factor exp(-gamma r12) via STG-6G,
  !> physicist layout Gao(I,J,K,L) = <IJ|f12|KL> (e1=I,K ; e2=J,L), normalized.
  subroutine build_geminal_ao(nsh, shl_l, shl_np, shl_e, shl_c, shl_r, gamma, nao, S, Gao)
    integer,  intent(in)  :: nsh, shl_l(nsh), shl_np(nsh)
    real(dp), intent(in)  :: shl_e(:,:), shl_c(:,:), shl_r(3,nsh), gamma, S(:,:)
    integer,  intent(in)  :: nao
    real(dp), intent(out) :: Gao(:,:,:,:)
    integer, allocatable :: ao_sh(:), ao_l(:,:)
    real(dp), allocatable :: nrm(:)
    real(dp) :: cc(6), omg(6)
    integer :: maxao, I, J, K, L, ng, nn, tmp
    maxao = 0
    do I = 1, nsh
      maxao = maxao + ao_ncart(shl_l(I))
    end do
    allocate(ao_sh(maxao), ao_l(3,maxao), nrm(maxao))
    call ao_enum(nsh, shl_l, nn, ao_sh, ao_l)
    do I = 1, nao
      nrm(I) = 1.0_dp/sqrt(S(I,I))   ! S already normalized -> nrm=1, but keep general
    end do
    call stg6(gamma, cc, omg, ng)
    do I=1,nao; do J=1,nao; do K=1,nao; do L=1,nao
      Gao(I,J,K,L) = cg(I,J,K,L)
    end do; end do; end do; end do
    deallocate(ao_sh, ao_l, nrm)
    tmp = nn

  contains
    real(dp) function cg(I,J,K,L) result(v)
      integer, intent(in) :: I,J,K,L
      integer :: si,sj,sk,sl, pi,pj,pk,pl, g
      real(dp) :: ci,cj,ck,cl, ei,ej,ek,el, RI(3),RJ(3),RK(3),RL(3)
      si=ao_sh(I); sj=ao_sh(J); sk=ao_sh(K); sl=ao_sh(L)
      RI=shl_r(:,si); RJ=shl_r(:,sj); RK=shl_r(:,sk); RL=shl_r(:,sl)
      v = 0.0_dp
      ! physicist <IJ|f12|KL>: e1=(I,K), e2=(J,L) -> geminal_cart(I,K,J,L)
      do pi=1,shl_np(si); ei=shl_e(pi,si); ci=shl_c(pi,si)*cart_norm(ao_l(:,I),ei)
        do pj=1,shl_np(sj); ej=shl_e(pj,sj); cj=shl_c(pj,sj)*cart_norm(ao_l(:,J),ej)
          do pk=1,shl_np(sk); ek=shl_e(pk,sk); ck=shl_c(pk,sk)*cart_norm(ao_l(:,K),ek)
            do pl=1,shl_np(sl); el=shl_e(pl,sl); cl=shl_c(pl,sl)*cart_norm(ao_l(:,L),el)
              do g=1,ng
                v = v + cc(g)*ci*cj*ck*cl* &
                    geminal_cart(ao_l(:,I),RI,ei, ao_l(:,K),RK,ek, &
                                 ao_l(:,J),RJ,ej, ao_l(:,L),RL,el, omg(g))
              end do
            end do
          end do
        end do
      end do
      v = v*nrm(I)*nrm(J)*nrm(K)*nrm(L)
    end function cg
  end subroutine build_geminal_ao

  !> STG-6G fit of exp(-gamma r12) -> sum_k c_k exp(-omega_k r12^2).
  pure subroutine stg6(gamma, c, omg, ng)
    real(dp), intent(in)  :: gamma
    real(dp), intent(out) :: c(:), omg(:)
    integer,  intent(out) :: ng
    real(dp), parameter :: CC(6) = [0.3144_dp,0.3037_dp,0.1681_dp,0.09811_dp,0.06024_dp,0.03726_dp]
    real(dp), parameter :: EE(6) = [0.2209_dp,1.004_dp,3.622_dp,12.16_dp,45.87_dp,254.4_dp]
    integer :: k
    ng = 6
    do k = 1, 6
      c(k) = CC(k); omg(k) = EE(k)*gamma*gamma
    end do
  end subroutine stg6

end module ptc_ao
