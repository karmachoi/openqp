!> Fast 2-electron (Ms=0) full-CI engine with direct (alpha,beta) indexing, so
!> the FCI Hamiltonian, S^2, and the geminal transcorrelation cost O(norb^4)
!> instead of the O(norb^6) of the general determinant engine -- this is what lets
!> the pTC-MRSF-CIS demo use diffuse/large bases (Rydberg states S0..S3).
!>
!> Determinant |p,q> = alpha electron in MO p, beta in MO q; index = (p-1)*norb+q.
!> H[(p',q'),(p,q)] = d_{q'q} h(p',p) + d_{p'p} h(q',q) + (p'p|q'q)  [+ enuc on diag]
!> S^2[(p,q),(p,q)] = 1 (p/=q), 0 (p=q);  S^2[(p,q),(q,p)] = -1 (p/=q).
module ptc_fci2e
  use precision, only: dp
  implicit none
  private
  public :: idx2e, build_H2e, build_S2_2e, build_geminal_T2_2e, cas22_2e

contains

  pure integer function idx2e(p, q, norb) result(k)   ! 1-based MO p (alpha), q (beta)
    integer, intent(in) :: p, q, norb
    k = (p-1)*norb + q
  end function idx2e

  subroutine build_H2e(h1, eri_c, enuc, norb, H)
    integer,  intent(in)  :: norb
    real(dp), intent(in)  :: h1(norb,norb), eri_c(norb,norb,norb,norb), enuc
    real(dp), intent(out) :: H(norb*norb, norb*norb)
    integer  :: p, q, pp, qq, r, c
    H = 0.0_dp
    !$omp parallel do collapse(2) default(shared) private(p,q,pp,qq,r,c)
    do p = 1, norb
      do q = 1, norb
        c = idx2e(p,q,norb)
        do pp = 1, norb
          do qq = 1, norb
            r = idx2e(pp,qq,norb)
            H(r,c) = eri_c(pp,p,qq,q)                       ! (p'p|q'q) chemist
            if (qq == q) H(r,c) = H(r,c) + h1(pp,p)         ! alpha 1-body
            if (pp == p) H(r,c) = H(r,c) + h1(qq,q)         ! beta  1-body
          end do
        end do
        H(c,c) = H(c,c) + enuc
      end do
    end do
  end subroutine build_H2e

  subroutine build_S2_2e(norb, S2)
    integer,  intent(in)  :: norb
    real(dp), intent(out) :: S2(norb*norb, norb*norb)
    integer  :: p, q
    S2 = 0.0_dp
    do p = 1, norb
      do q = 1, norb
        if (p /= q) then
          S2(idx2e(p,q,norb), idx2e(p,q,norb)) = 1.0_dp
          S2(idx2e(q,p,norb), idx2e(p,q,norb)) = -1.0_dp
        end if
      end do
    end do
  end subroutine build_S2_2e

  !> (2,2) frontier compact: determinants with both MOs in {1,2}.
  subroutine cas22_2e(norb, cas, nc)
    integer,              intent(in)  :: norb
    integer, allocatable, intent(out) :: cas(:)
    integer,              intent(out) :: nc
    integer :: p, q, cnt
    nc = 4; allocate(cas(4)); cnt = 0
    do p = 1, 2
      do q = 1, 2
        cnt = cnt + 1; cas(cnt) = idx2e(p,q,norb)
      end do
    end do
  end subroutine cas22_2e

  !> Geminal transcorrelation operator T (doubles) in the (alpha,beta) basis.
  !> Opposite-spin double |i,j> -> |a,b> (alpha i->a, beta j->b), amplitude
  !> -1/(2 gamma) * Gmo(i,j,a,b) [the 1/2 opposite-spin cusp; same-spin doubles
  !> do not exist for a 1-alpha/1-beta determinant]. i,j in iact, a,b in iext.
  subroutine build_geminal_T2_2e(Gmo, norb, iact, nact, iext, next, gamma, T2)
    integer,  intent(in)  :: norb, iact(:), nact, iext(:), next
    real(dp), intent(in)  :: Gmo(norb,norb,norb,norb), gamma
    real(dp), intent(out) :: T2(norb*norb, norb*norb)
    integer  :: ii, jj, aa, bb, i, j, a, b
    real(dp) :: amp
    T2 = 0.0_dp
    do ii = 1, nact
      do jj = 1, nact
        do aa = 1, next
          do bb = 1, next
            i = iact(ii); j = iact(jj); a = iext(aa); b = iext(bb)
            amp = -0.5_dp/gamma * Gmo(i,j,a,b)
            ! two opposite-spin paths (matching the general spin-orbital engine):
            !  (alpha i->a, beta j->b): |i,j> -> |a,b>
            T2(idx2e(a,b,norb), idx2e(i,j,norb)) = T2(idx2e(a,b,norb), idx2e(i,j,norb)) + amp
            !  (alpha j->b, beta i->a): |j,i> -> |b,a>
            T2(idx2e(b,a,norb), idx2e(j,i,norb)) = T2(idx2e(b,a,norb), idx2e(j,i,norb)) + amp
          end do
        end do
      end do
    end do
  end subroutine build_geminal_T2_2e

end module ptc_fci2e
