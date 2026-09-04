! Slice-level BLAS-3 contractions shared by the MRSF-TDDFT semilocal XC
! Hessian integrators.
!
! Every quantity below is an exact algebraic rearrangement of the former
! per-grid-point, per-AO-pair accumulations.  For a slice of m quadrature
! points with AO value matrix Phi(n,m) and AO gradient matrices dPhi_c(n,m),
! a symmetrized density S = D + D^T gives
!
!   rho(p)        = 1/2 sum_mu Phi(mu,p) (S Phi)(mu,p),
!   grad rho_c(p) = sum_mu dPhi_c(mu,p) (S Phi)(mu,p),
!
! and the fixed-grid nuclear derivatives of rho and grad rho only need
! (S Phi) and (S dPhi_c) restricted to the AOs of each atom.  The nuclear
! derivative of an AO Fock matrix is assembled from per-point potentials
! a_k(p), b_k(p) (fixed AO pair) and v(p), c(p) (moving AO pair) as
!
!   dF_k = H_k + H_k^T,   H_k(mu,nu) = sum_p Phi(mu,p) X_k(nu,p) + moving terms,
!
! where X_k = a_k/2 Phi + b_k . dPhi.  The moving-AO and moving-grid terms
! of coordinate (A,a) touch only the AOs centred on atom A plus the owner
! atom of the slice, which is exploited by scattering the per-point moving
! potentials into the same stacked right-hand factor.  All sums over points
! are dgemm calls; only the symmetrization H + H^T is left to the caller.
!
! Work arrays are explicit-shape dummies so that one per-thread allocation
! sized for the full AO count can be reused with the pruned AO count of a
! slice through sequence association.
module mod_dft_gridint_mrsf_xc_slice_gemm

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_der1,xc_der2_contr,xc_der3_contr

  implicit none
  private

  integer, parameter, public :: hmap(3,3)=reshape([1,4,6, 4,2,5, 6,5,3],[3,3])

  public :: slice_stack_values
  public :: slice_stack_fixed
  public :: slice_stack_second
  public :: slice_fock_derivative_accumulate
  public :: slice_chunk_size
  public :: build_unweighted_kernels
  public :: symmetrize_half_accumulator
  public :: gather_stack
  public :: scatter_half_accumulator
  public :: build_symmetrized_stack

contains

!> Number of grid points processed per chunk so that the largest per-thread
!> work array (n x nwide x m doubles) stays near the requested byte budget.
  pure integer function slice_chunk_size(n,nwide,max_points,budget_bytes)
    integer, intent(in) :: n,nwide,max_points
    integer(8), intent(in) :: budget_bytes
    integer(8) :: per_point
    per_point=8_8*int(n,8)*int(max(1,nwide),8)
    slice_chunk_size=int(max(1_8,budget_bytes/max(1_8,per_point)))
    slice_chunk_size=max(8,min(slice_chunk_size,256))
    slice_chunk_size=max(1,min(slice_chunk_size,max_points))
  end function slice_chunk_size

!> stack(mu,f,nu)=scale*(D_f(mu,nu)+D_f(nu,mu))*norm(mu)*norm(nu) for the
!> f-th matrix of a batch stored as density(:,:,f).
  subroutine build_symmetrized_stack(nbf,nf,f,density,norm,scale,stack)
    integer, intent(in) :: nbf,nf,f
    real(fp), intent(in) :: density(nbf,nbf),norm(nbf),scale
    real(fp), intent(inout) :: stack(nbf,nf,nbf)
    integer :: mu,nu
    do nu=1,nbf
      do mu=1,nbf
        stack(mu,f,nu)=scale*(density(mu,nu)+density(nu,mu))*norm(mu)*norm(nu)
      end do
    end do
  end subroutine build_symmetrized_stack

!> Restrict a stack to the pruned AO list of a slice.
  subroutine gather_stack(nbf,n,nf,index,stack,stack_p)
    integer, intent(in) :: nbf,n,nf,index(n)
    real(fp), intent(in) :: stack(nbf,nf,nbf)
    real(fp), intent(out) :: stack_p(n,nf,n)
    integer :: mu,nu
    do nu=1,n
      do mu=1,n
        stack_p(mu,:,nu)=stack(index(mu),:,index(nu))
      end do
    end do
  end subroutine gather_stack

!> Add a pruned-slice half accumulator into the full-AO accumulator.
  subroutine scatter_half_accumulator(nbf,n,nwide,index,acc_p,acc)
    integer, intent(in) :: nbf,n,nwide,index(n)
    real(fp), intent(in) :: acc_p(n,n,nwide)
    real(fp), intent(inout) :: acc(nbf,nbf,nwide)
    integer :: k,mu,nu
    do k=1,nwide
      do nu=1,n
        do mu=1,n
          acc(index(mu),index(nu),k)=acc(index(mu),index(nu),k)+acc_p(mu,nu,k)
        end do
      end do
    end do
  end subroutine scatter_half_accumulator

!> Contract a stack of nf symmetrized density matrices with the AO values of
!> one chunk of grid points.
!>
!> stack(mu,f,nu)=S_f(mu,nu).  For points p0..p0+m-1:
!>   t(mu,f,p)=sum_nu S_f(mu,nu) Phi(nu,p)
!>   value(f,p)=1/2 sum_mu Phi(mu,p) t(mu,f,p)
!>   gradient(c,f,p)=sum_mu dPhi_c(mu,p) t(mu,f,p)
  subroutine slice_stack_values(n,m,nf,stack,aov,aog1,p0,t,value,gradient)
    integer, intent(in) :: n,m,nf,p0
    real(fp), intent(in) :: stack(n,nf,n)
    real(fp), intent(in), contiguous :: aov(:,:),aog1(:,:,:)
    real(fp), intent(inout) :: t(n,nf,m)
    real(fp), intent(out) :: value(nf,m),gradient(3,nf,m)

    integer :: c,f,lda,p,q

    if(n<=0 .or. m<=0 .or. nf<=0) return
    lda=size(aov,1)
    call dgemm('N','N',n*nf,m,n,1.0_fp,stack,n*nf,aov(1,p0),lda,0.0_fp, &
      t,n*nf)
    do p=1,m
      q=p0+p-1
      do f=1,nf
        value(f,p)=0.5_fp*dot_product(aov(1:n,q),t(1:n,f,p))
        do c=1,3
          gradient(c,f,p)=dot_product(aog1(1:n,q,c),t(1:n,f,p))
        end do
      end do
    end do
  end subroutine slice_stack_values

!> Fixed-grid first nuclear derivatives of rho_f and grad rho_f for the
!> same chunk, given t from slice_stack_values:
!>   u(mu,f,p,c)=sum_nu S_f(mu,nu) dPhi_c(nu,p)
!>   fixed_d(a,A,f,p)=-sum_{mu on A} dPhi_a(mu,p) t(mu,f,p)
!>   fixed_g(c,a,A,f,p)=-sum_{mu on A} [d2Phi_ac(mu,p) t(mu,f,p)
!>                                     +dPhi_a(mu,p) u(mu,f,p,c)]
  subroutine slice_stack_fixed(n,m,nf,nat,stack,aog1,aog2,p0,ao_atom,t,u, &
      fixed_d,fixed_g)
    integer, intent(in) :: n,m,nf,nat,p0,ao_atom(n)
    real(fp), intent(in) :: stack(n,nf,n)
    real(fp), intent(in), contiguous :: aog1(:,:,:),aog2(:,:,:)
    real(fp), intent(in) :: t(n,nf,m)
    real(fp), intent(inout) :: u(n,nf,m,3)
    real(fp), intent(out) :: fixed_d(3,nat,nf,m),fixed_g(3,3,nat,nf,m)

    integer :: atom,c,cart,f,ldg,mu,p,q,space
    real(fp) :: tv,g1(3)

    fixed_d=0.0_fp
    fixed_g=0.0_fp
    if(n<=0 .or. m<=0 .or. nf<=0) return
    ldg=size(aog1,1)
    do c=1,3
      call dgemm('N','N',n*nf,m,n,1.0_fp,stack,n*nf,aog1(1,p0,c),ldg, &
        0.0_fp,u(1,1,1,c),n*nf)
    end do
    do p=1,m
      q=p0+p-1
      do f=1,nf
        do mu=1,n
          atom=ao_atom(mu)
          tv=t(mu,f,p)
          g1=aog1(mu,q,:)
          do cart=1,3
            fixed_d(cart,atom,f,p)=fixed_d(cart,atom,f,p)-g1(cart)*tv
            do space=1,3
              fixed_g(space,cart,atom,f,p)=fixed_g(space,cart,atom,f,p)- &
                aog2(mu,q,hmap(cart,space))*tv-g1(cart)*u(mu,f,p,space)
            end do
          end do
        end do
      end do
    end do
  end subroutine slice_stack_fixed

!> Fixed-grid second nuclear derivatives of rho_f and grad rho_f for a
!> chunk, given t=S Phi and u_c=S dPhi_c from slice_stack_fixed.  With
!> U^{b,B}=S_{:,B} dPhi_b|_B and W^{bc,B}=S_{:,B} d2Phi_bc|_B restricted to
!> the AO columns of atom B,
!>   d2rho(a,b,A,B)=sum_{mu on A} dPhi_a(mu) U^{b,B}(mu)
!>                  +delta_AB sum_{mu on A} d2Phi_ab(mu) t(mu)
!>   d2grho(c,a,b,A,B)=sum_{mu on A} [d2Phi_ac(mu) U^{b,B}(mu)
!>                     +dPhi_a(mu) W^{bc,B}(mu)]
!>                  +delta_AB sum_{mu on A} [d3Phi_abc(mu) t(mu)
!>                     +d2Phi_ab(mu) u_c(mu)],
!> which is the ordered-pair form of gga_density_nuclear_point_batch.  The
!> AOs of every atom must be contiguous (atom_first/atom_last).
  subroutine slice_stack_second(n,m,nf,nat,stack,aog1,aog2,aog3,p0,ao_atom, &
      atom_first,atom_last,t,u,ub,wb,fixed_d2,fixed_g2)
    integer, intent(in) :: n,m,nf,nat,p0,ao_atom(n),atom_first(nat), &
      atom_last(nat)
    real(fp), intent(in) :: stack(n,nf,n)
    real(fp), intent(in), contiguous :: aog1(:,:,:),aog2(:,:,:),aog3(:,:,:)
    real(fp), intent(in) :: t(n,nf,m),u(n,nf,m,3)
    real(fp), intent(inout) :: ub(n,nf,m,3,nat),wb(n,nf,m,6,nat)
    real(fp), intent(out) :: fixed_d2(3,3,nat,nat,nf,m), &
      fixed_g2(3,3,3,nat,nat,nf,m)

    integer, parameter :: tmap(3,3,3)=reshape([ &
      1,4,5, 4,6,10, 5,10,8, &
      4,6,10, 6,2,7, 10,7,9, &
      5,10,8, 10,7,9, 8,9,3],[3,3,3])
    integer :: a,b,c,atom,batom,f,ldg,ldh,mu,nb,p,q
    real(fp) :: g1(3),g2(6),tv,uv(3),ubv(3),wbv(6)

    fixed_d2=0.0_fp
    fixed_g2=0.0_fp
    if(n<=0 .or. m<=0 .or. nf<=0) return
    ldg=size(aog1,1)
    ldh=size(aog2,1)
    do batom=1,nat
      nb=atom_last(batom)-atom_first(batom)+1
      if(nb<=0) then
        ub(:,:,:,:,batom)=0.0_fp
        wb(:,:,:,:,batom)=0.0_fp
        cycle
      end if
      do b=1,3
        call dgemm('N','N',n*nf,m,nb,1.0_fp,stack(1,1,atom_first(batom)), &
          n*nf,aog1(atom_first(batom),p0,b),ldg,0.0_fp,ub(1,1,1,b,batom),n*nf)
      end do
      do b=1,6
        call dgemm('N','N',n*nf,m,nb,1.0_fp,stack(1,1,atom_first(batom)), &
          n*nf,aog2(atom_first(batom),p0,b),ldh,0.0_fp,wb(1,1,1,b,batom),n*nf)
      end do
    end do
    do p=1,m
      q=p0+p-1
      do f=1,nf
        do mu=1,n
          atom=ao_atom(mu)
          tv=t(mu,f,p)
          uv=u(mu,f,p,:)
          g1=aog1(mu,q,:)
          g2=aog2(mu,q,:)
          ! same-AO terms (both derivatives on phi_mu)
          do b=1,3
            do a=1,3
              fixed_d2(a,b,atom,atom,f,p)=fixed_d2(a,b,atom,atom,f,p)+ &
                g2(hmap(a,b))*tv
              do c=1,3
                fixed_g2(c,a,b,atom,atom,f,p)=fixed_g2(c,a,b,atom,atom,f,p)+ &
                  aog3(mu,q,tmap(a,b,c))*tv+g2(hmap(a,b))*uv(c)
              end do
            end do
          end do
          ! cross terms (one derivative on phi_mu, one on phi_nu of atom B)
          do batom=1,nat
            ubv=ub(mu,f,p,:,batom)
            wbv=wb(mu,f,p,:,batom)
            do b=1,3
              do a=1,3
                fixed_d2(a,b,atom,batom,f,p)=fixed_d2(a,b,atom,batom,f,p)+ &
                  g1(a)*ubv(b)
                do c=1,3
                  fixed_g2(c,a,b,atom,batom,f,p)= &
                    fixed_g2(c,a,b,atom,batom,f,p)+ &
                    g2(hmap(a,c))*ubv(b)+g1(a)*wbv(hmap(b,c))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine slice_stack_second

!> Accumulate the unsymmetrized half H_k of the nuclear derivative of an AO
!> operator for every (coordinate k, probe, spin) from per-point potentials.
!>
!> a(s,k,pr,p), b(:,s,k,pr,p): fixed-AO scalar and gradient potentials
!>   (weights included), giving sum_p [a phi_mu phi_nu + b.grad(phi_mu phi_nu)].
!> v(s,pr,p), c(:,s,pr,p): unweighted potentials multiplying the derivative
!>   of the AO pair when the grid point moves with its owner atom and the AO
!>   centres move with their atoms; fw(p) is the finite quadrature weight.
!> acc(:,:,k,pr,s) receives H_k; the caller forms H_k+H_k^T at the end.
  subroutine slice_fock_derivative_accumulate(n,m,ncart,nprobe,aov,aog1, &
      aog2,p0,ao_atom,owner,fw,a,b,v,c,xs,z,psiw,g,acc)
    integer, intent(in) :: n,m,ncart,nprobe,p0,ao_atom(n),owner
    real(fp), intent(in), contiguous :: aov(:,:),aog1(:,:,:),aog2(:,:,:)
    real(fp), intent(in) :: fw(m),a(2,ncart,nprobe,m),b(3,2,ncart,nprobe,m)
    real(fp), intent(in) :: v(2,nprobe,m),c(3,2,nprobe,m)
    real(fp), intent(inout) :: xs(n,ncart,nprobe,m),z(n,m),psiw(n,3,m), &
      g(n,n,3)
    real(fp), intent(inout) :: acc(n,n,ncart,nprobe,2)

    integer :: cart,k,kc,kown,lda,nu,p,pr,q,spin
    real(fp) :: y1,cc(3),vv

    if(n<=0 .or. m<=0 .or. nprobe<=0) return
    lda=size(aov,1)
    do p=1,m
      q=p0+p-1
      do cart=1,3
        psiw(1:n,cart,p)=fw(p)*aog1(1:n,q,cart)
      end do
    end do
    do spin=1,2
      do p=1,m
        q=p0+p-1
        do pr=1,nprobe
          do k=1,ncart
            xs(1:n,k,pr,p)=(0.5_fp*a(spin,k,pr,p))*aov(1:n,q)+ &
              b(1,spin,k,pr,p)*aog1(1:n,q,1)+ &
              b(2,spin,k,pr,p)*aog1(1:n,q,2)+ &
              b(3,spin,k,pr,p)*aog1(1:n,q,3)
          end do
          vv=fw(p)*v(spin,pr,p)
          cc=fw(p)*c(:,spin,pr,p)
          do cart=1,3
            kown=3*(owner-1)+cart
            do nu=1,n
              y1=vv*aog1(nu,q,cart)+cc(1)*aog2(nu,q,hmap(cart,1))+ &
                cc(2)*aog2(nu,q,hmap(cart,2))+cc(3)*aog2(nu,q,hmap(cart,3))
              kc=3*(ao_atom(nu)-1)+cart
              xs(nu,kown,pr,p)=xs(nu,kown,pr,p)+y1
              xs(nu,kc,pr,p)=xs(nu,kc,pr,p)-y1
            end do
          end do
        end do
      end do
      call dgemm('N','T',n,n*ncart*nprobe,m,1.0_fp,aov(1,p0),lda,xs, &
        n*ncart*nprobe,1.0_fp,acc(1,1,1,1,spin),n)
      do pr=1,nprobe
        do p=1,m
          q=p0+p-1
          z(1:n,p)=c(1,spin,pr,p)*aog1(1:n,q,1)+c(2,spin,pr,p)*aog1(1:n,q,2)+ &
            c(3,spin,pr,p)*aog1(1:n,q,3)
        end do
        call dgemm('N','T',n,3*n,m,1.0_fp,z,n,psiw,3*n,0.0_fp,g,n)
        do cart=1,3
          kown=3*(owner-1)+cart
          acc(:,:,kown,pr,spin)=acc(:,:,kown,pr,spin)+g(:,:,cart)
          do nu=1,n
            kc=3*(ao_atom(nu)-1)+cart
            acc(:,nu,kc,pr,spin)=acc(:,nu,kc,pr,spin)-g(:,nu,cart)
          end do
        end do
      end do
    end do
  end subroutine slice_fock_derivative_accumulate

!> Replace every half accumulator H by H+H^T.
  subroutine symmetrize_half_accumulator(n,nwide,acc)
    integer, intent(in) :: n,nwide
    real(fp), intent(inout) :: acc(n,n,nwide)
    integer :: k,mu,nu
    real(fp) :: x
    do k=1,nwide
      do nu=1,n
        do mu=1,nu-1
          x=acc(mu,nu,k)+acc(nu,mu,k)
          acc(mu,nu,k)=x
          acc(nu,mu,k)=x
        end do
        acc(nu,nu,k)=2.0_fp*acc(nu,nu,k)
      end do
    end do
  end subroutine symmetrize_half_accumulator

!> Unweighted first, second, and third functional derivatives with respect
!> to (rho_a,rho_b,sigma_aa,sigma_bb,sigma_ab) at one grid point.
  subroutine build_unweighted_kernels(xce,ipt,nvar,finite_weight,first, &
      second,third)
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: ipt,nvar
    real(fp), intent(in) :: finite_weight
    real(fp), intent(out) :: first(:),second(:,:),third(:,:,:)
    real(fp) :: dr(2),ds(3),dt(2),fr(2),fs(3),ft(2),ffs(3)
    real(fp) :: gr(2),gs(3),gt(2)
    real(fp) :: square(5,5),mixed(5)
    integer :: j,k

    call xc_der1(xce,.true.,ipt,dr,ds,dt)
    call join_direction(nvar,dr,ds,first)
    first=first/finite_weight
    second=0.0_fp
    do j=1,nvar
      call split_basis(nvar,j,dr,ds)
      dt=0.0_fp
      call xc_der2_contr(xce,.true.,ipt,dr,ds,dt,fr,fs,ft)
      call join_direction(nvar,fr,fs,second(:,j))
    end do
    second=0.5_fp*(second+transpose(second))/finite_weight
    square=0.0_fp
    mixed=0.0_fp
    do j=1,nvar
      call split_basis(nvar,j,dr,ds)
      dt=0.0_fp
      call xc_der3_contr(xce,ipt,dr,ds,dt,[0.0_fp,0.0_fp,0.0_fp], &
        ffs,gr,gs,gt)
      call join_direction(nvar,gr,gs,square(1:nvar,j))
    end do
    third=0.0_fp
    do j=1,nvar
      third(:,j,j)=square(1:nvar,j)/finite_weight
      do k=j+1,nvar
        call split_pair(nvar,j,k,dr,ds)
        dt=0.0_fp
        call xc_der3_contr(xce,ipt,dr,ds,dt,[0.0_fp,0.0_fp,0.0_fp], &
          ffs,gr,gs,gt)
        call join_direction(nvar,gr,gs,mixed(1:nvar))
        mixed(1:nvar)=0.5_fp*(mixed(1:nvar)-square(1:nvar,j)- &
          square(1:nvar,k))/finite_weight
        third(:,j,k)=mixed(1:nvar)
        third(:,k,j)=mixed(1:nvar)
      end do
    end do
  end subroutine build_unweighted_kernels

  pure subroutine split_basis(nvar,index,dr,ds)
    integer, intent(in) :: nvar,index
    real(fp), intent(out) :: dr(2),ds(3)
    dr=0.0_fp
    ds=0.0_fp
    if(index<=2) then
      dr(index)=1.0_fp
    else if(nvar==5) then
      ds(index-2)=1.0_fp
    end if
  end subroutine split_basis

  pure subroutine split_pair(nvar,index_a,index_b,dr,ds)
    integer, intent(in) :: nvar,index_a,index_b
    real(fp), intent(out) :: dr(2),ds(3)
    real(fp) :: dra(2),drb(2),dsa(3),dsb(3)
    call split_basis(nvar,index_a,dra,dsa)
    call split_basis(nvar,index_b,drb,dsb)
    dr=dra+drb
    ds=dsa+dsb
  end subroutine split_pair

  pure subroutine join_direction(nvar,dr,ds,result)
    integer, intent(in) :: nvar
    real(fp), intent(in) :: dr(2),ds(3)
    real(fp), intent(out) :: result(:)
    result=0.0_fp
    result(1:2)=dr
    if(nvar==5) result(3:5)=ds
  end subroutine join_direction

end module mod_dft_gridint_mrsf_xc_slice_gemm
