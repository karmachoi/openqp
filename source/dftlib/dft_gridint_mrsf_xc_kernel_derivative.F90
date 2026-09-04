! Analytic nuclear derivative of the spin-polarized semilocal XC-kernel
! action needed by the differentiated MRSF orbital-adjoint equations.
!
! The quadrature is organised per grid slice: the density, probe, and
! response fields of every point are obtained from stacked dgemm calls, the
! per-point potentials of all (coordinate, probe, spin) triples are formed
! from the first three functional derivatives, and the AO-matrix
! contributions of the slice are accumulated by one stacked dgemm per spin
! (mod_dft_gridint_mrsf_xc_slice_gemm).  The integrand is unchanged from the
! former per-AO-pair implementation; only the summation order differs.
module mod_dft_gridint_mrsf_xc_kernel_derivative

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc
  use mod_dft_partition_hessian, only: &
    partition_weight_nuclear_first_derivatives
  use mod_dft_gridint_tdgga_response, only: gga_add_owner_motion_first
  use mod_dft_gridint_mrsf_xc_hessian_point, only: &
    mrsf_xc_kernel_fock_coefficients
  use mod_dft_gridint_mrsf_xc_slice_gemm, only: slice_stack_values, &
    slice_stack_fixed,slice_fock_derivative_accumulate,slice_chunk_size, &
    build_unweighted_kernels,symmetrize_half_accumulator,gather_stack, &
    scatter_half_accumulator,build_symmetrized_stack

  implicit none
  private

  integer(8), parameter :: chunk_budget_bytes=24_8*1024_8*1024_8

  type :: kernel_slice_workspace_t
    ! chunk-level fields
    real(fp), allocatable :: t(:,:,:),u(:,:,:,:),value(:,:),gradient(:,:,:)
    real(fp), allocatable :: fixed_d(:,:,:,:),fixed_g(:,:,:,:,:)
    real(fp), allocatable :: rt(:,:,:),rvalue(:,:),rgradient(:,:,:)
    ! chunk-level potentials and accumulation scratch
    real(fp), allocatable :: a(:,:,:,:),b(:,:,:,:,:),v(:,:,:),c(:,:,:,:),fw(:)
    real(fp), allocatable :: xs(:,:,:,:),z(:,:),psiw(:,:,:),g(:,:,:)
    ! pruned-slice gathers (allocated on first use)
    real(fp), allocatable :: stack_p(:,:,:),rstack_p(:,:,:),acc_p(:,:,:,:,:)
    integer, allocatable :: ao_atom_p(:)
    ! point-level chain-rule scratch
    real(fp), allocatable :: dg_rho(:,:),dp_rho(:,:)
    real(fp), allocatable :: dg_grad(:,:,:),dp_grad(:,:,:)
    real(fp), allocatable :: total_d(:,:,:,:),total_g(:,:,:,:,:)
    real(fp), allocatable :: weights(:),dweights(:,:,:)
    real(fp), allocatable :: ground(:),probe(:),dground(:,:),dprobe(:,:)
    real(fp), allocatable :: first(:),second(:,:),third(:,:,:)
    real(fp), allocatable :: v_r(:),dv_r(:,:),dweight_flat(:)
    real(fp), allocatable :: coefficient(:,:),dcoefficient(:,:,:)
    integer :: mchunk=0
  contains
    procedure :: init => kernel_workspace_init
    procedure :: clean => kernel_workspace_clean
  end type kernel_slice_workspace_t

  type, extends(xc_consumer_t) :: kernel_derivative_consumer_t
    !> Half accumulators H (mu,nu,coordinate,probe,spin,thread).
    real(fp), allocatable :: derivative(:,:,:,:,:,:)
    !> Symmetrized reference and probe densities (mu,f,nu),
    !> f=spin+2*(field-1), field 1 = reference, field 1+p = probe p.
    real(fp), allocatable :: stack(:,:,:)
    !> Symmetrized nuclear responses (mu,j,nu),
    !> j=k+ncart*(spin-1)+2*ncart*(field-1).
    real(fp), allocatable :: rstack(:,:,:)
    type(kernel_slice_workspace_t), allocatable :: workspace(:)
    real(fp), allocatable :: atom_xyz(:,:),surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    real(fp), allocatable :: worker_error(:)
    integer :: part_fun_type=0
    integer :: nprobe=0
    logical :: has_surface_shift=.false.,is_gga=.false.
  contains
    procedure :: parallel_start => kernel_start
    procedure :: parallel_stop => kernel_stop
    procedure :: update => kernel_update
    procedure :: postUpdate => kernel_post
    procedure :: clean => kernel_clean
  end type kernel_derivative_consumer_t

  public :: mrsf_xc_kernel_fock_total_derivative
  public :: mrsf_xc_kernel_fock_total_derivative_batch

contains

!-----------------------------------------------------------------------------

  subroutine kernel_workspace_init(self,nbf,nat,nprobe,is_gga,max_points)
    class(kernel_slice_workspace_t), intent(inout) :: self
    integer, intent(in) :: nbf,nat,nprobe,max_points
    logical, intent(in) :: is_gga
    integer :: m,ncart,nf,nj,nvar

    call self%clean()
    ncart=3*nat
    nf=2*(nprobe+1)
    nj=ncart*nf
    nvar=merge(5,2,is_gga)
    m=slice_chunk_size(nbf,max(nj,ncart*nprobe,3*nf),max_points, &
      chunk_budget_bytes)
    self%mchunk=m
    allocate(self%t(nbf,nf,m),self%u(nbf,nf,m,3),self%value(nf,m), &
      self%gradient(3,nf,m),self%fixed_d(3,nat,nf,m), &
      self%fixed_g(3,3,nat,nf,m),self%rt(nbf,nj,m),self%rvalue(nj,m), &
      self%rgradient(3,nj,m),self%a(2,ncart,nprobe,m), &
      self%b(3,2,ncart,nprobe,m),self%v(2,nprobe,m),self%c(3,2,nprobe,m), &
      self%fw(m),self%xs(nbf,ncart,nprobe,m),self%z(nbf,m), &
      self%psiw(nbf,3,m),self%g(nbf,nbf,3))
    allocate(self%dg_rho(2,ncart),self%dp_rho(2,ncart), &
      self%dg_grad(3,2,ncart),self%dp_grad(3,2,ncart), &
      self%total_d(3,nat,2,nprobe+1),self%total_g(3,3,nat,2,nprobe+1), &
      self%weights(nat),self%dweights(3,nat,nat),self%ground(nvar), &
      self%probe(nvar),self%dground(nvar,ncart),self%dprobe(nvar,ncart), &
      self%first(nvar),self%second(nvar,nvar),self%third(nvar,nvar,nvar), &
      self%v_r(2),self%dv_r(2,ncart),self%dweight_flat(ncart), &
      self%coefficient(3,2),self%dcoefficient(3,2,ncart))
  end subroutine kernel_workspace_init

  subroutine kernel_workspace_clean(self)
    class(kernel_slice_workspace_t), intent(inout) :: self
    if(allocated(self%t)) deallocate(self%t,self%u,self%value, &
      self%gradient,self%fixed_d,self%fixed_g,self%rt,self%rvalue, &
      self%rgradient,self%a,self%b,self%v,self%c,self%fw,self%xs,self%z, &
      self%psiw,self%g)
    if(allocated(self%stack_p)) deallocate(self%stack_p)
    if(allocated(self%rstack_p)) deallocate(self%rstack_p)
    if(allocated(self%acc_p)) deallocate(self%acc_p)
    if(allocated(self%ao_atom_p)) deallocate(self%ao_atom_p)
    if(allocated(self%dg_rho)) deallocate(self%dg_rho,self%dp_rho, &
      self%dg_grad,self%dp_grad,self%total_d,self%total_g,self%weights, &
      self%dweights,self%ground,self%probe,self%dground,self%dprobe, &
      self%first,self%second,self%third,self%v_r,self%dv_r, &
      self%dweight_flat,self%coefficient,self%dcoefficient)
    self%mchunk=0
  end subroutine kernel_workspace_clean

!> Differentiate K_xc[D](Q) for physical alpha/beta reference and probe
!> densities D and Q.  dD and dQ are their AO-coefficient responses.  Moving
!> AO, atom-centred grid, partition-weight, fxc, and kxc terms are included.
!> No displaced geometry or electronic-state expansion is used.
  subroutine mrsf_xc_kernel_fock_total_derivative(basis,mol_grid,da,db, &
      qa,qb,dda,ddb,dqa,dqb,derivative_a,derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: da(:,:),db(:,:),qa(:,:),qb(:,:)
    real(fp), intent(in) :: dda(:,:,:),ddb(:,:,:),dqa(:,:,:),dqb(:,:,:)
    real(fp), intent(out) :: derivative_a(:,:,:),derivative_b(:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    real(fp), allocatable :: qa_batch(:,:,:),qb_batch(:,:,:)
    real(fp), allocatable :: dqa_batch(:,:,:,:),dqb_batch(:,:,:,:)
    real(fp), allocatable :: derivative_a_batch(:,:,:,:), &
      derivative_b_batch(:,:,:,:)
    integer :: nbf,ncart

    nbf=basis%nbf
    ncart=3*infos%mol_prop%natom
    allocate(qa_batch(nbf,nbf,1),qb_batch(nbf,nbf,1), &
      dqa_batch(nbf,nbf,ncart,1),dqb_batch(nbf,nbf,ncart,1), &
      derivative_a_batch(nbf,nbf,ncart,1), &
      derivative_b_batch(nbf,nbf,ncart,1),source=0.0_fp)
    qa_batch(:,:,1)=qa; qb_batch(:,:,1)=qb
    dqa_batch(:,:,:,1)=dqa; dqb_batch(:,:,:,1)=dqb
    call mrsf_xc_kernel_fock_total_derivative_batch(basis,mol_grid,da,db, &
      qa_batch,qb_batch,dda,ddb,dqa_batch,dqb_batch,derivative_a_batch, &
      derivative_b_batch,infos,status,threshold)
    derivative_a=derivative_a_batch(:,:,:,1)
    derivative_b=derivative_b_batch(:,:,:,1)
    deallocate(qa_batch,qb_batch,dqa_batch,dqb_batch, &
      derivative_a_batch,derivative_b_batch)
  end subroutine mrsf_xc_kernel_fock_total_derivative

!> Exact multi-probe derivative of K_xc[D](Q_p).  The physical reference D
!> and its nuclear response are common to every probe.  All probes share one
!> molecular-grid traversal and the same AO values; no grid, density, or
!> functional approximation is introduced.
  subroutine mrsf_xc_kernel_fock_total_derivative_batch(basis,mol_grid,da,db, &
      qa,qb,dda,ddb,dqa,dqb,derivative_a,derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: da(:,:),db(:,:),qa(:,:,:),qb(:,:,:)
    real(fp), intent(in) :: dda(:,:,:),ddb(:,:,:),dqa(:,:,:,:),dqb(:,:,:,:)
    real(fp), intent(out) :: derivative_a(:,:,:,:),derivative_b(:,:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    type(kernel_derivative_consumer_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: da_normalized(:,:),db_normalized(:,:)
    real(fp) :: grid_threshold
    integer :: first,i,j,k,last,natom,nbf,ncart,nf,nj,nprobe,probe,shell

    nbf=basis%nbf
    natom=infos%mol_prop%natom
    ncart=3*natom
    nprobe=size(qa,3)
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. natom<=0 .or. nprobe<=0 .or. &
       any(shape(da)/=[nbf,nbf]) .or. any(shape(db)/=[nbf,nbf]) .or. &
       any(shape(qa)/=[nbf,nbf,nprobe]) .or. &
       any(shape(qb)/=[nbf,nbf,nprobe]) .or. &
       any(shape(dda)/=[nbf,nbf,ncart]) .or. &
       any(shape(ddb)/=[nbf,nbf,ncart]) .or. &
       any(shape(dqa)/=[nbf,nbf,ncart,nprobe]) .or. &
       any(shape(dqb)/=[nbf,nbf,ncart,nprobe]) .or. &
       any(shape(derivative_a)/=[nbf,nbf,ncart,nprobe]) .or. &
       any(shape(derivative_b)/=[nbf,nbf,ncart,nprobe])) then
      status=-1
      return
    end if
    if(infos%functional%needTau .or. infos%functional%needLapl .or. &
       infos%dft%cam_flag) then
      status=-2
      return
    end if
    grid_threshold=infos%dft%grid_density_cutoff
    if(present(threshold)) grid_threshold=threshold
    if(grid_threshold<0.0_fp) then
      status=-3
      return
    end if

    nf=2*(nprobe+1)
    nj=ncart*nf
    allocate(da_normalized(nbf,nbf),db_normalized(nbf,nbf), &
      dat%stack(nbf,nf,nbf),dat%rstack(nbf,nj,nbf),dat%ao_atom(nbf))
    do i=1,nbf
      da_normalized(:,i)=da(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      db_normalized(:,i)=db(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
    end do
    call build_symmetrized_stack(nbf,nf,1,da,basis%bfnrm,1.0_fp,dat%stack)
    call build_symmetrized_stack(nbf,nf,2,db,basis%bfnrm,1.0_fp,dat%stack)
    do probe=1,nprobe
      call build_symmetrized_stack(nbf,nf,2*probe+1,qa(:,:,probe), &
        basis%bfnrm,1.0_fp,dat%stack)
      call build_symmetrized_stack(nbf,nf,2*probe+2,qb(:,:,probe), &
        basis%bfnrm,1.0_fp,dat%stack)
    end do
    do k=1,ncart
      call build_symmetrized_stack(nbf,nj,k,dda(:,:,k),basis%bfnrm,1.0_fp, &
        dat%rstack)
      call build_symmetrized_stack(nbf,nj,k+ncart,ddb(:,:,k),basis%bfnrm, &
        1.0_fp,dat%rstack)
      do probe=1,nprobe
        j=k+2*ncart*probe
        call build_symmetrized_stack(nbf,nj,j,dqa(:,:,k,probe), &
          basis%bfnrm,1.0_fp,dat%rstack)
        call build_symmetrized_stack(nbf,nj,j+ncart,dqb(:,:,k,probe), &
          basis%bfnrm,1.0_fp,dat%rstack)
      end do
    end do
    dat%ao_atom=0
    do shell=1,basis%nshell
      first=basis%ao_offset(shell)
      last=first+basis%naos(shell)-1
      dat%ao_atom(first:last)=basis%origin(shell)
    end do
    if(any(dat%ao_atom<1) .or. any(dat%ao_atom>natom)) then
      status=-4
      call dat%clean()
      deallocate(da_normalized,db_normalized)
      return
    end if
    dat%atom_xyz=infos%atoms%xyz
    dat%dummy_atom=mol_grid%dummyAtom
    dat%surface_shift=mol_grid%surfaceShift
    dat%part_fun_type=mol_grid%partFunType
    dat%nprobe=nprobe
    dat%has_surface_shift=mol_grid%hasSurfaceShift
    dat%is_gga=infos%functional%needGrd

    opts%isGGA=dat%is_gga
    opts%needTau=.false.
    opts%functional=>infos%functional
    opts%hasBeta=.true.
    opts%isWFVecs=.false.
    opts%numAOs=nbf
    opts%maxPts=mol_grid%maxSlicePts
    opts%limPts=mol_grid%maxNRadTimesNAng
    opts%numAtoms=natom
    opts%maxAngMom=basis%mxam
    opts%nDer=merge(1,2,dat%is_gga)
    opts%nXCDer=3
    opts%numOccAlpha=infos%mol_prop%nelec_A
    opts%numOccBeta=infos%mol_prop%nelec_B
    opts%wfAlpha=>da_normalized
    opts%wfBeta=>db_normalized
    opts%dft_threshold=grid_threshold
    ! AOs of prescreened-out shells are zeroed on every slice, but the
    ! compressed-AO slice layout is not used: with it the butadiene/SG-2
    ! Hessian deviated from the uncompressed result by 7e-6 Eh/bohr^2
    ! (2026-09-04), so the derivative quadratures keep the full AO layout.
    opts%ao_threshold=infos%dft%grid_ao_threshold
    opts%ao_sparsity_ratio=0.0_fp
    opts%molGrid=>mol_grid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(opts,dat,basis)
    if(allocated(dat%worker_error)) then
      if(dat%worker_error(1)>0.0_fp) status=-5
    end if
    if(status==0) then
      do probe=1,nprobe
        do k=1,ncart
          do i=1,nbf
            derivative_a(:,i,k,probe)=dat%derivative(:,i,k,probe,1,1)* &
              basis%bfnrm(:)*basis%bfnrm(i)
            derivative_b(:,i,k,probe)=dat%derivative(:,i,k,probe,2,1)* &
              basis%bfnrm(:)*basis%bfnrm(i)
          end do
        end do
      end do
    end if
    call dat%clean()
    deallocate(da_normalized,db_normalized)
  end subroutine mrsf_xc_kernel_fock_total_derivative_batch

!-----------------------------------------------------------------------------

  subroutine kernel_start(self,xce,nthreads)
    class(kernel_derivative_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: ncart,thread
    ncart=3*xce%numAtoms
    allocate(self%derivative(xce%numAOs,xce%numAOs,ncart,self%nprobe,2, &
      nthreads),self%worker_error(nthreads),source=0.0_fp)
    allocate(self%workspace(nthreads))
    do thread=1,nthreads
      call self%workspace(thread)%init(xce%numAOs,xce%numAtoms, &
        self%nprobe,self%is_gga,max(1,xce%maxPts))
    end do
  end subroutine kernel_start

  subroutine kernel_stop(self)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    integer :: nwide,spin
    if(size(self%derivative,6)>1) then
      self%derivative(:,:,:,:,:,1)=sum(self%derivative,dim=6)
      self%worker_error(1)=sum(self%worker_error)
    end if
    nwide=size(self%derivative,3)*size(self%derivative,4)*2
    call symmetrize_half_accumulator(size(self%derivative,1),nwide, &
      self%derivative(:,:,:,:,:,1))
    do spin=1,2
      call self%pe%allreduce(self%derivative(:,:,:,:,spin,1), &
        size(self%derivative(:,:,:,:,spin,1)))
    end do
    call self%pe%allreduce(self%worker_error(1),1)
  end subroutine kernel_stop

  subroutine kernel_post(self,xce,mythread)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
  end subroutine kernel_post

  subroutine kernel_clean(self)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    integer :: thread
    if(allocated(self%derivative)) deallocate(self%derivative)
    if(allocated(self%stack)) deallocate(self%stack)
    if(allocated(self%rstack)) deallocate(self%rstack)
    if(allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if(allocated(self%surface_shift)) deallocate(self%surface_shift)
    if(allocated(self%ao_atom)) deallocate(self%ao_atom)
    if(allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    if(allocated(self%worker_error)) deallocate(self%worker_error)
    if(allocated(self%workspace)) then
      do thread=1,size(self%workspace)
        call self%workspace(thread)%clean()
      end do
      deallocate(self%workspace)
    end if
  end subroutine kernel_clean

!-----------------------------------------------------------------------------

  subroutine kernel_update(self,xce,mythread)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    integer :: n,nbf,ncart,nf,nj,nwide,status

    nbf=xce%numAOs
    n=xce%numAOs_p
    ncart=3*xce%numAtoms
    nf=2*(self%nprobe+1)
    nj=ncart*nf
    nwide=ncart*self%nprobe*2
    if(n<=0 .or. xce%numPts<=0) return
    associate(ws=>self%workspace(mythread))
    if(xce%skip_p) then
      call kernel_slice(self,xce,mythread,n,self%stack,self%rstack, &
        self%ao_atom,self%derivative(:,:,:,:,:,mythread),status)
    else
      if(.not.allocated(ws%stack_p)) then
        allocate(ws%stack_p(nbf,nf,nbf),ws%rstack_p(nbf,nj,nbf), &
          ws%acc_p(nbf,nbf,ncart,self%nprobe,2),ws%ao_atom_p(nbf))
      end if
      call gather_stack(nbf,n,nf,xce%indices_p(1:n),self%stack,ws%stack_p)
      call gather_stack(nbf,n,nj,xce%indices_p(1:n),self%rstack,ws%rstack_p)
      ws%ao_atom_p(1:n)=self%ao_atom(xce%indices_p(1:n))
      ws%acc_p=0.0_fp
      call kernel_slice(self,xce,mythread,n,ws%stack_p,ws%rstack_p, &
        ws%ao_atom_p,ws%acc_p,status)
      call scatter_half_accumulator(nbf,n,nwide,xce%indices_p(1:n), &
        ws%acc_p,self%derivative(:,:,:,:,:,mythread))
    end if
    end associate
    if(status/=0) self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
  end subroutine kernel_update

!> Process one slice with the pruned AO count n.  stack/rstack/ao_atom/acc
!> are already restricted to the slice's AO list.
  subroutine kernel_slice(self,xce,mythread,n,stack,rstack,ao_atom,acc,status)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread,n
    real(fp), intent(in) :: stack(*),rstack(*)
    integer, intent(in) :: ao_atom(n)
    real(fp), intent(inout) :: acc(*)
    integer, intent(out) :: status

    integer :: m,nat,ncart,nf,nj,npts,p,p0,local_status

    nat=xce%numAtoms
    ncart=3*nat
    nf=2*(self%nprobe+1)
    nj=ncart*nf
    npts=xce%numPts
    status=0
    associate(ws=>self%workspace(mythread))
    do p0=1,npts,ws%mchunk
      m=min(ws%mchunk,npts-p0+1)
      call slice_stack_values(n,m,nf,stack,xce%aoV,xce%aoG1,p0,ws%t, &
        ws%value,ws%gradient)
      call slice_stack_fixed(n,m,nf,nat,stack,xce%aoG1,xce%aoG2,p0,ao_atom, &
        ws%t,ws%u,ws%fixed_d,ws%fixed_g)
      call slice_stack_values(n,m,nj,rstack,xce%aoV,xce%aoG1,p0,ws%rt, &
        ws%rvalue,ws%rgradient)
      do p=1,m
        call kernel_point_potentials(self,xce,p0+p-1,p,ws,local_status)
        if(local_status/=0) then
          status=local_status
          return
        end if
      end do
      call slice_fock_derivative_accumulate(n,m,ncart,self%nprobe,xce%aoV, &
        xce%aoG1,xce%aoG2,p0,ao_atom,xce%currAtom,ws%fw,ws%a,ws%b,ws%v, &
        ws%c,ws%xs,ws%z,ws%psiw,ws%g,acc)
    end do
    end associate
  end subroutine kernel_slice

!-----------------------------------------------------------------------------

!> Per-point potentials of every (coordinate, probe, spin) triple.  The
!> chain rule is unchanged: total density derivatives (fixed grid plus owner
!> motion plus AO-coefficient response) enter the second and third
!> functional derivatives through mrsf_xc_kernel_fock_coefficients, and the
!> partition-weight derivative multiplies the undifferentiated potential.
  subroutine kernel_point_potentials(self,xce,ipt,p,ws,status)
    class(kernel_derivative_consumer_t), intent(in) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: ipt,p
    type(kernel_slice_workspace_t), intent(inout) :: ws
    integer, intent(out) :: status

    real(fp) :: rho(2),prho(2),grad_rho(3,2),grad_probe(3,2)
    real(fp) :: finite_weight,quadrature_scale,scale
    integer :: atom,cart,coordinate,field,jg,jp,k,local_status
    integer :: nat,ncart,nprobe,nvar,owner,probe_index,spin

    nat=xce%numAtoms
    ncart=3*nat
    nprobe=self%nprobe
    nvar=merge(5,2,self%is_gga)
    owner=xce%currAtom
    status=0
    finite_weight=xce%xyzw(ipt,4)
    ws%fw(p)=finite_weight
    ws%a(:,:,:,p)=0.0_fp
    ws%b(:,:,:,:,p)=0.0_fp
    ws%v(:,:,p)=0.0_fp
    ws%c(:,:,:,p)=0.0_fp
    if(abs(finite_weight)<=tiny(1.0_fp)) then
      ws%fw(p)=0.0_fp
      return
    end if
    call partition_weight_nuclear_first_derivatives(self%atom_xyz, &
      xce%xyzw(ipt,1:3),owner,self%dummy_atom,self%part_fun_type, &
      self%has_surface_shift,self%surface_shift,ws%weights,ws%dweights, &
      local_status)
    if(local_status/=0 .or. ws%weights(owner)<=sqrt(tiny(1.0_fp))) then
      status=-2
      return
    end if
    quadrature_scale=finite_weight/ws%weights(owner)
    do atom=1,nat
      do cart=1,3
        coordinate=3*(atom-1)+cart
        ws%dweight_flat(coordinate)=ws%dweights(cart,atom,owner)
      end do
    end do

    call build_unweighted_kernels(xce,ipt,nvar,finite_weight,ws%first, &
      ws%second,ws%third)
    do field=1,nprobe+1
      do spin=1,2
        call gga_add_owner_motion_first(owner, &
          ws%fixed_d(:,:,spin+2*(field-1),p), &
          ws%fixed_g(:,:,:,spin+2*(field-1),p), &
          ws%total_d(:,:,spin,field),ws%total_g(:,:,:,spin,field))
      end do
    end do
    rho=ws%value(1:2,p)
    grad_rho=ws%gradient(:,1:2,p)
    do spin=1,2
      do coordinate=1,ncart
        atom=(coordinate-1)/3+1
        cart=coordinate-3*(atom-1)
        jg=coordinate+ncart*(spin-1)
        ws%dg_rho(spin,coordinate)=ws%total_d(cart,atom,spin,1)+ &
          ws%rvalue(jg,p)
        ws%dg_grad(:,spin,coordinate)=ws%total_g(:,cart,atom,spin,1)+ &
          ws%rgradient(:,jg,p)
      end do
    end do
    do probe_index=1,nprobe
      prho=ws%value(2*probe_index+1:2*probe_index+2,p)
      grad_probe=ws%gradient(:,2*probe_index+1:2*probe_index+2,p)
      do spin=1,2
        do coordinate=1,ncart
          atom=(coordinate-1)/3+1
          cart=coordinate-3*(atom-1)
          jp=coordinate+ncart*(spin-1)+2*ncart*probe_index
          ws%dp_rho(spin,coordinate)= &
            ws%total_d(cart,atom,spin,probe_index+1)+ws%rvalue(jp,p)
          ws%dp_grad(:,spin,coordinate)= &
            ws%total_g(:,cart,atom,spin,probe_index+1)+ws%rgradient(:,jp,p)
        end do
      end do
      call build_density_variables(self%is_gga,rho,grad_rho,prho,grad_probe, &
        ws%dg_rho,ws%dg_grad,ws%dp_rho,ws%dp_grad,ws%ground,ws%probe, &
        ws%dground,ws%dprobe)
      call mrsf_xc_kernel_fock_coefficients(ws%first,ws%second,ws%third, &
        ws%probe,ws%dground,ws%dprobe,self%is_gga,grad_rho,grad_probe, &
        ws%dg_grad,ws%dp_grad,ws%v_r,ws%coefficient,ws%dv_r, &
        ws%dcoefficient,local_status)
      if(local_status/=0) then
        status=-5
        return
      end if
      do spin=1,2
        do k=1,ncart
          scale=quadrature_scale*ws%dweight_flat(k)
          ws%a(spin,k,probe_index,p)=scale*ws%v_r(spin)+ &
            finite_weight*ws%dv_r(spin,k)
          ws%b(:,spin,k,probe_index,p)=scale*ws%coefficient(:,spin)+ &
            finite_weight*ws%dcoefficient(:,spin,k)
        end do
        ws%v(spin,probe_index,p)=ws%v_r(spin)
        ws%c(:,spin,probe_index,p)=ws%coefficient(:,spin)
      end do
    end do
  end subroutine kernel_point_potentials

!-----------------------------------------------------------------------------

  pure subroutine build_density_variables(is_gga,rho,grho,prho,pgrho, &
      drho,dgrho,dprho,dpgrho,ground,probe,dground,dprobe)
    logical, intent(in) :: is_gga
    real(fp), intent(in) :: rho(2),grho(3,2),prho(2),pgrho(3,2)
    real(fp), intent(in) :: drho(:,:),dgrho(:,:,:),dprho(:,:),dpgrho(:,:,:)
    real(fp), intent(out) :: ground(:),probe(:),dground(:,:),dprobe(:,:)
    integer :: coordinate

    ground=0.0_fp
    probe=0.0_fp
    dground=0.0_fp
    dprobe=0.0_fp
    ground(1:2)=rho
    probe(1:2)=prho
    dground(1:2,:)=drho
    dprobe(1:2,:)=dprho
    if(.not.is_gga) return
    ground(3)=dot_product(grho(:,1),grho(:,1))
    ground(4)=dot_product(grho(:,2),grho(:,2))
    ground(5)=dot_product(grho(:,1),grho(:,2))
    probe(3)=2.0_fp*dot_product(grho(:,1),pgrho(:,1))
    probe(4)=2.0_fp*dot_product(grho(:,2),pgrho(:,2))
    probe(5)=dot_product(grho(:,1),pgrho(:,2))+ &
      dot_product(grho(:,2),pgrho(:,1))
    do coordinate=1,size(drho,2)
      dground(3,coordinate)=2.0_fp*dot_product(grho(:,1), &
        dgrho(:,1,coordinate))
      dground(4,coordinate)=2.0_fp*dot_product(grho(:,2), &
        dgrho(:,2,coordinate))
      dground(5,coordinate)=dot_product(dgrho(:,1,coordinate),grho(:,2))+ &
        dot_product(grho(:,1),dgrho(:,2,coordinate))
      dprobe(3,coordinate)=2.0_fp*( &
        dot_product(dgrho(:,1,coordinate),pgrho(:,1))+ &
        dot_product(grho(:,1),dpgrho(:,1,coordinate)))
      dprobe(4,coordinate)=2.0_fp*( &
        dot_product(dgrho(:,2,coordinate),pgrho(:,2))+ &
        dot_product(grho(:,2),dpgrho(:,2,coordinate)))
      dprobe(5,coordinate)= &
        dot_product(dgrho(:,1,coordinate),pgrho(:,2))+ &
        dot_product(grho(:,1),dpgrho(:,2,coordinate))+ &
        dot_product(dgrho(:,2,coordinate),pgrho(:,1))+ &
        dot_product(grho(:,2),dpgrho(:,1,coordinate))
    end do
  end subroutine build_density_variables

end module mod_dft_gridint_mrsf_xc_kernel_derivative
