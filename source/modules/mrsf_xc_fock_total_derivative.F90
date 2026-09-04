! Analytic total nuclear derivative of the spin-resolved semilocal XC Fock
! matrices needed by MRSF-TDDFT and ROKS response equations.
module mrsf_xc_fock_total_derivative_mod

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr
  use mod_dft_gridint_fxc, only: utddft_fxc
  use mod_dft_partition_hessian, only: &
    partition_weight_nuclear_first_derivatives
  use mod_dft_gridint_tdgga_response, only: gga_add_owner_motion_first
  use mod_dft_gridint_mrsf_xc_slice_gemm, only: slice_stack_values, &
    slice_stack_fixed,slice_fock_derivative_accumulate,slice_chunk_size, &
    symmetrize_half_accumulator,gather_stack,scatter_half_accumulator, &
    build_symmetrized_stack

  implicit none
  private

  integer(8), parameter :: chunk_budget_bytes=24_8*1024_8*1024_8

  !> Per-thread slice workspace.  Chunk-level arrays are sized for the full
  !> AO count and reused with the pruned count through sequence association.
  type :: fock_derivative_workspace_t
    real(fp), allocatable :: t(:,:,:),u(:,:,:,:),value(:,:),gradient(:,:,:)
    real(fp), allocatable :: fixed_d(:,:,:,:),fixed_g(:,:,:,:,:)
    real(fp), allocatable :: a(:,:,:,:),b(:,:,:,:,:),v(:,:,:),c(:,:,:,:),fw(:)
    real(fp), allocatable :: xs(:,:,:,:),z(:,:),psiw(:,:,:),g(:,:,:)
    real(fp), allocatable :: stack_p(:,:,:),acc_p(:,:,:,:,:)
    integer, allocatable :: ao_atom_p(:)
    real(fp), allocatable :: total_d(:,:,:,:),total_g(:,:,:,:,:)
    real(fp), allocatable :: weights(:),dweights(:,:,:),dweight_flat(:)
    real(fp), allocatable :: drho(:,:),dgrad_rho(:,:,:),dvr(:,:),dvs(:,:)
    integer :: mchunk=0
  contains
    procedure :: init => workspace_init
    procedure :: clean => workspace_clean
  end type fock_derivative_workspace_t

  type, extends(xc_consumer_t) :: xc_fock_derivative_consumer_t
    !> Half accumulators H (mu,nu,coordinate,1,spin,thread).
    real(fp), allocatable :: derivative(:,:,:,:,:,:)
    !> Symmetrized reference densities (mu,spin,nu).
    real(fp), allocatable :: stack(:,:,:)
    type(fock_derivative_workspace_t), allocatable :: workspace(:)
    real(fp), allocatable :: atom_xyz(:,:),surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    real(fp), allocatable :: worker_error(:)
    integer :: part_fun_type=0
    logical :: has_surface_shift=.false.
    logical :: is_gga=.false.
  contains
    procedure :: parallel_start => derivative_parallel_start
    procedure :: parallel_stop => derivative_parallel_stop
    procedure :: update => derivative_update
    procedure :: postUpdate => derivative_post_update
    procedure :: clean => derivative_clean
  end type xc_fock_derivative_consumer_t

  public :: mrsf_xc_fock_total_derivative
  public :: mrsf_xc_fock_total_derivative_from_dmo

contains

!> Spin-resolved semilocal XC Fock total nuclear derivative.
!>
!> Methodological starting point: Hiroya Nakata's analytical TDHF/TDDFT
!> Hessian formulation.  The first derivative is separated exactly as
!>
!>   dVxc_s = (dVxc_s)_moving-AO/grid + sum_t fxc_st[dP_t].
!>
!> The first term is evaluated by analytic atom-centred grid integrals below;
!> the second reuses the production unrestricted `utddft_fxc` kernel action.
!> Nuclear finite differences and displaced geometries are not used.
!>
!> pa/pb and dpa/dpb are spin density matrices in the AO coefficient
!> representation.  The output coordinate order is x1,y1,z1,x2,y2,z2,... .
!> This density-only interface makes no determinant or Slater construction.
!> Global hybrids are accepted because their semilocal XC part has the same
!> derivative; the exact-exchange derivative is assembled by its own Fock
!> path.  Meta-GGA and range-separated/CAM cases fail closed.
  subroutine mrsf_xc_fock_total_derivative(basis,mol_grid,pa,pb,dpa,dpb, &
      derivative_a,derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: pa(:,:),pb(:,:),dpa(:,:,:),dpb(:,:,:)
    real(fp), intent(out) :: derivative_a(:,:,:),derivative_b(:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    type(xc_fock_derivative_consumer_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: pa_normalized(:,:),pb_normalized(:,:)
    real(fp), allocatable :: dpa_work(:,:,:),dpb_work(:,:,:)
    real(fp), allocatable :: kernel_a(:,:,:),kernel_b(:,:,:)
    real(fp) :: grid_threshold
    integer :: first,i,k,last,natom,nbf,ncart,shell

    nbf=basis%nbf
    natom=infos%mol_prop%natom
    ncart=3*natom
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. natom<=0 .or. any(shape(pa)/=[nbf,nbf]) .or. &
       any(shape(pb)/=[nbf,nbf]) .or. any(shape(dpa)/=[nbf,nbf,ncart]) .or. &
       any(shape(dpb)/=[nbf,nbf,ncart]) .or. &
       any(shape(derivative_a)/=[nbf,nbf,ncart]) .or. &
       any(shape(derivative_b)/=[nbf,nbf,ncart])) then
      status=-1
      return
    end if
    if(infos%functional%needTau .or. infos%functional%needLapl) then
      status=-2
      return
    end if
    if(infos%dft%cam_flag) then
      status=-3
      return
    end if

    grid_threshold=infos%dft%grid_density_cutoff
    if(present(threshold)) grid_threshold=threshold
    if(grid_threshold<0.0_fp) then
      status=-4
      return
    end if

    allocate(pa_normalized(nbf,nbf),pb_normalized(nbf,nbf), &
      dpa_work(nbf,nbf,ncart),dpb_work(nbf,nbf,ncart), &
      kernel_a(nbf,nbf,ncart),kernel_b(nbf,nbf,ncart), &
      dat%ao_atom(nbf),dat%stack(nbf,2,nbf))
    do i=1,nbf
      pa_normalized(:,i)=pa(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      pb_normalized(:,i)=pb(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
    end do
    call build_symmetrized_stack(nbf,2,1,pa,basis%bfnrm,1.0_fp,dat%stack)
    call build_symmetrized_stack(nbf,2,2,pb,basis%bfnrm,1.0_fp,dat%stack)
    dat%ao_atom=0
    do shell=1,basis%nshell
      first=basis%ao_offset(shell)
      last=first+basis%naos(shell)-1
      dat%ao_atom(first:last)=basis%origin(shell)
    end do
    if(any(dat%ao_atom<1) .or. any(dat%ao_atom>natom)) then
      status=-5
      call dat%clean()
      deallocate(pa_normalized,pb_normalized,dpa_work,dpb_work,kernel_a,kernel_b)
      return
    end if

    dat%atom_xyz=infos%atoms%xyz
    dat%dummy_atom=mol_grid%dummyAtom
    dat%surface_shift=mol_grid%surfaceShift
    dat%part_fun_type=mol_grid%partFunType
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
    ! run_xc adds one derivative for GGA.  Both branches therefore collocate
    ! AO values through G2, as required by d grad(phi_mu phi_nu)/dR.
    opts%nDer=merge(1,2,dat%is_gga)
    opts%nXCDer=2
    opts%numOccAlpha=infos%mol_prop%nelec_A
    opts%numOccBeta=infos%mol_prop%nelec_B
    opts%wfAlpha=>pa_normalized
    opts%wfBeta=>pb_normalized
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
      if(dat%worker_error(1)>0.0_fp) status=-6
    end if
    if(status==0) then
      do k=1,ncart
        do i=1,nbf
          derivative_a(:,i,k)=dat%derivative(:,i,k,1,1,1) &
            *basis%bfnrm(:)*basis%bfnrm(i)
          derivative_b(:,i,k)=dat%derivative(:,i,k,1,2,1) &
            *basis%bfnrm(:)*basis%bfnrm(i)
        end do
      end do
    end if
    call dat%clean()

    if(status==0) then
      ! utddft_fxc scales its density arguments temporarily, so use private
      ! copies even though it restores them before returning.
      dpa_work=dpa
      dpb_work=dpb
      kernel_a=0.0_fp
      kernel_b=0.0_fp
      call utddft_fxc(basis=basis,molGrid=mol_grid,isVecs=.false., &
        wfa=pa,wfb=pb,fxa=kernel_a,fxb=kernel_b,dxa=dpa_work,dxb=dpb_work, &
        nMtx=ncart,threshold=grid_threshold,infos=infos)
      derivative_a=derivative_a+kernel_a
      derivative_b=derivative_b+kernel_b
    end if
    deallocate(pa_normalized,pb_normalized,dpa_work,dpb_work,kernel_a,kernel_b)
  end subroutine mrsf_xc_fock_total_derivative

!-------------------------------------------------------------------------------

!> Orbital-response convenience interface.  Occupation vectors permit ROKS
!> and other spin-resolved ensembles without introducing configuration-state
!> or determinant objects.
  subroutine mrsf_xc_fock_total_derivative_from_dmo(basis,mol_grid,pa,pb, &
      mo_a,mo_b,dmo_a,dmo_b,occupation_a,occupation_b,derivative_a, &
      derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: pa(:,:),pb(:,:),mo_a(:,:),mo_b(:,:)
    real(fp), intent(in) :: dmo_a(:,:,:),dmo_b(:,:,:)
    real(fp), intent(in) :: occupation_a(:),occupation_b(:)
    real(fp), intent(out) :: derivative_a(:,:,:),derivative_b(:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    real(fp), allocatable :: dpa(:,:,:),dpb(:,:,:)
    integer :: i,j,k,ncart,nbf,orbital

    nbf=size(pa,1)
    ncart=size(dmo_a,3)
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. size(pa,2)/=nbf .or. any(shape(pb)/=[nbf,nbf]) .or. &
       size(mo_a,1)/=nbf .or. size(mo_b,1)/=nbf .or. &
       size(dmo_a,1)/=nbf .or. size(dmo_b,1)/=nbf .or. &
       size(dmo_a,2)/=size(mo_a,2) .or. size(dmo_b,2)/=size(mo_b,2) .or. &
       size(dmo_b,3)/=ncart .or. size(occupation_a)/=size(mo_a,2) .or. &
       size(occupation_b)/=size(mo_b,2)) then
      status=-1
      return
    end if
    allocate(dpa(nbf,nbf,ncart),dpb(nbf,nbf,ncart),source=0.0_fp)
    do k=1,ncart
      do orbital=1,size(mo_a,2)
        if(abs(occupation_a(orbital))<=tiny(1.0_fp)) cycle
        do j=1,nbf
          do i=1,nbf
            dpa(i,j,k)=dpa(i,j,k)+occupation_a(orbital)*( &
              dmo_a(i,orbital,k)*mo_a(j,orbital) &
              +mo_a(i,orbital)*dmo_a(j,orbital,k))
          end do
        end do
      end do
      do orbital=1,size(mo_b,2)
        if(abs(occupation_b(orbital))<=tiny(1.0_fp)) cycle
        do j=1,nbf
          do i=1,nbf
            dpb(i,j,k)=dpb(i,j,k)+occupation_b(orbital)*( &
              dmo_b(i,orbital,k)*mo_b(j,orbital) &
              +mo_b(i,orbital)*dmo_b(j,orbital,k))
          end do
        end do
      end do
    end do
    call mrsf_xc_fock_total_derivative(basis,mol_grid,pa,pb,dpa,dpb, &
      derivative_a,derivative_b,infos,status,threshold)
    deallocate(dpa,dpb)
  end subroutine mrsf_xc_fock_total_derivative_from_dmo

!-------------------------------------------------------------------------------

  subroutine workspace_init(self,nbf,nat,max_points)
    class(fock_derivative_workspace_t), intent(inout) :: self
    integer, intent(in) :: nbf,nat,max_points
    integer :: m,ncart

    call self%clean()
    ncart=3*nat
    m=slice_chunk_size(nbf,max(ncart,6),max_points,chunk_budget_bytes)
    self%mchunk=m
    allocate(self%t(nbf,2,m),self%u(nbf,2,m,3),self%value(2,m), &
      self%gradient(3,2,m),self%fixed_d(3,nat,2,m),self%fixed_g(3,3,nat,2,m), &
      self%a(2,ncart,1,m),self%b(3,2,ncart,1,m),self%v(2,1,m), &
      self%c(3,2,1,m),self%fw(m),self%xs(nbf,ncart,1,m),self%z(nbf,m), &
      self%psiw(nbf,3,m),self%g(nbf,nbf,3))
    allocate(self%total_d(3,nat,2,1),self%total_g(3,3,nat,2,1), &
      self%weights(nat),self%dweights(3,nat,nat),self%dweight_flat(ncart), &
      self%drho(2,ncart),self%dgrad_rho(3,2,ncart),self%dvr(2,ncart), &
      self%dvs(3,ncart))
  end subroutine workspace_init

  subroutine workspace_clean(self)
    class(fock_derivative_workspace_t), intent(inout) :: self
    if(allocated(self%t)) deallocate(self%t,self%u,self%value, &
      self%gradient,self%fixed_d,self%fixed_g,self%a,self%b,self%v,self%c, &
      self%fw,self%xs,self%z,self%psiw,self%g)
    if(allocated(self%stack_p)) deallocate(self%stack_p)
    if(allocated(self%acc_p)) deallocate(self%acc_p)
    if(allocated(self%ao_atom_p)) deallocate(self%ao_atom_p)
    if(allocated(self%total_d)) deallocate(self%total_d,self%total_g, &
      self%weights,self%dweights,self%dweight_flat,self%drho, &
      self%dgrad_rho,self%dvr,self%dvs)
    self%mchunk=0
  end subroutine workspace_clean

!-------------------------------------------------------------------------------

  subroutine derivative_parallel_start(self,xce,nthreads)
    class(xc_fock_derivative_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: ncart,thread
    ncart=3*xce%numAtoms
    allocate(self%derivative(xce%numAOs,xce%numAOs,ncart,1,2,nthreads), &
      self%worker_error(nthreads),source=0.0_fp)
    allocate(self%workspace(nthreads))
    do thread=1,nthreads
      call self%workspace(thread)%init(xce%numAOs,xce%numAtoms, &
        max(1,xce%maxPts))
    end do
  end subroutine derivative_parallel_start

  subroutine derivative_parallel_stop(self)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    integer :: spin
    if(size(self%derivative,6)>1) then
      self%derivative(:,:,:,:,:,1)=sum(self%derivative,dim=6)
      self%worker_error(1)=sum(self%worker_error)
    end if
    call symmetrize_half_accumulator(size(self%derivative,1), &
      2*size(self%derivative,3),self%derivative(:,:,:,:,:,1))
    do spin=1,2
      call self%pe%allreduce(self%derivative(:,:,:,:,spin,1), &
        size(self%derivative(:,:,:,:,spin,1)))
    end do
    call self%pe%allreduce(self%worker_error(1),1)
  end subroutine derivative_parallel_stop

  subroutine derivative_post_update(self,xce,mythread)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    ! Contributions are accumulated directly by update.
  end subroutine derivative_post_update

  subroutine derivative_clean(self)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    integer :: thread
    if(allocated(self%derivative)) deallocate(self%derivative)
    if(allocated(self%stack)) deallocate(self%stack)
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
  end subroutine derivative_clean

!-------------------------------------------------------------------------------

!> One grid slice.  The fixed-grid density derivatives, the per-point
!> potentials of all coordinates, and the AO-matrix accumulation are the
!> slice-level dgemm contractions of mod_dft_gridint_mrsf_xc_slice_gemm.
  subroutine derivative_update(self,xce,mythread)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    integer :: n,nbf,ncart,status

    nbf=xce%numAOs
    n=xce%numAOs_p
    ncart=3*xce%numAtoms
    if(n<=0 .or. xce%numPts<=0) return
    if(xce%currAtom<1 .or. xce%currAtom>xce%numAtoms) then
      self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
      return
    end if
    associate(ws=>self%workspace(mythread))
    if(xce%skip_p) then
      call derivative_slice(self,xce,mythread,n,self%stack,self%ao_atom, &
        self%derivative(:,:,:,:,:,mythread),status)
    else
      if(.not.allocated(ws%stack_p)) &
        allocate(ws%stack_p(nbf,2,nbf),ws%acc_p(nbf,nbf,ncart,1,2), &
          ws%ao_atom_p(nbf))
      call gather_stack(nbf,n,2,xce%indices_p(1:n),self%stack,ws%stack_p)
      ws%ao_atom_p(1:n)=self%ao_atom(xce%indices_p(1:n))
      ws%acc_p=0.0_fp
      call derivative_slice(self,xce,mythread,n,ws%stack_p,ws%ao_atom_p, &
        ws%acc_p,status)
      call scatter_half_accumulator(nbf,n,2*ncart,xce%indices_p(1:n), &
        ws%acc_p,self%derivative(:,:,:,:,:,mythread))
    end if
    end associate
    if(status/=0) self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
  end subroutine derivative_update

  subroutine derivative_slice(self,xce,mythread,n,stack,ao_atom,acc,status)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread,n
    real(fp), intent(in) :: stack(*)
    integer, intent(in) :: ao_atom(n)
    real(fp), intent(inout) :: acc(*)
    integer, intent(out) :: status

    integer :: m,nat,ncart,npts,p,p0,local_status

    nat=xce%numAtoms
    ncart=3*nat
    npts=xce%numPts
    status=0
    associate(ws=>self%workspace(mythread))
    do p0=1,npts,ws%mchunk
      m=min(ws%mchunk,npts-p0+1)
      call slice_stack_values(n,m,2,stack,xce%aoV,xce%aoG1,p0,ws%t, &
        ws%value,ws%gradient)
      call slice_stack_fixed(n,m,2,nat,stack,xce%aoG1,xce%aoG2,p0,ao_atom, &
        ws%t,ws%u,ws%fixed_d,ws%fixed_g)
      do p=1,m
        call point_potentials(self,xce,p0+p-1,p,ws,local_status)
        if(local_status/=0) then
          status=local_status
          return
        end if
      end do
      call slice_fock_derivative_accumulate(n,m,ncart,1,xce%aoV,xce%aoG1, &
        xce%aoG2,p0,ao_atom,xce%currAtom,ws%fw,ws%a,ws%b,ws%v,ws%c,ws%xs, &
        ws%z,ws%psiw,ws%g,acc)
    end do
    end associate
  end subroutine derivative_slice

!> Per-point potentials of every coordinate: the partition-weight derivative
!> multiplies the undifferentiated potential, and the moving-grid density
!> derivative (fixed grid plus owner motion) enters the second functional
!> derivative.  The AO-coefficient response is handled by utddft_fxc.
  subroutine point_potentials(self,xce,ipt,p,ws,status)
    class(xc_fock_derivative_consumer_t), intent(in) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: ipt,p
    type(fock_derivative_workspace_t), intent(inout) :: ws
    integer, intent(out) :: status

    real(fp) :: dr(2),ds(3),dt(2),fr(2),fs(3),ft(2)
    real(fp) :: vr(2),vs(3),vt(2),grad_rho(3,2)
    real(fp) :: coefficient(3,2),dcoefficient(3,2)
    real(fp) :: finite_weight,quadrature_scale,scale
    integer :: atom,cart,coordinate,nat,ncart,owner,partition_status,spin

    nat=xce%numAtoms
    ncart=3*nat
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
      partition_status)
    if(partition_status/=0 .or. ws%weights(owner)<=sqrt(tiny(1.0_fp))) then
      status=-1
      return
    end if
    quadrature_scale=finite_weight/ws%weights(owner)
    do atom=1,nat
      do cart=1,3
        coordinate=3*(atom-1)+cart
        ws%dweight_flat(coordinate)=ws%dweights(cart,atom,owner)
      end do
    end do

    do spin=1,2
      call gga_add_owner_motion_first(owner,ws%fixed_d(:,:,spin,p), &
        ws%fixed_g(:,:,:,spin,p),ws%total_d(:,:,spin,1), &
        ws%total_g(:,:,:,spin,1))
    end do
    do atom=1,nat
      do cart=1,3
        coordinate=3*(atom-1)+cart
        ws%drho(:,coordinate)=ws%total_d(cart,atom,:,1)
        ws%dgrad_rho(:,:,coordinate)=ws%total_g(:,cart,atom,:,1)
      end do
    end do

    call xc_der1(xce,.true.,ipt,vr,vs,vt)
    vr=vr/finite_weight
    vs=vs/finite_weight
    grad_rho=0.0_fp
    if(self%is_gga) then
      grad_rho(:,1)=xce%xclib%drho(1:3,ipt)
      grad_rho(:,2)=xce%xclib%drho(4:6,ipt)
    end if
    do coordinate=1,ncart
      dr=ws%drho(:,coordinate)
      ds=0.0_fp
      if(self%is_gga) then
        ds(1)=2.0_fp*dot_product(grad_rho(:,1),ws%dgrad_rho(:,1,coordinate))
        ds(2)=2.0_fp*dot_product(grad_rho(:,2),ws%dgrad_rho(:,2,coordinate))
        ds(3)=dot_product(ws%dgrad_rho(:,1,coordinate),grad_rho(:,2)) &
          +dot_product(grad_rho(:,1),ws%dgrad_rho(:,2,coordinate))
      end if
      dt=0.0_fp
      call xc_der2_contr(xce,.true.,ipt,dr,ds,dt,fr,fs,ft)
      ws%dvr(:,coordinate)=fr/finite_weight
      ws%dvs(:,coordinate)=fs/finite_weight
    end do

    coefficient=0.0_fp
    if(self%is_gga) then
      coefficient(:,1)=2.0_fp*vs(1)*grad_rho(:,1)+vs(3)*grad_rho(:,2)
      coefficient(:,2)=2.0_fp*vs(2)*grad_rho(:,2)+vs(3)*grad_rho(:,1)
    end if
    do coordinate=1,ncart
      dcoefficient=0.0_fp
      if(self%is_gga) then
        dcoefficient(:,1)=2.0_fp*ws%dvs(1,coordinate)*grad_rho(:,1) &
          +2.0_fp*vs(1)*ws%dgrad_rho(:,1,coordinate) &
          +ws%dvs(3,coordinate)*grad_rho(:,2) &
          +vs(3)*ws%dgrad_rho(:,2,coordinate)
        dcoefficient(:,2)=2.0_fp*ws%dvs(2,coordinate)*grad_rho(:,2) &
          +2.0_fp*vs(2)*ws%dgrad_rho(:,2,coordinate) &
          +ws%dvs(3,coordinate)*grad_rho(:,1) &
          +vs(3)*ws%dgrad_rho(:,1,coordinate)
      end if
      scale=quadrature_scale*ws%dweight_flat(coordinate)
      do spin=1,2
        ws%a(spin,coordinate,1,p)=scale*vr(spin)+ &
          finite_weight*ws%dvr(spin,coordinate)
        ws%b(:,spin,coordinate,1,p)=scale*coefficient(:,spin)+ &
          finite_weight*dcoefficient(:,spin)
      end do
    end do
    do spin=1,2
      ws%v(spin,1,p)=vr(spin)
      ws%c(:,spin,1,p)=coefficient(:,spin)
    end do
  end subroutine point_potentials

end module mrsf_xc_fock_total_derivative_mod
