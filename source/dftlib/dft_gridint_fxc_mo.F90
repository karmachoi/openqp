! Spin-polarized semilocal XC-kernel response evaluated directly in the
! occupied-virtual MO block.
!
! For a symmetric AO response density D_s = C_v X_s C_o^T + C_o X_s^T C_v^T
! of spin s, utddft_fxc forms the full AO matrix F_s = fxc[D] and the
! response solvers then transform its virtual-occupied block back to the MO
! basis.  With MO values on the grid, phi_v(p) = C_v^T Phi(p) and
! phi_o(p) = C_o^T Phi(p),
!
!   rho_s(p)      = 2 sum_ai X_ai phi_a(p) phi_i(p),
!   (C_v^T F C_o)_ai = sum_p [ f_r phi_a phi_i + c . (grad phi_a phi_i
!                                                     + phi_a grad phi_i) ],
!
! so both the density evaluation and the Fock contraction are dgemm calls
! of dimension nvir x nocc x npts per trial matrix instead of nbf x nbf x
! npts.  The kernel coefficients f_r and c are identical to those of
! utddft_fxc (mod_dft_gridint_fxc, UUpdate).
module mod_dft_gridint_fxc_mo

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr,OQP_FUNTYP_LDA

  implicit none
  private

  type :: fxc_mo_workspace_t
    real(fp), allocatable :: cmo(:,:,:)      ! gathered MO coefficients (n,nbf,2)
    real(fp), allocatable :: mov(:,:,:)      ! MO values (nbf,maxPts,2)
    real(fp), allocatable :: mog(:,:,:,:)    ! MO gradients (nbf,maxPts,3,2)
    real(fp), allocatable :: y(:,:),ytilde(:,:),ow(:,:),vw(:,:)
    real(fp), allocatable :: rho(:,:,:),grad(:,:,:,:)   ! (2,maxPts,nmtx),(3,2,maxPts,nmtx)
    real(fp), allocatable :: v(:,:,:),c(:,:,:,:)         ! (2,maxPts,nmtx),(3,2,maxPts,nmtx)
  end type fxc_mo_workspace_t

  type, extends(xc_consumer_t) :: fxc_mo_consumer_t
    real(fp), pointer :: mo(:,:,:) => null()          ! (nbf,nbf,2) normalized
    real(fp), pointer :: xa(:,:,:) => null()          ! (nvira,nocca,nmtx)
    real(fp), pointer :: xb(:,:,:) => null()          ! (nvirb,noccb,nmtx)
    real(fp), allocatable :: fa(:,:,:,:),fb(:,:,:,:)  ! (nvir,nocc,nmtx,nthreads)
    type(fxc_mo_workspace_t), allocatable :: workspace(:)
    integer :: nmtx=0,nocca=0,noccb=0
    logical :: same_mo=.true.
  contains
    procedure :: parallel_start => fxc_mo_start
    procedure :: parallel_stop => fxc_mo_stop
    procedure :: update => fxc_mo_update
    procedure :: postUpdate => fxc_mo_post
    procedure :: clean => fxc_mo_clean
  end type fxc_mo_consumer_t

  public :: utddft_fxc_mo

contains

!> Occupied-virtual MO block of the unrestricted XC-kernel response.
!>
!> mo_a/mo_b are the alpha/beta MO coefficient matrices with nocca/noccb
!> occupied columns; xa(nvira,nocca,nmtx) and xb(nvirb,noccb,nmtx) define the
!> symmetric AO response densities D_s = C_v X C_o^T + C_o X^T C_v^T.  On
!> return fa(a,i,k)=(C_v^T fxc[D_k] C_o)_ai for the alpha spin and likewise
!> fb, i.e. exactly the block that utddft_fxc followed by an AO-to-MO
!> transformation would give.  The results are added to fa/fb.
  subroutine utddft_fxc_mo(basis,molGrid,mo_a,mo_b,nocca,noccb,xa,xb,fa,fb, &
      nmtx,threshold,infos)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: molGrid
    real(fp), intent(in) :: mo_a(:,:),mo_b(:,:)
    integer, intent(in) :: nocca,noccb,nmtx
    real(fp), intent(in), target :: xa(:,:,:),xb(:,:,:)
    real(fp), intent(inout) :: fa(:,:,:),fb(:,:,:)
    real(fp), intent(in) :: threshold
    type(information), target, intent(in) :: infos

    type(fxc_mo_consumer_t) :: dat
    type(xc_options_t) :: xc_opts
    real(fp), allocatable, target :: mo(:,:,:)
    integer :: i,nbf

    nbf=basis%nbf
    if(nmtx<=0) return
    allocate(mo(nbf,nbf,2))
    do i=1,nbf
      mo(:,i,1)=mo_a(:,i)*basis%bfnrm(:)
      mo(:,i,2)=mo_b(:,i)*basis%bfnrm(:)
    end do
    dat%same_mo=all(mo_a==mo_b)
    dat%mo=>mo
    dat%xa=>xa
    dat%xb=>xb
    dat%nmtx=nmtx
    dat%nocca=nocca
    dat%noccb=noccb

    xc_opts%isGGA=infos%functional%needGrd
    xc_opts%needTau=infos%functional%needTau
    xc_opts%functional=>infos%functional
    xc_opts%hasBeta=.true.
    xc_opts%isWFVecs=.true.
    xc_opts%numAOs=nbf
    xc_opts%maxPts=molGrid%maxSlicePts
    xc_opts%limPts=molGrid%maxNRadTimesNAng
    xc_opts%numAtoms=infos%mol_prop%natom
    xc_opts%maxAngMom=basis%mxam
    xc_opts%nDer=0
    xc_opts%nXCDer=2
    xc_opts%numOccAlpha=infos%mol_prop%nelec_A
    xc_opts%numOccBeta=infos%mol_prop%nelec_B
    xc_opts%wfAlpha=>mo(:,:,1)
    xc_opts%wfBeta=>mo(:,:,2)
    xc_opts%dft_threshold=threshold
    xc_opts%molGrid=>molGrid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(xc_opts,dat,basis)
    fa=fa+dat%fa(:,:,:,1)
    fb=fb+dat%fb(:,:,:,1)
    call dat%clean()
    deallocate(mo)
  end subroutine utddft_fxc_mo

!-------------------------------------------------------------------------------

  subroutine fxc_mo_start(self,xce,nthreads)
    class(fxc_mo_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: nbf,nocca,noccb,nvira,nvirb,mp,thread
    nbf=xce%numAOs
    nocca=self%nocca
    noccb=self%noccb
    nvira=nbf-nocca
    nvirb=nbf-noccb
    mp=max(1,xce%maxPts)
    allocate(self%fa(nvira,nocca,self%nmtx,nthreads), &
      self%fb(nvirb,noccb,self%nmtx,nthreads),source=0.0_fp)
    allocate(self%workspace(nthreads))
    do thread=1,nthreads
      associate(ws=>self%workspace(thread))
      allocate(ws%cmo(nbf,nbf,2),ws%mov(nbf,mp,2),ws%mog(nbf,mp,3,2), &
        ws%y(nbf,mp),ws%ytilde(nbf,mp),ws%ow(nbf,mp),ws%vw(nbf,mp), &
        ws%rho(2,mp,self%nmtx),ws%grad(3,2,mp,self%nmtx), &
        ws%v(2,mp,self%nmtx),ws%c(3,2,mp,self%nmtx))
      end associate
    end do
  end subroutine fxc_mo_start

  subroutine fxc_mo_stop(self)
    class(fxc_mo_consumer_t), intent(inout) :: self
    if(size(self%fa,4)>1) then
      self%fa(:,:,:,1)=sum(self%fa,dim=4)
      self%fb(:,:,:,1)=sum(self%fb,dim=4)
    end if
    call self%pe%allreduce(self%fa(:,:,:,1),size(self%fa(:,:,:,1)))
    call self%pe%allreduce(self%fb(:,:,:,1),size(self%fb(:,:,:,1)))
  end subroutine fxc_mo_stop

  subroutine fxc_mo_post(self,xce,mythread)
    class(fxc_mo_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
  end subroutine fxc_mo_post

  subroutine fxc_mo_clean(self)
    class(fxc_mo_consumer_t), intent(inout) :: self
    integer :: thread
    if(allocated(self%fa)) deallocate(self%fa)
    if(allocated(self%fb)) deallocate(self%fb)
    if(allocated(self%workspace)) then
      do thread=1,size(self%workspace)
        associate(ws=>self%workspace(thread))
        if(allocated(ws%cmo)) deallocate(ws%cmo,ws%mov,ws%mog,ws%y, &
          ws%ytilde,ws%ow,ws%vw,ws%rho,ws%grad,ws%v,ws%c)
        end associate
      end do
      deallocate(self%workspace)
    end if
    nullify(self%mo,self%xa,self%xb)
  end subroutine fxc_mo_clean

!-------------------------------------------------------------------------------

  subroutine fxc_mo_update(self,xce,mythread)
    class(fxc_mo_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(fp) :: rhoab(2),sigma(3),tauab(2),d_r(2),d_s(3),d_t(2)
    real(fp) :: f_r(2),f_s(3),f_t(2),dsaa,dsab,dsbb,dsba,gr(3,2)
    integer :: c,i,k,n,nbf,nocc,nspin,nvir,npts,spin
    logical :: gga

    n=xce%numAOs_p
    nbf=xce%numAOs
    npts=xce%numPts
    if(n<=0 .or. npts<=0) return
    gga=xce%funTyp/=OQP_FUNTYP_LDA
    nspin=merge(1,2,self%same_mo)
    associate(ws=>self%workspace(mythread),aoV=>xce%aoV,aoG1=>xce%aoG1, &
      dRho=>xce%xclib%dRho)

    ! MO values and gradients on the slice.
    do spin=1,nspin
      if(xce%skip_p) then
        ws%cmo(:,:,spin)=self%mo(:,:,spin)
      else
        ws%cmo(1:n,:,spin)=self%mo(xce%indices_p(1:n),:,spin)
      end if
      call dgemm('T','N',nbf,npts,n,1.0_fp,ws%cmo(1,1,spin),nbf,aoV,n, &
        0.0_fp,ws%mov(1,1,spin),nbf)
      if(gga) then
        do c=1,3
          call dgemm('T','N',nbf,npts,n,1.0_fp,ws%cmo(1,1,spin),nbf, &
            aoG1(1,1,c),n,0.0_fp,ws%mog(1,1,c,spin),nbf)
        end do
      end if
    end do

    ! Response densities rho_k^s(p) and gradients for every trial matrix.
    do k=1,self%nmtx
      do spin=1,2
        call spin_dimensions(self,nbf,spin,nocc,nvir)
        call response_density(self,ws,xce,k,spin,merge(spin,1,nspin==2), &
          nbf,npts,nocc,nvir,gga)
      end do
    end do

    ! Kernel coefficients, identical to UUpdate of mod_dft_gridint_fxc.
    do k=1,self%nmtx
      do i=1,npts
        rhoab=ws%rho(:,i,k)
        sigma=0.0_fp
        tauab=0.0_fp
        if(gga) then
          gr=ws%grad(:,:,i,k)
          dsaa=dot_product(gr(:,1),dRho(1:3,i))
          dsab=dot_product(gr(:,1),dRho(4:6,i))
          dsbb=dot_product(gr(:,2),dRho(4:6,i))
          dsba=dot_product(gr(:,2),dRho(1:3,i))
          sigma=[2.0_fp*dsaa,2.0_fp*dsbb,dsab+dsba]
        end if
        call xc_der1(xce,.true.,i,d_r,d_s,d_t)
        call xc_der2_contr(xce,.true.,i,rhoab,sigma,tauab,f_r,f_s,f_t)
        if(gga) then
          if(maxval(abs([dsaa,dsbb,dsab,dsba]))<xce%threshold) f_s=0.0_fp
        end if
        ws%v(:,i,k)=f_r
        ws%c(:,:,i,k)=0.0_fp
        if(gga) then
          ws%c(:,1,i,k)=2.0_fp*f_s(1)*dRho(1:3,i)+f_s(3)*dRho(4:6,i)+ &
            2.0_fp*d_s(1)*gr(:,1)+d_s(3)*gr(:,2)
          ws%c(:,2,i,k)=2.0_fp*f_s(2)*dRho(4:6,i)+f_s(3)*dRho(1:3,i)+ &
            2.0_fp*d_s(2)*gr(:,2)+d_s(3)*gr(:,1)
        end if
      end do
    end do

    ! Occupied-virtual blocks of the kernel response.
    do k=1,self%nmtx
      do spin=1,2
        call spin_dimensions(self,nbf,spin,nocc,nvir)
        call fock_block(self,ws,mythread,k,spin,merge(spin,1,nspin==2), &
          nbf,npts,nocc,nvir,gga)
      end do
    end do
    end associate
  end subroutine fxc_mo_update

  subroutine spin_dimensions(self,nbf,spin,nocc,nvir)
    class(fxc_mo_consumer_t), intent(in) :: self
    integer, intent(in) :: nbf,spin
    integer, intent(out) :: nocc,nvir
    if(spin==1) then
      nocc=self%nocca
    else
      nocc=self%noccb
    end if
    nvir=nbf-nocc
  end subroutine spin_dimensions

!> rho(p)=2 sum_a phi_a(p) y_a(p), y=X phi_o, and
!> grad rho=2 sum_a grad phi_a y_a + 2 sum_i ytilde_i grad phi_i, ytilde=X^T phi_v.
  subroutine response_density(self,ws,xce,k,spin,mo_spin,nbf,npts,nocc,nvir, &
      gga)
    class(fxc_mo_consumer_t), intent(in) :: self
    type(fxc_mo_workspace_t), intent(inout) :: ws
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: k,spin,mo_spin,nbf,npts,nocc,nvir
    logical, intent(in) :: gga
    integer :: c,p
    real(fp) :: s

    if(nocc<=0 .or. nvir<=0) then
      ws%rho(spin,1:npts,k)=0.0_fp
      ws%grad(:,spin,1:npts,k)=0.0_fp
      return
    end if
    if(spin==1) then
      call dgemm('N','N',nvir,npts,nocc,1.0_fp,self%xa(1,1,k),nvir, &
        ws%mov(1,1,mo_spin),nbf,0.0_fp,ws%y,nbf)
      if(gga) call dgemm('T','N',nocc,npts,nvir,1.0_fp,self%xa(1,1,k),nvir, &
        ws%mov(nocc+1,1,mo_spin),nbf,0.0_fp,ws%ytilde,nbf)
    else
      call dgemm('N','N',nvir,npts,nocc,1.0_fp,self%xb(1,1,k),nvir, &
        ws%mov(1,1,mo_spin),nbf,0.0_fp,ws%y,nbf)
      if(gga) call dgemm('T','N',nocc,npts,nvir,1.0_fp,self%xb(1,1,k),nvir, &
        ws%mov(nocc+1,1,mo_spin),nbf,0.0_fp,ws%ytilde,nbf)
    end if
    do p=1,npts
      ws%rho(spin,p,k)=2.0_fp*dot_product(ws%mov(nocc+1:nbf,p,mo_spin), &
        ws%y(1:nvir,p))
    end do
    if(.not.gga) then
      ws%grad(:,spin,1:npts,k)=0.0_fp
      return
    end if
    do p=1,npts
      do c=1,3
        s=dot_product(ws%mog(nocc+1:nbf,p,c,mo_spin),ws%y(1:nvir,p))+ &
          dot_product(ws%ytilde(1:nocc,p),ws%mog(1:nocc,p,c,mo_spin))
        ws%grad(c,spin,p,k)=2.0_fp*s
      end do
    end do
  end subroutine response_density

!> F_ai += sum_p phi_a(p) ow_i(p) + vw_a(p) phi_i(p) with
!> ow_i=v phi_i + c.grad phi_i and vw_a=c.grad phi_a.
  subroutine fock_block(self,ws,mythread,k,spin,mo_spin,nbf,npts,nocc,nvir,gga)
    class(fxc_mo_consumer_t), intent(inout) :: self
    type(fxc_mo_workspace_t), intent(inout) :: ws
    integer, intent(in) :: mythread,k,spin,mo_spin,nbf,npts,nocc,nvir
    logical, intent(in) :: gga
    integer :: p
    real(fp) :: cc(3)

    if(nocc<=0 .or. nvir<=0) return
    do p=1,npts
      ws%ow(1:nocc,p)=ws%v(spin,p,k)*ws%mov(1:nocc,p,mo_spin)
      if(gga) then
        cc=ws%c(:,spin,p,k)
        ws%ow(1:nocc,p)=ws%ow(1:nocc,p)+cc(1)*ws%mog(1:nocc,p,1,mo_spin)+ &
          cc(2)*ws%mog(1:nocc,p,2,mo_spin)+cc(3)*ws%mog(1:nocc,p,3,mo_spin)
        ws%vw(1:nvir,p)=cc(1)*ws%mog(nocc+1:nbf,p,1,mo_spin)+ &
          cc(2)*ws%mog(nocc+1:nbf,p,2,mo_spin)+cc(3)*ws%mog(nocc+1:nbf,p,3,mo_spin)
      end if
    end do
    if(spin==1) then
      call dgemm('N','T',nvir,nocc,npts,1.0_fp,ws%mov(nocc+1,1,mo_spin),nbf, &
        ws%ow,nbf,1.0_fp,self%fa(1,1,k,mythread),nvir)
      if(gga) call dgemm('N','T',nvir,nocc,npts,1.0_fp,ws%vw,nbf, &
        ws%mov(1,1,mo_spin),nbf,1.0_fp,self%fa(1,1,k,mythread),nvir)
    else
      call dgemm('N','T',nvir,nocc,npts,1.0_fp,ws%mov(nocc+1,1,mo_spin),nbf, &
        ws%ow,nbf,1.0_fp,self%fb(1,1,k,mythread),nvir)
      if(gga) call dgemm('N','T',nvir,nocc,npts,1.0_fp,ws%vw,nbf, &
        ws%mov(1,1,mo_spin),nbf,1.0_fp,self%fb(1,1,k,mythread),nvir)
    end if
  end subroutine fock_block

end module mod_dft_gridint_fxc_mo
