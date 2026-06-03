module mod_dft_gridint_grad

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use mod_dft_gridint, only: OQP_FUNTYP_LDA, OQP_FUNTYP_MGGA
  use mod_dft_gridint, only: compAtGradRho, compAtGradDRho, compAtGradTau
  use mod_dft_partfunc, only: partition_function

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_grad_t
    real(kind=fp), allocatable :: bfgrad(:,:,:)
    real(kind=fp), allocatable :: tmp_(:,:)
    real(kind=fp), allocatable :: d1dsx(:,:,:) !< Temporary storage for dE/d\sigma
    !> @name Quadrature-weight (Becke/SSF) derivative contribution
    !> @{
    !> EXPERIMENTAL, opt-in via dftgrid.weight_derivatives. Adds the partition
    !> weight-derivative term to the nuclear gradient. Only consistent with the
    !> SSF-type fuzzy cell (dft_bfc_algo=0, no surface shifting).
    logical :: wtDeriv = .false.
    integer :: nat = 0
    real(kind=fp), allocatable :: atxyz(:,:)    !< atomic coordinates (3,nat)
    real(kind=fp), allocatable :: rijInv(:,:)   !< inverse interatomic distances
    logical, allocatable :: dummyAtom(:)        !< .true. for point-charge/dummy atoms
    type(partition_function) :: partfunc        !< partition function (eval + deriv)
    real(kind=fp), allocatable :: wtgrad(:,:,:) !< weight-deriv gradient (3,nat,nthreads)
    !> @}
  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: resetGradPointers
    procedure :: update
    procedure :: postUpdate
    procedure :: compWtGrad
    procedure :: clean
  end type

!-------------------------------------------------------------------------------

  private
  public derexc_blk

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nthreads)
    implicit none
    class(xc_consumer_grad_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    if (allocated(self%bfGrad)) deallocate(self%bfGrad)
    if (allocated(self%d1dsx)) deallocate(self%d1dsx)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    if (allocated(self%wtgrad)) deallocate(self%wtgrad)
    allocate( self%bfGrad(xce%numAOs, 3, nthreads) &
            , self%d1dsx(xce%maxPts, 3, nthreads) &
            , self%tmp_(xce%numAOs*3, nthreads) &
            , source=0.0d0)
    if (self%wtDeriv) &
      allocate(self%wtgrad(3, self%nat, nthreads), source=0.0d0)
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_grad_t), intent(inout) :: self
    if (ubound(self%bfGrad,3) /= 1) then
      self%bfGrad(:,:,lbound(self%bfGrad,3)) = sum(self%bfGrad, dim=3)
    end if

    call self%pe%allreduce(self%bfGrad(:,:,1), &
              size(self%bfGrad(:,:,1)))

    if (self%wtDeriv .and. allocated(self%wtgrad)) then
      if (ubound(self%wtgrad,3) /= 1) &
        self%wtgrad(:,:,1) = sum(self%wtgrad, dim=3)
      call self%pe%allreduce(self%wtgrad(:,:,1), size(self%wtgrad(:,:,1)))
    end if
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_grad_t), intent(inout) :: self
    if (allocated(self%bfGrad)) deallocate(self%bfGrad)
    if (allocated(self%d1dsx)) deallocate(self%d1dsx)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    if (allocated(self%wtgrad)) deallocate(self%wtgrad)
    if (allocated(self%atxyz)) deallocate(self%atxyz)
    if (allocated(self%rijInv)) deallocate(self%rijInv)
    if (allocated(self%dummyAtom)) deallocate(self%dummyAtom)
  end subroutine

!-------------------------------------------------------------------------------
!> @brief Adjust internal memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
 subroutine resetGradPointers(self, xce, tmp, myThread)
    class(xc_consumer_grad_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: tmp(:,:)
    integer, intent(in) :: myThread

!   pruned AOs or no pruned AOs
    associate ( numAOs => xce%numAOs_p &  ! number of pruned AOs
      )
      tmp(1:numAOs, 1:3) => self%tmp_(1:numAOs*3, myThread)
    end associate

 end subroutine

!-------------------------------------------------------------------------------

 subroutine update(self, xce, mythread)

    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread, i
    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread)

    associate ( bfgrad  => self%bfgrad(:,:,mythread) &
              , d1dsx   => self%d1dsx(:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoG2    => xce%aoG2 &
              , moVA    => xce%moVA &
              , moVB    => xce%moVB &
              , moG1A   => xce%moG1A &
              , moG1B   => xce%moG1B &
              , hasBeta => xce%hasBeta &
              , numPts  => xce%numPts &
              , xc      => xce%XCLib &
              , drho    => xce%xclib%drho  &
              , ids     => xce%XCLib%ids &
              , d1ds    => xce%XCLib%d1ds &
              , d1dr    => xce%XCLib%d1dr &
              , d1dt    => xce%XCLib%d1dt &
      )
      tmpGrad = 0.0d0

!     LDA gradient
      call compAtGradRho(tmpGrad, d1dr(1,:), moVA, aoG1, numPts)

!     GGA gradient
      if (xce%funTyp /= OQP_FUNTYP_LDA) then
          do i = 1, numPts
            d1dsx(i,1:3) = 2*d1ds(ids%ga,i)*drho(1:3,i)+d1ds(ids%gc,i)*drho(4:6,i)
          end do

          call compAtGradDRho(tmpGrad, d1dsx, moVA, moG1A, aoG1, aoG2, numPts)
      end if

!     Meta-GGA gradient
      if (xce%funTyp == OQP_FUNTYP_MGGA) &
        call compAtGradTau(tmpGrad, d1dt(1,:), moG1A, aoG2, numPts)


      if (hasBeta) then
        call compAtGradRho(tmpGrad, d1dr(2,:), moVB, aoG1, numPts)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
            do i = 1, numPts
              d1dsx(i,1:3) = 2*d1ds(ids%gb,i)*drho(4:6,i)+d1ds(ids%gc,i)*drho(1:3,i)
            end do
            call compAtGradDRho(tmpGrad, d1dsx, moVB, moG1B, aoG1, aoG2, numPts)
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) &
          call compAtGradTau(tmpGrad, d1dt(2,:), moG1B, aoG2, numPts)
      end if

    end associate

!   Optional Becke/SSF quadrature-weight derivative contribution
    if (self%wtDeriv) call self%compWtGrad(xce, mythread)
 end subroutine

 subroutine postUpdate(self, xce, mythread)

    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread)

    associate ( numAOs  => xce%numAOs_p &  ! number of pruned AOs
              , indices => xce%indices_p &
      )

      if (xce%skip_p) then
        self%bfGrad(:,:,myThread) = self%bfGrad(:,:,myThread) + tmpGrad
      else
        self%bfGrad(indices(1:numAOs), :, mythread) = &
          self%bfGrad(indices(1:numAOs), :, mythread) + tmpGrad(1:numAOs, :)
      end if

   end associate

 end subroutine

!-------------------------------------------------------------------------------

!> @brief Becke/SSF quadrature-weight derivative contribution to the gradient
!>
!> @details For a fuzzy-cell partition the total weight of a grid point owned
!>  by atom A is  w = w_ar * W_A(r),  where w_ar is the (rigid) radial*angular
!>  weight and  W_A = P_A / sum_k P_k  is the normalized cell function, with
!>  P_k = prod_{l/=k} p(mu_kl),  mu_kl = (r_k - r_l)/R_kl,  r_k = |r - R_k|.
!>  The XC energy is  E = sum_g w_ar(g) W_A(r_g) eps(r_g)  with eps the XC
!>  energy density (XCLib%exc). The missing nuclear-gradient term is
!>      dE^w/dR_C = sum_g eps(g) w_ar(g) dW_A/dR_C .
!>  Holding the quadrature point fixed in space, dW_A/dR_C is evaluated for
!>  every center C/=A and, by translational invariance (sum_C dW_A/dR_C = 0),
!>  the owner atom A receives minus the sum of those contributions. This pairs
!>  with the integrand (basis-function) term already accumulated in bfGrad.
!>
!> @note EXPERIMENTAL / opt-in (dftgrid.weight_derivatives). Consistent with
!>  the SSF-type fuzzy cell (dft_bfc_algo=0); Becke surface shifting (aij) is
!>  not included here. Must be validated against finite differences before use
!>  in production. Off by default.
!> @param[in]    xce       XC engine (provides exc, xyzw, numPts, curAtom)
!> @param[in]    mythread  OpenMP thread index
 subroutine compWtGrad(self, xce, mythread)
    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread

    integer :: p, a, c, k, l, nat
    real(kind=fp) :: w, eps, z, wa, w_ar, prefac, mu, pf, qf
    real(kind=fp) :: r(3), dmu(3), nhat(3), dpa(3), dz(3), gc(3)
    real(kind=fp), allocatable :: ri(:), cells(:), uhat(:,:)
    real(kind=fp), parameter :: tiny_ = 1.0e-12_fp

    nat = self%nat
    a = xce%curAtom
    if (a < 1 .or. a > nat) return
    if (self%dummyAtom(a)) return

    allocate(ri(nat), cells(nat), uhat(3,nat))

    associate (exc => xce%XCLib%exc, npts => xce%numPts, xyzw => xce%xyzw)

    do p = 1, npts
      w = xyzw(p,4)
      if (w == 0.0_fp) cycle
      eps = exc(p)
      r = xyzw(p,1:3)

!     Point-to-atom distances and unit vectors d|r-R_k|/dr = (r-R_k)/r_k
      do k = 1, nat
        if (self%dummyAtom(k)) then
          ri(k) = 0.0_fp
          uhat(:,k) = 0.0_fp
          cycle
        end if
        uhat(:,k) = r - self%atxyz(:,k)
        ri(k) = norm2(uhat(:,k))
        if (ri(k) > tiny_) uhat(:,k) = uhat(:,k)/ri(k)
      end do

!     Unnormalized cell functions P_k (same convention as do_bfc)
      cells = 1.0_fp
      where (self%dummyAtom) cells = 0.0_fp
      do k = 2, nat
        if (self%dummyAtom(k)) cycle
        do l = 1, k-1
          if (self%dummyAtom(l)) cycle
          mu = (ri(k)-ri(l))*self%rijInv(l,k)
          pf = self%partfunc%eval(mu)
          cells(k) = cells(k)*abs(pf)
          cells(l) = cells(l)*abs(1.0_fp-pf)
        end do
      end do

      z = sum(cells)
      if (z <= tiny_) cycle
      wa = cells(a)/z
      if (wa <= tiny_) cycle
      w_ar = w/wa
      prefac = eps*w_ar

!     dW_A/dR_C = (dP_A/dR_C - W_A dZ/dR_C)/Z for each center C /= A
      do c = 1, nat
        if (c == a) cycle
        if (self%dummyAtom(c)) cycle

!       dP_A/dR_C : only the pair (A,C) depends on R_C
        dpa = 0.0_fp
        mu = (ri(a)-ri(c))*self%rijInv(a,c)
        pf = self%partfunc%eval(mu)
        if (abs(pf) > tiny_) then
          qf = self%partfunc%deriv(mu)/pf
          nhat = (self%atxyz(:,a)-self%atxyz(:,c))*self%rijInv(a,c)
          dmu = uhat(:,c)*self%rijInv(a,c) + (mu*self%rijInv(a,c))*nhat
          dpa = cells(a)*qf*dmu
        end if

!       dZ/dR_C = sum_k dP_k/dR_C ; nonzero only for pairs containing C
        dz = 0.0_fp
!       k = C : pairs (C,l)
        do l = 1, nat
          if (l == c .or. self%dummyAtom(l)) cycle
          mu = (ri(c)-ri(l))*self%rijInv(c,l)
          pf = self%partfunc%eval(mu)
          if (abs(pf) <= tiny_) cycle
          qf = self%partfunc%deriv(mu)/pf
          nhat = (self%atxyz(:,c)-self%atxyz(:,l))*self%rijInv(c,l)
          dmu = -uhat(:,c)*self%rijInv(c,l) - (mu*self%rijInv(c,l))*nhat
          dz = dz + cells(c)*qf*dmu
        end do
!       k /= C : pairs (k,C)
        do k = 1, nat
          if (k == c .or. self%dummyAtom(k)) cycle
          mu = (ri(k)-ri(c))*self%rijInv(k,c)
          pf = self%partfunc%eval(mu)
          if (abs(pf) <= tiny_) cycle
          qf = self%partfunc%deriv(mu)/pf
          nhat = (self%atxyz(:,k)-self%atxyz(:,c))*self%rijInv(k,c)
          dmu = uhat(:,c)*self%rijInv(k,c) + (mu*self%rijInv(k,c))*nhat
          dz = dz + cells(k)*qf*dmu
        end do

        gc = (dpa - wa*dz)/z
        self%wtgrad(:,c,mythread) = self%wtgrad(:,c,mythread) + prefac*gc
        self%wtgrad(:,a,mythread) = self%wtgrad(:,a,mythread) - prefac*gc
      end do

    end do

    end associate

    deallocate(ri, cells, uhat)
 end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute grid XC contribution to the nuclear gradient
!> @note  By default weight derivatives are not applied here (good enough for
!>  fine, unpruned grids). For pruned/coarse grids the Becke/SSF weight
!>  derivative term can be enabled via dftgrid.weight_derivatives, which adds
!>  compWtGrad() above (EXPERIMENTAL; validate against finite differences).
!> @param[in]    da        density matrix, alpha-spin
!> @param[in]    db        density matrix, beta-spin
!> @param[inout] dedft     nuclear gradient
!> @param[out]   totele    electronic denisty integral
!> @param[out]   totkin    kinetic energy integral
!> @param[in]    mxAngMom  max. needed ang. mom. value (incl. derivatives)
!> @param[in]    nbf        basis set size
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    urohf     .TRUE. if open-shell calculation
!> @author Vladimir Mironov
  subroutine derexc_blk(basis, molGrid, da, db, dedft, &
                        totele, totkin, &
                        mxAngMom, nbf, dft_threshold, urohf, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid

    type(basis_set) :: basis
    logical, intent(IN) :: urohf
    integer, intent(IN) :: mxAngMom, nbf
    real(KIND=fp), intent(INOUT) :: totele, totkin
    real(KIND=fp), intent(INOUT) :: da(nbf, *), db(nbf, *), dedft(:, :)
    real(kind=fp), intent(in) :: dft_threshold

    type(xc_consumer_grad_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: j

    integer :: nat
    integer :: ia, ja
    real(kind=fp) :: dist

    real(KIND=fp), target, allocatable :: da2(:, :), db2(:, :)


    nat = infos%mol_prop%natom

    allocate (da2(nbf, nbf))
    do j = 1, nbf
      da2(:, j) = da(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
    end do
    if (urohf) then
      allocate (db2(nbf, nbf))
      do j = 1, nbf
        db2(:, j) = db(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
      end do
    end if

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = urohf
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = mxAngMom
    xc_opts%nDer = 1
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => da2
    xc_opts%wfBeta => db2
    xc_opts%dft_threshold = dft_threshold
    xc_opts%molGrid => molGrid

!   Optional: set up Becke/SSF quadrature-weight derivative contribution
    dat%wtDeriv = infos%dft%dft_wt_der
    if (dat%wtDeriv) then
      dat%nat = nat
      allocate(dat%atxyz(3, nat))
      dat%atxyz(1:3, 1:nat) = basis%atoms%xyz(1:3, 1:nat)
      allocate(dat%dummyAtom(nat), source=.false.)
      if (allocated(molGrid%dummyAtom)) then
        if (size(molGrid%dummyAtom) >= nat) &
          dat%dummyAtom(1:nat) = molGrid%dummyAtom(1:nat)
      end if
      allocate(dat%rijInv(nat, nat), source=0.0_fp)
      do ia = 1, nat
        do ja = 1, ia-1
          dist = norm2(dat%atxyz(:,ia) - dat%atxyz(:,ja))
          if (dist > 0.0_fp) then
            dat%rijInv(ia, ja) = 1.0_fp/dist
            dat%rijInv(ja, ia) = dat%rijInv(ia, ja)
          end if
        end do
      end do
      call dat%partfunc%set(infos%dft%dft_partfun)
    end if

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    totele = dat%N_elec
    totkin = dat%E_kin

    do j = 1, basis%nshell
      associate (atom => basis%origin(j), &
                 offset => basis%ao_offset(j), &
                 naos => basis%naos(j))
        dedft(1, atom) = dedft(1, atom)-sum(dat%bfGrad(offset:offset+naos-1, 1, 1))
        dedft(2, atom) = dedft(2, atom)-sum(dat%bfGrad(offset:offset+naos-1, 2, 1))
        dedft(3, atom) = dedft(3, atom)-sum(dat%bfGrad(offset:offset+naos-1, 3, 1))
      end associate
    end do

!   Add the Becke/SSF quadrature-weight derivative term (if enabled)
    if (dat%wtDeriv .and. allocated(dat%wtgrad)) then
      do ia = 1, nat
        dedft(1:3, ia) = dedft(1:3, ia) + dat%wtgrad(1:3, ia, 1)
      end do
    end if

    deallocate (da2)
    if (urohf) deallocate (db2)

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

end module mod_dft_gridint_grad
