!> @file qmrsf_icpt2_engine.F90
!> @brief QMRSF-icPT2 determinant-CI perturber engine (brute force) + downfold driver.
!>
!> @details
!>   Given the (frozen-core-dressed) MO integrals over the window of active+virtual
!>   orbitals, this builds the full determinant space, partitions it into the CAS
!>   model space P (all electrons confined to the first `nact` orbitals = the
!>   quintet SOMOs) and the external space Q, assembles H_PP / H_QP / diag(H)_Q by
!>   the Slater-Condon rules, diagonalizes the CAS, internally contracts the
!>   couplings, and applies the validated des-Cloizeaux multistate downfold
!>   (qmrsf_icpt2_downfold_mod) with Epstein-Nesbet denominators.
!>
!>   This is the BRUTE-FORCE perturber generation: the external space is built
!>   explicitly, so it is correct for small windows (the validation systems H4,
!>   small CBD) and is the bridge to the production contracted engine. It is a
!>   faithful port of the validated NumPy prototype and the standalone
!>   tools/qmrsf_pathways_proto/fortran/qmrsf_icpt2_full.f90 (matched 5.7e-14).
!>
!>   The window is frozen-core: the truly inactive doubly-occupied orbitals are
!>   folded into h via the core mean field (qmrsf_ao2mo) before this routine, so
!>   `norb` here is (nbf - ncore) and the four SOMOs are orbitals 1..nact.
module qmrsf_icpt2_engine_mod
  use precision, only: dp
  use eigen, only: diag_symm_full
  use qmrsf_icpt2_downfold_mod, only: icpt2_eff_hamiltonian
  implicit none
  private
  public :: qmrsf_icpt2_dress

contains

  !> @brief Build det-CI external space and return the icPT2-dressed energies.
  !> @param[in]  h(norb,norb)            window MO 1e integrals (frozen-core dressed)
  !> @param[in]  eri(norb^4)             window MO 2e integrals, CHEMIST (pq|rs)
  !> @param[in]  eps(norb)               window MO (ROHF) orbital energies (Dyall H0)
  !> @param[in]  norb,na,nb              window orbitals and alpha/beta electron counts
  !> @param[in]  nact                    CAS active orbitals (=4); P = electrons in 1..nact
  !> @param[in]  nPd                     number of dressed roots to return
  !> @param[out] eP(nPd)                 bare CAS roots (lowest nPd)
  !> @param[out] edr_en(nPd)             icPT2 Epstein-Nesbet dressed energies
  !> @param[out] edr_dy(nPd)             icPT2 Dyall dressed energies
  !> @param[out] ndet,nP,nQ             determinant counts
  !> @param[out] herm                    H_eff Hermiticity residual
  !>
  !> @details Dyall H0: for the frozen-core window (no inactive/core orbitals inside
  !> the window) the reference has no occupied non-active orbitals, so the only
  !> external excitations are active->virtual and the Dyall zeroth-order denominator
  !> reduces to  d_q = -(sum of MO energies of the virtual orbitals occupied in q)
  !> (the active part is treated exactly and cancels against the CAS reference). This
  !> is exactly the NumPy prototype's dyall_denoms restricted to the frozen-core case.
  subroutine qmrsf_icpt2_dress(h, eri, eps, norb, na, nb, nact, nPd, eP, &
                               edr_en, edr_dy, ndet, nP, nQ, herm)
    integer,  intent(in)  :: norb, na, nb, nact, nPd
    real(dp), intent(in)  :: h(norb,norb), eri(norb,norb,norb,norb), eps(norb)
    real(dp), intent(out) :: eP(nPd), edr_en(nPd), edr_dy(nPd)
    integer,  intent(out) :: ndet, nP, nQ
    real(dp), intent(out) :: herm

    integer :: nso, nelec, nalp, nbet, i, j, k, l, q, ca, cb, p, ierr
    real(dp), allocatable :: H1(:,:), g(:,:,:,:)
    integer,  allocatable :: astr(:,:), bstr(:,:), dets(:,:), Pid(:), Qid(:)
    real(dp), allocatable :: HPP(:,:), HQP(:,:), Hqq(:), cP(:,:), ePall(:)
    real(dp), allocatable :: coup(:,:), invd(:,:), sumv(:)
    real(dp) :: d, hdum
    logical :: inP

    nso = 2*norb; nelec = na + nb
    allocate(H1(nso,nso), g(nso,nso,nso,nso))
    call build_spinorb(h, eri, norb, nso, H1, g)

    nalp = ncomb(norb, na); nbet = ncomb(norb, nb)
    ndet = nalp*nbet
    allocate(astr(na,nalp), bstr(nb,nbet), dets(nelec,ndet))
    call all_combs(norb, na, astr, nalp)
    call all_combs(norb, nb, bstr, nbet)
    l = 0
    do ca = 1, nalp
      do cb = 1, nbet
        l = l + 1
        do i = 1, na; dets(i,l)    = astr(i,ca);          end do
        do i = 1, nb; dets(na+i,l) = bstr(i,cb) + norb;   end do
        call isort(dets(:,l), nelec)
      end do
    end do

    ! P = CAS: no occupied spin-orbital has spatial index > nact (all electrons in 1..nact)
    allocate(Pid(ndet), Qid(ndet)); nP = 0; nQ = 0
    do l = 1, ndet
      inP = .true.
      do i = 1, nelec
        p = mod(dets(i,l)-1, norb) + 1
        if (p > nact) then; inP = .false.; exit; end if
      end do
      if (inP) then; nP = nP + 1; Pid(nP) = l
      else;          nQ = nQ + 1; Qid(nQ) = l; end if
    end do

    allocate(HPP(nP,nP), HQP(nQ,nP), Hqq(nQ))
    do i = 1, nP; do j = 1, nP
      HPP(i,j) = melem(dets(:,Pid(i)), dets(:,Pid(j)), H1, g, nelec, nso)
    end do; end do
    do i = 1, nQ; do j = 1, nP
      HQP(i,j) = melem(dets(:,Qid(i)), dets(:,Pid(j)), H1, g, nelec, nso)
    end do; end do
    do i = 1, nQ
      Hqq(i) = melem(dets(:,Qid(i)), dets(:,Qid(i)), H1, g, nelec, nso)
    end do

    allocate(ePall(nP), cP(nP,nP)); cP = HPP
    call diag_symm_full(0, nP, cP, nP, ePall, ierr)
    eP = ePall(1:nPd)

    if (nQ == 0) then
      edr_en = eP; edr_dy = eP; herm = 0.0_dp
      return
    end if

    allocate(coup(nQ,nPd), invd(nQ,nPd), sumv(nQ))
    coup = matmul(HQP, cP(:,1:nPd))

    ! Dyall sum of occupied-virtual MO energies per Q determinant
    do q = 1, nQ
      sumv(q) = 0.0_dp
      do i = 1, nelec
        p = mod(dets(i,Qid(q))-1, norb) + 1
        if (p > nact) sumv(q) = sumv(q) + eps(p)
      end do
    end do

    ! Epstein-Nesbet
    do k = 1, nPd
      do q = 1, nQ
        d = eP(k) - Hqq(q)
        if (abs(d) < 1.0e-6_dp) d = sign(1.0e-6_dp, d) + 1.0e-30_dp
        invd(q,k) = 1.0_dp / d
      end do
    end do
    call icpt2_eff_hamiltonian(nPd, nQ, eP, coup, invd, edr_en, herm)

    ! Dyall (frozen-core window: d_q = -sum_v, root-independent)
    do k = 1, nPd
      do q = 1, nQ
        d = -sumv(q)
        if (abs(d) < 1.0e-6_dp) d = sign(1.0e-6_dp, d) + 1.0e-30_dp
        invd(q,k) = 1.0_dp / d
      end do
    end do
    call icpt2_eff_hamiltonian(nPd, nQ, eP, coup, invd, edr_dy, hdum)
  end subroutine qmrsf_icpt2_dress

!-------------------------------------------------------------------------------

  subroutine build_spinorb(h, eri, norb, nso, H1, g)
    integer,  intent(in)  :: norb, nso
    real(dp), intent(in)  :: h(norb,norb), eri(norb,norb,norb,norb)
    real(dp), intent(out) :: H1(nso,nso), g(nso,nso,nso,nso)
    integer :: P,Q,R,S, sp,sq,sr,ss, ap,aq,ar,as
    real(dp) :: a, b
    H1 = 0.0_dp
    do P = 1, nso; do Q = 1, nso
      if ((P-1)/norb == (Q-1)/norb) H1(P,Q) = h(mod(P-1,norb)+1, mod(Q-1,norb)+1)
    end do; end do
    g = 0.0_dp
    do P = 1, nso; sp=(P-1)/norb; ap=mod(P-1,norb)+1
    do Q = 1, nso; sq=(Q-1)/norb; aq=mod(Q-1,norb)+1
    do R = 1, nso; sr=(R-1)/norb; ar=mod(R-1,norb)+1
    do S = 1, nso; ss=(S-1)/norb; as=mod(S-1,norb)+1
      a = 0.0_dp; b = 0.0_dp
      if (sp==sr .and. sq==ss) a = eri(ap,ar,aq,as)
      if (sp==ss .and. sq==sr) b = eri(ap,as,aq,ar)
      g(P,Q,R,S) = a - b
    end do; end do; end do; end do
  end subroutine build_spinorb

  real(dp) function melem(D1, D2, H1, g, nelec, nso)
    integer,  intent(in) :: nelec, nso
    integer,  intent(in) :: D1(nelec), D2(nelec)
    real(dp), intent(in) :: H1(nso,nso), g(nso,nso,nso,nso)
    integer :: holes(nelec), parts(nelec), common(nelec), nh, np, nc
    integer :: occ(nelec), nocc, i, k, idx, cnt
    integer :: p1,p2,ho1,ho2, pp, hh
    real(dp) :: sgn, val, e
    nh=0; np=0; nc=0
    do i=1,nelec
      if (.not. any(D1==D2(i))) then; nh=nh+1; holes(nh)=D2(i); end if
    end do
    do i=1,nelec
      if (.not. any(D2==D1(i))) then; np=np+1; parts(np)=D1(i)
      else;                            nc=nc+1; common(nc)=D1(i); end if
    end do
    if (nh > 2) then; melem = 0.0_dp; return; end if
    occ = D2; nocc = nelec; sgn = 1.0_dp
    do k=1,nh
      idx = 0
      do i=1,nocc
        if (occ(i)==holes(k)) then; idx=i; exit; end if
      end do
      if (mod(idx-1,2)==1) sgn = -sgn
      do i=idx,nocc-1; occ(i)=occ(i+1); end do
      nocc = nocc-1
    end do
    do k=np,1,-1
      cnt = 0
      do i=1,nocc
        if (occ(i) < parts(k)) cnt = cnt+1
      end do
      if (mod(cnt,2)==1) sgn = -sgn
      do i=nocc,cnt+1,-1; occ(i+1)=occ(i); end do
      occ(cnt+1)=parts(k); nocc=nocc+1
    end do
    if (nh==0) then
      e = 0.0_dp
      do i=1,nelec; e = e + H1(D1(i),D1(i)); end do
      do i=1,nelec-1
        do k=i+1,nelec; e = e + g(D1(i),D1(k),D1(i),D1(k)); end do
      end do
      melem = e
    else if (nh==1) then
      pp = parts(1); hh = holes(1)
      val = H1(pp,hh)
      do i=1,nc; val = val + g(pp,common(i),hh,common(i)); end do
      melem = sgn*val
    else
      p1=parts(1); p2=parts(2); ho1=holes(1); ho2=holes(2)
      melem = sgn*g(p1,p2,ho1,ho2)
    end if
  end function melem

  integer function ncomb(m,k)
    integer, intent(in) :: m,k
    integer :: i
    ncomb = 1
    do i=1,k; ncomb = ncomb*(m-k+i)/i; end do
  end function ncomb

  subroutine all_combs(m, k, out, ncol)
    integer, intent(in)  :: m, k, ncol
    integer, intent(out) :: out(k,ncol)
    integer :: cc(k), j, col
    do j=1,k; cc(j)=j; end do
    col = 0
    do
      col = col + 1; out(:,col) = cc
      if (.not. nextc(cc,k,m)) exit
    end do
  end subroutine all_combs

  logical function nextc(cc,k,m)
    integer, intent(inout) :: cc(k)
    integer, intent(in) :: k,m
    integer :: i,j
    i=k
    do while (i>=1)
      if (cc(i) /= m-k+i) exit
      i=i-1
    end do
    if (i<1) then; nextc=.false.; return; end if
    cc(i)=cc(i)+1
    do j=i+1,k; cc(j)=cc(j-1)+1; end do
    nextc=.true.
  end function nextc

  subroutine isort(a,n)
    integer, intent(inout) :: a(n)
    integer, intent(in) :: n
    integer :: i,j,t
    do i=2,n
      t=a(i); j=i-1
      do while (j>=1)
        if (a(j)<=t) exit
        a(j+1)=a(j); j=j-1
      end do
      a(j+1)=t
    end do
  end subroutine isort

end module qmrsf_icpt2_engine_mod
