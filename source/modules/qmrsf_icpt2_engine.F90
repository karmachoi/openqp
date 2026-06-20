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
  use qmrsf_icpt2_downfold_mod, only: icpt2_eff_hamiltonian, icpt2_safe_inv
  implicit none
  private
  public :: qmrsf_icpt2_dress             !< brute-force (full det space; small windows)
  public :: qmrsf_icpt2_dress_contracted  !< contracted (no FCI list; production)
  public :: qmrsf_icpt2_count_perturbers  !< feasibility estimate for the contracted path

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

!-------------------------------------------------------------------------------
! CONTRACTED external-Q engine: never builds the FCI determinant list. Enumerates
! perturbers by per-spin (active-residual x virtual) blocks, contracts each
! coupling against the CAS eigenvectors (Slater-Condon, rank-pruned), and STREAMS
! the des-Cloizeaux dressing into H_eff. Port of the validated standalone
! tools/qmrsf_pathways_proto/fortran/qmrsf_icpt2_contracted.f90 (PASS <1e-13).
!-------------------------------------------------------------------------------

  !> @brief Total external-Q perturber count for a frozen-core window (feasibility).
  integer(8) function qmrsf_icpt2_count_perturbers(norb, na, nb, nact) result(nQ)
    integer, intent(in) :: norb, na, nb, nact
    integer :: nvirt
    nvirt = norb - nact
    nQ = spin_choices(na, nact, nvirt) * spin_choices(nb, nact, nvirt) &
         - int(ncomb(nact,na),8)*int(ncomb(nact,nb),8)
  contains
    integer(8) function spin_choices(ns, na_, nv_) result(c)
      integer, intent(in) :: ns, na_, nv_
      integer :: nv
      c = 0_8
      do nv = 0, ns
        if (ns-nv <= na_ .and. nv <= nv_) c = c + int(ncomb(na_,ns-nv),8)*int(ncomb(nv_,nv),8)
      end do
    end function
  end function qmrsf_icpt2_count_perturbers

  subroutine qmrsf_icpt2_dress_contracted(h, eri, eps, shift, norb, na, nb, nact, nPd, &
                                          eP, edr_en, edr_dy, nQ, herm)
    integer,  intent(in)  :: norb, na, nb, nact, nPd
    real(dp), intent(in)  :: h(norb,norb), eri(norb,norb,norb,norb), eps(norb)
    real(dp), intent(in)  :: shift           !< EN imaginary level shift (Eh); 0 = off
    real(dp), intent(out) :: eP(nPd), edr_en(nPd), edr_dy(nPd)
    integer,  intent(out) :: nQ
    real(dp), intent(out) :: herm

    integer :: nso, nelec, nvirt, nP, na_c, nb_c, i, j, k, l, t, p, ierr
    integer :: nblkA, nblkB, ia, ib, ndiff, cc
    real(dp), allocatable :: H1(:,:), g(:,:,:,:), HPP(:,:), cP(:,:), ePall(:)
    real(dp), allocatable :: Heff_en(:,:), Heff_dy(:,:), melv(:), cvec(:), inv_en(:)
    integer,  allocatable :: Pdets(:,:), acomb(:,:), bcomb(:,:)
    integer,  allocatable :: blkA(:,:), blkB(:,:), nvA(:), nvB(:), qd(:)
    real(dp) :: hqq, sv, dyd, hd

    nso = 2*norb; nelec = na+nb; nvirt = norb - nact
    allocate(H1(nso,nso), g(nso,nso,nso,nso))
    call build_spinorb(h, eri, norb, nso, H1, g)

    ! CAS(4,4) P determinants (na alpha + nb beta in orbitals 1..nact) + eigenpairs
    na_c = ncomb(nact,na); nb_c = ncomb(nact,nb); nP = na_c*nb_c
    allocate(acomb(na,na_c), bcomb(nb,nb_c), Pdets(nelec,nP))
    call all_combs(nact, na, acomb, na_c)
    call all_combs(nact, nb, bcomb, nb_c)
    cc = 0
    do i = 1, na_c; do j = 1, nb_c
      cc = cc + 1
      do t = 1, na; Pdets(t,cc)    = acomb(t,i);        end do
      do t = 1, nb; Pdets(na+t,cc) = bcomb(t,j) + norb; end do
      call isort(Pdets(:,cc), nelec)
    end do; end do

    allocate(HPP(nP,nP), cP(nP,nP), ePall(nP))
    do i = 1, nP; do j = 1, nP; HPP(i,j) = melem(Pdets(:,i),Pdets(:,j),H1,g,nelec,nso); end do; end do
    cP = HPP
    call diag_symm_full(0, nP, cP, nP, ePall, ierr)
    eP = ePall(1:nPd)

    call gen_spin_blocks(na, 0,    nact, nvirt, blkA, nvA, nblkA)
    call gen_spin_blocks(nb, norb, nact, nvirt, blkB, nvB, nblkB)

    allocate(Heff_en(nPd,nPd), Heff_dy(nPd,nPd), melv(nP), cvec(nPd), inv_en(nPd), qd(nelec))
    Heff_en = 0.0_dp; Heff_dy = 0.0_dp; nQ = 0
    do ia = 1, nblkA
      do ib = 1, nblkB
        if (nvA(ia)+nvB(ib) == 0) cycle
        nQ = nQ + 1
        do t = 1, na; qd(t)    = blkA(t,ia); end do
        do t = 1, nb; qd(na+t) = blkB(t,ib); end do
        call isort(qd, nelec)
        do j = 1, nP
          ndiff = 0
          do t = 1, nelec; if (.not. any(Pdets(:,j)==qd(t))) ndiff = ndiff + 1; end do
          if (ndiff > 2) then; melv(j) = 0.0_dp
          else;                melv(j) = melem(qd, Pdets(:,j), H1, g, nelec, nso); end if
        end do
        cvec = matmul(melv, cP(:,1:nPd))
        hqq = melem(qd, qd, H1, g, nelec, nso)
        sv = 0.0_dp
        do t = 1, nelec; p = mod(qd(t)-1,norb)+1; if (p > nact) sv = sv + eps(p); end do
        do k = 1, nPd; inv_en(k) = en_inv(eP(k)-hqq, shift); end do
        dyd = icpt2_safe_inv(-sv)
        do k = 1, nPd
          do l = 1, nPd
            Heff_en(k,l) = Heff_en(k,l) + 0.5_dp*cvec(k)*cvec(l)*(inv_en(k)+inv_en(l))
            Heff_dy(k,l) = Heff_dy(k,l) + cvec(k)*cvec(l)*dyd
          end do
        end do
      end do
    end do
    do k = 1, nPd; Heff_en(k,k) = Heff_en(k,k) + eP(k); Heff_dy(k,k) = Heff_dy(k,k) + eP(k); end do

    herm = 0.0_dp
    do k = 1, nPd; do l = 1, nPd; herm = max(herm, abs(Heff_en(k,l)-Heff_en(l,k))); end do; end do
    Heff_en = 0.5_dp*(Heff_en + transpose(Heff_en))
    Heff_dy = 0.5_dp*(Heff_dy + transpose(Heff_dy))
    call diag_symm_full(0, nPd, Heff_en, nPd, edr_en, ierr)
    call diag_symm_full(0, nPd, Heff_dy, nPd, edr_dy, ierr)
    hd = herm
  end subroutine qmrsf_icpt2_dress_contracted

  !> Epstein-Nesbet inverse denominator with optional imaginary level shift
  !> (Forsberg-Malmqvist): Re[1/(d + i*beta)] = d/(d^2 + beta^2). beta=0 -> 1/d
  !> (intruder-guarded). beta>0 smoothly damps near-singular intruder denominators.
  real(dp) function en_inv(d, beta) result(r)
    real(dp), intent(in) :: d, beta
    if (beta > 0.0_dp) then
      r = d / (d*d + beta*beta)
    else
      r = icpt2_safe_inv(d)
    end if
  end function en_inv

  !> per-spin perturber blocks: place n_sig electrons as (n_sig-nv) active + nv virtual.
  subroutine gen_spin_blocks(n_sig, base, nact, nvirt, blk_occ, blk_nv, nblk)
    integer, intent(in)  :: n_sig, base, nact, nvirt
    integer, allocatable, intent(out) :: blk_occ(:,:), blk_nv(:)
    integer, intent(out) :: nblk
    integer :: nv, na_act, nac, nvc, ia, iv, t, maxblk
    integer, allocatable :: acomb(:,:), vcomb(:,:)
    maxblk = 0
    do nv = 0, n_sig
      na_act = n_sig - nv
      if (na_act > nact .or. nv > nvirt) cycle
      maxblk = maxblk + ncomb(nact,na_act)*ncomb(nvirt,nv)
    end do
    allocate(blk_occ(max(n_sig,1),maxblk), blk_nv(maxblk))
    nblk = 0
    do nv = 0, n_sig
      na_act = n_sig - nv
      if (na_act > nact .or. nv > nvirt) cycle
      nac = ncomb(nact,na_act); nvc = ncomb(nvirt,nv)
      allocate(acomb(na_act,nac), vcomb(nv,nvc))
      call all_combs(nact, na_act, acomb, nac)
      call all_combs(nvirt, nv, vcomb, nvc)
      do ia = 1, nac
        do iv = 1, nvc
          nblk = nblk + 1
          do t = 1, na_act; blk_occ(t,nblk)        = base + acomb(t,ia);        end do
          do t = 1, nv;     blk_occ(na_act+t,nblk) = base + nact + vcomb(t,iv); end do
          blk_nv(nblk) = nv
        end do
      end do
      deallocate(acomb, vcomb)
    end do
  end subroutine gen_spin_blocks

end module qmrsf_icpt2_engine_mod
