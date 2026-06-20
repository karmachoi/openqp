! QMRSF-icPT2 CONTRACTED external-Q engine (standalone, build-testable).
! Faithful Fortran port of the validated NumPy contracted prototype
! (tools/qmrsf_pathways_proto/qmrsf_icpt2_contracted_proto.py). It never builds the
! FCI determinant list: it enumerates external-Q perturbers INTRINSICALLY by per-spin
! (active-residual subset) x (virtual subset) blocks, contracts each coupling against
! the 36 CAS eigenvectors via Slater-Condon (pruned by excitation rank <=2), and
! STREAMS the des-Cloizeaux dressing directly into H_eff (no nQ x nPd coupling stored).
!
! Reads stageB/qmrsf_icpt2_full_live.dat (the live H4/6-31G window) and reproduces the
! dumped EN and Dyall dressed spectra. Frozen-core window: NO core orbitals inside it,
! 4 active electrons (na=nb=2), virtuals = orbitals nact+1..norb.
module icpt2_con_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: norb, nso, nelec
contains
  pure integer function spat(P)
    integer, intent(in) :: P
    spat = mod(P-1, norb) + 1
  end function
  pure integer function spn(P)
    integer, intent(in) :: P
    spn = (P-1) / norb
  end function

  subroutine build_spinorb(h_mo, eri, H1, g)
    real(dp), intent(in)  :: h_mo(norb,norb), eri(norb,norb,norb,norb)
    real(dp), intent(out) :: H1(nso,nso), g(nso,nso,nso,nso)
    integer :: P,Q,R,S
    real(dp) :: a,b
    H1 = 0.0_dp
    do P=1,nso; do Q=1,nso
      if (spn(P)==spn(Q)) H1(P,Q) = h_mo(spat(P),spat(Q))
    end do; end do
    g = 0.0_dp
    do P=1,nso; do Q=1,nso; do R=1,nso; do S=1,nso
      a=0.0_dp; b=0.0_dp
      if (spn(P)==spn(R) .and. spn(Q)==spn(S)) a = eri(spat(P),spat(R),spat(Q),spat(S))
      if (spn(P)==spn(S) .and. spn(Q)==spn(R)) b = eri(spat(P),spat(S),spat(Q),spat(R))
      g(P,Q,R,S) = a - b
    end do; end do; end do; end do
  end subroutine

  real(dp) function melem(D1, D2, H1, g)
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
      idx=0
      do i=1,nocc; if (occ(i)==holes(k)) then; idx=i; exit; end if; end do
      if (mod(idx-1,2)==1) sgn=-sgn
      do i=idx,nocc-1; occ(i)=occ(i+1); end do
      nocc=nocc-1
    end do
    do k=np,1,-1
      cnt=0
      do i=1,nocc; if (occ(i)<parts(k)) cnt=cnt+1; end do
      if (mod(cnt,2)==1) sgn=-sgn
      do i=nocc,cnt+1,-1; occ(i+1)=occ(i); end do
      occ(cnt+1)=parts(k); nocc=nocc+1
    end do
    if (nh==0) then
      e=0.0_dp
      do i=1,nelec; e=e+H1(D1(i),D1(i)); end do
      do i=1,nelec-1; do k=i+1,nelec; e=e+g(D1(i),D1(k),D1(i),D1(k)); end do; end do
      melem=e
    else if (nh==1) then
      pp=parts(1); hh=holes(1)
      val=H1(pp,hh)
      do i=1,nc; val=val+g(pp,common(i),hh,common(i)); end do
      melem=sgn*val
    else
      p1=parts(1); p2=parts(2); ho1=holes(1); ho2=holes(2)
      melem=sgn*g(p1,p2,ho1,ho2)
    end if
  end function

  integer function ncomb(m,k)
    integer, intent(in) :: m,k
    integer :: i
    ncomb=1
    do i=1,k; ncomb=ncomb*(m-k+i)/i; end do
  end function
  subroutine all_combs(m, k, out, ncol)
    integer, intent(in)  :: m, k, ncol
    integer, intent(out) :: out(k,ncol)
    integer :: cc(k), j, col
    do j=1,k; cc(j)=j; end do
    col=0
    do
      col=col+1; if (k>0) out(:,col)=cc
      if (.not. nextc(cc,k,m)) exit
    end do
  end subroutine
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
  end function
  subroutine isort(a,n)
    integer, intent(inout) :: a(n)
    integer, intent(in) :: n
    integer :: i,j,t
    do i=2,n
      t=a(i); j=i-1
      do while (j>=1); if (a(j)<=t) exit; a(j+1)=a(j); j=j-1; end do
      a(j+1)=t
    end do
  end subroutine

  !> per-spin perturber blocks: place n_sig electrons as (n_sig-nv) active + nv virtual.
  !> base = 0 (alpha) or norb (beta). Returns blk_occ(n_sig,nblk) and blk_nv(nblk).
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
  end subroutine

  elemental real(dp) function ginv(d)
    real(dp), intent(in) :: d
    real(dp) :: dd
    dd = d
    if (abs(dd) < 1.0e-6_dp) dd = sign(1.0e-6_dp,dd) + 1.0e-30_dp
    ginv = 1.0_dp/dd
  end function
end module icpt2_con_mod

program test_icpt2_contracted
  use icpt2_con_mod
  implicit none
  integer :: nPd, nact, nvirt, na, nb, i, j, k, l, p, q, r, s, info, lwork, ierr
  integer :: nP, nblkA, nblkB, ia, ib, t, nQ, ndiff
  real(dp), allocatable :: h_mo(:,:), eri(:,:,:,:), H1(:,:), g(:,:,:,:), eps(:)
  real(dp), allocatable :: HPP(:,:), cP(:,:), ePall(:), eP(:), refEN(:), refDY(:)
  real(dp), allocatable :: Heff_en(:,:), Heff_dy(:,:), edEN(:), edDY(:), work(:)
  real(dp), allocatable :: melv(:), cvec(:), inv_en(:)
  integer,  allocatable :: Pdets(:,:), blkA(:,:), blkB(:,:), nvA(:), nvB(:), qd(:)
  real(dp) :: ecore, hqq, sv, dyd, dEN, dDY, herm
  integer,  allocatable :: acomb(:,:), bcomb(:,:)
  integer :: cc, na_c, nb_c

  open(10, file="qmrsf_icpt2_full_live.dat", status="old", action="read")
  read(10,*) norb, nPd
  allocate(h_mo(norb,norb), eri(norb,norb,norb,norb), eps(norb))
  do i=1,norb; read(10,*) (h_mo(i,j), j=1,norb); end do
  do i=1,norb; do j=1,norb; do k=1,norb; read(10,*) (eri(i,j,k,l), l=1,norb); end do; end do; end do
  read(10,*) ecore
  read(10,*) (eps(i), i=1,norb)
  allocate(eP(nPd), refEN(nPd), refDY(nPd))
  read(10,*) (eP(i), i=1,nPd)      ! live bare CAS (overwritten below by our own)
  read(10,*) (refEN(i), i=1,nPd)
  read(10,*) (refDY(i), i=1,nPd)
  close(10)

  nso = 2*norb; na = 2; nb = 2; nelec = na+nb; nact = 4; nvirt = norb - nact
  allocate(H1(nso,nso), g(nso,nso,nso,nso))
  call build_spinorb(h_mo, eri, H1, g)

  ! CAS(4,4) P determinants: na alpha + nb beta confined to orbitals 1..nact
  na_c = ncomb(nact,na); nb_c = ncomb(nact,nb); nP = na_c*nb_c
  allocate(acomb(na,na_c), bcomb(nb,nb_c), Pdets(nelec,nP))
  call all_combs(nact, na, acomb, na_c)
  call all_combs(nact, nb, bcomb, nb_c)
  cc = 0
  do ia=1,na_c; do ib=1,nb_c
    cc = cc+1
    do t=1,na; Pdets(t,cc)    = acomb(t,ia);        end do
    do t=1,nb; Pdets(na+t,cc) = bcomb(t,ib) + norb; end do
    call isort(Pdets(:,cc), nelec)
  end do; end do

  allocate(HPP(nP,nP), cP(nP,nP), ePall(nP))
  do i=1,nP; do j=1,nP; HPP(i,j)=melem(Pdets(:,i),Pdets(:,j),H1,g); end do; end do
  cP = HPP
  lwork = 64*nP; allocate(work(lwork))
  call dsyev('V','U', nP, cP, nP, ePall, work, lwork, info)
  if (info/=0) stop 'dsyev HPP'
  eP = ePall(1:nPd)

  ! per-spin perturber blocks
  call gen_spin_blocks(na, 0,    nact, nvirt, blkA, nvA, nblkA)
  call gen_spin_blocks(nb, norb, nact, nvirt, blkB, nvB, nblkB)

  allocate(Heff_en(nPd,nPd), Heff_dy(nPd,nPd))
  allocate(melv(nP), cvec(nPd), inv_en(nPd), qd(nelec))
  Heff_en = 0.0_dp; Heff_dy = 0.0_dp; nQ = 0
  do ia=1,nblkA
    do ib=1,nblkB
      if (nvA(ia)+nvB(ib) == 0) cycle          ! all-active = P space, skip
      nQ = nQ + 1
      do t=1,na; qd(t)    = blkA(t,ia); end do
      do t=1,nb; qd(na+t) = blkB(t,ib); end do
      call isort(qd, nelec)
      ! contracted coupling: cvec(k) = sum_j melem(q,Pdet_j) cP(j,k), pruned by rank
      do j=1,nP
        ndiff = 0
        do t=1,nelec; if (.not. any(Pdets(:,j)==qd(t))) ndiff=ndiff+1; end do
        if (ndiff > 2) then; melv(j)=0.0_dp; else; melv(j)=melem(qd,Pdets(:,j),H1,g); end if
      end do
      cvec = matmul(melv, cP(:,1:nPd))
      hqq = melem(qd, qd, H1, g)
      sv = 0.0_dp
      do t=1,nelec; p=mod(qd(t)-1,norb)+1; if (p>nact) sv = sv + eps(p); end do
      ! stream into H_eff (EN: state-specific denom; Dyall: d=-sv, state-independent)
      do k=1,nPd; inv_en(k) = ginv(eP(k)-hqq); end do
      dyd = ginv(-sv)
      do k=1,nPd
        do l=1,nPd
          Heff_en(k,l) = Heff_en(k,l) + 0.5_dp*cvec(k)*cvec(l)*(inv_en(k)+inv_en(l))
          Heff_dy(k,l) = Heff_dy(k,l) + cvec(k)*cvec(l)*dyd
        end do
      end do
    end do
  end do
  do k=1,nPd; Heff_en(k,k)=Heff_en(k,k)+eP(k); Heff_dy(k,k)=Heff_dy(k,k)+eP(k); end do

  herm = 0.0_dp
  do k=1,nPd; do l=1,nPd; herm=max(herm,abs(Heff_en(k,l)-Heff_en(l,k))); end do; end do
  Heff_en = 0.5_dp*(Heff_en+transpose(Heff_en))
  Heff_dy = 0.5_dp*(Heff_dy+transpose(Heff_dy))
  allocate(edEN(nPd), edDY(nPd))
  deallocate(work); lwork=64*nPd; allocate(work(lwork))
  call dsyev('N','U', nPd, Heff_en, nPd, edEN, work, lwork, info); if (info/=0) stop 'dsyev en'
  call dsyev('N','U', nPd, Heff_dy, nPd, edDY, work, lwork, info); if (info/=0) stop 'dsyev dy'

  dEN = maxval(abs(edEN - refEN))
  dDY = maxval(abs(edDY - refDY))
  print '(a)', "==== QMRSF-icPT2 CONTRACTED engine (Fortran) vs live dump ===="
  print '(a,i0,a,i0,a,i0)', "  norb=",norb,"  nvirt=",nvirt,"  CAS nP=",nP
  print '(a,i0,a,i0,a,i0)', "  perturber blocks: alpha=",nblkA,"  beta=",nblkB,"  nQ=",nQ
  print '(a,es10.2)', "  H_eff Hermiticity  = ", herm
  print '(a,es10.2)', "  max|EN - ref|      = ", dEN
  print '(a,es10.2)', "  max|Dyall - ref|   = ", dDY
  print '(a,f16.8,a,f16.8)', "  ground EN total = ", edEN(1)+ecore, "   Dyall = ", edDY(1)+ecore
  if (dEN < 1.0d-9 .and. dDY < 1.0d-9) then
     print '(a)', "  RESULT: PASS  (contracted engine reproduces brute-force EN+Dyall to <1e-9)"
  else
     print '(a)', "  RESULT: FAIL"
  end if
end program test_icpt2_contracted
