! QMRSF-icPT2 full pipeline (standalone, build-testable): MO integrals -> spin-orbital
! tensors -> full determinant CI -> P/Q partition -> internally-contracted external-Q
! EN downfold (des Cloizeaux multistate). Faithful Fortran port of the validated NumPy
! prototype (qmrsf_icpt2_multistate.py / qmrsf_icpt2_ppp_proto.py). This is the BRUTE-FORCE
! perturber generation (builds the full external space explicitly) -- correct for small
! systems and the validation bridge for the production contracted engine; it consumes the
! same downfold algebra already validated in qmrsf_icpt2_downfold.f90.
!
! Reads icpt2_full_ref.dat (norb,na,nb,ncore,nact,nPdress; h_mo; chemist eri_mo; ref EN).
module icpt2_full_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: norb, nso, nelec, na, nb
contains

  pure integer function spat(P)        ! 1-based spatial index of spin-orbital P
    integer, intent(in) :: P
    spat = mod(P-1, norb) + 1
  end function
  pure integer function spn(P)         ! 0 alpha (P<=norb), 1 beta
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
      a = 0.0_dp; b = 0.0_dp
      if (spn(P)==spn(R) .and. spn(Q)==spn(S)) a = eri(spat(P),spat(R),spat(Q),spat(S))
      if (spn(P)==spn(S) .and. spn(Q)==spn(R)) b = eri(spat(P),spat(S),spat(Q),spat(R))
      g(P,Q,R,S) = a - b
    end do; end do; end do; end do
  end subroutine

  ! next combination of k indices from 1..m in lexicographic order; .false. when done
  logical function next_comb(c, k, m)
    integer, intent(inout) :: c(k)
    integer, intent(in) :: k, m
    integer :: i, j
    i = k
    do while (i >= 1)
      if (c(i) /= m - k + i) exit
      i = i - 1
    end do
    if (i < 1) then; next_comb = .false.; return; end if
    c(i) = c(i) + 1
    do j = i+1, k
      c(j) = c(j-1) + 1
    end do
    next_comb = .true.
  end function

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
    ! phase: annihilate holes (ascending) then create parts (descending) on D2
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
  end function

end module icpt2_full_mod

program test_icpt2_full
  use icpt2_full_mod
  implicit none
  integer :: ncore, nact, nPd, i, j, k, l, q, p, info, lwork, ierr
  integer :: nalp, nbet, ndet, nP, nQ, ca, cb, vlo
  real(dp), allocatable :: h_mo(:,:), eri(:,:,:,:), H1(:,:), g(:,:,:,:)
  integer,  allocatable :: dets(:,:), astr(:,:), bstr(:,:)
  integer,  allocatable :: Pid(:), Qid(:)
  real(dp), allocatable :: Hfull(:,:), HPP(:,:), HQP(:,:), Hqq(:)
  real(dp), allocatable :: ePall(:), cP(:,:), eP(:), coup(:,:), invd(:,:)
  real(dp), allocatable :: refEN(:), edEN(:), work(:), Heff(:,:)
  integer,  allocatable :: c(:)
  logical :: inP
  real(dp) :: s, d, dEN, herm

  open(10, file="icpt2_full_ref.dat", status="old", action="read")
  read(10,*) norb, na, nb, ncore, nact, nPd
  nso = 2*norb; nelec = na + nb
  allocate(h_mo(norb,norb), eri(norb,norb,norb,norb))
  do i=1,norb; read(10,*) (h_mo(i,j), j=1,norb); end do
  do i=1,norb; do j=1,norb; do k=1,norb
     read(10,*) (eri(i,j,k,l), l=1,norb)
  end do; end do; end do
  allocate(refEN(nPd)); read(10,*) (refEN(i), i=1,nPd)
  close(10)

  allocate(H1(nso,nso), g(nso,nso,nso,nso))
  call build_spinorb(h_mo, eri, H1, g)

  ! enumerate alpha/beta strings and merge into spin-orbital determinants
  nalp = ncomb(norb, na); nbet = ncomb(norb, nb)
  ndet = nalp*nbet
  allocate(astr(na,nalp), bstr(nb,nbet), dets(nelec,ndet))
  call all_combs(norb, na, astr, nalp)
  call all_combs(norb, nb, bstr, nbet)
  l = 0
  do ca=1,nalp
    do cb=1,nbet
      l = l+1
      do i=1,na; dets(i,l) = astr(i,ca); end do
      do i=1,nb; dets(na+i,l) = bstr(i,cb) + norb; end do
      call isort(dets(:,l), nelec)
    end do
  end do

  ! partition: P = core doubly occupied AND no virtual occupied
  vlo = ncore + nact            ! virtual spatial orbitals = vlo+1 .. norb
  allocate(Pid(ndet), Qid(ndet)); nP=0; nQ=0
  do l=1,ndet
    inP = .true.
    do p=1,ncore
      if (.not. (any(dets(:,l)==p) .and. any(dets(:,l)==p+norb))) inP=.false.
    end do
    do p=vlo+1,norb
      if (any(dets(:,l)==p) .or. any(dets(:,l)==p+norb)) inP=.false.
    end do
    if (inP) then; nP=nP+1; Pid(nP)=l; else; nQ=nQ+1; Qid(nQ)=l; end if
  end do

  ! H_PP, H_QP, diag(H)_Q
  allocate(HPP(nP,nP), HQP(nQ,nP), Hqq(nQ))
  do i=1,nP; do j=1,nP; HPP(i,j)=melem(dets(:,Pid(i)),dets(:,Pid(j)),H1,g); end do; end do
  do i=1,nQ; do j=1,nP; HQP(i,j)=melem(dets(:,Qid(i)),dets(:,Pid(j)),H1,g); end do; end do
  do i=1,nQ; Hqq(i)=melem(dets(:,Qid(i)),dets(:,Qid(i)),H1,g); end do

  ! diagonalize H_PP -> CAS roots/vectors
  allocate(ePall(nP), cP(nP,nP)); cP = HPP
  lwork = 64*nP; allocate(work(lwork))
  call dsyev('V','U', nP, cP, nP, ePall, work, lwork, info)
  if (info/=0) stop 'dsyev HPP'
  allocate(eP(nPd)); eP = ePall(1:nPd)

  ! contracted couplings + EN denominators
  allocate(coup(nQ,nPd), invd(nQ,nPd))
  coup = matmul(HQP, cP(:,1:nPd))
  do k=1,nPd
    do q=1,nQ
      d = eP(k) - Hqq(q)
      if (abs(d) < 1.0e-6_dp) d = sign(1.0e-6_dp,d) + 1.0e-30_dp
      invd(q,k) = 1.0_dp/d
    end do
  end do

  ! des Cloizeaux multistate effective Hamiltonian + spectrum
  allocate(Heff(nPd,nPd), edEN(nPd)); Heff = 0.0_dp
  do k=1,nPd; Heff(k,k)=eP(k); end do
  do k=1,nPd; do l=1,nPd
    s = 0.0_dp
    do q=1,nQ; s = s + coup(q,k)*coup(q,l)*0.5_dp*(invd(q,k)+invd(q,l)); end do
    Heff(k,l) = Heff(k,l) + s
  end do; end do
  herm = 0.0_dp
  do k=1,nPd; do l=1,nPd; herm=max(herm,abs(Heff(k,l)-Heff(l,k))); end do; end do
  Heff = 0.5_dp*(Heff+transpose(Heff))
  deallocate(work); lwork=64*nPd; allocate(work(lwork))
  call dsyev('N','U', nPd, Heff, nPd, edEN, work, lwork, info)
  if (info/=0) stop 'dsyev Heff'

  dEN = maxval(abs(edEN - refEN))
  print '(a)', "==== QMRSF-icPT2 FULL pipeline (Fortran) vs NumPy prototype ===="
  print '(a,i0,a,i0,a,i0,a,i0)', "  norb=",norb,"  ndet=",ndet,"  nP=",nP,"  nQ=",nQ
  print '(a,4f14.8)', "  CAS roots eP      = ", eP
  print '(a,4f14.8)', "  dressed EN (Fort) = ", edEN
  print '(a,4f14.8)', "  dressed EN (ref)  = ", refEN
  print '(a,es10.2)', "  max|EN - ref|     = ", dEN
  print '(a,es10.2)', "  Heff Hermiticity  = ", herm
  if (dEN < 1.0d-9) then
     print '(a)', "  RESULT: PASS  (integrals -> det-CI -> downfold matches NumPy to <1e-9)"
  else
     print '(a)', "  RESULT: FAIL"
  end if

contains
  integer function ncomb(m,k)
    integer, intent(in) :: m,k
    integer :: i
    ncomb = 1
    do i=1,k; ncomb = ncomb*(m-k+i)/i; end do
  end function
  subroutine all_combs(m, k, out, ncol)
    integer, intent(in) :: m, k, ncol
    integer, intent(out) :: out(k,ncol)
    integer :: cc(k), j, col
    do j=1,k; cc(j)=j; end do
    col = 0
    do
      col = col+1; out(:,col) = cc
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
      do while (j>=1)
        if (a(j)<=t) exit
        a(j+1)=a(j); j=j-1
      end do
      a(j+1)=t
    end do
  end subroutine
end program test_icpt2_full
