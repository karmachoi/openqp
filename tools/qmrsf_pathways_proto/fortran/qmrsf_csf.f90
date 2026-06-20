! QMRSF CSF spin-adaptation (Fortran, build-testable) -- stage 1.
! Takes the DSF-generated CAS(4,4) Ms=0 determinant Hamiltonian (Slater-Condon, same core as
! qmrsf_backbone_core.f90) and builds the CSF spin-adaptation U: 0OS determinants are singlets;
! 2OS pairs are Hadamard-coupled (singlet = antisym, triplet = sym); the 4OS sextet is spin-coupled
! numerically (diagonalizing S^2 within the 6-determinant block = a genealogical-CG-equivalent basis).
! Validation: U^T H U is block-diagonal in S, the blocks have dims 20/15/1, <S^2> = 0/2/6 exactly,
! and the union of the spin-block spectra reproduces the full 36-state spectrum -- i.e. the CSF
! construction yields the spin-pure blocks (no S^2 projection of the production operator needed).
module qmrsf_csf_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: NACT = 4, NSO = 8, NDET = 36
contains
  subroutine build_spinorb(h_act, eri_act, H1, g)
    real(dp), intent(in)  :: h_act(NACT,NACT), eri_act(NACT,NACT,NACT,NACT)
    real(dp), intent(out) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
    integer :: P,Q,R,S, spat(NSO), spin(NSO), i
    real(dp) :: a,b
    do i=1,NSO
      if (i<=NACT) then; spat(i)=i; spin(i)=0; else; spat(i)=i-NACT; spin(i)=1; end if
    end do
    H1=0.0_dp
    do P=1,NSO; do Q=1,NSO; if (spin(P)==spin(Q)) H1(P,Q)=h_act(spat(P),spat(Q)); end do; end do
    g=0.0_dp
    do P=1,NSO; do Q=1,NSO; do R=1,NSO; do S=1,NSO
      a=0.0_dp; b=0.0_dp
      if (spin(P)==spin(R).and.spin(Q)==spin(S)) a=eri_act(spat(P),spat(R),spat(Q),spat(S))
      if (spin(P)==spin(S).and.spin(Q)==spin(R)) b=eri_act(spat(P),spat(S),spat(Q),spat(R))
      g(P,Q,R,S)=a-b
    end do;end do;end do;end do
  end subroutine

  subroutine gen_dets(dets)
    integer, intent(out) :: dets(4,NDET)
    integer :: a1,a2,b1,b2,k,t(4),i,j,tmp
    k=0
    do a1=1,NACT-1; do a2=a1+1,NACT; do b1=1,NACT-1; do b2=b1+1,NACT
      k=k+1; t=(/a1,a2,b1+NACT,b2+NACT/)
      do i=1,3; do j=i+1,4; if (t(j)<t(i)) then; tmp=t(i);t(i)=t(j);t(j)=tmp; end if; end do; end do
      dets(:,k)=t
    end do;end do;end do;end do
  end subroutine

  real(dp) function melem(D1,D2,H1,g)
    integer,intent(in)::D1(4),D2(4)
    real(dp),intent(in)::H1(NSO,NSO),g(NSO,NSO,NSO,NSO)
    integer::holes(4),parts(4),common(4),nh,np,nc,occ(4),nocc,i,idx,k,p1,p2,ho1,ho2,Pp,Hh,Qc
    real(dp)::sgn,val,e
    nh=0;np=0;nc=0
    do i=1,4; if (.not.any(D1==D2(i))) then; nh=nh+1; holes(nh)=D2(i); end if; end do
    do i=1,4
      if (.not.any(D2==D1(i))) then; np=np+1; parts(np)=D1(i); end if
      if (any(D2==D1(i))) then; nc=nc+1; common(nc)=D1(i); end if
    end do
    if (nh>2) then; melem=0.0_dp; return; end if
    occ=D2; nocc=4; sgn=1.0_dp
    do k=1,nh
      idx=0; do i=1,nocc; if (occ(i)==holes(k)) then; idx=i; exit; end if; end do
      if (mod(idx-1,2)==1) sgn=-sgn
      do i=idx,nocc-1; occ(i)=occ(i+1); end do; nocc=nocc-1
    end do
    do k=np,1,-1
      idx=1; do i=1,nocc; if (occ(i)<parts(k)) idx=idx+1; end do
      if (mod(idx-1,2)==1) sgn=-sgn
      do i=nocc,idx,-1; occ(i+1)=occ(i); end do; occ(idx)=parts(k); nocc=nocc+1
    end do
    if (nh==0) then
      e=0.0_dp; do i=1,4; e=e+H1(D1(i),D1(i)); end do
      do i=1,3; do k=i+1,4; e=e+g(D1(i),D1(k),D1(i),D1(k)); end do; end do
      melem=e
    else if (nh==1) then
      Pp=parts(1);Hh=holes(1); val=H1(Pp,Hh)
      do i=1,nc; Qc=common(i); val=val+g(Pp,Qc,Hh,Qc); end do
      melem=sgn*val
    else
      p1=parts(1);p2=parts(2);ho1=holes(1);ho2=holes(2); melem=sgn*g(p1,p2,ho1,ho2)
    end if
  end function

  ! a^dag_cre a_ann on sorted spin-orbital det D -> Dnew, phase (0 if annihilated)
  subroutine apply_ex(D,cre,ann,Dnew,phase)
    integer,intent(in)::D(4),cre,ann
    integer,intent(out)::Dnew(4)
    real(dp),intent(out)::phase
    integer::occ(5),nocc,i,idx
    Dnew=0; phase=0.0_dp
    if (.not.any(D==ann)) return
    if (cre/=ann .and. any(D==cre)) return
    occ(1:4)=D; nocc=4; phase=1.0_dp
    idx=0; do i=1,nocc; if (occ(i)==ann) then; idx=i; exit; end if; end do
    if (mod(idx-1,2)==1) phase=-phase
    do i=idx,nocc-1; occ(i)=occ(i+1); end do; nocc=nocc-1
    idx=1; do i=1,nocc; if (occ(i)<cre) idx=idx+1; end do
    if (mod(idx-1,2)==1) phase=-phase
    do i=nocc,idx,-1; occ(i+1)=occ(i); end do; occ(idx)=cre; nocc=nocc+1
    Dnew=occ(1:4)
  end subroutine

  ! S^2 = S_- S_+  (Ms=0 sector) in the determinant basis
  subroutine build_S2(dets,S2)
    integer,intent(in)::dets(4,NDET)
    real(dp),intent(out)::S2(NDET,NDET)
    integer::j,p,q,Di(4),Dk(4),idx,m
    real(dp)::ph1,ph2
    S2=0.0_dp
    do j=1,NDET
      do p=1,NACT                      ! S_+ = sum_p a^dag_{p,alpha} a_{p,beta}
        call apply_ex(dets(:,j), p, p+NACT, Di, ph1)
        if (ph1==0.0_dp) cycle
        do q=1,NACT                    ! S_- = sum_q a^dag_{q,beta} a_{q,alpha}
          call apply_ex(Di, q+NACT, q, Dk, ph2)
          if (ph2==0.0_dp) cycle
          idx=0; do m=1,NDET; if (all(dets(:,m)==Dk)) then; idx=m; exit; end if; end do
          if (idx>0) S2(idx,j)=S2(idx,j)+ph1*ph2
        end do
      end do
    end do
  end subroutine

  ! spatial occupation vector (0/1/2 per active orbital) of a determinant
  subroutine occvec(D,occ)
    integer,intent(in)::D(4)
    integer,intent(out)::occ(NACT)
    integer::i,p
    occ=0
    do i=1,4; if (D(i)<=NACT) then; p=D(i); else; p=D(i)-NACT; end if; occ(p)=occ(p)+1; end do
  end subroutine
end module qmrsf_csf_mod


program qmrsf_csf
  use qmrsf_csf_mod
  implicit none
  real(dp) :: h_act(NACT,NACT), eri_act(NACT,NACT,NACT,NACT)
  real(dp) :: H1(NSO,NSO), g(NSO,NSO,NSO,NSO)
  integer  :: dets(4,NDET), p,q,r,s,i,j,k,nr,nev,info
  real(dp) :: Hmat(NDET,NDET), S2(NDET,NDET), U(NDET,NDET), Hc(NDET,NDET)
  real(dp) :: evfull(NDET), work(2048), s2lab(NDET), tmpv(NDET)
  integer  :: occ_i(NACT), occ_j(NACT), gmem(NDET), ng, used(NDET)
  real(dp) :: sub(6,6), subev(6), subvec(6,6), col(NDET)
  integer  :: ncsf, n0,n1,n2, gidx(NDET)
  real(dp) :: blk(NDET,NDET), blkev(NDET), allev(NDET); integer::nall
  real(dp) :: offmax, herm2, dmax, ssq

  open(10,file="qmrsf_cas_ref.dat",status="old",action="read")
  read(10,*) nr
  do p=1,NACT; read(10,*) (h_act(p,q),q=1,NACT); end do
  do p=1,NACT;do q=1,NACT;do r=1,NACT; read(10,*) (eri_act(p,q,r,s),s=1,NACT); end do;end do;end do
  read(10,*) nev; read(10,*) (evfull(i),i=1,nev); close(10)

  call build_spinorb(h_act,eri_act,H1,g)
  call gen_dets(dets)
  do i=1,NDET; do j=1,NDET; Hmat(i,j)=melem(dets(:,i),dets(:,j),H1,g); end do; end do
  call build_S2(dets,S2)
  herm2=0.0_dp; do i=1,NDET;do j=1,NDET; herm2=max(herm2,abs(S2(i,j)-S2(j,i))); end do;end do

  ! ---- CSF construction ----
  U=0.0_dp; ncsf=0; used=0
  do i=1,NDET
    if (used(i)==1) cycle
    call occvec(dets(:,i),occ_i)
    ! gather all dets with the same spatial occupation
    ng=0
    do j=1,NDET
      call occvec(dets(:,j),occ_j)
      if (all(occ_j==occ_i)) then; ng=ng+1; gmem(ng)=j; used(j)=1; end if
    end do
    if (ng==1) then                                  ! 0OS: singlet determinant
      ncsf=ncsf+1; U(gmem(1),ncsf)=1.0_dp
    else if (ng==2) then                             ! 2OS: Hadamard pair
      ncsf=ncsf+1; U(gmem(1),ncsf)= 1.0_dp/sqrt(2.0_dp); U(gmem(2),ncsf)=-1.0_dp/sqrt(2.0_dp)
      ncsf=ncsf+1; U(gmem(1),ncsf)= 1.0_dp/sqrt(2.0_dp); U(gmem(2),ncsf)= 1.0_dp/sqrt(2.0_dp)
    else                                             ! 4OS sextet: spin-couple via S^2 in the block
      do p=1,ng; do q=1,ng; sub(p,q)=S2(gmem(p),gmem(q)); end do; end do
      subvec=sub; call dsyev('V','U',ng,subvec,6,subev,work,2048,info)
      do p=1,ng
        ncsf=ncsf+1
        do q=1,ng; U(gmem(q),ncsf)=subvec(q,p); end do
      end do
    end if
  end do

  ! ---- label each CSF by <S^2> and validate block structure ----
  Hc=matmul(transpose(U),matmul(Hmat,U))
  n0=0;n1=0;n2=0
  do k=1,NDET
    col=U(:,k); tmpv=matmul(S2,col); ssq=dot_product(col,tmpv); s2lab(k)=ssq
    if (abs(ssq)<1d-8) then; n0=n0+1; gidx(k)=0
    else if (abs(ssq-2.0_dp)<1d-8) then; n1=n1+1; gidx(k)=1
    else if (abs(ssq-6.0_dp)<1d-8) then; n2=n2+1; gidx(k)=2
    else; gidx(k)=-1; end if
  end do
  ! off-block-diagonal coupling between different-S CSFs (must vanish)
  offmax=0.0_dp
  do i=1,NDET; do j=1,NDET; if (gidx(i)/=gidx(j)) offmax=max(offmax,abs(Hc(i,j))); end do; end do
  ! union of spin-block spectra vs full spectrum
  nall=0
  do s=0,2
    k=0
    do i=1,NDET; if (gidx(i)==s) then; k=k+1; gmem(k)=i; end if; end do
    if (k==0) cycle
    do i=1,k; do j=1,k; blk(i,j)=Hc(gmem(i),gmem(j)); end do; end do
    call dsyev('N','U',k,blk,NDET,blkev,work,2048,info)
    do i=1,k; nall=nall+1; allev(nall)=blkev(i); end do
  end do
  call dsort(allev,nall); call dsort(evfull,nev)
  dmax=0.0_dp; do i=1,nev; dmax=max(dmax,abs(allev(i)-evfull(i))); end do

  print '(a)', "==== QMRSF CSF spin-adaptation (Fortran) ===="
  print '(a,es12.3)', "  S^2 Hermiticity            = ", herm2
  print '(a,i0,a,i0,a,i0,a,i0)', "  CSF counts: singlet=",n0,"  triplet=",n1,"  quintet=",n2,"  total=",ncsf
  print '(a,es12.3)', "  max off-S-block |U^T H U|   = ", offmax
  print '(a,es12.3)', "  max|block spectra union - full FCI| = ", dmax
  if (n0==20 .and. n1==15 .and. n2==1 .and. offmax<1d-9 .and. dmax<1d-9) then
     print '(a)', "  RESULT: PASS  (CSF construction -> spin-pure 20/15/1 blocks, reproduce full spectrum)"
  else
     print '(a)', "  RESULT: FAIL"
  end if
contains
  subroutine dsort(a,n)
    real(dp),intent(inout)::a(*); integer,intent(in)::n
    integer::i,j; real(dp)::t
    do i=1,n-1; do j=1,n-i; if (a(j)>a(j+1)) then; t=a(j);a(j)=a(j+1);a(j+1)=t; end if; end do; end do
  end subroutine
end program
