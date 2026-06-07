// M5: libintRot <-> OpenQP C interface (host backend).
//
// Computes one contracted Cartesian shell-quartet block (la,lb,lc,ld) over RAW
// monomial Gaussians using the spherical-resolution engine (eri_assemble), in
// OpenQP's Cartesian component order (constants.F90 CART_X/Y/Z) and OpenQP's
// integral memory layout  ints(nd,nc,nb,na)  (Fortran column-major):
//     out[ nd + nbf4*(nc + nbf3*(nb + nbf2*na)) ] = (a_na b_nb | c_nc d_nd).
//
// OpenQP applies primitive normalization through `cc` (leading component) and
// per-component normalization afterwards via normalize_ints(shells_pnrm2), which
// matches what rys/libint feed it -- so this routine emits the raw-monomial
// block with OpenQP's `cc` and nothing else.  Compiled with a plain C++ compiler
// (host path of port.h); links into liboqp.
#include "sph_eri/assemble_cuda.cuh"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <complex>

// CBLAS complex GEMM, resolved from the MKL/BLAS already linked into liboqp.
// Declared here to avoid a cblas.h / mkl.h build dependency in the OpenQP tree.
// liboqp links ILP64 OpenBLAS, whose CBLAS dimension args are 64-bit (blasint=long);
// the CBLAS_ORDER/TRANSPOSE enums stay 32-bit int.  Declaring M..ldc as long matches
// the ILP64 ABI (chc4 used MKL-ILP64 with identical widths).
extern "C" void cblas_zgemm(int Order,int TransA,int TransB,long M,long N,long K,
    const void* alpha,const void* A,long lda,const void* B,long ldb,
    const void* beta,void* C,long ldc);
enum { ROT_CblasRowMajor=101, ROT_CblasNoTrans=111, ROT_CblasTrans=112 };

namespace {
sph_eri::RadialTable g_rt;
sph_eri::GauntTable  g_gt;
bool g_ready = false;

// OpenQP CART_X/Y/Z(:,0:6) -- component (lx,ly,lz) per shell, OpenQP order.
const int CX[7][28] = {
 {0},
 {1,0,0},
 {2,0,0,1,1,0},
 {3,0,0,2,2,1,0,1,0,1},
 {4,0,0,3,3,1,0,1,0,2,2,0,2,1,1},
 {5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1},
 {6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2}};
const int CY[7][28] = {
 {0},
 {0,1,0},
 {0,2,0,1,0,1},
 {0,3,0,1,0,2,2,0,1,1},
 {0,4,0,1,0,3,3,0,1,2,0,2,1,2,1},
 {0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2},
 {0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2}};
const int CZ[7][28] = {
 {0},
 {0,0,1},
 {0,0,2,0,1,1},
 {0,0,3,0,1,0,1,2,2,1},
 {0,0,4,0,1,0,1,3,3,0,2,2,1,1,2},
 {0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2},
 {0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2}};
inline int ncartf(int l){ return (l+1)*(l+2)/2; }

// GEMM-cast of the per-component Hermite contraction (the 95-99% hot loop).
//   out[cbra,cket] += cabcd*pref * Re{ Pbra . RTmat . Pket^T }
// where  Pbra[cbra][t] = Eab(cbra,t) i^{|t|},  Pket[cket][u] = Ecd(cket,u) i^{-|u|},
//        RTmat[t][u]   = RT[t+u]  (the MD angular sum, complex).
// Two complex GEMMs (MKL zgemm) replace nc^4 scalar component contractions; this is
// the only place MKL accelerates the rot recompute.  Accumulates into out (OpenQP layout).
typedef std::complex<double> cd_t;
inline void gemm_contract(int la,int lb,int lc,int ld,int n1,int n2,int n3,int n4,
    int Wb,int Wk,int lb1,int ld1,int L1,int Lmax,double p,double q,double cabcd,
    const double* EBx,const double* EBy,const double* EBz,
    const double* EKx,const double* EKy,const double* EKz,
    const double* rtre,const double* rtim,double* out){
  const int nBraC=n1*n2, nKetC=n3*n4;
  const int Tb=la+lb+1, Ub=lc+ld+1;
  const int nT=Tb*Tb*Tb, nU=Ub*Ub*Ub;
  static thread_local std::vector<cd_t> Pbra,Pket,RTm,Mm,OC;
  Pbra.assign((size_t)nBraC*nT, cd_t(0,0));
  Pket.assign((size_t)nKetC*nU, cd_t(0,0));
  RTm.resize((size_t)nT*nU); Mm.resize((size_t)nBraC*nU); OC.resize((size_t)nBraC*nKetC);
  // Pbra[cbra][t] = Eab(t) * i^{dt}
  for(int i=0;i<n1;++i){ int ax=CX[la][i],ay=CY[la][i],az=CZ[la][i];
   for(int j=0;j<n2;++j){ int bx=CX[lb][j],by=CY[lb][j],bz=CZ[lb][j]; int cbra=i*n2+j;
    const double* ex=&EBx[(ax*lb1+bx)*Wb]; const double* ey=&EBy[(ay*lb1+by)*Wb]; const double* ez=&EBz[(az*lb1+bz)*Wb];
    int axb=ax+bx,ayb=ay+by,azb=az+bz;
    for(int tx=0;tx<=axb;++tx)for(int ty=0;ty<=ayb;++ty)for(int tz=0;tz<=azb;++tz){
      double pr,pi; sph_eri::ipow(tx+ty+tz,pr,pi); double e=ex[tx]*ey[ty]*ez[tz];
      Pbra[(size_t)cbra*nT + (tx*Tb+ty)*Tb+tz] = cd_t(e*pr,e*pi);
    }}}
  // Pket[cket][u] = Ecd(u) * i^{-du}
  for(int k=0;k<n3;++k){ int cx=CX[lc][k],cy=CY[lc][k],cz=CZ[lc][k];
   for(int m=0;m<n4;++m){ int dx=CX[ld][m],dy=CY[ld][m],dz=CZ[ld][m]; int cket=k*n4+m;
    const double* fx=&EKx[(cx*ld1+dx)*Wk]; const double* fy=&EKy[(cy*ld1+dy)*Wk]; const double* fz=&EKz[(cz*ld1+dz)*Wk];
    int cxd=cx+dx,cyd=cy+dy,czd=cz+dz;
    for(int ux=0;ux<=cxd;++ux)for(int uy=0;uy<=cyd;++uy)for(int uz=0;uz<=czd;++uz){
      double pr,pi; sph_eri::ipow(-(ux+uy+uz),pr,pi); double f=fx[ux]*fy[uy]*fz[uz];
      Pket[(size_t)cket*nU + (ux*Ub+uy)*Ub+uz] = cd_t(f*pr,f*pi);
    }}}
  // RTmat[t][u] = RT[t+u]  (only sum<=Lmax is meaningful; rest never multiplied by nonzero)
  for(int tx=0;tx<Tb;++tx)for(int ty=0;ty<Tb;++ty)for(int tz=0;tz<Tb;++tz){
    int tf=(tx*Tb+ty)*Tb+tz;
    for(int ux=0;ux<Ub;++ux)for(int uy=0;uy<Ub;++uy)for(int uz=0;uz<Ub;++uz){
      int uf=(ux*Ub+uy)*Ub+uz; int sx=tx+ux,sy=ty+uy,sz=tz+uz;
      if(sx+sy+sz<=Lmax){ size_t Ti=((size_t)sx*L1+sy)*L1+sz; RTm[(size_t)tf*nU+uf]=cd_t(rtre[Ti],rtim[Ti]); }
      else RTm[(size_t)tf*nU+uf]=cd_t(0,0);
    }}
  const cd_t one(1,0), zero(0,0);
  // Mm = Pbra (nBraC x nT) * RTm (nT x nU)
  cblas_zgemm(ROT_CblasRowMajor,ROT_CblasNoTrans,ROT_CblasNoTrans, nBraC,nU,nT,
              &one, Pbra.data(),nT, RTm.data(),nU, &zero, Mm.data(),nU);
  // OC = Mm (nBraC x nU) * Pket^T (nU x nKetC)
  cblas_zgemm(ROT_CblasRowMajor,ROT_CblasNoTrans,ROT_CblasTrans, nBraC,nKetC,nU,
              &one, Mm.data(),nU, Pket.data(),nU, &zero, OC.data(),nKetC);
  const double s=cabcd*(2.0/M_PI)*std::pow(M_PI/p,1.5)*std::pow(M_PI/q,1.5);
  for(int cb=0;cb<nBraC;++cb)for(int ck=0;ck<nKetC;++ck)
    out[(size_t)cb*nKetC+ck] += s*OC[(size_t)cb*nKetC+ck].real();
}
} // namespace

extern "C" {

// Diagnostics/verification: OQP_ROT_SCALE multiplies every rot integral (default
// 1.0); a value != 1 must shift the SCF energy iff the rot path is truly in use.
// g_calls counts shell quartets that went through this backend.
static double g_scale = 1.0;
static long   g_calls = 0;
static bool   g_gemm  = true;    // MKL zgemm contraction (default); OQP_ROT_GEMM=0 -> scalar

// Load the precomputed radial + Gaunt tables (call once at SCF start).
int librot_init(const char* radial_path, const char* gaunt_path){
  g_rt = sph_eri::radial_load(radial_path);
  g_gt = sph_eri::gaunt_load(gaunt_path);
  const char* sc = std::getenv("OQP_ROT_SCALE");
  if (sc) g_scale = std::atof(sc);
  const char* gm = std::getenv("OQP_ROT_GEMM");
  if (gm && (gm[0]=='0'||gm[0]=='n'||gm[0]=='N'||gm[0]=='f'||gm[0]=='F')) g_gemm = false;
  // The zgemm contraction is parallelised over quartets by OpenQP's OMP, so each
  // GEMM should be sequential -- otherwise MKL nests under OMP and oversubscribes
  // (~25% slower: 11.1 s vs 8.9 s on H2O/cc-pVTZ).  Calling MKL's thread-control
  // API from here segfaults under the LD_PRELOAD-only MKL setup, so we set the env
  // (best effort; for a guaranteed effect export MKL_NUM_THREADS=1 before running).
  if (g_gemm && !std::getenv("MKL_NUM_THREADS")) setenv("MKL_NUM_THREADS","1",1);
  g_ready = true;
  std::fprintf(stderr, "[libintRot] rot backend initialised (scale=%.6g, contraction=%s, MKL_NUM_THREADS=%s)\n",
               g_scale, g_gemm?"MKL-zgemm":"scalar", std::getenv("MKL_NUM_THREADS")?std::getenv("MKL_NUM_THREADS"):"(unset)");
  return 0;
}

long librot_calls(void){ return g_calls; }

int librot_ready(void){ return g_ready ? 1 : 0; }

// One contracted Cartesian shell-quartet block.  out must hold nbf1*nbf2*nbf3*nbf4
// doubles (= product of (l+1)(l+2)/2).  Returns 0 on success, <0 on bad input.
int librot_eri_cart(int la,int lb,int lc,int ld,
    const double* A,const double* B,const double* C,const double* D,
    const double* ea,const double* ca,int Ka,
    const double* eb,const double* cb,int Kb,
    const double* ec,const double* cc,int Kc,
    const double* ed,const double* cd,int Kd,
    double* out){
  if(!g_ready) return -1;
  if(la<0||lb<0||lc<0||ld<0||la>6||lb>6||lc>6||ld>6) return -2;
  #pragma omp atomic
  ++g_calls;
  // Use the zgemm path unless its (nT x nU) buffers would be very large (high L);
  // then fall back to the scalar contraction.
  const size_t _nT=(size_t)(la+lb+1)*(la+lb+1)*(la+lb+1);
  const size_t _nU=(size_t)(lc+ld+1)*(lc+ld+1)*(lc+ld+1);
  const bool use_gemm = g_gemm && (_nT*_nU <= 4000000u);
  int n1=ncartf(la), n2=ncartf(lb), n3=ncartf(lc), n4=ncartf(ld);
  size_t nout=(size_t)n1*n2*n3*n4;
  for(size_t i=0;i<nout;++i) out[i]=0.0;

  // MD angular-sum hoist (the contracted-shell hot path used by the J/K engine):
  // the geometry-only g_l/Y_lm and the angular sum S(T)=sum_{l,m} i^l A_lm g_l Y_lm
  // are IDENTICAL for every Cartesian component of a primitive quartet, so build
  // them ONCE per primitive quartet and reuse across all n1*n2*n3*n4 components --
  // instead of calling the from-scratch eri_assemble (special functions + Gaunt
  // sum) per component, which recomputed all of this n1*n2*n3*n4 times.
  const int Lmax = la+lb+lc+ld, L1 = Lmax+1, Yw = 2*Lmax+1;
  const size_t nTcell = (size_t)L1*L1*L1;
  const int Wb = la+lb+1, Wk = lc+ld+1;            // per-pair Hermite widths
  const int lb1=lb+1, ld1=ld+1;
  const size_t nBax=(size_t)(la+1)*lb1*Wb;         // one bra Hermite axis table
  const size_t nKax=(size_t)(lc+1)*ld1*Wk;         // one ket Hermite axis table
  static thread_local std::vector<double> gbuf, ybuf, rtre, rtim, EB, EK;
  gbuf.resize((size_t)L1*L1); ybuf.resize((size_t)L1*Yw);
  rtre.resize(nTcell); rtim.resize(nTcell);

  // Hermite E is per-PRIMITIVE-PAIR data: the bra table depends only on (pa,pb)
  // and the ket table only on (pc,pd) -- NOT on the full quartet.  Build each once
  // per pair (bra: Ka*Kb, ket: Kc*Kd) and reuse, instead of rebuilding it Kc*Kd /
  // Ka*Kb times inside the primitive-quartet loop.  dispB/dispK are geometry-only,
  // so they are computed once here, not per primitive.
  const double dispB[3]={A[0]-B[0],A[1]-B[1],A[2]-B[2]};
  const double dispK[3]={C[0]-D[0],C[1]-D[1],C[2]-D[2]};
  EB.resize((size_t)Ka*Kb*3*nBax);
  EK.resize((size_t)Kc*Kd*3*nKax);
  for(int pa=0;pa<Ka;++pa)for(int pb=0;pb<Kb;++pb){
    double* base=&EB[((size_t)(pa*Kb+pb))*3*nBax];
    sph_eri::precompute_hermite<double>(la,lb,dispB,ea[pa],eb[pb], base,base+nBax,base+2*nBax);
  }
  for(int pc=0;pc<Kc;++pc)for(int pd=0;pd<Kd;++pd){
    double* base=&EK[((size_t)(pc*Kd+pd))*3*nKax];
    sph_eri::precompute_hermite<double>(lc,ld,dispK,ec[pc],ed[pd], base,base+nKax,base+2*nKax);
  }
  const sph_eri::GauntTable& gt = g_gt;

  for(int pa=0;pa<Ka;++pa)for(int pb=0;pb<Kb;++pb)for(int pc=0;pc<Kc;++pc)for(int pd=0;pd<Kd;++pd){
    const double cabcd = ca[pa]*cb[pb]*cc[pc]*cd[pd];
    const double p = ea[pa]+eb[pb], q = ec[pc]+ed[pd];
    // (1) g_l / Y_lm once per primitive quartet
    sph_eri::precompute_gY<double>(A,ea[pa], B,eb[pb], C,ec[pc], D,ed[pd], g_rt, Lmax,
                                   gbuf.data(), ybuf.data());
    // (2) angular sum S(T) once per primitive quartet.  eri_contract_RT_pre only ever
    //     reads cells with Tx+Ty+Tz<=Lmax (bra Hermite degree<=la+lb, ket<=lc+ld), so
    //     iterate just that tetrahedron -- ~1/6 of the L1^3 box -- instead of the full
    //     box with most cells zeroed.  Cells outside it are never read.
    for(int Tx=0; Tx<=Lmax; ++Tx)
     for(int Ty=0, Tym=Lmax-Tx; Ty<=Tym; ++Ty)
      for(int Tz=0, Tzm=Lmax-Tx-Ty; Tz<=Tzm; ++Tz){
        size_t Ti=((size_t)Tx*L1+Ty)*L1+Tz; int nT=Tx+Ty+Tz;
        double sre=0.0, sim=0.0;
        if(Tx<=gt.maxdeg && Ty<=gt.maxdeg && Tz<=gt.maxdeg){
          int key=sph_eri::gaunt_encode(Tx,Ty,Tz,gt.maxdeg);
          int s=gt.d_start[key], cnt=gt.d_count[key];
          for(int e=0;e<cnt;++e){
            int l=gt.d_l[s+e], mm=gt.d_m[s+e]; double Almv=gt.d_A[s+e];
            double ilr,ili; sph_eri::ipow(l,ilr,ili);
            double s0 = Almv * gbuf[l*L1+nT] * ybuf[l*Yw+(mm+Lmax)];
            sre += s0*ilr; sim += s0*ili;
          }
        }
        rtre[Ti]=sre; rtim[Ti]=sim;
      }
    // (3) contraction against the precomputed per-pair Hermite E.  Two paths:
    //     scalar per-component (default), or the MKL zgemm cast (OQP_ROT_GEMM=1).
    const double* EBx=&EB[((size_t)(pa*Kb+pb))*3*nBax]; const double* EBy=EBx+nBax; const double* EBz=EBx+2*nBax;
    const double* EKx=&EK[((size_t)(pc*Kd+pd))*3*nKax]; const double* EKy=EKx+nKax; const double* EKz=EKx+2*nKax;
    if(use_gemm){
      gemm_contract(la,lb,lc,ld,n1,n2,n3,n4,Wb,Wk,lb1,ld1,L1,Lmax,p,q,cabcd,
                    EBx,EBy,EBz, EKx,EKy,EKz, rtre.data(),rtim.data(), out);
    } else {
    for(int i=0;i<n1;++i){ int ax=CX[la][i],ay=CY[la][i],az=CZ[la][i];
     for(int j=0;j<n2;++j){ int bx=CX[lb][j],by=CY[lb][j],bz=CZ[lb][j];
      const double* ex=&EBx[(ax*lb1+bx)*Wb]; const double* ey=&EBy[(ay*lb1+by)*Wb]; const double* ez=&EBz[(az*lb1+bz)*Wb];
      int axb=ax+bx, ayb=ay+by, azb=az+bz;
      for(int k=0;k<n3;++k){ int cx=CX[lc][k],cy=CY[lc][k],cz=CZ[lc][k];
       for(int m=0;m<n4;++m){ int dx=CX[ld][m],dy=CY[ld][m],dz=CZ[ld][m];
        const double* fx=&EKx[(cx*ld1+dx)*Wk]; const double* fy=&EKy[(cy*ld1+dy)*Wk]; const double* fz=&EKz[(cz*ld1+dz)*Wk];
        double v = sph_eri::eri_contract_RT_pre<double>(
                     axb,ayb,azb, cx+dx,cy+dy,cz+dz, ex,ey,ez, fx,fy,fz,
                     rtre.data(), rtim.data(), Lmax, p, q);
        // OpenQP layout: ints(nd,nc,nb,na), nd fastest
        out[ ((((size_t)i*n2 + j)*n3 + k)*n4 + m) ] += cabcd * v;
       }}}}
    }
  }
  if(g_scale!=1.0) for(size_t i=0;i<nout;++i) out[i]*=g_scale;
  return 0;
}

} // extern "C"
