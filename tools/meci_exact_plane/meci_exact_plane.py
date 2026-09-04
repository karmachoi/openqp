"""Projected steepest descent on the seam with the exact branching plane (x from grad difference, y from analytic NAC).
Bearpark-Robb composite gradient: F = 2 dE x_hat + P_perp(mean gradient); fixed-trust-radius steps."""
import sys, numpy as np
from oqp.pyoqp import Runner
from oqp.library.libscipy import MECIOpt
from oqp.library.nac_analytic import analytic_nac
inp=sys.argv[1]; maxit=int(sys.argv[2]) if len(sys.argv)>2 else 150; trust=float(sys.argv[3]) if len(sys.argv)>3 else 0.02  # bohr max displacement
geomfile=sys.argv[4] if len(sys.argv)>4 else None
r=Runner(input_file=inp, log=inp.replace('.inp','.log')); mol=r.mol; opt=MECIOpt(mol); opt.meci_search='ubp'
I,J=opt.istate,opt.jstate; log=open('sd_trace.txt','a')
def composite(coords):
    opt.itr+=1
    energies,grads=opt.evaluate(coords) if hasattr(opt,'evaluate') else (None,None)
    return energies,grads
# use MECIOpt.one_step's machinery through a work_func hook that stores energies/grads
store={}
def hook(coordinates,energies,grads):
    store['e']=np.array(energies,dtype=float); store['g']=[np.array(g,dtype=float).reshape(-1) for g in grads]; return 0.0, np.zeros_like(coordinates)
opt.work_func=hook
x=np.array(mol.get_system(),dtype=float).reshape(-1)
if geomfile: x=np.loadtxt(geomfile).reshape(-1); log.write(f'restart from {geomfile}\n')
prev_e=None
st={}
for it in range(1,maxit+1):
    opt.one_step(x)
    e=store['e']; gi,gj=store['g'][I],store['g'][J]
    dE=e[J]-e[I]; gdiff=(gj-gi)/2; gmean=(gj+gi)/2; xh=gdiff/np.linalg.norm(gdiff)
    nacv,_=analytic_nac(mol); h=np.array(nacv[I-1,J-1]).reshape(-1); y=h-np.dot(h,xh)*xh; yn=np.linalg.norm(y); yh=y/yn if yn>1e-12 else y
    Pg=gmean-np.dot(gmean,xh)*xh-np.dot(gmean,yh)*yh
    F=2.0*dE*xh+Pg
    fmax=np.abs(F).max(); frms=np.sqrt(np.mean(F**2))
    gx=np.linalg.norm(gdiff); hh=np.linalg.norm(h); sx=float(gmean@xh); sy=float(gmean@yh)
    dgh=np.sqrt((gx**2+hh**2)/2); Dgh=(gx**2-hh**2)/(gx**2+hh**2); P=(sx/gx)**2+(sy/hh)**2; Bi=(abs(sx)/gx)**(2/3)+(abs(sy)/hh)**(2/3)
    log.write(f"   invariants: |g|={gx:.5f} |h|={hh:.5f} s_x={sx:.5f} s_y={sy:.5f} delta_gh={dgh:.5f} Delta_gh={Dgh:.4f} P={P:.4f} B={Bi:.4f} E(S1)={e[I]:.10f} E(S2)={e[J]:.10f} E(S0)={e[1] if len(e)>1 else float('nan'):.10f}\n")
    np.savez(f'point_{it}.npz',x=x,gdiff=gdiff,h=h,gmean=gmean,gi=gi,gj=gj,e=e)
    log.write(f"{it:4d} Ei={e[I]:.10f} Ej={e[J]:.10f} gap={dE:.3e} |g|={np.linalg.norm(gdiff):.4e} |h|={np.linalg.norm(h):.4e} |Pg|max={np.abs(Pg).max():.3e} Fmax={fmax:.3e} Frms={frms:.3e}\n"); log.flush()
    if fmax<3e-5 and abs(dE)<5e-6:
        s_=gmean; gx=np.linalg.norm(gdiff); hh=np.linalg.norm(h); sx=float(s_@xh); sy=float(s_@yh)
        dgh=np.sqrt((gx**2+hh**2)/2); Dgh=(gx**2-hh**2)/(gx**2+hh**2); P=(sx/gx)**2+(sy/hh)**2; Bi=(abs(sx)/gx)**(2/3)+(abs(sy)/hh)**(2/3)
        log.write(f"CONVERGED |g|={gx:.5f} |h|={hh:.5f} s_x={sx:.5f} s_y={sy:.5f} delta_gh={dgh:.5f} Delta_gh={Dgh:.4f} P={P:.4f} B={Bi:.4f}\n"); log.flush()
        np.savetxt('converged_geom_bohr.txt',x.reshape(-1,3)); break
    n=x.size
    if 'H' not in st:
        st['H']=np.eye(n)*2.0   # inverse Hessian guess (bohr^2/Eh)
    else:
        sv=x-st['x']; yv=F-st['F']; sy=float(sv@yv)
        if sy>1e-10:
            Hk=st['H']; rho=1.0/sy; I_=np.eye(n)
            st['H']=(I_-rho*np.outer(sv,yv))@Hk@(I_-rho*np.outer(yv,sv))+rho*np.outer(sv,sv)
    st['x']=x.copy(); st['F']=F.copy()
    step=-st['H']@F
    if float(step@F)>0: step=-F; st['H']=np.eye(n)*2.0   # not a descent direction: reset
    smax=np.abs(step).max(); scale=min(1.0, trust/smax) if smax>0 else 0.0
    log.write(f'     |step|max={smax:.3e} scale={scale:.3f}\n')
    x=x+scale*step
    np.savetxt('sd_last_geom.txt',x.reshape(-1,3))
print("SD_DONE")
