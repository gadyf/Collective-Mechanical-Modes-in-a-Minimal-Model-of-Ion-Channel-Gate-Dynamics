"""Independent re-implementation of the 4-rod/10-spring toy model (geometry copied from the notebook).
Checks rank of the spring compatibility matrix, eigenvalues, zero modes and SI scaling. Run: python check_model.py"""
import numpy as np
np.set_printoptions(precision=4, suppress=True, linewidth=160)
rod=[(-1,0,np.pi/2-np.pi/9,3),(1,0,np.pi/2+np.pi/9,3),(-1.1,0.5,np.pi/2-np.pi/20,2),(1.1,0.5,np.pi/2+np.pi/20,2)]
conn=[("r1","r2"),("r1","r2"),("r1","r3"),("r2","r4"),("r1","r3"),("r2","r4"),("r1","f2"),("r3","f1"),("r2","f4"),("r4","f3")]
loc=[[1,1],[-1,-1],[-0.5,-1.0],[-0.5,-1.0],[1.4,0.85],[1.4,0.85],[-1.5,0],[1,0],[-1.5,0],[1,0]]
fp=[(-3,2),(-3,-2),(3,2),(3,-2)]
kk=np.array([10,20,1000,1000,10,10,10,10,10,10.])
m=np.array([3,3,2,2.]); I=m*np.array([r[3] for r in rod])**2/12
q0=np.array([v for r in rod for v in r[:3]])
W=np.array([w for i in range(4) for w in (m[i],m[i],I[i])])
def pt(q,end,d):
    if end[0]=='f': return np.array(fp[int(end[1:])-1],float)
    i=int(end[1:])-1; x,y,t=q[3*i:3*i+3]; return np.array([x+d*np.cos(t),y+d*np.sin(t)])
def lengths(q): return np.array([np.linalg.norm(pt(q,a,loc[j][0])-pt(q,b,loc[j][1])) for j,(a,b) in enumerate(conn)])
L0=lengths(q0)
def U(q): return 0.5*np.sum(kk*(lengths(q)-L0)**2)
# compatibility (rigidity) matrix B = dl/dq  (10x12)
h=1e-7; B=np.zeros((10,12))
for c in range(12):
    e=np.zeros(12); e[c]=h; B[:,c]=(lengths(q0+e)-lengths(q0-e))/(2*h)
sv=np.linalg.svd(B,compute_uv=False); print("singular values of B (10x12):",sv); print("rank B =",np.sum(sv>1e-8*sv[0]))
Hs=B.T@np.diag(kk)@B             # material stiffness (unstressed => full Hessian)
Mi=np.diag(1/np.sqrt(W)); D=Mi@Hs@Mi
ev,V=np.linalg.eigh(D); print("eigenvalues of mass-weighted Hessian (code units):",ev)
print("ratio lowest two / largest:",ev[:2]/ev[-1])
# compare to notebook values
nb=[1.28785085249803e-14,2.29210996660496e-12,0.0448555198659439,0.350260025297473,2.91728075148219,9.65103558280504,12.2995617754985,33.4739892102471,36.0306679260537,41.1579461851789,2231.72059267258,2236.54022345842]
print("notebook eig:",np.array(nb))
# scaling to SI: omega_SI = omega_code*sqrt(fk/fm)
fm=1.44e-24
for fk,lab in [(0.1,"k = 1, 2, 100 N/m (used)"),(1.0,"k = 10, 20, 1000 N/m (model values read as N/m)")]:
    print(lab, "omega SI:", np.array([np.sqrt(max(x,0))*np.sqrt(fk/fm) for x in ev]))
# residual force at equilibrium
g=np.array([(U(q0+np.eye(12)[c]*1e-6)-U(q0-np.eye(12)[c]*1e-6))/2e-6 for c in range(12)]); print("grad U at q0 (should be 0):",np.abs(g).max())
# Null space of B in physical coords (mechanisms)
_,s,Vt=np.linalg.svd(B); N=Vt[-2:].T   # 12x2
names=[f"{c}{i}" for i in range(1,5) for c in ("x","y","th")]
# symmetry adapt: mirror x->-x maps rod1<->rod2, rod3<->rod4, x->-x, y->y, th->pi-th (dth->-dth)
P=np.zeros((12,12))
for a,b in [(0,1),(1,0),(2,3),(3,2)]:
    P[3*b,3*a]=-1; P[3*b+1,3*a+1]=1; P[3*b+2,3*a+2]=-1
Pn=N.T@P@N; w,u=np.linalg.eigh(Pn); Ns=N@u
for j in range(2):
    v=Ns[:,j]/np.abs(Ns[:,j]).max()
    print(f"\nzero mode {j+1}: mirror parity {w[j]:+.3f}")
    print("  ", {n:round(x,3) for n,x in zip(names,v)})
    # effect on channel geometry: width at top/mid/bottom and y-shift
    def ends(qv,i,sgn): x,y,t=qv[3*i:3*i+3]; Lh=rod[i][3]/2; return np.array([x+sgn*Lh*np.cos(t), y+sgn*Lh*np.sin(t)])
    eps=1e-6; qa=q0+eps*v
    for nm,s_ in [("top",1),("bottom",-1)]:
        d=(ends(qa,1,s_)-ends(qa,0,s_))-(ends(q0,1,s_)-ends(q0,0,s_)); print(f"   d(width,{nm})={d[0]/eps:+.3f}  d(yshift,{nm})={d[1]/eps:+.3f}")
    dm=(qa[3:5]-qa[0:2])-(q0[3:5]-q0[0:2]); print(f"   d(width,mid)={dm[0]/eps:+.3f}  d(yshift,mid)={dm[1]/eps:+.3f}")
    # second-order energy along mechanism (finite amplitude): does U grow as a^4?
    for a in [1e-3,1e-2,1e-1]:
        print(f"   U(a={a}) = {U(q0+a*v):.3e}")
# Which rods participate: per-rod norm
for j in range(2):
    v=Ns[:,j]; print("rod participation zero mode",j+1,[round(np.linalg.norm(v[3*i:3*i+3]),3) for i in range(4)])
import os; np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)),"nullspace.npy"), Ns)
