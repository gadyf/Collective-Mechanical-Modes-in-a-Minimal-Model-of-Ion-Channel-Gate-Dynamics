"""Tests whether the zero modes are finite mechanisms (relaxed energy along the null space). Run check_model.py first."""
import os
HERE=os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.optimize import minimize
exec(open(os.path.join(HERE,"check_model.py")).read().split("# Null space")[0].replace('print(','(lambda *a,**k:None)('))
_,s,Vt=np.linalg.svd(B); N=np.load(os.path.join(HERE,"nullspace.npy"))
Perp=Vt[:10].T
for j in range(2):
    v=N[:,j]/np.linalg.norm(N[:,j]); other=N[:,1-j]/np.linalg.norm(N[:,1-j])
    print(f"--- null direction {j+1}")
    for a in [0.01,0.05,0.1,0.2,0.4]:
        # constrain projection on v to a; allow the 10 rigid dirs AND the other null dir to relax
        basis=np.column_stack([Perp,other])
        f=lambda c: U(q0+a*v+basis@c)
        r=minimize(f,np.zeros(11),method='BFGS',options={'gtol':1e-13,'maxiter':5000})
        dq=a*v+basis@r.x
        print(f"a={a}: Umin={r.fun:.3e}  |relax|={np.linalg.norm(r.x):.3f}  dy_mid={(dq[4]-dq[1]):+.3f}  dwidth_mid={(dq[3]-dq[0]):+.4f}  max|spring strain|={np.abs(lengths(q0+dq)-L0).max():.2e}")
