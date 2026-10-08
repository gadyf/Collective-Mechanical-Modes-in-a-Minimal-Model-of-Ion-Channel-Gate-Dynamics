"""SI-unit analysis of the 4-rod / 10-spring K2P gate model (geometry identical to the notebook, cell 12).

Produces, for the parameter set used for Table 1 (toy k x 0.1 -> 1, 2, 100 N/m; masses 4.32e-24 / 2.88e-24 kg;
lengths 3 / 2 nm):
  1. normal-mode frequencies (check against Table 1) and mirror parity of each mode
  2. equipartition (310 K) rms amplitudes of channel width / asymmetric shift (top, mid, bottom) and rod rotations
  3. damping: (a) single-sphere Stokes estimate of Q, (b) slender-rod friction matrix -> modal Q and
     exact damped eigenvalues (state-space), for several viscosities
  4. sensitivity of frequencies / mode identity to the spring constants
Run:  python analysis_tables.py  (writes analysis_output.txt next to this script)
Paper mode numbers = ascending order including the two zero modes (1, 2), as in the manuscript.
"""
import os, sys
import numpy as np
from scipy.linalg import eigh, eig
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, "analysis_output.txt"), "w", encoding="utf-8")
def p(*a):
    s = " ".join(str(x) for x in a); print(s); out.write(s + "\n")

# ---------------- model (toy geometry, lengths in nm) ----------------
NM = 1e-9
rod = [(-1, 0, np.pi/2 - np.pi/9, 3), (1, 0, np.pi/2 + np.pi/9, 3),
       (-1.1, 0.5, np.pi/2 - np.pi/20, 2), (1.1, 0.5, np.pi/2 + np.pi/20, 2)]
conn = [("r1","r2"),("r1","r2"),("r1","r3"),("r2","r4"),("r1","r3"),("r2","r4"),
        ("r1","f2"),("r3","f1"),("r2","f4"),("r4","f3")]
loc = [[1,1],[-1,-1],[-0.5,-1.0],[-0.5,-1.0],[1.4,0.85],[1.4,0.85],[-1.5,0],[1,0],[-1.5,0],[1,0]]
fp = [(-3,2),(-3,-2),(3,2),(3,-2)]
k_toy = np.array([10,20,1000,1000,10,10,10,10,10,10.])
FK = 0.1                                   # model unit of spring constant, N/m (model value 10 -> 1 N/m)
k_SI = k_toy * FK                          # N/m  -> 1, 2, 100, 100, 1, ...
m_SI = np.array([4.32e-24, 4.32e-24, 2.88e-24, 2.88e-24])   # kg (model masses 3, 2 x 1.44e-24 kg)
L_SI = np.array([r[3] for r in rod]) * NM
I_SI = m_SI * L_SI**2 / 12
W = np.array([w for i in range(4) for w in (m_SI[i], m_SI[i], I_SI[i])])
q0 = np.array([v for i, r in enumerate(rod) for v in (r[0]*NM, r[1]*NM, r[2])])
kB, T = 1.380649e-23, 310.0
kT = kB*T
names_paper_spring = {0:"k1",1:"k2",2:"k9",3:"k10",4:"k7",5:"k8",6:"k5",7:"k3",8:"k6",9:"k4"}  # code index -> Fig.1B label (from notebook cell 17 label positions)

def pt(q, end, d):
    if end[0] == 'f': return np.array(fp[int(end[1:])-1], float)*NM
    i = int(end[1:])-1; x, y, t = q[3*i:3*i+3]; d = d*NM
    return np.array([x + d*np.cos(t), y + d*np.sin(t)])
def lengths(q): return np.array([np.linalg.norm(pt(q,a,loc[j][0]) - pt(q,b,loc[j][1])) for j,(a,b) in enumerate(conn)])
def Bmat(q):
    B = np.zeros((10,12))
    for c in range(12):
        h = 1e-6*NM if c % 3 != 2 else 1e-6
        e = np.zeros(12); e[c] = h; B[:,c] = (lengths(q+e) - lengths(q-e))/(2*h)
    return B
B0 = Bmat(q0)
def modes(k):
    H = B0.T @ np.diag(k) @ B0
    lam, Phi = eigh(H, np.diag(W))          # Phi^T M Phi = 1
    lam[lam < 0] = 0
    return np.sqrt(lam), Phi, H
w, Phi, H = modes(k_SI)

# mirror operator (x -> -x): rod1<->rod2, rod3<->rod4, dx->-dx, dy->dy, dtheta->-dtheta
P = np.zeros((12,12))
for a, b in [(0,1),(1,0),(2,3),(3,2)]:
    P[3*b,3*a] = -1; P[3*b+1,3*a+1] = 1; P[3*b+2,3*a+2] = -1
def parity(v): return float(v @ P @ v / (v @ v))

# geometric observables (linear in dq): width and y-shift at top/mid/bottom of inner rods (Eqs. 8-9)
def end_jac(i, s):
    x, y, t = q0[3*i:3*i+3]; h = s*L_SI[i]/2
    J = np.zeros((2,12)); J[0,3*i] = 1; J[1,3*i+1] = 1; J[0,3*i+2] = -h*np.sin(t); J[1,3*i+2] = h*np.cos(t)
    return J
obs = {}
for nm, s in [("top",1),("mid",0),("bot",-1)]:
    J = end_jac(0, s), end_jac(1, s)
    obs["dx_"+nm] = J[0][0] - J[1][0]      # x1 - x2 as in Eq. 8 (sign irrelevant for rms)
    obs["dy_"+nm] = J[0][1] - J[1][1]      # y1 - y2 as in Eq. 9

# ---------------- 1. frequencies ----------------
table1 = [2.99e4,3.99e5,5.58e10,1.56e11,4.50e11,8.19e11,9.24e11,1.52e12,1.58e12,1.69e12,1.24e13,1.25e13]
p("="*100); p("1. FREQUENCIES (SI, k = 1/2/100 N/m), mirror parity (+1 symmetric = width-type, -1 antisymmetric = shift-type)")
p(f"{'mode':>4} {'omega (rad/s)':>14} {'Table 1':>10} {'f (GHz)':>9} {'T (ps)':>9} {'parity':>7}")
for n in range(12):
    T_ps = 2*np.pi/w[n]*1e12 if w[n] > 1 else float('inf')
    p(f"{n+1:>4} {w[n]:>14.3e} {table1[n]:>10.2e} {w[n]/2/np.pi/1e9:>9.2f} {T_ps:>9.2f} {parity(Phi[:,n]):>+7.2f}")
Bs = B0 / np.sqrt(W)[None,:]            # mass-weighted compatibility matrix (consistent units)
sv = np.linalg.svd(Bs, compute_uv=False)
p("Modes 1-2: exact zero modes; rank of mass-weighted compatibility matrix = %d (singular values: %s)." % (np.sum(sv > 1e-8*sv[0]), ", ".join(f"{x:.2e}" for x in sv)))
p("Mode 2 printed value is round-off (notebook gave 3.99e5, this run %.1e): no physical frequency." % w[1])
p("If the model values 10/20/1000 were read as N/m instead, every omega above would be multiplied by sqrt(10) = 3.162.")

# ---------------- 2. thermal amplitudes ----------------
p("="*100); p("2. EQUIPARTITION rms AMPLITUDES at 310 K (exact, from eigenvectors). Valid with or without damping.")
p("   rms of observable g in mode n = |g . phi_n| * sqrt(kT) / omega_n   (phi normalised phi^T M phi = 1)")
hdr = f"{'mode':>4} " + " ".join(f"{k:>7}" for k in obs) + "  " + " ".join(f"th{i+1}(deg)" for i in range(4))
p("   lengths in Angstrom"); p(hdr)
amp = {}
for n in range(2, 12):
    a = {k: abs(g @ Phi[:,n])*np.sqrt(kT)/w[n]*1e10 for k, g in obs.items()}
    th = [abs(Phi[3*i+2,n])*np.sqrt(kT)/w[n]*180/np.pi for i in range(4)]
    amp[n] = a
    p(f"{n+1:>4} " + " ".join(f"{a[k]:>7.3f}" for k in obs) + "  " + " ".join(f"{x:>9.3f}" for x in th))
tot = {k: np.sqrt(sum(amp[n][k]**2 for n in amp))*1e0 for k in obs}
p("all " + " ".join(f"{tot[k]:>7.3f}" for k in obs) + "   (modes 3-12 combined; zero modes excluded - undefined in harmonic model)")
# largest single-rod centre-of-mass displacement per mode
p("   max rod-centre rms displacement per mode (A) and upper bound sqrt(kT/(m_light w^2)) (A):")
for n in range(2,12):
    d = max(np.hypot(Phi[3*i,n],Phi[3*i+1,n]) for i in range(4))*np.sqrt(kT)/w[n]*1e10
    p(f"   mode {n+1:>2}: {d:6.3f}   bound {np.sqrt(kT/(2.88e-24*w[n]**2))*1e10:6.3f}")

# ---------------- 3. damping ----------------
p("="*100); p("3. DAMPING")
p("3a. Single-sphere estimate Q = m*omega/gamma, Stokes sphere gamma = 6*pi*eta*a")
for lab, m_, gam in [("m=4.32e-24 kg, gamma=2e-11 kg/s", 4.32e-24, 2e-11),
                     ("a=0.5 nm, eta=0.7 mPa s, m=2.88e-24 kg", 2.88e-24, 6*np.pi*0.7e-3*0.5e-9)]:
    p("   " + lab + ": " + ", ".join(f"m{n+1}:{m_*w[n]/gam:.2f}" for n in range(2,12)))

D_ROD = 1.0e-9   # effective helix diameter (m)
def rod_friction(eta):
    G = np.zeros((12,12))
    for i in range(4):
        L = L_SI[i]; pr = L/D_ROD; t = q0[3*i+2]; u = np.array([np.cos(t), np.sin(t)])
        # Tirado & Garcia de la Torre (1979, 1980) short-cylinder corrections
        g_par = 2*np.pi*eta*L/(np.log(pr) - 0.207 + 0.980/pr - 0.133/pr**2)
        g_perp = 4*np.pi*eta*L/(np.log(pr) + 0.839 + 0.185/pr + 0.233/pr**2)
        g_rot = np.pi*eta*L**3/(3*(np.log(pr) - 0.662 + 0.917/pr - 0.050/pr**2))
        G[3*i:3*i+2,3*i:3*i+2] = g_par*np.outer(u,u) + g_perp*(np.eye(2) - np.outer(u,u))
        G[3*i+2,3*i+2] = g_rot
    return G
p("3b. Slender-rod friction matrix (Tirado-Garcia de la Torre cylinders, d = 1 nm), modal Q_n = omega_n / (phi_n^T Gamma phi_n)")
p("    and exact damped dynamics M q'' + Gamma q' + H q = 0 (24x24 state space): 'osc' = complex pair (underdamped), 'relax' = real roots")
etas = [(0.7e-3,"water 0.7 mPa s"),(7e-3,"10x water"),(7e-2,"100x water (protein/lipid-like)"),(0.7,"1000x water (lipid)")]
Mm = np.diag(W)
for eta, lab in etas:
    G = rod_friction(eta)
    Qn = [w[n]/(Phi[:,n] @ G @ Phi[:,n]) for n in range(2,12)]
    A = np.block([[np.zeros((12,12)), np.eye(12)], [-np.linalg.solve(Mm,H), -np.linalg.solve(Mm,G)]])
    ev = eig(A, right=False)
    osc = ev[np.abs(ev.imag) > 1e-6*np.abs(ev).max()]
    osc_w = np.sort(np.unique(np.round(np.abs(osc.imag), -6)))
    p(f"   eta = {eta:g} Pa s ({lab}):")
    p("      modal Q: " + ", ".join(f"m{n+3}:{q:.2f}" for n, q in enumerate(Qn)))
    p("      exact: %d oscillatory pairs; damped frequencies (rad/s): %s" % (len(osc)//2, ", ".join(f"{x:.2e}" for x in osc_w)))
    relax = ev[np.abs(ev.imag) <= 1e-6*np.abs(ev).max()].real
    relax = np.sort(-relax[relax < -1])
    if len(relax): p("      relaxation times (ps) of overdamped roots: " + ", ".join(f"{1e12/x:.2g}" for x in relax))
# critical viscosity per mode (modal approx.): Q = 0.5
p("   viscosity at which each mode becomes critically damped (modal Q = 0.5), Pa s:")
G1 = rod_friction(1.0)
p("      " + ", ".join(f"m{n+1}:{w[n]/(0.5*(Phi[:,n] @ G1 @ Phi[:,n])):.2g}" for n in range(2,12)))

# ---------------- 4. sensitivity ----------------
p("="*100); p("4. SENSITIVITY TO SPRING CONSTANTS")
p("   Exact: k -> s*k scales all omega by sqrt(s) and leaves mode shapes unchanged. Only ratios matter.")
ref_w, ref_Phi = w, Phi
Mh = np.sqrt(W)
ref_V = (Mh[:,None]*ref_Phi)[:, 2:]   # mass-weighted unit vectors of modes 3-12
def match(k):
    w2, Ph2, _ = modes(k)
    V2 = (Mh[:,None]*Ph2)[:, 2:]
    O = np.abs(ref_V.T @ V2)
    r, c = linear_sum_assignment(-O)
    return w2[2:][c], O[r, c], np.sum(w2 < 1e-6*w2.max())
groups = {"weak springs k1,k3-k8 (1 N/m)": [0,4,5,6,7,8,9], "k2 (2 N/m)": [1], "stiff k9,k10 (100 N/m)": [2,3]}
for g, idx in groups.items():
    for f in (0.5, 1.5):
        k = k_SI.copy(); k[idx] *= f
        w2, ov, nz = match(k)
        p(f"   {g} x{f}: zero modes={nz}; freq ratio m3..m12: " + " ".join(f"{x:.2f}" for x in w2/ref_w[2:]) +
          "; min overlap " + f"{ov.min():.2f} (mode {np.argmin(ov)+3})")
rng = np.random.default_rng(1)
pairs = [(0,),(1,),(2,3),(4,5),(6,8),(7,9)]          # mirror-partner springs (code indices)
groups_deg = {"8-10": [5,6,7], "11-12": [8,9]}       # near-degenerate groups (indices into modes 3..12)
def run(sym, NS=2000):
    ratios, overl, nzs, sub = [], [], [], {g: [] for g in groups_deg}
    for _ in range(NS):
        if sym:
            k = k_SI.copy()
            for pr in pairs: k[list(pr)] *= rng.uniform(0.5, 1.5)
        else:
            k = k_SI*rng.uniform(0.5, 1.5, 10)
        w2, Ph2, _ = modes(k); V2 = (Mh[:,None]*Ph2)[:, 2:]
        O = np.abs(ref_V.T @ V2); r, c = linear_sum_assignment(-O)
        ratios.append(w2[2:][c]/ref_w[2:]); overl.append(O[r,c]); nzs.append(int(np.sum(w2 < 1e-6*w2.max())))
        for g, idx in groups_deg.items():
            # subspace overlap: smallest singular value of ref_V[:,idx]^T V2[:,same positions in sorted order]
            sv_ = np.linalg.svd(ref_V[:,idx].T @ V2[:,idx], compute_uv=False); sub[g].append(sv_.min())
    ratios, overl = np.array(ratios), np.array(overl)
    p(f"   {'mode':>4} {'freq ratio median [5%,95%]':>30} {'% overlap>0.9':>14} {'% overlap>0.8':>14}")
    for j in range(10):
        p(f"   {j+3:>4} {np.median(ratios[:,j]):>14.2f} [{np.percentile(ratios[:,j],5):.2f},{np.percentile(ratios[:,j],95):.2f}] "
          f"{100*np.mean(overl[:,j]>0.9):>12.0f}% {100*np.mean(overl[:,j]>0.8):>13.0f}%")
    for g in groups_deg:
        a = np.array(sub[g]); p(f"   subspace of modes {g}: min principal overlap median {np.median(a):.2f}, >0.9 in {100*np.mean(a>0.9):.0f}% of samples")
    p(f"   zero modes in every sample: {set(nzs)}")
p("   (A) SYMMETRY-PRESERVING: mirror-partner springs scaled together, each pair x U(0.5,1.5), 2000 samples")
run(True)
p("   (B) SYMMETRY-BREAKING: each of the 10 springs independently x U(0.5,1.5), 2000 samples")
run(False)
out.close()
