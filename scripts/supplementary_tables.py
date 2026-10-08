"""Numbers for Supplementary Tables S1-S4 of the revised K2P paper.

Writes  ../supplementary/tables_data.json   (read by the Word generator make_supplementary_docx.js)
and     supplementary_tables_output.txt      (human-readable copy of every number).

Model and methods are identical to analysis_tables.py (same geometry, SI parameters, observables, friction
matrix and random seed), so Tables S2-S4 reproduce analysis_output.txt. Additional here:
  * Table S1: geometry and parameters in SI, springs numbered as in Fig. 1B;
  * Table S4: rigid-limit check (springs 9 and 10 x10 and x100).
Run from the repository root (the folder that contains scripts/):  python scripts/supplementary_tables.py
"""
import os, json
import numpy as np
from scipy.linalg import eigh, eig
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(os.path.dirname(HERE), "supplementary")
txt = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); txt.append(s)

# ---------------- model (notebook cell 12; toy units -> SI as in analysis_tables.py) ----------------
NM = 1e-9
rod = [(-1, 0, np.pi/2 - np.pi/9, 3), (1, 0, np.pi/2 + np.pi/9, 3),
       (-1.1, 0.5, np.pi/2 - np.pi/20, 2), (1.1, 0.5, np.pi/2 + np.pi/20, 2)]
conn = [("r1","r2"),("r1","r2"),("r1","r3"),("r2","r4"),("r1","r3"),("r2","r4"),
        ("r1","f2"),("r3","f1"),("r2","f4"),("r4","f3")]
loc = [[1,1],[-1,-1],[-0.5,-1.0],[-0.5,-1.0],[1.4,0.85],[1.4,0.85],[-1.5,0],[1,0],[-1.5,0],[1,0]]
fp = [(-3,2),(-3,-2),(3,2),(3,-2)]
k_toy = np.array([10,20,1000,1000,10,10,10,10,10,10.])
FK, FM = 0.1, 1.44e-24                     # toy -> SI: spring constant (N/m), mass (kg); length unit = 1 nm
k_SI = k_toy*FK
m_SI = np.array([3, 3, 2, 2.])*FM
L_SI = np.array([r[3] for r in rod])*NM
I_SI = m_SI*L_SI**2/12
W = np.array([w for i in range(4) for w in (m_SI[i], m_SI[i], I_SI[i])])
q0 = np.array([v for r in rod for v in (r[0]*NM, r[1]*NM, r[2])])
kB, T = 1.380649e-23, 310.0; kT = kB*T
DA = 1.66053906660e-27                     # kg per dalton
# Fig. 1B spring labels (notebook cell 17 label positions) for code indices 0..9
fig_label = {0:1, 1:2, 2:9, 3:10, 4:7, 5:8, 6:5, 7:3, 8:6, 9:4}

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
    lam, Phi = eigh(H, np.diag(W)); lam[lam < 0] = 0
    return np.sqrt(lam), Phi, H
w, Phi, H = modes(k_SI)
P = np.zeros((12,12))
for a, b in [(0,1),(1,0),(2,3),(3,2)]:
    P[3*b,3*a] = -1; P[3*b+1,3*a+1] = 1; P[3*b+2,3*a+2] = -1
def parity(v): return float(v @ P @ v/(v @ v))

data = {}
# ================= Table S1 =================
p("="*90); p("TABLE S1  geometry and parameters (SI)")
rods = []
roles = ["inner (pore wall), left", "inner (pore wall), right", "outer, left", "outer, right"]
for i, r in enumerate(rod):
    rods.append(dict(rod=i+1, role=roles[i], x=r[0], y=r[1], theta_deg=np.degrees(r[2]), L_nm=r[3],
                     m_kg=m_SI[i], m_kDa=m_SI[i]/DA/1e3, I_kgm2=I_SI[i]))
    p(f"rod {i+1}: centre ({r[0]:+.2f}, {r[1]:+.2f}) nm, theta {np.degrees(r[2]):.1f} deg, L {r[3]} nm, "
      f"m {m_SI[i]:.3e} kg = {m_SI[i]/DA/1e3:.2f} kDa, I {I_SI[i]:.3e} kg m^2")
# Fixed points labelled as in the paper (Fig. 1): A1 = upper pair, A2 = lower pair; left/right added for uniqueness.
# (code f1..f4 = (-3,2), (-3,-2), (3,2), (3,-2))
ANCH = {1: "A1 (left)", 2: "A2 (left)", 3: "A1 (right)", 4: "A2 (right)"}
anch = [dict(name=ANCH[j+1], code=f"f{j+1}", x=f[0], y=f[1]) for j, f in enumerate(fp)]
for a in anch: p(f"anchoring point {a['name']}: ({a['x']:+.1f}, {a['y']:+.1f}) nm")
L0 = lengths(q0)
def endname(e): return f"rod {e[1:]}" if e[0] == 'r' else f"anchoring point {ANCH[int(e[1:])]}"
springs = []
for j, (a, b) in enumerate(conn):
    pa, pb = pt(q0, a, loc[j][0])/NM, pt(q0, b, loc[j][1])/NM
    springs.append(dict(spring=fig_label[j], code_index=j, end1=endname(a), s1=loc[j][0] if a[0]=='r' else None,
                        p1=[pa[0], pa[1]], end2=endname(b), s2=loc[j][1] if b[0]=='r' else None, p2=[pb[0], pb[1]],
                        L0_nm=L0[j]/NM, k=k_SI[j]))
springs.sort(key=lambda s: s["spring"])
for s in springs:
    p(f"spring {s['spring']:>2} (code {s['code_index']}): {s['end1']} s={s['s1']} at ({s['p1'][0]:+.3f},{s['p1'][1]:+.3f}) -- "
      f"{s['end2']} s={s['s2']} at ({s['p2'][0]:+.3f},{s['p2'][1]:+.3f}); natural length {s['L0_nm']:.3f} nm; k {s['k']:g} N/m")
# residual force at equilibrium with natural length = equilibrium length (unstressed network)
p("natural lengths = equilibrium lengths -> network unstressed (zero force at equilibrium by construction)")
data["S1"] = dict(rods=rods, anchors=anch, springs=springs)

# ================= Table S2 =================
def end_jac(i, s):
    x, y, t = q0[3*i:3*i+3]; h = s*L_SI[i]/2
    J = np.zeros((2,12)); J[0,3*i] = 1; J[1,3*i+1] = 1; J[0,3*i+2] = -h*np.sin(t); J[1,3*i+2] = h*np.cos(t)
    return J
obs = {}
for nm, s in [("top",1),("mid",0),("bot",-1)]:
    J1, J2 = end_jac(0, s), end_jac(1, s)
    obs["dx_"+nm] = J2[0] - J1[0]                     # width change, right minus left (Eq. 8)
    obs["dy_"+nm] = J2[1] - J1[1]                     # asymmetric shift (Eq. 9)
e = np.eye(12)
obs["sway"] = (e[0] + e[3])/2                          # common sideways displacement of the walls (centres)
obs["piston"] = (e[1] + e[4])/2                        # common axial displacement of the walls (centres)
p("="*90); p("TABLE S2  thermal rms amplitudes at 310 K (A; angles in deg)")
cols = list(obs) + ["th_in", "th_out"]
p("mode " + " ".join(f"{c:>7}" for c in cols))
S2 = []; acc = {c: 0.0 for c in cols}
for n in range(2, 12):
    row = {c: abs(g @ Phi[:,n])*np.sqrt(kT)/w[n]*1e10 for c, g in obs.items()}
    row["th_in"] = abs(Phi[2,n])*np.sqrt(kT)/w[n]*180/np.pi
    row["th_out"] = abs(Phi[8,n])*np.sqrt(kT)/w[n]*180/np.pi
    for c in cols: acc[c] += row[c]**2
    row["mode"] = n+1; S2.append(row)
    p(f"{n+1:>4} " + " ".join(f"{row[c]:>7.3f}" for c in cols))
tot = {c: np.sqrt(acc[c]) for c in cols}
p(" all " + " ".join(f"{tot[c]:>7.3f}" for c in cols))
data["S2"] = dict(rows=S2, total=tot)

# ================= Table S3 =================
D_ROD = 1.0e-9
def rod_friction(eta):
    G = np.zeros((12,12))
    for i in range(4):
        L = L_SI[i]; pr = L/D_ROD; t = q0[3*i+2]; u = np.array([np.cos(t), np.sin(t)])
        g_par = 2*np.pi*eta*L/(np.log(pr) - 0.207 + 0.980/pr - 0.133/pr**2)
        g_perp = 4*np.pi*eta*L/(np.log(pr) + 0.839 + 0.185/pr + 0.233/pr**2)
        g_rot = np.pi*eta*L**3/(3*(np.log(pr) - 0.662 + 0.917/pr - 0.050/pr**2))
        G[3*i:3*i+2,3*i:3*i+2] = g_par*np.outer(u,u) + g_perp*(np.eye(2) - np.outer(u,u))
        G[3*i+2,3*i+2] = g_rot
    return G
p("="*90); p("TABLE S3  damping")
# Modal approximation: zeta'' + gamma_n zeta' + omega_n^2 zeta = 0, gamma_n = v_n^T Gamma~ v_n = omega_n/Q_n.
#   Overdamped (Q_n < 1/2): slow root s = -(gamma_n - sqrt(gamma_n^2 - 4 omega_n^2))/2, tau_n = 1/|s|  (~ 1/(Q_n omega_n) for Q_n << 1).
# Exact: eigenvalues s of the 24x24 matrix [[0, 1], [-M^-1 H, -M^-1 Gamma]] (friction coupling between modes kept).
#   Real s < 0 -> relaxation time 1/|s|. Each real eigenvector (x, s x) is projected on the normal modes:
#   weight_n = (v_n . sqrt(W) x)^2 / |sqrt(W) x|^2 (sums to 1). Each mode n is assigned the slow root with the largest weight_n.
etas = [0.7e-3, 7e-3, 7e-2]
Q = {}; exact = {}
Mm = np.diag(W); V = np.sqrt(W)[:, None]*Phi           # orthonormal mass-weighted eigenvectors
def tau_single(g, wn):
    return 2/(g - np.sqrt(g*g - 4*wn*wn)) if g > 2*wn else None
for eta in etas:
    G = rod_friction(eta)
    gam = [Phi[:,n] @ G @ Phi[:,n] for n in range(12)]
    Q[eta] = [w[n]/gam[n] for n in range(2,12)]
    A = np.block([[np.zeros((12,12)), np.eye(12)], [-np.linalg.solve(Mm,H), -np.linalg.solve(Mm,G)]])
    ev, R = eig(A)
    big = np.abs(ev).max()
    osc = ev[np.abs(ev.imag) > 1e-6*big]
    roots = []                                           # (tau_ps, weights) for real, nonzero roots
    for i, sv in enumerate(ev):
        if abs(sv.imag) > 1e-6*big or abs(sv) < 1e-9*big: continue   # oscillatory, or numerical zero root of modes 1-2
        x = np.sqrt(W)*np.real(R[:12, i]); x /= np.linalg.norm(x)
        roots.append((1e12/abs(sv.real), (V.T @ x)**2))
    roots.sort(key=lambda r: -r[0])                      # slowest first
    slow = roots[:8]                                     # the slow branch (modes 3-10); fast roots are < 0.4 ps
    assign = {}
    for n in range(2, 10):
        t, wt = max(slow, key=lambda r: r[1][n])
        assign[n+1] = dict(tau_ps=t, weight=float(wt[n]))
    exact[eta] = dict(n_osc_pairs=len(osc)//2, osc_w=sorted(set(np.round(np.abs(osc.imag), -9).tolist())),
                      slow_tau_ps=[r[0] for r in slow],
                      slow_weights=[[round(float(x), 3) for x in r[1]] for r in slow], assign=assign,
                      tau_modal_ps={n+1: (1e12*tau_single(gam[n], w[n]) if tau_single(gam[n], w[n]) else None) for n in range(2,12)})
    p(f"eta {eta*1e3:g} mPa s: modal Q " + ", ".join(f"m{n+3}:{q:.3f}" for n, q in enumerate(Q[eta])))
    p(f"   exact: {len(osc)//2} oscillatory pairs {exact[eta]['osc_w']}; eight slowest relaxation roots:")
    for t, wt in slow:
        p(f"      tau {t:9.4g} ps; mode weights m1..m12: " + " ".join(f"{x:.2f}" for x in wt))
    for n in range(3, 11):
        tm = exact[eta]["tau_modal_ps"][n]
        p(f"   mode {n:>2}: modal (single-mode root) {tm:.4g} ps, 1/(Q w) {1e12*gam[n-1]/w[n-1]**2:.4g} ps; "
          f"exact root with largest weight {assign[n]['tau_ps']:.4g} ps (weight {assign[n]['weight']:.2f})")
G1 = rod_friction(1.0)
eta_crit = [w[n]/(0.5*(Phi[:,n] @ G1 @ Phi[:,n])) for n in range(2,12)]
Q_sphere = [4.32e-24*w[n]/2e-11 for n in range(2,12)]
S3 = []
ew = exact[0.7e-3]
for j, n in enumerate(range(2,12)):
    a_ = ew["assign"].get(n+1)
    row = dict(mode=n+1, omega=w[n], Q_sphere=Q_sphere[j], Q_water=Q[0.7e-3][j], Q_10=Q[7e-3][j], Q_100=Q[7e-2][j],
               eta_crit_mPas=eta_crit[j]*1e3, tau_modal_ps=ew["tau_modal_ps"][n+1],
               tau_exact_ps=a_["tau_ps"] if a_ else None, exact_weight=a_["weight"] if a_ else None)
    S3.append(row)
    p(f"mode {n+1:>2}: omega {w[n]:.3e}; Q sphere {row['Q_sphere']:.2f}; Q water/10x/100x {row['Q_water']:.3f}/{row['Q_10']:.4f}/{row['Q_100']:.5f}; "
      f"eta_crit {row['eta_crit_mPas']:.2g} mPa s; tau modal {row['tau_modal_ps']}; tau exact {row['tau_exact_ps']} (weight {row['exact_weight']})")
data["S3"] = dict(rows=S3, exact={f"{e*1e3:g}": exact[e] for e in etas})

# ================= Table S4 =================
p("="*90); p("TABLE S4  sensitivity to spring constants")
Mh = np.sqrt(W)
ref_V = (Mh[:,None]*Phi)[:, 2:]
def match(k):
    w2, Ph2, _ = modes(k); V2 = (Mh[:,None]*Ph2)[:, 2:]
    O = np.abs(ref_V.T @ V2); r, c = linear_sum_assignment(-O)
    return w2[2:][c], O[r, c], int(np.sum(w2 < 1e-6*w2.max()))
groups = [("Springs 1, 3-8 (1 N/m)", [0,4,5,6,7,8,9], (0.5, 1.5)), ("Spring 2 (2 N/m)", [1], (0.5, 1.5)),
          ("Springs 9, 10 (100 N/m)", [2,3], (0.5, 1.5, 10, 100))]
S4a = []
for g, idx, facs in groups:
    for f in facs:
        k = k_SI.copy(); k[idx] *= f
        w2, ov, nz = match(k); rat = w2/w[2:]
        S4a.append(dict(group=g, factor=f, ratios=rat.tolist(), min_overlap=float(ov.min()), min_overlap_mode=int(np.argmin(ov)+3),
                        zero_modes=nz, max_change_3_10=float(np.abs(rat[:8]-1).max())))
        p(f"{g} x{f:g}: zero modes {nz}; ratios m3..m12 " + " ".join(f"{x:.4f}" for x in rat)
          + f"; min overlap {ov.min():.5f} (mode {np.argmin(ov)+3}); max |change| modes 3-10 {100*np.abs(rat[:8]-1).max():.3f}%")
rng = np.random.default_rng(1)
pairs = [(0,),(1,),(2,3),(4,5),(6,8),(7,9)]
groups_deg = {"8-10": [5,6,7], "11-12": [8,9]}
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
        ratios.append(w2[2:][c]/w[2:]); overl.append(O[r,c]); nzs.append(int(np.sum(w2 < 1e-6*w2.max())))
        for g, idx in groups_deg.items():
            sub[g].append(np.linalg.svd(ref_V[:,idx].T @ V2[:,idx], compute_uv=False).min())
    ratios, overl = np.array(ratios), np.array(overl)
    rows = []
    for j in range(10):
        rows.append(dict(mode=j+3, median=float(np.median(ratios[:,j])), p5=float(np.percentile(ratios[:,j],5)),
                         p95=float(np.percentile(ratios[:,j],95)), pct09=float(100*np.mean(overl[:,j]>0.9)),
                         pct08=float(100*np.mean(overl[:,j]>0.8))))
        r_ = rows[-1]; p(f"   mode {j+3:>2}: ratio {r_['median']:.2f} [{r_['p5']:.2f}, {r_['p95']:.2f}]; overlap>0.9 {r_['pct09']:.0f}%; >0.8 {r_['pct08']:.0f}%")
    subs = {g: float(100*np.mean(np.array(v) > 0.9)) for g, v in sub.items()}
    p(f"   subspaces preserved (>0.9): {subs}; zero modes in every sample: {sorted(set(nzs))}")
    return dict(rows=rows, subspace_pct=subs, zero_modes=sorted(set(nzs)))
p("(A) symmetry-preserving, 2000 samples"); A_ = run(True)
p("(B) each spring independently, 2000 samples"); B_ = run(False)
data["S4"] = dict(groups=S4a, sym=A_, indep=B_)

os.makedirs(OUTDIR, exist_ok=True)
json.dump(data, open(os.path.join(OUTDIR, "tables_data.json"), "w"), indent=1, default=float)
open(os.path.join(HERE, "supplementary_tables_output.txt"), "w", encoding="utf-8").write("\n".join(txt) + "\n")
