"""Independent check of Eq. (12) and of the fluctuation-dissipation relation (Eq. 13) by direct simulation.

Simulates the Langevin equation  W q'' = C q - Gamma q' + xi,  <xi xi^T> = 2 kT Gamma delta(t-t')
(physical coordinates, SI, k = 1/2/100 N/m, rod friction as in analysis_tables.py) with the EXACT
discrete-time propagator of this linear stochastic equation (Van Loan method), so no time-step error.

Checks:
  1. rms channel width (Dx top/mid/bottom) from the simulation  vs  equipartition / Eq. (12)
     (0.826 / 0.401 / 0.601 A), in water and in 10x water -> must agree and be independent of friction.
  2. mean kinetic energy per coordinate = kT/2 (equipartition of velocities, tests the noise strength).
  3. zero modes: the asymmetric shift Dy (which contains the zero modes) is NOT stationary -
     its variance keeps growing (free diffusion), while Dx is stationary.
Run: python langevin_check.py   (writes langevin_check_output.txt; takes a few minutes)
"""
import os, numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "analysis_tables.py"), encoding="utf-8").read()
exec(src.split("# ---------------- 1. frequencies")[0].replace('out = open(os.path.join(HERE, "analysis_output.txt"), "w", encoding="utf-8")', 'out = open(os.devnull, "w")'))
D_ROD = 1.0e-9
def rod_friction(eta):
    G = np.zeros((12, 12))
    for i in range(4):
        L = L_SI[i]; pr = L / D_ROD; t = q0[3*i+2]; u = np.array([np.cos(t), np.sin(t)])
        g_par = 2*np.pi*eta*L/(np.log(pr) - 0.207 + 0.980/pr - 0.133/pr**2)
        g_perp = 4*np.pi*eta*L/(np.log(pr) + 0.839 + 0.185/pr + 0.233/pr**2)
        g_rot = np.pi*eta*L**3/(3*(np.log(pr) - 0.662 + 0.917/pr - 0.050/pr**2))
        G[3*i:3*i+2, 3*i:3*i+2] = g_par*np.outer(u, u) + g_perp*(np.eye(2) - np.outer(u, u))
        G[3*i+2, 3*i+2] = g_rot
    return G

fout = open(os.path.join(HERE, "langevin_check_output.txt"), "w", encoding="utf-8")
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); fout.write(s + "\n")

H = B0.T @ np.diag(k_SI) @ B0          # stiffness, C = -H (SI)
# ---- well-conditioned units: mass-weighted coordinates (paper's q~), time in ps, energy in kT ----
# q~ = sqrt(W) q ; x = q~ / (sqrt(kT)*TU)  ->  x'' = -Ct x - Gt x' + noise,  <noise noise> = 2 Gt delta
TU = 1e-12                                         # time unit: 1 ps
sW = np.sqrt(W)
Ht = (H / np.outer(sW, sW)) * TU**2                # = -C~ in ps^-2 (dimensionless numbers)
G_dx = {k: obs[k] for k in ("dx_top", "dx_mid", "dx_bot")}
G_dy = obs["dy_mid"]
def obs_x(g):                                      # observable g.q in metres from x:  g.q = sum_k g_k/sqrt(W_k) * sqrt(kT)*TU * x_k
    return (g / sW) * np.sqrt(kT) * TU
w, Phi, _ = modes(k_SI)
pred = {k: np.sqrt(sum((g @ Phi[:, n])**2 * kT / w[n]**2 for n in range(2, 12))) for k, g in G_dx.items()}
pr("Conditioning check: eigenvalues of H~ (ps^-2):", np.array2string(np.linalg.eigvalsh(Ht), precision=3))

rng = np.random.default_rng(7)
NTRAJ, DT = 64, 0.05                               # dt in ps
for eta, label, t_burn, t_run in [(0.7e-3, "water", 8000.0, 20000.0), (7e-3, "10x water", 60000.0, 60000.0)]:
    Gt = (rod_friction(eta) / np.outer(sW, sW)) * TU     # Gamma~ in ps^-1
    A = np.block([[np.zeros((12, 12)), np.eye(12)], [-Ht, -Gt]])
    Dn = np.zeros((24, 24)); Dn[12:, 12:] = 2*Gt          # noise covariance rate (kT = 1)
    F = expm(np.block([[-A, Dn], [np.zeros((24, 24)), A.T]]) * DT)   # Van Loan
    Phd = F[24:, 24:].T; Qd = Phd @ F[:24, 24:]; Qd = (Qd + Qd.T)/2
    ev, U = np.linalg.eigh(Qd); Lq = U * np.sqrt(np.clip(ev, 0, None))
    z = np.zeros((24, NTRAJ))
    nb, nr = int(t_burn/DT), int(t_run/DT)
    for _ in range(nb):
        z = Phd @ z + Lq @ rng.standard_normal((24, NTRAJ))
    acc = {k: 0.0 for k in G_dx}; ke = np.zeros(12); dy_var_t = []
    gdx = {k: obs_x(g) for k, g in G_dx.items()}; gdy = obs_x(G_dy)
    dy0 = gdy @ z[:12]
    for s in range(nr):
        z = Phd @ z + Lq @ rng.standard_normal((24, NTRAJ))
        for k, g in gdx.items():
            acc[k] += np.mean((g @ z[:12])**2)
        ke += np.mean(z[12:]**2, axis=1) / 2               # kinetic energy per coordinate in kT
        if s % (nr//6) == 0:
            dy_var_t.append((s*DT/1000, np.var(gdy @ z[:12] - dy0)))
    pr(f"=== {label} (eta = {eta} Pa s), {NTRAJ} trajectories, burn-in {t_burn/1000:.0f} ns, sampling {t_run/1000:.0f} ns, exact propagator dt = {DT} ps")
    for k in G_dx:
        sim = np.sqrt(acc[k]/nr)
        pr(f"  rms {k}: simulation {sim*1e10:.3f} A   Eq.(12) {pred[k]*1e10:.3f} A   ratio {sim/pred[k]:.3f}")
    pr(f"  mean kinetic energy per coordinate / (kT/2): min {np.min(ke/nr/0.5):.3f}, max {np.max(ke/nr/0.5):.3f}")
    pr("  zero-mode diffusion: variance of Dy_mid displacement since start of sampling (A^2) vs time (ns): " +
       ", ".join(f"{t:.0f} ns: {v*1e20:.3g}" for t, v in dy_var_t))
fout.close()
