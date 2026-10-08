"""Full decomposition of every normal mode (SI model, k = 1/2/100 N/m), to characterise the motion seen in the animations.

For each mode n (paper numbering, ascending omega):
  * mirror parity  (+1: mirror-symmetric, -1: mirror-antisymmetric)
  * kinetic-energy (= mass-weighted eigenvector) fractions in: inner rods (x, y, theta), outer rods (x, y, theta)
  * thermal (310 K) rms of inner-rod collective coordinates:
        width   Dx = x2 - x1          (opening/closing, symmetric)
        sway    Sx = (x1 + x2)/2      (whole channel moves sideways, antisymmetric)
        piston  Sy = (y1 + y2)/2      (both walls move axially together, symmetric)
        shift   Dy = y2 - y1          (walls slide axially against each other, antisymmetric)
        tilt    th1, th2 (deg): counter-tilt (mirror, symmetric) or co-tilt (same sense, antisymmetric)
    and outer-rod tilt (deg)
  * thermal rms of the channel observables of the paper (Eqs. 8-9) at top/mid/bottom of the inner rods
Run: python mode_character.py   (prints and writes mode_character_output.txt)
"""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "analysis_tables.py"), encoding="utf-8").read()
# reuse model definitions only (everything up to the frequency printout)
exec(src.split("# ---------------- 1. frequencies")[0].replace('out = open(os.path.join(HERE, "analysis_output.txt"), "w", encoding="utf-8")', 'out = open(os.devnull, "w")'))
fout = open(os.path.join(HERE, "mode_character_output.txt"), "w", encoding="utf-8")
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); fout.write(s + "\n")

sq = np.sqrt(kT)
A = 1e10; DEG = 180/np.pi
pr("Mode decomposition, k = 1/2/100 N/m, T = 310 K. Rods: 1 inner-left, 2 inner-right (long, channel walls), 3 outer-left, 4 outer-right.")
pr("KE fractions: share of the mode's kinetic energy in each group (sums to 1).")
pr(f"{'mode':>4} {'parity':>6} | {'in_x':>5} {'in_y':>5} {'in_th':>5} | {'out_x':>5} {'out_y':>5} {'out_th':>6} || "
   f"{'width':>6} {'sway':>6} {'piston':>6} {'shift':>6} {'th_in':>6} {'th_out':>6} || {'dx_top':>6} {'dx_bot':>6} {'dy_top':>6} {'dy_bot':>6}")
pr(" "*13 + "(fractions)" + " "*32 + "(thermal rms: lengths in Angstrom, angles in degrees)")
for n in range(12):
    phi = Phi[:, n]; u2 = W*phi**2; u2 /= u2.sum()
    par = parity(phi)
    f = lambda idx: u2[idx].sum()
    fr = [f([0,3]), f([1,4]), f([2,5]), f([6,9]), f([7,10]), f([8,11])]
    if w[n] < 1e-6*w.max():
        pr(f"{n+1:>4} {par:>+6.0f} | " + " ".join(f"{x:>5.2f}" for x in fr[:3]) + " | " + " ".join(f"{x:>5.2f}" for x in fr[3:]) +
           " || zero mode: thermal amplitude undefined in the harmonic model")
        continue
    s = sq/w[n]
    width = abs(phi[3]-phi[0])*s*A; sway = abs(phi[0]+phi[3])/2*s*A
    piston = abs(phi[1]+phi[4])/2*s*A; shift = abs(phi[4]-phi[1])*s*A
    th_in = abs(phi[2])*s*DEG; th_out = abs(phi[8])*s*DEG
    o = {k: abs(g @ phi)*s*A for k, g in obs.items()}
    pr(f"{n+1:>4} {par:>+6.0f} | " + " ".join(f"{x:>5.2f}" for x in fr[:3]) + " | " + " ".join(f"{x:>5.2f}" for x in fr[3:]) + " || " +
       " ".join(f"{x:>6.3f}" for x in (width, sway, piston, shift, th_in, th_out)) + " || " +
       " ".join(f"{o[k]:>6.3f}" for k in ("dx_top","dx_bot","dy_top","dy_bot")))
pr("")
pr("Sign relations (exact, from mirror symmetry):")
pr("  parity +1: dx2 = -dx1, dy2 = dy1, dth2 = -dth1  -> width, piston, counter-tilt;  sway = shift = 0")
pr("  parity -1: dx2 =  dx1, dy2 = -dy1, dth2 =  dth1 -> sway, shift, co-tilt;         width = piston = 0")
# relative phase of inner vs outer tilt and of inner translation vs tilt (for descriptive labels)
pr("")
pr("END-POINT VIEW (what the animation shows). Per rod: rms displacement of rod ends split into")
pr("  T = centre-of-mass translation |d_cm|,  R = rotational end displacement (L/2)|dtheta|.")
pr("  'anim' = peak amplitude in the notebook animation units (u/sqrt(W_toy)*0.2, toy length units);")
pr("  'thermal' = rms at 310 K in Angstrom.  Inner = rods 1,2 (long, channel walls); outer = rods 3,4 (short).")
pr(f"{'mode':>4} | {'anim: inner T':>13} {'inner R':>8} {'outer T':>8} {'outer R':>8} | {'thermal: inner T':>16} {'inner R':>8} {'outer T':>8} {'outer R':>8} | inner T dir (deg from x)")
Wtoy = np.array([w_ for i in range(4) for w_ in ((3,3,2,2)[i],(3,3,2,2)[i],(3,3,2,2)[i]*np.array([r[3] for r in rod])[i]**2/12)])
for n in range(12):
    phi = Phi[:, n]; u = phi*np.sqrt(W); u /= np.linalg.norm(u)
    da = u/np.sqrt(Wtoy)*0.2                      # animation displacement (toy units / rad)
    if w[n] < 1e-6*w.max():
        th = None
    else:
        th = phi*np.sqrt(kT)/w[n]                 # thermal rms (m / rad)
    def tr(v, i, Lh): return np.hypot(v[3*i], v[3*i+1]), Lh*abs(v[3*i+2])
    aT_in, aR_in = tr(da, 0, rod[0][3]/2); aT_out, aR_out = tr(da, 2, rod[2][3]/2)
    ang = np.degrees(np.arctan2(da[1], da[0])) % 180
    if th is None:
        pr(f"{n+1:>4} | {aT_in:>13.3f} {aR_in:>8.3f} {aT_out:>8.3f} {aR_out:>8.3f} | {'zero mode':>16} {'':>8} {'':>8} {'':>8} | {ang:5.0f}")
    else:
        tT_in, tR_in = tr(th, 0, L_SI[0]/2); tT_out, tR_out = tr(th, 2, L_SI[2]/2)
        pr(f"{n+1:>4} | {aT_in:>13.3f} {aR_in:>8.3f} {aT_out:>8.3f} {aR_out:>8.3f} | {tT_in*1e10:>16.3f} {tR_in*1e10:>8.3f} {tT_out*1e10:>8.3f} {tR_out*1e10:>8.3f} | {ang:5.0f}")
pr("(inner rods are at 70/110 deg to x: a direction of ~70 or ~110 deg = motion along the rod axis; ~0/180 = horizontal)")
fout.close()

# ---------------- pivot / axial-perpendicular decomposition ----------------
fout = open(os.path.join(HERE, "mode_character_output.txt"), "a", encoding="utf-8")
pr("")
pr("AXIAL / PERPENDICULAR / PIVOT VIEW (rods 1 and 3; rods 2 and 4 follow by mirror symmetry).")
pr("  a = CM displacement along the rod axis (sliding), p = CM displacement perpendicular to the axis, r = (L/2)*dtheta (end displacement from rotation),")
pr("  pivot s* = point on the rod axis with zero perpendicular displacement, in units of L/2 measured from the centre (+ = upper end).")
pr("  |s*| < 1: rod visibly tilts about a point on itself;  |s*| >> 1: mostly sideways translation.  Values in animation units (toy length).")
pr(f"{'mode':>4} | {'inner a':>8} {'p':>7} {'r':>7} {'pivot':>7} | {'outer a':>8} {'p':>7} {'r':>7} {'pivot':>7}")
for n in range(12):
    phi = Phi[:, n]; u = phi*np.sqrt(W); u /= np.linalg.norm(u)
    da = u/np.sqrt(Wtoy)*0.2
    row = []
    for i in (0, 2):
        t = rod[i][2]; ax = np.array([np.cos(t), np.sin(t)]); nn = np.array([-np.sin(t), np.cos(t)])
        d = da[3*i:3*i+2]; dth = da[3*i+2]; Lh = rod[i][3]/2
        a = d @ ax; pp = d @ nn; r = Lh*dth
        piv = (-pp/dth)/Lh if abs(dth) > 1e-12 else np.inf
        row += [a, pp, r, piv]
    pr(f"{n+1:>4} | {row[0]:>+8.3f} {row[1]:>+7.3f} {row[2]:>+7.3f} {row[3]:>+7.2f} | {row[4]:>+8.3f} {row[5]:>+7.3f} {row[6]:>+7.3f} {row[7]:>+7.2f}")
fout.close()
