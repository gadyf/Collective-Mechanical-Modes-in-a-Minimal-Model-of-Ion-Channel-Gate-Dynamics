"""Check that each published animation supplementary/Animation_SNN_mode_N.gif (N = 3..12) shows normal mode N.
(These are the notebook GIFs rods_animation_mode_N_one_period.gif, renamed by supplementary_animations.py.)

Method (independent of the notebook): the rod positions are measured from the GIF pixels and compared
with the normal modes computed by an independent re-implementation of the model (as check_model.py).
  * Each rod is drawn in its own matplotlib colour (C0..C3). For a frame, the rod centre is the centroid
    of its pixels and the rod angle is the principal axis of those pixels (both in data coordinates).
  * The GIF shows q(t) = q0 + A cos(2 pi t/T); frame 0 is +A and frame 49 is -A (99 frames per period).
    The measured mode displacement is d = (q(frame 0) - q(frame 49))/2, which cancels drawing offsets.
  * d is mass-weighted (sqrt(W) d) and compared with every eigenvector v_n of C~ by |cos| of the angle.
Result: the mode with |cos| close to 1 is the mode shown in the GIF.
Run from the repository root (the folder that contains scripts/):  python scripts/verify_animation_modes.py
"""
import os
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
SUP = os.path.join(os.path.dirname(HERE), "supplementary")

rod = [(-1, 0, np.pi/2 - np.pi/9, 3), (1, 0, np.pi/2 + np.pi/9, 3),
       (-1.1, 0.5, np.pi/2 - np.pi/20, 2), (1.1, 0.5, np.pi/2 + np.pi/20, 2)]
conn = [("r1","r2"),("r1","r2"),("r1","r3"),("r2","r4"),("r1","r3"),("r2","r4"),
        ("r1","f2"),("r3","f1"),("r2","f4"),("r4","f3")]
loc = [[1,1],[-1,-1],[-0.5,-1.0],[-0.5,-1.0],[1.4,0.85],[1.4,0.85],[-1.5,0],[1,0],[-1.5,0],[1,0]]
fp = [(-3,2), (-3,-2), (3,2), (3,-2)]
kk = np.array([1, 2, 100, 100, 1, 1, 1, 1, 1, 1.])          # ratios only matter for the mode shapes
m = np.array([3, 3, 2, 2.]); I = m*np.array([r[3] for r in rod])**2/12
W = np.array([w for i in range(4) for w in (m[i], m[i], I[i])])
q0 = np.array([v for r in rod for v in r[:3]])
def pt(q, end, d):
    if end[0] == 'f': return np.array(fp[int(end[1:])-1], float)
    i = int(end[1:])-1; x, y, t = q[3*i:3*i+3]
    return np.array([x + d*np.cos(t), y + d*np.sin(t)])
def lengths(q): return np.array([np.linalg.norm(pt(q,a,loc[j][0]) - pt(q,b,loc[j][1])) for j,(a,b) in enumerate(conn)])
h = 1e-7
B = np.array([(lengths(q0+h*e) - lengths(q0-h*e))/(2*h) for e in np.eye(12)]).T
Ct = (B/np.sqrt(W)).T @ np.diag(kk) @ (B/np.sqrt(W))          # = -C~ (positive semidefinite)
lam, V = np.linalg.eigh(Ct)                                    # ascending: modes 1..12

# default matplotlib figure 640x480, axes [0.125, 0.11, 0.775, 0.77], limits -3..3
def to_data(px, py):
    return -3 + (px - 80)/(0.775*640)*6, -3 + ((480 - py) - 0.11*480)/(0.77*480)*6
COL = np.array([(31,119,180), (255,127,14), (44,160,44), (214,39,40)], float)

def measure(img):
    a = np.asarray(img.convert('RGB'), float)
    out = []
    for c in COL:
        mask = np.linalg.norm(a - c, axis=2) < 40
        py, px = np.nonzero(mask)
        x, y = to_data(px, py)
        xc, yc = x.mean(), y.mean()
        ev, evec = np.linalg.eigh(np.cov(np.vstack([x - xc, y - yc])))
        ux, uy = evec[:, -1]
        th = np.arctan2(uy, ux) % np.pi                       # rods are near-vertical: angle in (0, pi)
        out += [xc, yc, th]
    return np.array(out)

lines = ["Animation file -> best-matching computed mode (|cos|), second best", ""]
for n in range(3, 13):
    g = Image.open(os.path.join(SUP, f"Animation_S{n:02d}_mode_{n}.gif"))
    g.seek(0); qa = measure(g)
    g.seek(49); qb = measure(g)
    d = (qa - qb)/2
    dm = np.sqrt(W)*d; dm /= np.linalg.norm(dm)
    c = np.abs(V.T @ dm); order = np.argsort(c)[::-1]
    ok = "OK" if order[0] + 1 == n else "MISMATCH"
    lines.append(f"Animation S{n:<2d} -> mode {order[0]+1:2d} (|cos| = {c[order[0]]:.3f});"
                 f" next: mode {order[1]+1:2d} ({c[order[1]]:.3f})   {ok}")
txt = "\n".join(lines)
print(txt)
open(os.path.join(HERE, "verify_animation_modes_output.txt"), "w", encoding="utf-8").write(txt + "\n")
