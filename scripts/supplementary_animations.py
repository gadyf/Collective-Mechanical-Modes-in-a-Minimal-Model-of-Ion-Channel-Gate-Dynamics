"""Supplementary animations and Fig. S1 for the revised K2P paper.

Produces in supplementary/ (next to scripts/):
  * Animation_S03_mode_3.gif ... Animation_S12_mode_12.gif : copies of the existing notebook GIFs
    (rods_animation_mode_N_one_period.gif, N = 3..12), renamed by mode number (Animation S3-S12).
  * the same as .mp4 (frame-by-frame conversion, identical content).
  * Animation_S01_zero_mode_1.gif/.mp4 and Animation_S02_zero_mode_2.gif/.mp4 : the two zero modes
    (basis chosen to separate wall and outer-rod motion; see below), drawn like the notebook animations
    (cell 44), linear (harmonic) mode shape. Each is a ONE-WAY sweep from -AMPZ to +AMPZ (3 s, 0.5 s hold
    at each end), then the loop restarts. A zero mode has no restoring force and no frequency, so it is
    not shown as an oscillation and the speed has no physical meaning.
  * Fig_S1_zero_modes.png (600 dpi) and .svg : equilibrium and the +/-AMPZ configurations.

Requires: numpy, matplotlib, pillow; for MP4 also imageio + imageio-ffmpeg.
Run from the repository root (the folder that contains scripts/):  python scripts/supplementary_animations.py
"""
import os, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # repository root (contains scripts/, notebook/, supplementary/)
# folders searched for the notebook GIFs rods_animation_mode_N_one_period.gif (written by the notebook, cell 45)
GIF_DIRS = [os.path.join(ROOT, "notebook"), os.path.dirname(ROOT)]
OUT = os.path.join(ROOT, "supplementary")
os.makedirs(OUT, exist_ok=True)

# ---------------- model (identical to notebook cell 12 and scripts/check_model.py; toy units) -------------
rod = [(-1, 0, np.pi/2 - np.pi/9, 3), (1, 0, np.pi/2 + np.pi/9, 3),
       (-1.1, 0.5, np.pi/2 - np.pi/20, 2), (1.1, 0.5, np.pi/2 + np.pi/20, 2)]
conn = [("r1","r2"),("r1","r2"),("r1","r3"),("r2","r4"),("r1","r3"),("r2","r4"),
        ("r1","f2"),("r3","f1"),("r2","f4"),("r4","f3")]
loc = [[1,1],[-1,-1],[-0.5,-1.0],[-0.5,-1.0],[1.4,0.85],[1.4,0.85],[-1.5,0],[1,0],[-1.5,0],[1,0]]
fp = [(-3,2), (-3,-2), (3,2), (3,-2)]
m = np.array([3, 3, 2, 2.]); I = m*np.array([r[3] for r in rod])**2/12
W = np.array([w for i in range(4) for w in (m[i], m[i], I[i])])
q0 = np.array([v for r in rod for v in r[:3]])

def pt(q, end, d):
    if end[0] == 'f': return np.array(fp[int(end[1:])-1], float)
    i = int(end[1:])-1; x, y, t = q[3*i:3*i+3]
    return np.array([x + d*np.cos(t), y + d*np.sin(t)])
def lengths(q): return np.array([np.linalg.norm(pt(q,a,loc[j][0]) - pt(q,b,loc[j][1])) for j,(a,b) in enumerate(conn)])
h = 1e-7
B = np.array([(lengths(q0+h*e) - lengths(q0-h*e))/(2*h) for e in np.eye(12)]).T   # 10x12

# zero modes in mass-weighted coordinates: null space of B M^{-1/2}
Bm = B / np.sqrt(W)[None, :]
_, sv, Vt = np.linalg.svd(Bm)
N = Vt[-2:].T                                            # 12x2, orthonormal (mass-weighted)
# mirror operator (x -> -x): rod1<->rod2, rod3<->rod4; dx -> -dx, dy -> dy, dtheta -> -dtheta
P = np.zeros((12,12))
for a, b in [(0,1),(1,0),(2,3),(3,2)]:
    P[3*b,3*a] = -1; P[3*b+1,3*a+1] = 1; P[3*b+2,3*a+2] = -1
# Both zero modes are mirror-antisymmetric, so the mirror symmetry does not fix a unique basis of the
# degenerate pair. We choose the basis that separates the motion of the long inner rods (pore walls)
# from that of the short outer rods as far as possible: eigenvectors of the inner-rod projector
# restricted to the null space. Zero mode 1 = maximally wall-dominated, zero mode 2 = maximally outer-rod-dominated.
Pin = np.diag([1.0]*6 + [0.0]*6)
share, U = np.linalg.eigh(N.T @ Pin @ N)
Z = (N @ U)[:, ::-1]                                     # descending inner-rod share
for k in range(2):
    z = Z[:, k]
    if z[1] > 0: Z[:, k] = -z                            # sign convention: left wall moves down first
log = []
log.append(f"singular values of mass-weighted B: {np.array2string(sv, precision=3)}")
for k in range(2):
    z = Z[:, k]; dq = z/np.sqrt(W)
    log.append(f"zero mode {k+1}: parity {z@P@z:+.3f}; |B dq| = {np.linalg.norm(B@dq):.1e}; "
               f"share inner rods {np.linalg.norm(z[:6])**2:.2f}, outer {np.linalg.norm(z[6:])**2:.2f}")

# ---------------- drawing helpers (as notebook cell 44) ----------------
def rod_ends(i, dq):
    x, y, t, L = rod[i]
    xc, yc, tt = x + dq[3*i], y + dq[3*i+1], t + dq[3*i+2]
    return ((xc - L/2*np.cos(tt), xc + L/2*np.cos(tt)), (yc - L/2*np.sin(tt), yc + L/2*np.sin(tt)))

def label(ax, txt):
    ax.text(0.065, 0.9, txt, transform=ax.transAxes, fontsize=24, va='top', ha='left',
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='black', boxstyle='round'))

AMPZ = 0.7     # zero-mode amplitude (mass-weighted), 3.5x the 0.2 of the notebook GIFs, for visibility

def zdq(k, a):
    """Linear (harmonic) zero-mode displacement at mass-weighted amplitude a."""
    return Z[:, k]/np.sqrt(W)*a

def rod_end_disp(dq):
    return max(np.linalg.norm(np.array(rod_ends(i, dq)) - np.array(rod_ends(i, np.zeros(12))), axis=0).max()
               for i in range(4))

for k in range(2):
    dq = zdq(k, AMPZ)
    log.append(f"zero mode {k+1} drawn at amplitude {AMPZ}: max rod-end displacement {rod_end_disp(dq)*10:.2f} A; "
               f"mid dy(right-left) {(dq[4]-dq[1])*10:+.2f} A; mid width change {(dq[3]-dq[0])*10:+.3f} A")
NMOVE, NHOLD = 150, 25                                   # 3 s sweep + 0.5 s hold at each end (50 fps)

def make_zero_mode_animation(k):
    """One-way sweep along the zero mode, from -AMPZ to +AMPZ; then the loop restarts.
    A zero mode has no restoring force and no frequency, so it is not drawn as an oscillation."""
    amps = np.concatenate([np.full(NHOLD, -AMPZ), np.linspace(-AMPZ, AMPZ, NMOVE), np.full(NHOLD, AMPZ)])
    fig, ax = plt.subplots()
    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.grid(True)
    label(ax, f"({k+1})")
    for i in range(4):
        xs, ys = rod_ends(i, np.zeros(12))
        ax.plot(xs, ys, '-', color='0.6', lw=5, alpha=0.35, solid_capstyle='round', zorder=1)
    ax.text(0.5, 0.03, "zero mode: one-way displacement along the mode (no restoring force, no time scale)",
            transform=ax.transAxes, ha='center', fontsize=8,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    lines = [ax.plot([], [], 'o-', lw=2)[0] for _ in range(4)]
    def update(f):
        for i in range(4):
            xs, ys = rod_ends(i, zdq(k, amps[f])); lines[i].set_data(xs, ys)
        return lines
    ani = animation.FuncAnimation(fig, update, frames=len(amps), blit=True, interval=20)
    base = os.path.join(OUT, f"Animation_S{k+1:02d}_zero_mode_{k+1}")
    ani.save(base + ".gif", writer="pillow", fps=50)
    plt.close(fig)
    return base + ".gif"

def make_fig_s1():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
    colors = ['C0', 'C1', 'C2', 'C3']
    for k, ax in enumerate(axes):
        dplus, dminus = zdq(k, AMPZ), zdq(k, -AMPZ)
        for j, (a, b) in enumerate(conn):
            p1, p2 = pt(q0, a, loc[j][0]), pt(q0, b, loc[j][1])
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], '--', color='0.75', lw=1, zorder=1)
        for x, y in fp:
            ax.plot(x, y, 'o', color='0.6', ms=10, zorder=1)
        for i in range(4):
            xs, ys = rod_ends(i, np.zeros(12))
            ax.plot(xs, ys, '-', color='0.6', lw=6, alpha=0.45, solid_capstyle='round', zorder=2)
            for dq, ls in [(dplus, '-'), (dminus, ':')]:
                xs, ys = rod_ends(i, dq)
                ax.plot(xs, ys, ls, color=colors[i], lw=2.5, zorder=3)
            x0s, y0s = rod_ends(i, np.zeros(12)); x1s, y1s = rod_ends(i, dplus)
            for e in (0, 1):
                ax.annotate('', xy=(x1s[e], y1s[e]), xytext=(x0s[e], y0s[e]),
                            arrowprops=dict(arrowstyle='->', color=colors[i], lw=1.2), zorder=4)
        ax.set_aspect('equal'); ax.set_xlim(-3.4, 3.4); ax.set_ylim(-2.6, 2.6); ax.grid(True, alpha=0.4)
        ax.set_xlabel('x (nm)'); ax.set_ylabel('y (nm)')
        sh = np.linalg.norm(Z[:6, k])**2
        ax.set_title(f"Zero mode {k+1}: {sh*100:.0f}% of the motion in the long inner rods (pore walls)", fontsize=10.5)
        ax.text(0.5, 0.07, f"({k+1})", transform=ax.transAxes, fontsize=20, va='center', ha='center',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round'))
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color='0.6', lw=6, alpha=0.45, label='equilibrium'),
               Line2D([], [], color='k', lw=2.5, ls='-', label='displaced (+)'),
               Line2D([], [], color='k', lw=2.5, ls=':', label='displaced (−)'),
               Line2D([], [], color='0.75', lw=1, ls='--', label='springs'),
               Line2D([], [], color='0.6', marker='o', ls='', ms=8, label='fixed points')]
    fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(OUT, "Fig_S1_zero_modes.png"), dpi=600)
    fig.savefig(os.path.join(OUT, "Fig_S1_zero_modes.svg"))
    plt.close(fig)

def to_mp4(gif, mp4):
    try:
        import imageio.v2 as iio
    except ImportError:
        log.append("imageio not available: MP4 skipped"); return False
    from PIL import Image
    g = Image.open(gif); frames = []
    FPS = 50                                             # all animations are drawn at 20 ms per frame
    for f in range(g.n_frames):
        g.seek(f); fr = np.array(g.convert('RGB'))
        hh, ww = fr.shape[0] - fr.shape[0] % 2, fr.shape[1] - fr.shape[1] % 2   # even size for H.264
        rep = max(1, round((g.info.get('duration', 20) or 20)*FPS/1000))     # merged (held) frames
        frames += [fr[:hh, :ww]]*rep
    iio.mimwrite(mp4, frames, fps=FPS, codec='libx264', quality=8, macro_block_size=1)
    return True

if __name__ == "__main__":
    # 1. existing notebook GIFs, modes 3-12 -> Animation S3-S12
    for n in range(3, 13):
        dst = os.path.join(OUT, f"Animation_S{n:02d}_mode_{n}.gif")
        src = next((os.path.join(d, f"rods_animation_mode_{n}_one_period.gif") for d in GIF_DIRS
                    if os.path.exists(os.path.join(d, f"rods_animation_mode_{n}_one_period.gif"))), None)
        if src: shutil.copyfile(src, dst)
        elif not os.path.exists(dst):
            log.append(f"mode {n}: notebook GIF not found (run the notebook first); skipped"); continue
        to_mp4(dst, dst[:-4] + ".mp4")
    # 2. zero-mode animations S1, S2
    for k in range(2):
        gif = make_zero_mode_animation(k)
        to_mp4(gif, gif[:-4] + ".mp4")
    # 3. Fig. S1
    make_fig_s1()
    log.append("files: " + ", ".join(sorted(os.listdir(OUT))))
    open(os.path.join(OUT, "supplementary_animations_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
    print("\n".join(log))
