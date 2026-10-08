# Collective mechanical modes of a minimal ion-channel gate model

Code and data for the paper

> G. Frenkel, L. Mourokh, Y. Ben-Abu, **"Collective Mechanical Modes in a Minimal Model of Ion Channel Gate Dynamics"**, *Computational and Structural Biotechnology Journal* (under review, 2026). DOI: [to be added]

Archived version of this repository: [Zenodo DOI to be added]

## The model in brief

A two-dimensional model of the gate of a K2P potassium channel. It has four rigid rods and ten linear springs:
- rods 1 and 2 are the long inner rods (the pore walls);
- rods 3 and 4 are the short outer rods (peripheral elements);
- the springs couple the rods to each other and to four fixed anchoring points (A1: upper pair, A2: lower pair).

The Euler–Lagrange equations of the 12 coordinates (x, y and θ of each rod) are derived symbolically and linearised about the equilibrium configuration. They are then solved for the normal modes. All springs are at their natural length at equilibrium (an unstressed network). Ten springs therefore restrain at most ten of the twelve coordinate combinations, so the model has exactly two zero-frequency modes (modes 1 and 2) and ten oscillatory modes (modes 3–12).

The repository also contains the analysis added in the revision:
- thermal (equipartition) amplitudes of the pore width and asymmetric shift at 310 K;
- damping of every mode, using cylinder friction coefficients, and the exact damped solution;
- sensitivity of frequencies and mode shapes to the spring constants;
- an independent Langevin-dynamics check;
- all supplementary tables, figures and animations.

## Repository contents

```
notebook/
  k2p_rod_spring_model.ipynb     symbolic derivation (SymPy), normal modes, Table 1, Fig. 1B, Fig. 2/3 panels and animation of one mode per run
scripts/
  check_model.py                 independent re-implementation: rank of the compatibility matrix, eigenvalues, zero modes
  zero_mode_finite.py            are the zero modes finite mechanisms? (run check_model.py first)
  analysis_tables.py             frequencies in SI, thermal amplitudes, damping, sensitivity  -> analysis_output.txt
  mode_character.py              decomposition of each mode (sliding, tilting, pivot point)  -> mode_character_output.txt
  langevin_check.py              direct simulation of the Langevin equation (Eq. 13)          -> langevin_check_output.txt
  supplementary_tables.py        all numbers of Supplementary Tables S1-S4  -> supplementary/tables_data.json
  make_supplementary_docx.js     builds supplementary/Supplementary_Material.docx from tables_data.json
  supplementary_animations.py    Fig. S1 and Animations S1-S12 (GIF and MP4)
  verify_animation_modes.py      checks that each Animation S3-S12 shows the correct mode (from the GIF pixels)
  *_output.txt                   reference output of each script
supplementary/
  Supplementary_Material.docx/.pdf   Tables S1-S4, Fig. S1, captions
  Fig_S1_zero_modes.png/.svg         Fig. S1
  Animation_S01 ... Animation_S12    animations of modes 1-12 (GIF and identical MP4)
  tables_data.json                   numbers behind Tables S1-S4
  README.md                          file list and captions
```

## Installation

Python 3.10 or later:

```
pip install -r requirements.txt
```

For rebuilding the Word document only (optional), Node.js 18 or later:

```
npm install
```

## Reproducing the results

Run everything from the repository root.

| Paper item | Command | Output |
|---|---|---|
| Equations of motion, normal modes, Table 1 frequencies (SI), Fig. 1B; for the mode selected by `index_mode`: its Fig. 2 and Fig. 3 panels and its animation GIF | open `notebook/k2p_rod_spring_model.ipynb` in Jupyter and run all cells (set `index_mode` and rerun from that cell for other modes) | `figure_1b.svg`, `Time Evolution of coordinates - new mode N.svg` (Fig. 2 panel), `Figure 3 new_mode N.svg`, `rods_animation_mode_N_one_period.gif` in `notebook/` |
| Zero modes and rank (Methods, Eq. 10) | `python scripts/check_model.py` then `python scripts/zero_mode_finite.py` | printed |
| Table 1 (SI frequencies, amplitudes, Q), damping and sensitivity results | `python scripts/analysis_tables.py` | `scripts/analysis_output.txt` |
| Mode descriptions (Table 1, Results) | `python scripts/mode_character.py` | `scripts/mode_character_output.txt` |
| Langevin check of Eq. 12 (amplitudes independent of friction) | `python scripts/langevin_check.py` (several minutes) | `scripts/langevin_check_output.txt` |
| Supplementary Tables S1–S4 | `python scripts/supplementary_tables.py`, then `node scripts/make_supplementary_docx.js` | `supplementary/tables_data.json`, `supplementary/Supplementary_Material.docx` |
| Fig. S1, Animations S1–S12 | `python scripts/supplementary_animations.py` | `supplementary/` |
| Check of the animation–mode correspondence | `python scripts/verify_animation_modes.py` | `scripts/verify_animation_modes_output.txt` |

The panels of Figs. 2 and 3 in the paper were assembled from the per-mode figures written by the notebook. The random sampling in the sensitivity analysis uses a fixed seed (`numpy.random.default_rng(1)`), so the numbers are reproduced exactly. Tested with Python 3.10.6, NumPy 1.26.4, SciPy 1.12.0, Matplotlib 3.8.3, SymPy 1.14.0, Pillow 10.2.0, imageio 2.38.0 and imageio-ffmpeg 0.6.0.

## Units

The notebook works in model units chosen so that all numbers are of order 1:

| Quantity | Model unit | Model values | SI values used in the paper |
|---|---|---|---|
| length | 1 nm | rod lengths 3, 2 | 3 nm, 2 nm |
| mass | 1.44×10⁻²⁴ kg | 3, 2 | 4.32×10⁻²⁴ kg (2.6 kDa), 2.88×10⁻²⁴ kg (1.73 kDa) |
| spring constant | 0.1 N/m | 10, 20, 1000 | 1, 2, 100 N/m |

Frequencies convert as ω_SI = ω_model · √(0.1 N m⁻¹ / 1.44×10⁻²⁴ kg). The scripts work directly in SI units.

## Conventions

- **Mode numbers** are those of the paper: ascending frequency, including the two zero modes (1 and 2).
- **Spring numbers** are those of Fig. 1B. Note that the code lists the springs in a different order; Table S1 gives both positions of every spring.
- **Rods:** 1 inner left, 2 inner right, 3 outer left, 4 outer right. x is perpendicular to the pore axis and y is along it.

## Citation

If you use this code, please cite the paper above and the archived version of this repository (`CITATION.cff`).

## License

- **Code** (`notebook/`, `scripts/` and the configuration files): MIT License, see [`LICENSE`](LICENSE).
- **Figures, animations, supplementary documents and data** (`supplementary/`): Creative Commons Attribution 4.0 International (CC BY 4.0), see [`supplementary/LICENSE.md`](supplementary/LICENSE.md).
