// Builds supplementary/Supplementary_Material.docx (Tables S1-S4, Fig. S1, animation captions)
// from supplementary/tables_data.json (written by supplementary_tables.py) and Fig_S1_zero_modes.png.
// Requires the npm package "docx".  Run from the repository root (the folder that contains scripts/):  node scripts/make_supplementary_docx.js
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell, ImageRun,
  WidthType, ShadingType, AlignmentType, BorderStyle, Footer, PageNumber, PageBreak, VerticalAlign,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const SUP = path.join(ROOT, "supplementary");
const D = JSON.parse(fs.readFileSync(path.join(SUP, "tables_data.json"), "utf8"));

const FONT = "Calibri";
const PAGE_W = 11906, MARGIN = 1134, CONTENT_W = PAGE_W - 2 * MARGIN;   // A4, 2 cm margins

// ---------- text helpers: **bold**, _italic_, ^superscript^, ~subscript~ ----------
const TR = (t, base) => new TextRun({ text: t.split(US).join("_"), ...base });
const US = "\u0001";                          // placeholder for a literal underscore, written "\\_" in the strings below
function runs(text, base = {}) {             // markup may be nested, e.g. _g~k~_
  text = text.split("\\_").join(US);
  const out = [];
  const re = /(\*\*[^*]+\*\*|_[^_]+_|\^[^^]+\^|~[^~]+~)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) out.push(TR(text.slice(last, m.index), base));
    const t = m[0], bold = t.startsWith("**"), inner = t.slice(bold ? 2 : 1, bold ? -2 : -1);
    const add = bold ? { bold: true } : t[0] === "_" ? { italics: true } : t[0] === "^" ? { superScript: true } : { subScript: true };
    out.push(...runs(inner, { ...base, ...add }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(TR(text.slice(last), base));
  return out;
}
const P = (t, opt = {}) => new Paragraph({ children: runs(t), spacing: { after: 120, line: 276 }, ...opt });
const NOTE = (t) => new Paragraph({ children: runs(t, { size: 18 }), spacing: { after: 60, line: 252 } });
const H1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun(t)], spacing: { before: 240, after: 120 } });
const CAP = (label, t) => new Paragraph({ children: [new TextRun({ text: label + " ", bold: true }), ...runs(t)], spacing: { before: 120, after: 100, line: 264 } });
const SUB = (t) => new Paragraph({ children: runs(t, { bold: true, size: 20 }), spacing: { before: 160, after: 60 } });

const border = { style: BorderStyle.SINGLE, size: 4, color: "8C96A5" };
const borders = { top: border, bottom: border, left: border, right: border };
function table(headers, rows, widths, opt = {}) {
  const sum = widths.reduce((a, b) => a + b, 0);
  const size = opt.size || 17;
  const cell = (txt, i, j) => new TableCell({
    borders, width: { size: widths[j], type: WidthType.DXA }, verticalAlign: VerticalAlign.CENTER,
    shading: i === 0 ? { fill: "DCE6F2", type: ShadingType.CLEAR, color: "auto" }
      : (opt.boldLast && i === rows.length ? { fill: "F2F2F2", type: ShadingType.CLEAR, color: "auto" } : undefined),
    margins: { top: 40, bottom: 40, left: 70, right: 70 },
    children: String(txt).split("\n").map(line => new Paragraph({
      alignment: j === 0 && opt.leftFirst !== false ? AlignmentType.LEFT : (opt.left && opt.left.includes(j) ? AlignmentType.LEFT : AlignmentType.CENTER),
      children: runs(line, { size, bold: i === 0 || (opt.boldLast && i === rows.length) ? true : undefined }),
    })),
  });
  return new Table({
    width: { size: sum, type: WidthType.DXA }, columnWidths: widths,
    rows: [headers, ...rows].map((r, i) => new TableRow({ tableHeader: i === 0, cantSplit: true, children: r.map((c, j) => cell(c, i, j)) })),
  });
}

// ---------- number formatting ----------
const f = (x, d) => Number(x).toFixed(d);
function sci(x, d = 2) {                       // 5.58×10^10^
  const e = Math.floor(Math.log10(Math.abs(x))); const m = x / Math.pow(10, e);
  return `${m.toFixed(d)}×10^${e}^`;
}
function sig(x, n = 2) {                       // significant figures, plain decimal
  if (x === 0) return "0";
  const e = Math.floor(Math.log10(Math.abs(x)));
  if (e >= n) { const q = Math.pow(10, e - n + 1); return String(Math.round(x / q) * q); }
  return Number(x).toFixed(Math.max(0, n - 1 - e));
}
const amp = (x) => (x < 0.0005 ? "0" : f(x, 2 + (x < 0.01 ? 1 : 0)));
const pos = (p) => `(${f(p[0], 3)}, ${f(p[1], 3)})`;
const s1 = D.S1, s2 = D.S2, s3 = D.S3, s4 = D.S4;

const children = [];
// ================= title =================
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 120 }, children: [new TextRun({ text: "Supplementary Material", bold: true, size: 32 })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 120 }, children: [new TextRun({ text: "Collective Mechanical Modes in a Minimal Model of Ion Channel Gate Dynamics", bold: true, size: 26 })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 }, children: [new TextRun({ text: "Gady Frenkel, Lev Mourokh, Yuval Ben-Abu", size: 22 })] }));
children.push(SUB("Contents"));
[
  "Table S1. Equilibrium geometry and parameters of the model",
  "Figure S1. The two zero-frequency modes (modes 1 and 2)",
  "Table S2. Thermal rms amplitudes at 310 K, modes 3–12",
  "Table S3. Damping of the normal modes",
  "Table S4. Sensitivity of frequencies and mode shapes to the spring constants",
  "Animations S1–S12. Motion of the rods in modes 1–12 (separate files; captions at the end of this document)",
].forEach(t => children.push(new Paragraph({ children: runs(t), spacing: { after: 40 }, indent: { left: 360 } })));
children.push(P("Coordinates: _x_ is perpendicular to the pore axis and _y_ is along it, as in Fig. 1B; rod orientations _θ_ are measured from the _x_ axis. Rods 1 and 2 are the long inner rods (pore walls), rods 3 and 4 the short outer rods. Springs are numbered as in Fig. 1B. Mode numbers are those of Table 1 (modes 1 and 2 are the zero-frequency modes).", { spacing: { before: 160, after: 120 } }));

// ================= Table S1 =================
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Table S1. Equilibrium geometry and parameters of the model"));
children.push(SUB("(a) Rods"));
children.push(table(
  ["Rod", "Role", "Centre (_x_, _y_) (nm)", "Orientation _θ_ (°)", "Length _L_ (nm)", "Mass _m_", "Moment of inertia _I_ = _mL_^2^/12 (kg m^2^)"],
  s1.rods.map(r => [String(r.rod), r.role, `(${f(r.x, 2)}, ${f(r.y, 2)})`, f(r.theta_deg, 1), f(r.L_nm, 0),
    `${sci(r.m_kg)} kg\n(${f(r.m_kDa, 2)} kDa)`, sci(r.I_kgm2)]),
  [600, 1900, 1350, 1150, 950, 1600, 2088], { left: [1] }));
children.push(SUB("(b) Fixed anchoring points"));
children.push(table(["Point", "(_x_, _y_) (nm)", "Attached spring"],
  [...s1.anchors].sort((a, b) => a.name.localeCompare(b.name)).map(a => {
    const sp = s1.springs.find(s => s.end2.endsWith(a.name));
    return [a.name, `(${f(a.x, 1)}, ${f(a.y, 1)})`, sp ? `spring ${sp.spring} (${sp.end1})` : "–"];
  }), [1500, 2400, 3000]));
children.push(NOTE("A1: upper and A2: lower fixed anchoring points, labelled as in Fig. 1 (left and right added here to distinguish the mirror partners). All four are fixed in space."));
children.push(SUB("(c) Springs"));
const endTxt = (e, s, p) => e.startsWith("rod") ? `${e}, _s_ = ${s > 0 ? "+" : ""}${s}\n${pos(p)}` : `${e.replace("anchoring point ", "")}\n${pos(p)}`;
children.push(table(
  ["Spring", "End 1: element, _s_ (nm)\nposition (nm)", "End 2: element, _s_ (nm)\nposition (nm)", "Natural length (nm)", "_k_ (N m^−1^)"],
  s1.springs.map(s => [String(s.spring), endTxt(s.end1, s.s1, s.p1), endTxt(s.end2, s.s2, s.p2), f(s.L0_nm, 3), String(s.k)]),
  [900, 2700, 2700, 1700, 1638]));
children.push(NOTE("_s_: position of the attachment point along the rod axis, measured from the rod centre (positive towards the upper end of the rod in Fig. 1B). The natural length of every spring equals its length in the equilibrium configuration, so the network is unstressed (no spring force at equilibrium); this is why the model has exactly two zero modes (Methods, Eq. 10)."));
children.push(NOTE("Units used in the code: the notebook works in model units (length 1 nm, mass 1.44×10^−24^ kg, spring constant 0.1 N m^−1^), so that the model values (3, 2 for the masses; 10, 20, 1000 for the springs) correspond to the SI values above. Frequencies scale as √(_k_/_m_); Table 1 corresponds to the spring constants listed here. Other parameters: temperature _T_ = 310 K; for the damping estimates (Table S3) each rod is a cylinder of diameter 1 nm, in water of viscosity 0.7 mPa s at 310 K."));

// ================= Figure S1 =================
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Figure S1. The two zero-frequency modes"));
const img = fs.readFileSync(path.join(SUP, "Fig_S1_zero_modes.png"));
const IW = 640;
children.push(new Paragraph({ alignment: AlignmentType.CENTER, children: [new ImageRun({ type: "png", data: img, transformation: { width: IW, height: Math.round(IW * 3240 / 6600) }, altText: { title: "Figure S1", description: "Zero modes 1 and 2 of the rod-spring model", name: "FigS1" } })] }));
children.push(CAP("Figure S1.", "The two zero-frequency modes (modes 1 and 2) of the unstressed spring network. Grey: equilibrium configuration, with the springs (dashed) and the fixed anchoring points (grey circles). Coloured solid and dotted lines: the rods displaced in the two opposite directions along each mode; the arrows show the displacement of the rod ends for the solid configuration. Blue and orange: long inner rods (pore walls); green and red: short outer rods. Linear mode shapes, with displacements exaggerated for visibility (largest rod-end displacement 3.3 Å in mode 1 and 3.6 Å in mode 2). The two modes have the same (zero) frequency, so any combination of them is also a zero mode; both are mirror-antisymmetric. The basis shown separates the motion of the pore walls from that of the outer rods. In zero mode 1, all of the (mass-weighted) motion is in the pore walls, which slide along their axes in opposite directions, with a smaller sideways shift and tilt in the same sense; this changes the asymmetric shift Δ_y_ but not the pore width. In zero mode 2, 99% of the motion is in the outer rods, which slide axially in opposite directions and tilt; the pore width is unchanged and Δ_y_ changes by less than 5% of its change in mode 1. Neither motion has a restoring force in the harmonic model. See also Animations S1 and S2."));

// ================= Table S2 =================
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Table S2. Thermal rms amplitudes at 310 K, modes 3–12"));
const c2 = ["dx_top", "dx_mid", "dx_bot", "dy_top", "dy_mid", "dy_bot", "sway", "piston"];
children.push(table(
  ["Mode", "Δ_x_ top\n(Å)", "Δ_x_ middle\n(Å)", "Δ_x_ bottom\n(Å)", "Δ_y_ top\n(Å)", "Δ_y_ middle\n(Å)", "Δ_y_ bottom\n(Å)", "Sway\n(Å)", "Piston\n(Å)", "Wall tilt\n(°)", "Outer-rod tilt\n(°)"],
  [...s2.rows.map(r => [String(r.mode), ...c2.map(c => amp(r[c])), f(r.th_in, 2), f(r.th_out, 2)]),
   ["3–12", ...c2.map(c => f(s2.total[c], 2)), "", ""]],
  [700, 840, 900, 900, 840, 900, 900, 840, 840, 900, 1078], { boldLast: true }));
children.push(NOTE("Root-mean-square (rms) thermal amplitudes from equipartition, rms(_g_) = |Σ~k~ _g~k~ v~kn~_/√_W~k~_| (_k~B~T_)^1/2^/_ω~n~_ (Eq. 12), for _T_ = 310 K. They do not depend on damping. Δ_x_ (pore width) and Δ_y_ (asymmetric shift) are defined by Eqs. (8)–(9) at the upper ends (top), centres (middle) and lower ends (bottom) of the two long rods. Sway and piston are the common sideways and common axial displacements of the centres of the two long rods, (_x~1~_ + _x~2~_)/2 and (_y~1~_ + _y~2~_)/2. Wall tilt and outer-rod tilt are the rotations of each long and each short rod (equal for the left and right rods by symmetry)."));
children.push(NOTE("Mirror-symmetric modes (3, 4, 7, 8, 10, 12) change the width and the piston coordinate, with Δ_y_ and sway exactly zero; mirror-antisymmetric modes (5, 6, 9, 11) change Δ_y_ and the sway, with width and piston zero to first order. Last row: modes 3–12 combined (root sum of squares; the modes are statistically independent). The zero modes 1 and 2 are excluded, because their amplitude is not defined in the harmonic model; they change only Δ_y_ (Fig. S1)."));

// ================= Table S3 =================
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Table S3. Damping of the normal modes"));
const ex = s3.exact["0.7"];
children.push(table(
  ["Mode", "_ω~n~_\n(rad s^−1^)", "_Q_, single sphere", "_Q_, rods,\nwater", "_Q_, rods,\n10× water", "_Q_, rods,\n100× water", "_η_ for critical damping (mPa s)", "Relaxation time in water, modal (ps)", "Relaxation time in water, full solution (ps) [mode weight]"],
  s3.rows.map(r => [String(r.mode), sci(r.omega), sig(r.Q_sphere), sig(r.Q_water), sig(r.Q_10), sig(r.Q_100), sig(r.eta_crit_mPas),
    r.tau_modal_ps == null ? "oscillatory" : sig(r.tau_modal_ps, 2),
    r.tau_exact_ps == null ? "oscillatory" : `${sig(r.tau_exact_ps, 2)} [${f(r.exact_weight, 2)}]`]),
  [650, 1250, 950, 900, 950, 950, 1150, 1250, 1588]));
children.push(NOTE("Single sphere: _Q_ = _mω~n~_/_γ_ with _m_ = 4.32×10^−24^ kg and a Stokes friction coefficient _γ_ = 2×10^−11^ kg s^−1^, the order-of-magnitude estimate for one protein element. Rods: each rod is treated as a rigid cylinder of its length and a diameter of 1 nm, with translational friction parallel and perpendicular to its axis and rotational friction from the cylinder formulas of Tirado, Martínez and García de la Torre; the resulting 12×12 friction matrix Γ is projected on each mode, _Q~n~_ = _ω~n~_/(_v~n~_^T^ Γ̃ _v~n~_) (Eq. 14). Water: _η_ = 0.7 mPa s (310 K). A mode is underdamped (oscillates) only for _Q~n~_ > 1/2; the critical viscosity is the viscosity at which _Q~n~_ = 1/2."));
children.push(NOTE(`Relaxation times. Modal: the slow root of the single-mode equation _ζ̈~n~_ + _γ~n~ζ̇~n~_ + _ω~n~_^2^_ζ~n~_ = 0 with _γ~n~_ = _ω~n~_/_Q~n~_, i.e. _τ~n~_ = 2/[_γ~n~_ − (_γ~n~_^2^ − 4_ω~n~_^2^)^1/2^] ≈ 1/(_Q~n~ω~n~_) for _Q~n~_ ≪ 1. Full solution: the friction matrix is not diagonal in the basis of the normal modes, so it couples modes of the same mirror symmetry. The damped equations of motion (Eq. 13 without noise) were therefore also solved without approximation; their solutions are e^_st_^, with _s_ the eigenvalues of the 24×24 matrix [[0, 1], [−_W_^−1^_H_, −_W_^−1^Γ]] (_H_ = −_C_), and each real _s_ < 0 is a relaxation with time 1/|_s_|. Each such solution is in general a mixture of normal modes; the mode weight is the squared projection of its (mass-weighted, normalised) displacement on _v~n~_ (the weights of one solution sum to 1). For each mode the table gives the solution with the largest weight of that mode. Modes 3, 4, 5, 7 and 9 dominate their own solution (weight ≥ 0.86); mode 6 is mixed with mode 5 (0.55/0.37), and the two fastest of these solutions (2.3 and 2.0 ps) are mixtures of modes 7, 8 and 10, so mode 8 has no solution of its own. The relaxation times nevertheless span the same range, ≈2 ps to ≈1.1 ns. In water the exact solution has ${ex.n_osc_pairs} oscillatory pairs (modes 11 and 12; damped frequency ≈ ${sci(ex.osc_w[0], 2)} rad s^−1^); at 10× and 100× the viscosity of water no mode oscillates. These estimates are semi-quantitative: friction at picosecond time scales is frequency dependent, and the pore elements are only partly exposed to water (partial burial in protein and lipid increases the friction).`));

// ================= Table S4 =================
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Table S4. Sensitivity of frequencies and mode shapes to the spring constants"));
children.push(SUB("(a) Scaling of groups of springs: frequency ratio _ω~n~_′/_ω~n~_"));
const G = s4.groups;
const gHead = G.map(g => `${g.group.replace(/ \(.*\)/, "").replace("3-8", "3–8").replace("Springs", "Springs\n").replace("Spring ", "Spring\n")}\n×${g.factor}`);
const rowsA = [];
for (let j = 0; j < 10; j++) rowsA.push([String(j + 3), ...G.map(g => f(g.ratios[j], 3))]);
rowsA.push(["Min. overlap", ...G.map(g => `${f(g.min_overlap, 3)}\n(mode ${g.min_overlap_mode})`)]);
rowsA.push(["Zero modes", ...G.map(g => String(g.zero_modes))]);
children.push(table(["Mode", ...gHead], rowsA, [1238, ...Array(8).fill(1050)], { size: 16 }));
children.push(NOTE("Springs 1 and 3–8 (1 N m^−1^), spring 2 (2 N m^−1^) and springs 9 and 10 (100 N m^−1^) were scaled by the factor shown. Overlap: |_v~n~_ · _v~n~_′|, the scalar product of the normalised mass-weighted eigenvectors of the original and the modified model (1 = identical shape); modes were matched by maximal overlap. Scaling all spring constants by a common factor _s_ multiplies every frequency by √_s_ and leaves all mode shapes unchanged (exact). Rigid limit: making springs 9 and 10 ten or one hundred times stiffer changes the frequencies of modes 3–10 by at most " + f(100 * Math.max(...G.filter(g => g.factor >= 10).map(g => g.max_change_3_10)), 2) + "%, with unchanged shapes; only modes 11 and 12 change (by ≈√factor). Springs 9 and 10 therefore act as a nearly rigid connection."));
children.push(SUB("(b) Random variations: every spring constant multiplied by a random factor between 0.5 and 1.5 (2000 samples)"));
const rowsB = s4.sym.rows.map((r, j) => {
  const q = s4.indep.rows[j];
  return [String(r.mode), `${f(r.p5, 2)}–${f(r.p95, 2)} (${f(r.median, 2)})`, `${f(r.pct09, 0)}%`, `${f(q.p5, 2)}–${f(q.p95, 2)} (${f(q.median, 2)})`, `${f(q.pct09, 0)}%`];
});
children.push(table(["Mode", "Symmetry-preserving:\n_ω~n~_′/_ω~n~_ in 90% of samples (median)", "Symmetry-preserving:\nshape preserved (overlap > 0.9)", "Independent:\n_ω~n~_′/_ω~n~_ in 90% of samples (median)", "Independent:\nshape preserved (overlap > 0.9)"],
  rowsB, [900, 2300, 1900, 2300, 1900]));
const sp = (o) => Math.min(o.subspace_pct["8-10"], o.subspace_pct["11-12"]);
children.push(NOTE("Each spring constant was multiplied by a random factor drawn from the uniform distribution on [0.5, 1.5] (any value between 50% and 150% of the original, all equally likely), and the modes were recomputed; this was repeated for 2000 random parameter sets. Frequency column: the range between the 5th and 95th percentiles of _ω~n~_′/_ω~n~_, which contains 90% of the samples, is the measure of sensitivity (e.g. mode 3: in 90% of the parameter sets its frequency stays within −18% and +13% of the Table 1 value). The median (in parentheses) is close to 1 because the factors are centred on 1; it shows only that there is no systematic shift. For comparison, changing all springs together by ×0.5 or ×1.5 changes every frequency by ×0.71 or ×1.22. Shape column: percentage of samples in which the mode keeps an overlap > 0.9 with its original shape (Table S4a)."));
children.push(NOTE(`Symmetry-preserving: the two springs of each mirror pair (3/4, 5/6, 7/8, 9/10) were scaled by the same random factor, and springs 1 and 2 by their own factors. Independent: each of the ten springs was scaled by its own random factor, which breaks the mirror symmetry. When the symmetry is broken, the nearly degenerate pairs 9/10 and 11/12 mix (11 and 12 become left- or right-localised motions of the outer rods), while the subspaces they span are preserved: the subspace of modes 8–10 and that of modes 11–12 have all principal overlaps > 0.9 in ≥ ${f(Math.floor(Math.min(sp(s4.sym), sp(s4.indep))), 0)}% of the samples. The number of zero modes was two in every sample.`));

// ================= Animations =================
children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H1("Animations S1–S12"));
children.push(P("Separate files Animation\\_S01\\_zero\\_mode\\_1 … Animation\\_S12\\_mode\\_12 (animated GIF; identical MP4 versions are also provided). Animation S_n_ shows mode _n_; the number in the box in each animation is the mode number."));
children.push(CAP("Animations S1 and S2.", "Zero modes 1 and 2 (see Figure S1). The rods are moved one way along the mode, from one end of the displayed range to the other, and the loop then restarts. A zero mode has no restoring force and no frequency, so it is not shown as an oscillation, and the speed of the animation has no physical meaning. Grey: equilibrium configuration. Same amplitude as Figure S1."));
children.push(CAP("Animations S3–S12.", "Normal modes 3–12 (Animation S_n_ shows mode _n_), over one period of oscillation. Frequencies are given in Table 1. Blue and orange: long inner rods (pore walls); green and red: short outer rods. All modes are shown with the same normalised (mass-weighted) amplitude, exaggerated for visibility. The physical rms amplitudes at 310 K are given in Table 1 and Table S2: because the thermal amplitude scales as 1/_ω_, the soft modes 3 and 4 fluctuate much more than they appear to here. Axes in nm."));
children.push(P("Code and data: the Jupyter notebook and the Python scripts that produce all tables, figures and animations are available at https://github.com/gadyf/Collective-Mechanical-Modes-in-a-Minimal-Model-of-Ion-Channel-Gate-Dynamics and archived at Zenodo, doi:10.5281/zenodo.23687315.", { spacing: { before: 240 } }));

const doc = new Document({
  creator: "G. Frenkel", title: "Supplementary Material",
  styles: {
    default: { document: { run: { font: FONT, size: 21 } } },
    paragraphStyles: [{ id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { size: 26, bold: true, font: FONT, color: "1F3864" }, paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 0 } }],
  },
  sections: [{
    properties: { page: { size: { width: PAGE_W, height: 16838 }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: "S", size: 18 }), new TextRun({ children: [PageNumber.CURRENT], size: 18 })] })] }) },
    children,
  }],
});
Packer.toBuffer(doc).then(buf => {
  const out = process.env.SUPP_OUT || path.join(SUP, "Supplementary_Material.docx");   // SUPP_OUT: alternative output path
  fs.writeFileSync(out, buf); console.log("written", out);
});
