# NV Center Physics — Master Reference

**CryoLab Technical Document** | Created 2026-03-02 | v1.0
**Scope:** Ground-state spin physics, ODMR, magnetometry, fitting functions
**Audience:** Jacob — for fitting, analysis, and paper writing

> **Citation note:** All equations and values in this document are sourced from the
> papers in the `NV_Physics_References/` folder. Inline citations use [Author Year]
> format; full references with PDF filenames are in §10.

---

## Table of Contents
1. [Crystal Structure and the NV Defect](#1-crystal-structure-and-the-nv-defect)
2. [Electronic Structure](#2-electronic-structure)
3. [Ground-State Spin Hamiltonian](#3-ground-state-spin-hamiltonian)
4. [Zero-Field Splitting: D and E](#4-zero-field-splitting-d-and-e)
5. [Zeeman Interaction and Transition Frequencies](#5-zeeman-interaction-and-transition-frequencies)
6. [The Four NV Axes — Peak Pairing and Geometry](#6-the-four-nv-axes--peak-pairing-and-geometry)
7. [ODMR Transition Frequencies — The Fitting Function](#7-odmr-transition-frequencies--the-fitting-function)
8. [Linewidth and FWHM](#8-linewidth-and-fwhm)
9. [CW-ODMR Sensitivity](#9-cw-odmr-sensitivity)
10. [References](#10-references)
11. [Appendix A: Numerical Constants](#appendix-a-numerical-constants)
12. [Appendix B: Our Measured Values](#appendix-b-our-measured-values)

---

## 1. Crystal Structure and the NV Defect

Diamond is face-centered cubic (FCC) with a two-atom basis. Each carbon is
sp3-bonded to four nearest neighbors along the tetrahedral ⟨111⟩ directions.
The lattice constant is a = 3.567 A [Doherty 2013, §2.1].

The nitrogen-vacancy (NV) center is a **point defect**: one carbon atom replaced
by nitrogen (N) adjacent to a vacant lattice site (V). The N-V axis is along one
of the four crystallographically equivalent ⟨111⟩ directions. The defect has
**C_3v symmetry** — 3-fold rotational symmetry about the N-V axis, plus three
mirror planes containing that axis [Maze 2011, §II].

The relevant charge state for magnetometry is **NV^-** (negatively charged),
which has 6 electrons: 2 from the vacancy dangling bonds, 2 from the nitrogen,
and 1 captured from the lattice. This gives a spin-1 (S=1) ground state
[Doherty 2013, §3].

### Why NV^- and not NV^0?

NV^0 has spin-1/2 and no optically detected magnetic resonance (ODMR) contrast.
NV^- has spin-1 with a spin-dependent fluorescence mechanism that enables ODMR.
The charge state can switch under illumination; typical NV experiments ensure
stable NV^- population via appropriate laser power and wavelength (532 nm)
[Mittiga 2018].


## 2. Electronic Structure

The NV^- ground state is a spin triplet: **^3A_2** (orbital singlet, spin S=1).
The excited state is also a spin triplet: **^3E** (orbital doublet, spin S=1).

Key energy levels [Doherty 2013, §3; Maze 2011]:

```
                     ^3E (excited state)
                    /    \
                 m_s=0    m_s=±1
                    |         |
   532 nm excite    |    ISC to singlet
                    |    (spin-selective)
                    |         |
   637-800 nm       |    ^1A_1 → ^1E
   fluorescence     |    (non-radiative or 1042 nm IR)
                    |         |
                    \    /
                     ^3A_2 (ground state)
                    /    \
                 m_s=0    m_s=±1    ← separated by D ≈ 2.87 GHz
```

**The ODMR mechanism:**
1. Green laser (532 nm) excites both m_s = 0 and m_s = ±1 to ^3E
2. m_s = 0 decays radiatively back to m_s = 0 ground state (bright)
3. m_s = ±1 has a competing pathway through the singlet states (intersystem
   crossing, ISC) that is **non-radiative** (dark)
4. The singlet pathway preferentially decays to m_s = 0 (spin polarization)
5. Net effect: m_s = ±1 fluoresces ~30% less than m_s = 0

Therefore: when microwaves drive the m_s = 0 → m_s = ±1 transition,
fluorescence **drops**. This is the ODMR dip.


## 3. Ground-State Spin Hamiltonian

The full ground-state Hamiltonian for NV^- in an external magnetic field is
[Doherty 2013, Eq. 20; Maze 2011, Eq. 1; Rondin 2014, Eq. 1]:

```
H = H_ZFS + H_Zeeman + H_hyperfine + H_electric + H_strain
```

For ODMR magnetometry, the dominant terms are:

```
H = D·S_z^2 + E·(S_x^2 - S_y^2) + gamma_e·B·S
```

where:
- **S = (S_x, S_y, S_z)** are the spin-1 operators (3x3 matrices)
- **z-axis** is along the NV symmetry axis (the N→V direction)
- **D** is the axial zero-field splitting parameter
- **E** is the transverse zero-field splitting (strain/electric field)
- **gamma_e** is the NV gyromagnetic ratio = g_e * mu_B / h
- **B** = (B_x, B_y, B_z) is the external magnetic field in the NV frame

### Spin-1 Matrices

In the {|+1⟩, |0⟩, |-1⟩} basis (eigenstates of S_z):

```
         ⎡ 0  1  0 ⎤              ⎡ 0 -i  0 ⎤              ⎡ 1  0  0 ⎤
S_x = 1  ⎢ 1  0  1 ⎥    S_y = 1  ⎢ i  0 -i ⎥    S_z =    ⎢ 0  0  0 ⎥
      √2 ⎣ 0  1  0 ⎦         √2  ⎣ 0  i  0 ⎦              ⎣ 0  0 -1 ⎦
```

### Full Hamiltonian Matrix (in NV frame)

Writing B_∥ = B_z (along NV axis) and B_⊥ = sqrt(B_x^2 + B_y^2):

```
      ⎡ D + gamma·B_∥        gamma·B_⊥/√2         E          ⎤
H =   ⎢ gamma·B_⊥/√2              0           gamma·B_⊥/√2   ⎥
      ⎣      E             gamma·B_⊥/√2       D - gamma·B_∥   ⎦
```

Note: The off-diagonal E terms couple |+1⟩ and |-1⟩. The B_⊥ terms couple
|0⟩ to |±1⟩. When E = 0 and B_⊥ = 0, the eigenvalues are trivially
D + gamma·B_∥, 0, and D - gamma·B_∥.

### Full Hamiltonian Including Electric/Strain Fields

For completeness, the full ground-state Hamiltonian including Stark and
strain interactions is [Doherty 2013, Eq. 20]:

```
H_GS = D · [S_z^2 - (2/3)]
     + gamma_e · S · B
     + d_∥ · (E_z + delta_z) · [S_z^2 - (2/3)]
     + d_⊥ · (E_x + delta_x) · (S_y^2 - S_x^2)
     + d_⊥ · (E_y + delta_y) · (S_x·S_y + S_y·S_x)
```

where E is the electric field and delta is the strain field, with coupling
constants d_∥ = 0.35 Hz·cm/V and d_⊥ = 17 Hz·cm/V. The transverse terms
(d_⊥) are the physical origin of the E parameter in the simplified
Hamiltonian — they arise from strain inhomogeneity in the crystal.

In practice, we absorb the strain/electric field effects into the effective
D and E parameters:
- D_eff = D + d_∥ · (E_z + delta_z)  ← axial shift (tiny, negligible)
- E_eff = d_⊥ · sqrt((E_x + delta_x)^2 + (E_y + delta_y)^2)  ← transverse splitting

### What we neglect (and why it's okay)

| Term | Magnitude | Why negligible |
|------|-----------|----------------|
| Nuclear hyperfine (^14N) | A_∥ = -2.14 MHz, A_⊥ = -2.70 MHz [Felton 2009] | Below our linewidth (~10 MHz FWHM). Each m_s level splits into 3 (I=1), but unresolved in ensemble ODMR |
| Nuclear hyperfine (^13C) | ~13 MHz (nearest neighbor) | Natural abundance 1.1%. Contributes to inhomogeneous broadening, not peak positions |
| Electric field | d_∥ = 0.35 Hz/(V/cm), d_⊥ = 17 Hz/(V/cm) [Doherty 2013] | Negligible at lab fields. Strain term (below) is the relevant analog |
| Quadrupole (^14N) | P = -5.01 MHz [Felton 2009] | Only relevant for hyperfine-resolved experiments |

**For ensemble CW-ODMR magnetometry (our case), the 3-parameter Hamiltonian
H = D·S_z^2 + E·(S_x^2 - S_y^2) + gamma·B·S captures all resolvable physics.**


## 4. Zero-Field Splitting: D and E

### D — Axial Zero-Field Splitting

**Physical origin:** Spin-spin dipolar interaction between the two unpaired
electrons in the ^3A_2 state. The NV axis provides an anisotropic environment:
the two electrons' dipolar coupling is different along vs. perpendicular to
the N-V axis [Doherty 2013, §3.2].

D = (3/4) * (mu_0 / 4pi) * (g_e * mu_B)^2 * ⟨(1 - 3cos^2(theta)) / r^3⟩

where the average is over the two-electron wavefunction. This purely
quantum-mechanical quantity depends on the spatial extent of the defect
wavefunction within the diamond lattice.

**Values:**
- Room temperature: D = 2870.2 ± 0.1 MHz [Acosta 2010]
- At T = 0 K (extrapolated): D ≈ 2877.4 MHz [Acosta 2010]
- Our measurement at ~12K: **D_cryo ≈ 2878 MHz** (from center of dips at 10G)

### Temperature Dependence of D

D(T) decreases with increasing temperature. The physical origin is
~85% electron-phonon interaction and ~15% thermal lattice expansion
[Doherty 2014].

The original measurement [Acosta 2010] established:
```
dD/dT = -74.2(7) kHz/K   (near room temperature, 280-330K)
```

**Modern best-practice model: Cambria 2023 two-phonon Bose-Einstein:**

```
D(T) = D_0 + c_1 * n_1(T) + c_2 * n_2(T)
```

where n_i(T) = 1 / [exp(Delta_i / k_B T) - 1] is the Bose-Einstein
phonon occupation number.

| Parameter | Value | Description |
|-----------|-------|-------------|
| D_0 | ~2877.7 MHz | Zero-temperature ZFS |
| c_1 | -54.91 MHz | Coupling to phonon mode 1 |
| c_2 | -249.6 MHz | Coupling to phonon mode 2 |
| Delta_1 | 58.73 meV (681 K) | Energy of phonon mode 1 |
| Delta_2 | 145.5 meV (1689 K) | Energy of phonon mode 2 |

This model is valid from 15K to 500K without divergence (unlike polynomial
fits). It naturally produces T^4 scaling at low T consistent with the Debye
model. [Cambria 2023, Phys. Rev. B 108, L180102; arXiv: 2306.05318]

Earlier work: Chen et al. 2011 (APL 99, 161903) measured D(T) from 5.6K to
295K and fit with a fifth-order polynomial. Doherty et al. 2014 (Phys. Rev. B
90, 041201(R)) identified the thermal expansion vs electron-phonon decomposition.

**Effective thermal shift rates:**

| Temperature | dD/dT (kHz/K) | Phonon occupation |
|-------------|---------------|-------------------|
| 300 K       | -74           | Both modes active |
| 200 K       | -53           | Mode 2 freezing out |
| 100 K       | -17           | Only mode 1 contributes |
| 50 K        | ~-2           | Mode 1 nearly frozen |
| 12 K        | ~0 (< 0.1)   | Both modes exponentially frozen |

**At 12K:** n_1(12K) = 1/[exp(681/12)-1] ~ exp(-57) ≈ 0.
Both phonon modes are exponentially frozen out. D is a constant.

At our operating temperature of ~12K, a 0.1K fluctuation shifts D by
< 10 kHz — far below our 0.1 MHz bin size. This is **why we work at cryo**:
D is thermally stable, eliminating the dominant noise source in
room-temperature NV magnetometry.

### E — Transverse Zero-Field Splitting (Strain Parameter)

**Physical origin:** Breaks the C_3v symmetry. Any perturbation that
distinguishes x from y in the plane perpendicular to the NV axis:
- Crystal strain (most common in bulk diamond)
- Local electric fields (from nearby charges)
- Surface proximity (for shallow NVs)

E lifts the degeneracy of m_s = +1 and m_s = -1 at zero field. In an
ideal unstrained crystal, E = 0. In practice:

**Typical values:**
- High-quality CVD diamond: E < 1 MHz
- Natural diamond: E ~ 1-10 MHz
- Nanodiamonds: E ~ 5-20 MHz
- Single NV (can vary): E = 0 to ~10 MHz

**Effect on ODMR spectrum:**
- At zero field: two dips at D + E and D - E (instead of one degenerate dip at D)
- At finite field (B >> E/gamma): E causes a small perturbative shift, enters
  as sqrt((gamma·B_∥)^2 + E^2) in the transition frequencies
- For our fields (10-30G): gamma·B_∥ ~ 20-80 MHz >> E (typically), so E is a
  small correction

**How to measure E:** Zero-field ODMR. Two dips separated by 2E. The mean of
the two dip frequencies gives D. If E < linewidth, the two dips merge and E
is unresolved.


## 5. Zeeman Interaction and Transition Frequencies

### Decomposition of B in the NV frame

For a single NV orientation with axis n_hat:

```
B_∥ = B_vec · n_hat           (parallel component)
B_⊥ = |B_vec - B_∥ * n_hat|   (perpendicular component)
    = sqrt(|B|^2 - B_∥^2)
```

### Eigenvalues of the Ground-State Hamiltonian

The 3x3 Hamiltonian matrix has three eigenvalues. Exact analytical solutions
exist but are cumbersome. The two ODMR-active transitions are m_s = 0 → ±1.

**Exact transition frequencies** (from diagonalizing the 3x3 matrix):

For E = 0, the secular equation gives [Rondin 2014, Eq. 3]:

```
f_± = D + (gamma·B_⊥)^2 / (2D) ± gamma·B_∥ · sqrt(1 + (gamma·B_⊥)^2 / (D^2 - (gamma·B_∥)^2))
```

**Exact inverse formulas** [Rondin 2014]: Given measured f_upper and f_lower,
the field magnitude and angle can be reconstructed:

```
|B| = (1 / (sqrt(3) * gamma)) * sqrt(f_u^2 + f_l^2 - f_u*f_l - D^2)

cos^2(theta) = [-(f_u + f_l)^3 + 3*f_u^3 + 3*f_l^3] / [27*D*(gamma*B)^2]
               + 2*D^2 / [27*(gamma*B)^2] + 1/3
```

These are useful for single-NV measurements but for ensemble ODMR with
4 axes, we use the perturbative form below (fit all 8 peaks simultaneously).

This simplifies in the common regime where gamma·|B| << D (always true for
us: gamma·30G = 84 MHz << 2878 MHz):

### Perturbative result (gamma·|B| << D):

```
f_± = D + (gamma·B_⊥)^2 / (2D) ± sqrt((gamma·B_∥)^2 + E^2)
```

This is **the key fitting equation** [Doherty 2013, Eq. 22; Rondin 2014, §II.B].

Breaking it down:

| Term | Physical meaning |
|------|-----------------|
| D | Zero-field splitting (axial). Shifts both transitions equally up. |
| (gamma·B_⊥)^2 / (2D) | Second-order perpendicular Zeeman. Small upward shift of both transitions. Symmetric. |
| gamma·B_∥ | First-order parallel Zeeman. Splits the two transitions linearly. |
| E | Strain. Adds in quadrature with the parallel Zeeman. Important only when gamma·B_∥ is comparable to E. |

### Low-field limit (gamma·B_∥ >> E):

```
f_± ≈ D ± gamma·B_∥     (linear Zeeman)
```

Splitting: f_+ - f_- = 2·gamma·B_∥

### Zero-field limit (B = 0):

```
f_± = D ± E
```

Two dips at D+E and D-E. Splitting = 2E.

### Intermediate case — the complete expression:

For NV axis i with parallel projection B_∥,i and perpendicular projection B_⊥,i:

```
f_±,i = D + gamma^2 · B_⊥,i^2 / (2D) ± sqrt((gamma · B_∥,i)^2 + E^2)
```

This is the **master fitting function** for a single NV axis.


## 6. The Four NV Axes — Peak Pairing and Geometry

### The Four Orientations

In the diamond lattice, there are exactly four distinct ⟨111⟩ directions,
corresponding to the four NV axis orientations [Doherty 2013, §2.1]:

```
n_1 = [+1, +1, +1] / sqrt(3)
n_2 = [+1, -1, -1] / sqrt(3)
n_3 = [-1, +1, -1] / sqrt(3)
n_4 = [-1, -1, +1] / sqrt(3)
```

These are the four body diagonals of a cube. Each pair subtends the
tetrahedral angle: arccos(-1/3) ≈ 109.47 degrees.

### Projections onto B-field

Given an external field B = |B| * (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
in the crystal frame, the parallel projection onto each NV axis is:

```
B_∥,i = |B_vec · n_hat_i|
B_⊥,i = sqrt(|B|^2 - B_∥,i^2)
```

Explicitly:
```
B_∥,1 = |B| * |sin(theta)cos(phi) + sin(theta)sin(phi) + cos(theta)| / sqrt(3)
B_∥,2 = |B| * |sin(theta)cos(phi) - sin(theta)sin(phi) - cos(theta)| / sqrt(3)
B_∥,3 = |B| * |-sin(theta)cos(phi) + sin(theta)sin(phi) - cos(theta)| / sqrt(3)
B_∥,4 = |B| * |-sin(theta)cos(phi) - sin(theta)sin(phi) + cos(theta)| / sqrt(3)
```

**Key property:** The four projections are constrained. They are NOT independent.
Given |B|, theta, and phi, all four projections are determined.

### ODMR Spectrum: How Many Peaks?

Each NV axis i gives two transitions: f_+,i and f_-,i.
Four axes → up to **8 peaks** total (4 pairs).

**But in practice:**
- At low field, the four pairs of dips overlap because the splittings are
  comparable to the linewidth (~10 MHz FWHM in ensemble)
- The number of resolved dips depends on the angular separation of projections
  AND the linewidth

**Degenerate cases:**
- B along ⟨100⟩: All four projections are equal → 2 dips (all degenerate)
- B along ⟨110⟩: Two distinct projections → 4 dips (2 degenerate pairs)
- B along ⟨111⟩: One axis aligned (large splitting), three equivalent (small) → 4 dips
- General direction: Four distinct projections → 8 dips

### Peak Pairing Rules

Peaks come in pairs (f_+,i, f_-,i) symmetric about the center frequency:

```
center_i = D + gamma^2 · B_⊥,i^2 / (2D)    ≈ D (small correction)
splitting_i = 2 * sqrt((gamma · B_∥,i)^2 + E^2)
```

**Pairing constraints:**
1. Each pair is (approximately) symmetric about D
2. Pairs with larger splitting have their f_- dip lower and f_+ dip higher
3. The pair centers shift slightly upward with increasing B_⊥ (second-order effect)
4. Sum rule: In the ODMR spectrum, the average of all 8 dip frequencies ≈ D
   (to first order in B_⊥^2/D)

### When Do All 8 Peaks Resolve?

Condition: the minimum inter-dip gap between adjacent peaks from different NV
axes must exceed the FWHM.

For our geometry (magnet roughly along ⟨100⟩ with some tilt):
- At 10G: 2 visible dips (pairs overlap) — only the outermost pair (n4) is
  partially resolved from the cluster
- At 20G: Still partially overlapping — maybe 4-6 distinguishable features
- At 30G: **All 8 dips resolve** (minimum gap ≈ 10 MHz ≈ FWHM)
- At 35G+: Well-resolved with comfortable gaps

This is why we scan to high field: **resolving all 8 peaks uniquely determines
the B-field vector direction (theta, phi) in the crystal frame.**

### Vector Reconstruction from 4 Splittings

Given measured splittings S_i = f_+,i - f_-,i for each NV axis:

```
B_∥,i = sqrt(S_i^2/4 - E^2) / gamma
```

The parallel projections {B_∥,1, B_∥,2, B_∥,3, B_∥,4} determine B_vec.
Define signed projections p_i = s_i * B_∥,i where s_i = ±1 (sign ambiguity),
constrained by sum(p_i) having consistent sign pattern. Then:

```
Bx = (sqrt(3)/4) * (p_1 + p_2 - p_3 - p_4)
By = (sqrt(3)/4) * (p_1 - p_2 + p_3 - p_4)
Bz = (sqrt(3)/4) * (p_1 - p_2 - p_3 + p_4)
```

This has a 2-fold sign ambiguity (B and -B). Physical knowledge of the
magnet geometry (e.g., field points in -x direction) resolves this.


## 7. ODMR Transition Frequencies — The Fitting Function

### Complete Model for Ensemble ODMR

The observed ODMR spectrum is the sum of Lorentzian (or Gaussian, or Voigt)
dips from all NV axes:

```
P(f) = P_0 * [1 - sum_{i=1}^{4} sum_{s=+,-} C_{s,i} * L(f, f_{s,i}, w_i)]
```

where:
- **P_0** = baseline (off-resonance) power
- **C_{s,i}** = contrast of transition s on axis i (typically 1-10%)
- **L(f, f_0, w)** = lineshape function centered at f_0 with width w
- **f_{s,i}** = transition frequency for axis i, transition s (+ or -)

### Transition Frequencies (the physics):

For each NV axis i (i = 1,2,3,4):

```
f_±,i = D + gamma^2 · B_⊥,i^2 / (2D) ± sqrt((gamma · B_∥,i)^2 + E^2)
```

where:
```
B_∥,i = B_vec · n_hat_i
B_⊥,i = sqrt(|B|^2 - B_∥,i^2)
```

### Global Fit Parameters

**Physical parameters (5):**
| Parameter | Symbol | Expected range | Determined by |
|-----------|--------|---------------|---------------|
| Axial ZFS | D | 2877-2879 MHz | Zero-field ODMR |
| Transverse ZFS | E | 0-10 MHz | Zero-field ODMR |
| Field magnitude | \|B\| | 0-37 G | Splitting vs field |
| Polar angle | theta_B | 0-180 deg | High-field 8-peak pattern |
| Azimuthal angle | phi_B | 0-360 deg | High-field 8-peak pattern |

**Lineshape parameters (per axis, but often assumed equal across axes):**
| Parameter | Symbol | Expected range |
|-----------|--------|---------------|
| FWHM | w_i | 5-15 MHz (ensemble) |
| Contrast | C_i | 1-10% |

### Fitting Strategy (for this campaign)

**Step 1: Zero-field measurement** (Saturday)
- Measure ODMR with no applied field (Earth's field only, ~0.5G)
- Observe 2 dips at D ± E
- Extract D and E directly from dip positions
- D = (f_+ + f_-) / 2
- E = (f_+ - f_-) / 2

**Step 2: High-field measurement** (30G+ data from this campaign)
- At 30G, all 8 dips resolved
- D and E known from Step 1 (fixed)
- Fit theta_B and phi_B to match the 4 splitting patterns
- This is an over-determined system: 4 splittings, 2 unknowns

**Step 3: Intermediate fields** (10-28G data)
- D, E, theta_B, phi_B all known
- Verify Zeeman linearity: splitting should scale as 2*gamma*B_∥,i
- Fit individual dip positions to extract per-axis linewidths, contrasts
- Any systematic deviation from the model indicates higher-order effects
  or field inhomogeneity

**Step 4: Sensitivity characterization**
- At operating field (e.g., 10G), use the steepest dip slope + noise floor
- delta_B = sigma_noise / (dP/dB) = sigma_noise / (dP/df * df/dB)
- df/dB = gamma for the parallel projection of the most-aligned axis


## 8. Linewidth and FWHM

The observed ODMR linewidth has multiple contributions [Barry 2020, §III]:

### Intrinsic Broadening Mechanisms

**1. Inhomogeneous broadening (dominant in ensemble)**
- Different NV centers see slightly different local environments
- Sources: strain distribution, ^13C hyperfine, local charge variations
- Gives a Gaussian contribution to the lineshape
- Typical: 1-10 MHz for diamond with natural ^13C abundance (1.1%)
- Reduced in isotopically enriched ^12C diamond to ~0.1 MHz (single NV)

**2. ^14N Hyperfine**
- Each NV has a ^14N nucleus (I = 1)
- Splits each m_s transition into 3 lines separated by A_∥ ≈ 2.14 MHz
- In ensemble with typical linewidth > 5 MHz, these are unresolved
- Contributes ~4.3 MHz apparent broadening (spread of the 3 hyperfine lines)

**3. ^13C Hyperfine**
- Natural abundance ^13C nuclei near the NV
- Each nearby ^13C splits the line by ~13 MHz (nearest neighbor) down to kHz (distant)
- Ensemble: statistical distribution of ^13C configurations → Gaussian broadening
- Contribution: ~1-3 MHz for natural abundance, negligible for ^12C enriched

**4. Power broadening**
- Excessive MW power broadens the ODMR dip (Rabi physics)
- For CW-ODMR: FWHM_power = sqrt(FWHM_0^2 + Omega^2) where Omega is Rabi frequency
- At our +5 dBm: may contribute a few MHz. Trade-off with contrast (more power → more contrast but broader)

**5. Homogeneous broadening (T2-related)**
- Intrinsic T2 of NV^- ground state: ~1 ms in bulk diamond
- Contribution to linewidth: ~1/pi/T2 ≈ 0.3 kHz — negligible for CW-ODMR
- T2* (ensemble dephasing) is what matters: typically 1-10 us → ~30-300 kHz
  Still small compared to inhomogeneous broadening in ensemble

### Lineshape Model

For ensemble NV centers in diamond with natural ^13C abundance:

**Best model: Voigt profile** (convolution of Lorentzian + Gaussian)
- Lorentzian component: power broadening + homogeneous (small)
- Gaussian component: inhomogeneous strain, ^13C, ^14N hyperfine

**Simpler alternatives:**
- Single Lorentzian: acceptable fit for most ensemble ODMR
- Single Gaussian: sometimes better for strain-dominated broadening

For our data: start with Lorentzian, try Gaussian, and if neither fits
well in the wings, use Voigt.

### Our Observed Linewidth

From 10G campaign data:
- FWHM ≈ 10 MHz (per dip, from filtered grand average)
- This is consistent with ensemble NV in mm-scale diamond with
  natural ^13C abundance and moderate strain

The linewidth should NOT change significantly with applied field (it's
dominated by local environment, not the external field) — unless
field inhomogeneity across the diamond becomes significant.


## 9. CW-ODMR Sensitivity

### Shot-Noise Limited Sensitivity [Barry 2020, Eq. 17]

For CW-ODMR magnetic field sensing, the photon-shot-noise-limited
minimum detectable field change per unit bandwidth is:

```
eta_cw = P_lineshape / gamma_e * (Delta_f / (C * sqrt(R)))    [T / sqrt(Hz)]
```

where the lineshape prefactor for a **Lorentzian** dip is:

```
P_lineshape = 4 / (3 * sqrt(3)) ≈ 0.770
```

and:
- Delta_f = FWHM of the ODMR dip (Hz)
- C = ODMR contrast (dimensionless, e.g., 0.05 for 5%)
- R = photon detection rate (photons/s)
- gamma_e = 28.024 GHz/T (NV gyromagnetic ratio)

The sensitivity improves (smaller eta) with:
- Narrower linewidth (smaller Delta_f)
- Higher contrast (larger C)
- More photons (larger R, scales as 1/sqrt(R))

### Comparison of Magnetometry Protocols

| Protocol | Coherence limit | Bandwidth | Typical sensitivity |
|----------|----------------|-----------|-------------------|
| CW-ODMR | ~1/linewidth | DC (broadband) | ~nT/sqrt(Hz) |
| Ramsey (DC) | T2* | DC (broadband) | ~460 fT/sqrt(Hz) |
| Hahn Echo (AC) | T2 | Narrowband | ~210 fT/sqrt(Hz) |
| CPMG/XY8 (AC) | T2 (extended) | Tunable kHz-MHz | Best narrowband |

**Key distinction:** DC sensitivity is limited by T2*; AC sensitivity by T2.
Since T2 >> T2* (typically 1-2 orders of magnitude), pulsed AC protocols
achieve much better sensitivity. CW-ODMR (our method) is the simplest and
gives broadband DC sensitivity — adequate for our magnetometry goals.

**In terms of our experimental observables:**

```
eta = sigma_noise / max_slope    [T / sqrt(Hz)]
```

where:
- sigma_noise = RMS noise floor of the filtered ODMR spectrum (per sqrt(Hz))
- max_slope = maximum |dP/df| on the steepest part of the dip

This is the **slope-based sensitivity** — directly measured from data without
needing to know photon rates.

### Relation to Our Analysis

From the analysis script output:
```
delta_f = wing_noise / max_slope      (frequency noise floor, MHz)
delta_B = delta_f / gamma_NV          (field noise floor, Gauss)
```

This gives the noise-equivalent field per measurement (per grand average).
To get sensitivity per sqrt(Hz), divide by sqrt(T_measurement).

### Bin-Limited Regime

When the noise floor drops below the frequency resolution:

```
bin_effect = max_slope * bin_size
```

If sigma_noise < bin_effect, we are **bin-limited**: the measurement
precision is limited by our 0.1 MHz frequency step, not by noise.
This is always our target — it means we have more than enough signal
averaging.

### Sensitivity Scaling

```
eta(N) = eta(1) / sqrt(N)    (white noise regime)
```

Until drift onset (typically N ≈ 15-30 pairs in our setup), after which
eta plateaus. The optimal measurement strategy is to average within the
white noise regime and repeat independent measurements.


## 10. References

All PDFs are in `analysis/NV_Physics_References/`.

### [Doherty 2013]
M. W. Doherty, N. B. Manson, P. Delaney, F. Jelezko, J. Wrachtrup, and
L. C. L. Hollenberg, "The nitrogen-vacancy colour centre in diamond,"
*Physics Reports* **528**, 1-45 (2013).
arXiv: 1302.3288
**PDF:** `Doherty_2013_NV_Centre_Review_PhysRep.pdf`
**Used for:** §1-5. Comprehensive review of NV structure, Hamiltonian, symmetry.

### [Maze 2011]
J. R. Maze, A. Gali, E. Togan, Y. Chu, A. Trifonov, E. Kaxiras, and
M. D. Lukin, "Properties of nitrogen-vacancy centers in diamond: the
group theoretic approach," *New Journal of Physics* **13**, 025025 (2011).
arXiv: 1010.1338
**PDF:** `Maze_2011_NV_Properties_Group_Theory.pdf`
**Used for:** §3. Group-theoretic derivation of the Hamiltonian, spin-orbit, spin-spin.

### [Barry 2020]
J. F. Barry, J. M. Schloss, E. Bauch, M. J. Turner, C. A. Hart, L. M. Pham,
and R. L. Walsworth, "Sensitivity optimization for NV-diamond magnetometry,"
*Reviews of Modern Physics* **92**, 015004 (2020).
arXiv: 1903.08176
**PDF:** `Barry_2020_NV_Magnetometry_Sensitivity_RevModPhys.pdf`
**Used for:** §8-9. Definitive treatment of sensitivity, linewidth mechanisms, noise.

### [Rondin 2014]
L. Rondin, J.-P. Tetienne, T. Hingant, J.-F. Roch, P. Maletinsky, and
V. Jacques, "Magnetometry with nitrogen-vacancy defects in diamond,"
*Reports on Progress in Physics* **77**, 056503 (2014).
arXiv: 1311.5214
**PDF:** `Rondin_2014_NV_Magnetometry_Review.pdf`
**Used for:** §5-6. Zeeman splitting, transition frequencies, angular dependence.

### [Doherty 2012]
M. W. Doherty, F. Dolde, H. Fedder, F. Jelezko, J. Wrachtrup, N. B. Manson,
and L. C. L. Hollenberg, "Theory of the ground-state spin of the NV- center
in diamond," *Physical Review B* **85**, 205203 (2012).
arXiv: 1108.5546
**PDF:** `Doherty_2012_Ground_State_Spin_Theory.pdf`
**Used for:** §3-4. Detailed ground-state spin Hamiltonian, spin-spin interaction origin of D.

### [Acosta 2010]
V. M. Acosta, E. Bauch, M. P. Ledbetter, A. Waxman, L.-S. Bouchard, and
D. Budker, "Temperature dependence of the nitrogen-vacancy magnetic resonance
in diamond," *Physical Review Letters* **104**, 070801 (2010).
arXiv: 1007.0011
**PDF:** `Acosta_2010_Temperature_Dependence_ZFS.pdf`
**Used for:** §4. D(T) temperature dependence, cryogenic thermal stability.

### [Felton 2009]
S. Felton, A. M. Edmonds, M. E. Newton, P. M. Martineau, D. Fisher,
D. J. Twitchen, and J. M. Baker, "Hyperfine interaction in the ground state
of the negatively charged nitrogen vacancy center in diamond," *Physical
Review B* **79**, 075203 (2009).
arXiv: 0807.0768
**PDF:** `Felton_2009_Hyperfine_Interaction_NV.pdf`
**Used for:** §3, §8. ^14N and ^13C hyperfine constants, contribution to linewidth.

### [Mittiga 2018]
T. Mittiga, S. Hsieh, C. Zu, B. Kobrin, F. Machado, P. Bhatt, N. Y. Yao,
and M. S. Ziabari, "Imaging the local charge environment of nitrogen-vacancy
centers in diamond," *Physical Review Letters* **121**, 246402 (2018).
arXiv: 1809.01548
**PDF:** `Mittiga_2018_Local_Charge_Environment_NV.pdf`
**Used for:** §1. Charge state stability (NV^- vs NV^0), local charge environment.

### [Clevenson 2015]
H. Clevenson, M. E. Trusheim, C. Teale, T. Schr\"oder, D. Braje, and
D. Englund, "Broadband magnetometry and temperature sensing with a
light-trapping diamond waveguide," *Nature Physics* **11**, 393-397 (2015).
arXiv: 1406.3248
**PDF:** `Clevenson_2015_Broadband_Magnetometry_Diamond_Waveguide.pdf`
**Used for:** §9. Ensemble magnetometry sensitivity, light-trapping geometry.

### [Gruber 1997]
A. Gruber, A. Drabenstedt, C. Tietz, L. Fleury, J. Wrachtrup, and
C. von Borczyskowski, "Scanning confocal optical microscopy and magnetic
resonance on single defect centers," *Science* **276**, 2012-2014 (1997).
DOI: 10.1126/science.276.5321.2012
**PDF:** Not available (behind paywall). Jacob can obtain institutional access.
**Used for:** Historical — first single-NV ODMR demonstration.

### [Cambria 2023]
M. C. Cambria, A. Norambuena, H. T. Dinani, G. Thiering, A. Gardill,
I. Kemeny, Y. Li, V. Lordi, A. Gali, J. R. Maze, and S. Kolkowitz,
"Physically motivated analytical expression for the temperature dependence
of the zero-field splitting of the nitrogen-vacancy center in diamond,"
*Physical Review B* **108**, L180102 (2023).
arXiv: 2306.05318
**PDF:** Not downloaded — Jacob can obtain from arXiv.
**Used for:** §4. Modern D(T) model: two-phonon Bose-Einstein expression valid 15-500K.

### [Chen 2011]
X.-D. Chen, C.-H. Dong, F.-W. Sun, C.-L. Zou, J.-M. Cui, Z.-F. Han, and
G.-C. Guo, "Temperature dependent energy level shifts of nitrogen-vacancy
centers in diamond," *Applied Physics Letters* **99**, 161903 (2011).
**PDF:** Not downloaded — Jacob can check institutional access.
**Used for:** §4. Fifth-order polynomial D(T) from 5.6K to 295K.

### [Doherty 2014]
M. W. Doherty, V. M. Acosta, A. Jarmola, M. S. J. Barson, N. B. Manson,
D. Budker, and L. C. L. Hollenberg, "Temperature shifts of the resonances
of the NV- center in diamond," *Physical Review B* **90**, 041201(R) (2014).
arXiv: 1310.7303
**PDF:** Not downloaded — available on arXiv.
**Used for:** §4. Decomposition: thermal expansion (~15%) + electron-phonon (~85%).

---

## Appendix A: Numerical Constants

| Constant | Symbol | Value | Unit |
|----------|--------|-------|------|
| Electron g-factor | g_e | 2.0028 | dimensionless |
| Bohr magneton | mu_B | 9.2741 x 10^-24 | J/T |
| NV gyromagnetic ratio | gamma_NV | 2.8024 | MHz/G |
| (same in SI) | gamma_NV | 28.024 | GHz/T |
| Zero-field splitting (RT) | D_RT | 2870.2 | MHz |
| Zero-field splitting (cryo) | D_cryo | ~2878 | MHz |
| ^14N hyperfine A_∥ | A_∥ | -2.14 | MHz |
| ^14N hyperfine A_⊥ | A_⊥ | -2.70 | MHz |
| ^14N quadrupole | P | -5.01 | MHz |
| Planck constant | h | 6.6261 x 10^-34 | J·s |
| Diamond lattice constant | a | 3.567 | Angstrom |
| Tetrahedral angle | - | 109.47 | degrees |
| GSLAC field | B_GSLAC | ~1024 | Gauss |

### NV Axis Vectors (Crystal Frame)

```
n_1 = [+1, +1, +1] / sqrt(3) = [0.5774, 0.5774, 0.5774]
n_2 = [+1, -1, -1] / sqrt(3) = [0.5774, -0.5774, -0.5774]
n_3 = [-1, +1, -1] / sqrt(3) = [-0.5774, 0.5774, -0.5774]
n_4 = [-1, -1, +1] / sqrt(3) = [-0.5774, -0.5774, 0.5774]
```

Dot products between axes:
```
n_i · n_j = -1/3    (for i ≠ j)
arccos(-1/3) = 109.47 degrees
```

## Appendix B: Our Measured Values

From the Calibrated NV Cryo campaign series (2026-02-25 onwards):

| Quantity | Value | Source |
|----------|-------|--------|
| D_cryo | ~2878 MHz | Center of dips at 10G |
| Splitting rate (n4) | 4.26 MHz/G | Linear fit across 10.0-10.5G |
| Implied projection angle | ~40.5 deg | arccos(4.26 / 2*2.8024) |
| Ensemble FWHM | ~10 MHz | Filtered grand average at 10G |
| Contrast (10G) | ~6-7.5% | Per-pair, varies with pair |
| Contrast trend | decreasing with field | 7.55% at 10.1G → 6.88% at 13.0G |
| Baseline power | ~0.8 uW | 532nm, PM101 at 632nm setting |
| Temperature | 11-13 K | Lakeshore 335, sensor A |
| Operating dwell | 20 ms | Optimized (FOM = 9.72e-6) |
| Sensitivity (10G) | ~14 uT/sqrt(Hz) | From noise floor analysis |
| Drift onset | N ≈ 15-30 pairs | ~10-20 min wall clock |
| Bin-limited at | N ≈ 1-5 (with normalization) | After baseline normalization + Gaussian LP |

---

*This document will be expanded with Rabi oscillation physics, pulsed ODMR,
T1/T2 relaxometry, and dark matter sensitivity projections in future versions.*
