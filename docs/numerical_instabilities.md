# Numerical Instabilities in WAVEWATCH III (and third-generation spectral wave models generally)

WW3 solves the action balance equation on a 4-D grid (2 geographic + 2 spectral: frequency, direction), advancing it via **operator splitting**: separate sub-steps for geographic propagation, spectral (frequency/direction) propagation ("intra-spectral" propagation, i.e. refraction and current-induced frequency shifting), and source term integration. Almost every instability class below traces back to how one of these three sub-steps is discretized and how they interact through the splitting.

---

## 1. Geographic (spatial) advection — classic CFL violation

**Reason:** WW3's default explicit propagation schemes (first-order upwind, or higher-order UQ / ULTIMATE-QUICKEST / UNO2) are explicit in time. Stability requires the Courant number for the fastest-traveling spectral component (highest group velocity, i.e. lowest frequency / longest period in deep water, or shallowest depth) to stay below 1:

Δt ≤ C_CFL · Δx / c_g,max

**Pattern:** Blow-up (NaNs, unbounded energy growth) that appears first at the grid cells with the largest group velocity / smallest Δx (e.g., long-period swell in coarse, high-latitude, or curvilinear grid cells; or in unstructured-grid triangles with very small edges). On spherical grids it is worst near the poles where meridians converge and the effective Δx shrinks.

**Detection:** Sudden appearance of energy spikes/NaNs in Hs or spectral output, often localized initially, growing outward; log files reporting time-step violations if diagnostics are on; visually, "checkerboard" or single-cell blow-up patterns.

**Solution:**
- WW3 uses a **dynamically adjusted / sub-stepped time integration**: the "overall" time step is split into smaller spatial propagation sub-steps chosen automatically (or set by the user) to respect the CFL limit everywhere on the grid.
- Reduce the user-specified maximum global/propagation time step (`DTMAX`, spatial propagation `DTXY` in `ww3_grid.inp`).
- Use a coarser high-frequency cutoff or courser grid near problem areas, or locally refine/smooth the grid to avoid very small Δx cells.
- For unstructured grids (implicit / PDLIB solver), an implicit scheme removes the CFL restriction on Δx at the cost of solving a sparse linear system each step — this is the standard fix for very fine unstructured coastal meshes.

---

## 2. Spectral-space (intra-spectral) advection — refraction & current-induced frequency shift

**Reason:** The rate of directional turning (refraction, dθ/dt) and frequency shifting (dσ/dt) due to depth/current gradients can be locally very large — e.g., steep bathymetry, strong current shear, or wave blocking. This is itself an advection process in (σ,θ) space and has its own CFL-type constraint, separate from and generally stricter than the geographic one.

**Pattern:** Spurious energy piling up at the edges of the directional/frequency spectrum (energy "wrapping" or ringing near θ_max/θ_min or the highest resolved frequency), oscillatory noise concentrated in specific grid cells with strong bathymetric gradients or current jets, sometimes visible as a checkerboard pattern in directional spread output.

**Detection:** Inspect 2-D spectra at flagged points for unphysical concentration at bin edges; monitor the diagnostic refraction/frequency-shift time step chosen by the adaptive algorithm — if it's collapsing to very small values repeatedly at certain points, that's the signature.

**Solution:**
- WW3 automatically sub-steps this term too (separate `DTKTH` / refraction time step, historically recommended ≈ half the overall time step as a starting point).
- Locally smooth the bathymetry/current fields feeding the model (unresolved sharp depth gradients are numerically unrealistic anyway).
- Increase directional/frequency resolution in high-gradient regions.
- Cap or limit refraction source in extreme-gradient cells (used historically, though this is a workaround, not a physical fix).

---

## 3. Current blocking / opposing-current singularities

**Reason:** When waves propagate into a strongly opposing current, the group velocity relative to the ground can approach zero (wave blocking) or the action density formulation becomes singular as the intrinsic frequency solution for the dispersion relation loses uniqueness/stability. This is a genuine physical singularity that the discretization struggles with.

**Pattern:** Energy pileup at the blocking point, unbounded growth in specific spectral bins associated with the blocked component, followed by NaNs propagating outward from that cell.

**Detection:** Localize blow-up onset to grid cells with strong current gradients opposing wave direction; check current-field magnitude/gradient at flagged cells against the local group velocity.

**Solution:** No fully general fix exists — this is a known open weak point of the action balance / geometric-optics formulation. Practical mitigations: limit maximum opposing current magnitude used by the model (clip/smooth current fields), increase spectral resolution near blocking frequencies, reduce time step locally, or accept energy loss/dissipation of blocked components as a physical proxy (some dissipation parameterizations partially cover this via enhanced breaking near strong currents).

---

## 4. Source term integration instability / stiffness

**Reason:** The physics source terms (wind input S_in, nonlinear interactions S_nl, whitecapping/dissipation S_ds, bottom friction, depth-induced breaking S_db) can have very short intrinsic timescales relative to the propagation time step — especially strong wind forcing or shallow-water breaking. Naive explicit (forward Euler) integration of these stiff terms blows up or oscillates if Δt exceeds the local physical relaxation time.

**Pattern:** Oscillatory or explosive growth in wave energy at points with strong forcing (high winds, very shallow water) even though the propagation scheme itself is stable; sometimes shows as high-frequency temporal "ringing" in Hs time series at a single point rather than spatial spreading.

**Detection:** Point-time-series diagnostics showing oscillation/divergence uncorrelated with propagation CFL; occurs preferentially under high wind-speed events or in very shallow nested grids.

**Solution — this is the most thoroughly engineered part of WW3:**
- **Dynamic adaptive time stepping with a minimum source-term time step** (Tolman 2002a): the source-term sub-step is automatically shrunk in regions/times of rapid spectral change, down to a user-set minimum (often ~5–10 s).
- **Semi-/fully-implicit source term integration** — WW3's default is the fully implicit scheme of **Hargreaves and Annan (2000)**, combined with the adaptive step + limiter method of Tolman (2002a), which removes most of the stiffness problem outright.
- **Limiters on the maximum spectral change per time step** — inherited conceptually from the older WAM "dynamic-implicit" limiter approach; WW3 uses a non-convergent-style limiter that keeps growth bounded but becomes part of the effective solution (a known trade-off: fine for operational/engineering runs, less clean for research since the limiter's influence can't be separated from the "pure" physics).
- Practically: reduce `DTMIN` (minimum source term time step) and/or reduce the overall time step in cases with hurricane-force winds or very shallow, highly dissipative nested domains.

---

## 5. Depth-induced breaking stiffness (shallow water)

**Reason:** A special case of #4 — the depth-limited breaking source term (e.g., Battjes-Janssen-type formulations) becomes extremely stiff as depth → 0, since the breaking timescale shrinks rapidly with depth.

**Pattern:** Instability or unphysical energy collapse concentrated in the surf zone / very shallow nested grids; sensitive to grid resolution near the shoreline.

**Detection:** Isolate to shallow-water/surf-zone cells; test sensitivity by refining `DTMIN` there.

**Solution:** Smaller minimum source-term time step in shallow nests; implicit treatment of the breaking term (as above); avoid resolving unrealistically shallow (near-zero) depths without adequate wetting/drying treatment.

---

## 6. Nonlinear wave-wave interactions (DIA) — aliasing and unphysical energy transfer

**Reason:** WW3's default fast approximation for the quadruplet nonlinear interactions is the **Discrete Interaction Approximation (DIA)**, a coarse, computationally cheap simplification of the full Boltzmann integral. Because it only samples a small, fixed configuration of interacting wavenumber quadruplets, it can alias energy into unphysical parts of the spectrum, particularly when combined with a coarse directional/frequency grid.

**Pattern:** Spurious secondary peaks in the spectrum, unphysically rapid or slow spectral downshifting, and — in aggregate with other stiff terms — can seed local instabilities that then interact with S_in/S_ds feedback loops.

**Detection:** Compare DIA output against the exact nonlinear transfer (WRT / full Boltzmann integral) offline for representative spectra; look for spurious high-frequency bumps.

**Solution:** No general instability fix beyond the source-term time-stepping controls in #4 (DIA itself isn't usually the direct cause of blow-up, more of an accuracy/bias issue); using the more expensive exact method (e.g., WRT, or the newer neural-network-based approximations) reduces aliasing artifacts if affordable.

---

## 7. Garden Sprinkler Effect (GSE) — not strictly an instability, but a related and pervasive numerical pathology

**Reason:** Because frequency and direction are discretized, swell that should disperse continuously in space instead disperses as if it were a finite number of discrete "rays," each with a slightly different discrete speed and direction. Over long propagation distances this causes an initially continuous swell field to break up into separate discrete "sprinkler jets."

**Pattern:** Swell fields downstream from a storm show banded, comb-like or "sprinkler" patterns in Hs, especially visible when directional/frequency resolution is coarse and propagation distance is long. Higher-order, more accurate advection schemes (like the UQ/ULTIMATE-QUICKEST scheme WW3 uses) actually **highlight** GSE rather than mask it, because they have less inherent numerical diffusion.

**Detection:** Visual inspection of Hs maps downstream of a distant storm for streaky/discretized bands; sensitivity test — repeat with doubled directional resolution and see if the banding disappears.

**Solution:**
- Increase directional (and to a lesser extent frequency) spectral resolution — the most robust fix, but computationally expensive.
- Add an explicit **diffusion term in spectral space** aligned with the spectral direction θ (Booij & Holthuijsen 1987-style GSE correction terms for frequency and directional dispersion).
- **Averaging methods in θ-space** (Lavrenov & Onvlee 1995) as a cheaper alternative to full diffusion tensors.
- Tolman (2002b) "divergence" method: adds divergence to the advection field as a cheaper correction than a full diffusion tensor, implemented as an option in WW3.
- Note: GSE alleviation methods are generally only available for structured/regular grids in WW3; **unstructured grids currently have no GSE alleviation scheme**, which is a known limitation.

---

## 8. Time-step variability instability (operator-splitting / variable-Δt effects)

**Reason:** Even when each individual sub-step (propagation, refraction, source terms) is comfortably within its own stability limit, operator splitting between sub-processes with *different, adaptively varying* time steps can introduce instability that would not occur with a constant time step — a general property of split explicit schemes noted in the broader numerical PDE literature, not unique to WW3 but directly relevant given how aggressively WW3 varies its sub-steps.

**Pattern:** Instability appears intermittently, correlated with times/places where the adaptive time-stepping algorithm is rapidly changing its sub-step size (e.g., transitioning between calm and storm conditions), rather than being tied to a single obviously-too-large Δt.

**Detection:** Hard to catch analytically; empirically, check whether instability onset coincides with large changes in the diagnosed adaptive time step from one step to the next, rather than with the step size itself being large.

**Solution:** Limit how fast the adaptive algorithm is allowed to change the sub-step size between iterations (a smoother ramping), or fall back to a more conservative fixed minimum step during rapidly evolving conditions.

---

## 9. Multi-grid nesting instabilities

**Reason:** Boundary spectra passed from a coarse "parent" grid to a finer nested grid are only available at the parent's resolution and update frequency. Interpolation errors, resolution mismatches (in space, direction, or frequency), or update-frequency mismatches can inject spurious high-wavenumber noise or energy imbalance at the nest boundary.

**Pattern:** Noise or unphysical banding concentrated near nest boundaries; energy discontinuities visible where the fine grid meets the coarse grid; can grow if the fine grid has stiffer source terms (shallower water) than the coarse grid assumed.

**Detection:** Inspect Hs/spectra fields specifically along nest boundaries; compare energy balance across the boundary.

**Solution:** Increase boundary update frequency, use consistent (or finer) spectral resolution in the nest relative to the parent, use a buffer zone with gradually increasing resolution rather than an abrupt jump, and ensure the nest's own source-term time-stepping is fine enough for its shallower depths.

---

## 10. Negative or unphysical energy densities

**Reason:** Explicit/limiter-based dissipation or breaking corrections can, in pathological cases, overcorrect a spectral bin, driving its energy density below zero, which has no physical meaning and can propagate as an error into neighboring bins via the nonlinear interaction terms.

**Pattern:** Small negative values appear in raw spectral output at isolated frequency/direction bins, typically at the high-frequency tail or in heavily dissipating shallow cells; can seed slowly growing numerical noise if unclipped.

**Detection:** Direct inspection of raw 2-D spectral output (not just integrated parameters like Hs, which can mask small negative bins); automated post-processing checks for E(f,θ) < 0.

**Solution:** WW3 (like WAM/SWAN) applies **clipping/limiters** to keep spectral densities non-negative after each source-term update; if you see this in outputs, it usually indicates the source-term time step is too large for local conditions (tie back to fix #4).

---

## Summary table

| # | Instability | Where it shows up | Primary fix |
|---|---|---|---|
| 1 | Geographic CFL violation | Long swell, small Δx, poles, unstructured fine cells | Adaptive/sub-stepped propagation Δt; implicit solver for unstructured grids |
| 2 | Spectral-space CFL (refraction/freq shift) | Steep bathymetry/current gradients | Sub-stepped refraction Δt; smoother forcing fields |
| 3 | Current blocking singularity | Strong opposing currents | Clip/smooth currents; finer spectral res; no full fix exists |
| 4 | Source term stiffness | High wind, shallow water | Implicit (Hargreaves & Annan) scheme + adaptive min. time step + limiter (Tolman 2002a) |
| 5 | Depth-breaking stiffness | Surf zone / very shallow nests | Small `DTMIN`; implicit breaking treatment |
| 6 | DIA nonlinear aliasing | Coarse spectral grids | Finer resolution; exact nonlinear solver if affordable |
| 7 | Garden Sprinkler Effect | Long-range swell, coarse directional res | Finer directional res; GSE diffusion/averaging/divergence corrections (structured grids only) |
| 8 | Variable-Δt splitting instability | Rapidly changing adaptive steps | Smooth step-size changes; conservative minimum step |
| 9 | Nest boundary instability | Parent/child grid interfaces | Finer/matched nest resolution; frequent boundary updates |
| 10 | Negative energy densities | High-frequency tail, shallow dissipative cells | Non-negativity clipping; smaller source term Δt |

---

### Key references
- Tolman, H.L. (1992), *Effects of numerics on the physics in a third-generation wind-wave model*, J. Phys. Oceanogr. 22.
- Tolman, H.L. (2002a), adaptive time-stepping and limiter method (WW3 Tech Note / User Manual).
- Tolman, H.L. (2002b), *Alleviating the Garden Sprinkler Effect in wind wave models*, Ocean Modelling.
- Hargreaves, J.C. and Annan, J.D. (2000), fully implicit source term integration scheme.
- Booij, N. and Holthuijsen, L.H. (1987), *Propagation of ocean waves in discrete spectral wave models*, J. Comput. Phys.
- WW3 Development Group, *WAVEWATCH III Version X.XX User Manual and System Documentation* (NOAA/NCEP), section on numerics/time stepping.

*Caveat: several of these (current blocking, splitting-induced instability, nesting artifacts) don't have a single canonical "solved" fix in the literature — they're active areas of model development and are usually managed with conservative settings rather than eliminated outright.*