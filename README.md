# GBH Rayleigh Phase-Coherence Test on Planck SMICA

This repository contains a minimal implementation of the Rayleigh phase-coherence
test used in the paper:

> **"The Gödelian Boundary: A Structural Hypothesis on Minimal Physical Change"**  
> Jacob Guedalia (2025)

The goal is to illustrate that the Gödelian Boundary Hypothesis (GBH) is empirically
operationalizable. GBH predicts that a minimal pre-inflation event, if it leaves any
observable trace, must do so via residual **phase structure** rather than amplitude
in the CMB. The structural prediction is: **anomalous phase alignment on a thin ring
in multipole space**, detectable via a Rayleigh statistic on the spherical-harmonic
phases.

This code performs a **preliminary test** of that prediction using Planck SMICA
temperature data and a simple Gaussian Cℓ-matched null model.

---

## What the script does

The script `gbh_rayleigh_planck.py`:

1. **Downloads** the Planck SMICA CMB temperature map  
   (`COM_CMB_IQU-smica_2048_R3.00_full.fits`) from the Planck Legacy Archive.

2. **Computes** the angular power spectrum Cℓ from SMICA using `healpy.anafast`.

3. **Generates** spherical harmonic coefficients aℓm up to a configurable `LMAX`
   from the SMICA map.

4. For several **ℓ-bands** (e.g., [20,40], [40,80], …), it:
   - extracts the phases φℓm for m > 0 (independent modes),
   - computes the **Rayleigh resultant**  
     \( R = \left| \frac{1}{M} \sum e^{i\phi} \right| \),
   - computes the **test statistic** \( T = 2 M R^2 \).

5. **Monte Carlo null**:  
   It then generates `N_MC` Gaussian CMB simulations using the same Cℓ, computes
   the Rayleigh T statistic for each simulation and band, and builds a null
   distribution \(T_{\text{null}}\).

6. **Reports** per ℓ-band:
   - `T_real` for Planck SMICA  
   - mean ± std of the null T  
   - the fraction of null realizations exceeding the 3σ Rayleigh threshold
     (T ≥ 11.8 for χ²₂)  
   - an empirical p-value: P(T_null ≥ T_real)

This is a **conceptual check**, not a full CMB likelihood analysis. It ignores
galactic masks, beam and noise modeling, and multiple-testing corrections.
It is sufficient to demonstrate that the GBH prediction is testable and to show
that, in this preliminary implementation, Planck SMICA is consistent with the
Gaussian null model.

---

## Requirements

- Python 3.9+ (tested with 3.12)
- [`healpy`](https://healpy.readthedocs.io/)
- `numpy`

You can install them via:

```bash
pip install healpy numpy
