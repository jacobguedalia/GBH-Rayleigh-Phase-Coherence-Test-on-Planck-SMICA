# ============================================
# 0. INSTALL DEPENDENCIES
# ============================================

# !pip install healpy
# !apt-get install -y libcfitsio-dev

import os
import urllib.request
import numpy as np
import healpy as hp

# ============================================
# 1. DOWNLOAD PLANCK SMICA MAP (T map)
# ============================================

PLANCK_URL = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_3/"
    "all-sky-maps/maps/component-maps/cmb/"
    "COM_CMB_IQU-smica_2048_R3.00_full.fits"
)

PLANCK_FILE = "COM_CMB_IQU-smica_2048_R3.00_full.fits"

def download_planck():
    if not os.path.exists(PLANCK_FILE):
        print("Downloading Planck SMICA (≈250 MB)...")
        urllib.request.urlretrieve(PLANCK_URL, PLANCK_FILE)
        print("Download complete.")
    else:
        print("Planck file already exists.")

download_planck()

# ============================================
# 2. LOAD MAP & DEFINE PARAMETERS
# ============================================

print("Loading Planck T map...")
smica_T = hp.read_map(PLANCK_FILE, field=0)  # Temperature map
nside = hp.get_nside(smica_T)
print(f"NSIDE = {nside}")

# lmax for analysis
LMAX = 1000

# Define multipole bands to test
ELL_BANDS = [
    (20, 40),
    (40, 80),
    (80, 120),
    (120, 200),
    (200, 400),
    (400, 700),
]

# Number of simulations
N_MC = 100

# ============================================
# 3. COMPUTE C_ℓ FROM SMICA (NULL MODEL)
# ============================================

print("Computing C_ell from SMICA map...")
cl_smica = hp.anafast(smica_T, lmax=LMAX)

# ============================================
# 4. RAYLEIGH STATISTIC IN a_lm SPACE
# ============================================

def phases_from_alm(alm, lmin, lmax, lmax_global):
    phases = []
    for ell in range(lmin, lmax+1):
        for m in range(1, ell+1):
            idx = hp.Alm.getidx(lmax_global, ell, m)
            a = alm[idx]
            if a != 0j:
                phases.append(np.angle(a))
    return np.array(phases)

def rayleigh_T_from_phases(phases):
    M = len(phases)
    if M == 0:
        return np.nan, np.nan
    R = np.abs(np.mean(np.exp(1j * phases)))
    T = 2.0 * M * R**2
    return R, T

# ============================================
# 5. PLANCK (REAL) STATISTIC PER ℓ-BAND
# ============================================

print("Computing alm for SMICA up to LMAX...")
alm_smica = hp.map2alm(smica_T, lmax=LMAX, iter=0)

band_results = {}

for (lmin, lmax) in ELL_BANDS:
    phases_real = phases_from_alm(alm_smica, lmin, lmax, LMAX)
    R_real, T_real = rayleigh_T_from_phases(phases_real)
    band_results[(lmin,lmax)] = {
        'T_real': T_real,
        'phases_real_count': len(phases_real),
        'T_null_samples': []
    }

print("\nRayleigh T for Planck SMICA (real data):")
for band, res in band_results.items():
    lmin, lmax = band
    print(f"  ℓ ∈ [{lmin:4d}, {lmax:4d}]: T_real = {res['T_real']:.3f} (M={res['phases_real_count']})")

# ============================================
# 6. MONTE CARLO NULL REALIZATIONS
# ============================================

print(f"\nGenerating {N_MC} Monte Carlo null realizations...")
rng = np.random.default_rng(12345)

for i_mc in range(N_MC):
    if (i_mc+1) % 10 == 0:
        print(f"  MC {i_mc+1}/{N_MC}...")

    mc_map = hp.synfast(cl_smica, nside, lmax=LMAX, new=True, verbose=False)
    alm_mc = hp.map2alm(mc_map, lmax=LMAX, iter=0)

    for band in ELL_BANDS:
        lmin, lmax = band
        phases_null = phases_from_alm(alm_mc, lmin, lmax, LMAX)
        R_null, T_null = rayleigh_T_from_phases(phases_null)
        band_results[band]['T_null_samples'].append(T_null)

# ============================================
# 7. SUMMARY STATISTICS + P-VALUES
# ============================================

THRESH_3SIG = 11.8  # chi^2_2 3σ threshold

print("\n========== GBH PHASE-COHERENCE TEST (SMICA vs GAUSSIAN NULL) ==========")
print(f"LMAX = {LMAX}, N_MC = {N_MC}")
print(f"3σ threshold (chi^2_2): T >= {THRESH_3SIG:.1f}\n")

for band, res in band_results.items():
    lmin, lmax = band
    T_real = res['T_real']
    T_null = np.array(res['T_null_samples'])

    mean_null = np.nanmean(T_null)
    std_null  = np.nanstd(T_null)
    frac_null_3sig = np.mean(T_null >= THRESH_3SIG)

    p_emp = np.mean(T_null >= T_real) if np.isfinite(T_real) else np.nan

    print(f"ℓ ∈ [{lmin:4d}, {lmax:4d}]")
    print(f"  T_real                 = {T_real:.3f}")
    print(f"  Null mean (±σ)         = {mean_null:.3f} ± {std_null:.3f}")
    print(f"  frac_null(T>=11.8)     = {frac_null_3sig:.3f}")
    print(f"  Empirical p(T>=T_real) = {p_emp:.3f}")
    if T_real >= THRESH_3SIG:
        print("  -> Planck band exceeds nominal 3σ Rayleigh threshold.")
    print()
