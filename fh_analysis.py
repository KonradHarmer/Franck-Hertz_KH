"""Analysis tools for Franck-Hertz anode-current sweeps (mercury tube).

Every curve is handled as a 2 x N array: row 0 = accelerating voltage U2 (V),
row 1 = anode signal (V at the amplifier output; 1 V corresponds to 1 nA).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, medfilt, savgol_filter

DATA = Path(__file__).parent / "data"

# Mercury levels that can be excited from the 6^1S0 ground state (eV, nm)
HG_LINES = {"6³P₀": (4.67, 265.6), "6³P₁": (4.89, 253.7)}

U2_DIVIDER = 10  # the interface records U2 through a 1:10 divider


# --------------------------------------------------------------------------- I/O

def load_fh(path):
    """Load a sweep file. Columns: 0 = not used, 1 = U2 / 10 (V), 2 = anode signal (V)."""
    data = np.loadtxt(path, usecols=(1, 2)).T
    data[0] *= U2_DIVIDER
    return data


def _step(U2):
    """Mean U2 increment per sample (V)."""
    return (U2.max() - U2.min()) / (len(U2) - 1)


# ------------------------------------------------------------ signal processing

def despike(data, window_V=1.0, thresh=0.1):
    """Remove short glitches from the anode signal.

    Points deviating by more than `thresh` (V) from a running median of width
    `window_V` (in U2) are replaced by linear interpolation; the flanks of each
    spike are masked too.
    """
    U2, UA = data
    n = int(window_V / _step(U2)) | 1                    # odd window length in samples
    bad = np.abs(UA - medfilt(UA, n)) > thresh
    bad = binary_dilation(bad, iterations=n // 10)
    i = np.arange(len(UA))
    UA_clean = UA.copy()
    UA_clean[bad] = np.interp(i[bad], i[~bad], UA[~bad])
    return np.array([U2, UA_clean])


def find_extrema(data, smooth_V=0.5, min_sep_V=3.0, prominence=0.02, edge_V=1.0):
    """Locate the maxima and minima of a Franck-Hertz curve.

    The curve is despiked, then smoothed along the sample index with a
    Savitzky-Golay filter (U2 itself is noisy, so it is smoothed as well).
    Extrema closer than `edge_V` to the end of the sweep are discarded.

    Returns (maxima, minima, smoothed); each is a 2 x K array (U2, signal).
    """
    U2, UA = despike(data)
    dU = _step(U2)
    win = int(smooth_V / dU) | 1
    U2_s = savgol_filter(U2, win, 1)
    UA_s = savgol_filter(UA, win, 3)

    dist = int(min_sep_V / dU)
    imax, _ = find_peaks(UA_s, distance=dist, prominence=prominence)
    imin, _ = find_peaks(-UA_s, distance=dist, prominence=prominence)
    imax = imax[U2_s[imax] < U2_s.max() - edge_V]
    imin = imin[U2_s[imin] < U2_s.max() - edge_V]

    return (np.array([U2_s[imax], UA_s[imax]]),
            np.array([U2_s[imin], UA_s[imin]]),
            np.array([U2_s, UA_s]))


def contrast_snr(data, U_noise=3.0, max_gap_V=3.0, E_guess=4.9):
    """Peak contrast and signal-to-noise ratio for every maximum of a sweep.

    C   = (I_max - I_min) / (I_max + I_min), pairing each maximum with the next
          minimum (if it lies within `max_gap_V`)
    SNR = I_max² / I_noise², with I_noise the std of the signal for U2 < U_noise
          (before the current starts to rise)
    The order of a maximum is estimated as round(U2 / E_guess).
    """
    maxima, minima, _ = find_extrema(data)
    U2, UA = despike(data)
    I_noise = UA[U2 < U_noise].std()

    rows = []
    for U_max, I_max in maxima.T:
        nxt = minima[:, (minima[0] > U_max) & (minima[0] < U_max + max_gap_V)]
        if nxt.shape[1] == 0:
            continue
        I_min = nxt[1, 0]
        rows.append({"order": int(round(U_max / E_guess)),
                     "C": (I_max - I_min) / (I_max + I_min),
                     "SNR": I_max**2 / I_noise**2})
    return pd.DataFrame(rows)


def compare_settings(curves, index_names):
    """Contrast/SNR summary for a grid of operating settings.

    `curves` maps a settings tuple to a sweep. Returns (C per order, SNR per
    order, summary). The summary only uses orders present in every setting, so
    all settings are compared on the same peaks.
    """
    res = pd.concat([contrast_snr(d).assign(**dict(zip(index_names, k)))
                     for k, d in curves.items()])
    C = res.pivot_table(index=index_names, columns="order", values="C")
    SNR = res.pivot_table(index=index_names, columns="order", values="SNR")
    common = C.columns[C.notna().all() & SNR.notna().all()]
    lo, hi = common.min(), common.max()
    summary = pd.DataFrame({
        f"mean C (orders {lo}–{hi})": C[common].mean(axis=1),
        f"min SNR (orders {lo}–{hi})": SNR[common].min(axis=1),
        "max signal (V)": [curves[k][1].max() for k in C.index],   # clipping check, ADC limit 5 V
    })
    return C, SNR, summary


# ------------------------------------------------------------------ the model

def rapior(n, Ea, r):
    """Spacing between extremum n-1 and n (Rapior et al. 2006).

    ΔE(n) = [1 + (2n - 1) λ/L] E_a, with r = λ/L the mean free path of the
    electrons over the accelerating distance. E_a = ΔE(n = 1/2).
    """
    return (1 + (2 * n - 1) * r) * Ea


def extremum_positions(runs, which, k=5, smooth_V=0.5):
    """U2 positions of the first k maxima (which=0) or minima (which=1) of each run.

    Runs with fewer than k extrema of that kind are dropped.
    Returns (array of shape runs x k, indices of the runs kept).
    """
    U = [find_extrema(d, smooth_V=smooth_V)[which][0, :k] for d in runs]
    keep = [i for i, u in enumerate(U) if len(u) == k]
    return np.array([U[i] for i in keep]), keep


def fit_spacings(dU, full_cov=False):
    """Fit the Rapior model to per-run spacings dU (runs x spacings).

    The data points are the mean spacings over runs. With full_cov=False their
    errors are the standard errors of the means, treated as independent; with
    full_cov=True the full covariance matrix of the means is used (adjacent
    spacings share an extremum, so they are correlated).

    Returns a dict with the mean spacings y, their errors sy, the orders n, the
    parameters p = (E_a, λ/L), their covariance, χ²/ndf, the slope m = dΔE/dn
    and its error, and the error on E_a scaled by √max(χ²/ndf, 1) (PDG).
    """
    N, k = dU.shape
    n = np.arange(2, 2 + k)          # spacing between extremum n-1 and n; first extremum = 1
    y = dU.mean(axis=0)
    V = np.cov(dU, rowvar=False, ddof=1) / N
    sy = np.sqrt(np.diag(V))
    sigma = V if full_cov else sy
    p, cov = curve_fit(rapior, n, y, p0=(4.7, 0.01), sigma=sigma, absolute_sigma=True)

    resid = y - rapior(n, *p)
    chi2 = resid @ np.linalg.solve(V, resid) if full_cov else np.sum((resid / sy) ** 2)
    chi2r = chi2 / (k - 2)
    scale = np.sqrt(max(chi2r, 1))

    Ea, r = p
    m = 2 * Ea * r
    s_m = 2 * np.sqrt(r**2 * cov[0, 0] + Ea**2 * cov[1, 1] + 2 * Ea * r * cov[0, 1])
    return dict(n=n, y=y, sy=sy, p=p, cov=cov, chi2r=chi2r,
                Ea=Ea, sEa=np.sqrt(cov[0, 0]), sEa_scaled=np.sqrt(cov[0, 0]) * scale,
                r=r, sr=np.sqrt(cov[1, 1]) * scale, m=m, sm=s_m * scale,
                pulls=resid / sy)
