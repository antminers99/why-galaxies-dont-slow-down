#!/usr/bin/env python3
"""
Replication script for B'' (B-double-prime) galaxy rotation curve model.
Reproduces all key results from the N=45 quality subsample analysis.

Requirements: Python 3.7+, numpy (no other dependencies)
Usage: python replicate_from_scratch.py

Input:  N45_final_dataset.csv (must be in same directory)
Output: Prints all coefficients, metrics, and verification checks to stdout
"""

import csv
import os
import sys
import numpy as np
from pathlib import Path


def load_data(filepath):
    """Load the N=45 dataset from CSV."""
    data = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append(row)
    print(f"Loaded {len(data)} galaxies from {filepath}")
    return data


def ols_fit(Y, X):
    """
    Ordinary Least Squares regression.
    Y: (n,) target vector
    X: (n, p) predictor matrix (WITHOUT intercept column)
    Returns: beta (p+1,), residuals (n,), predictions (n,), stats dict
    """
    n = len(Y)
    Xa = np.column_stack([np.ones(n), X])
    p = Xa.shape[1]

    beta = np.linalg.lstsq(Xa, Y, rcond=None)[0]
    predictions = Xa @ beta
    residuals = Y - predictions

    rss = np.sum(residuals**2)
    tss = np.sum((Y - np.mean(Y))**2)
    r2 = 1 - rss / tss
    adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p)
    se = np.sqrt(rss / (n - p))

    XtX_inv = np.linalg.inv(Xa.T @ Xa)
    se_beta = np.sqrt(np.diag(XtX_inv) * rss / (n - p))
    t_stats = beta / se_beta

    return {
        'beta': beta,
        'se_beta': se_beta,
        't_stats': t_stats,
        'residuals': residuals,
        'predictions': predictions,
        'r2': r2,
        'adj_r2': adj_r2,
        'se': se,
        'rss': rss,
        'tss': tss,
        'n': n,
        'k': p
    }


def loo_cv(Y, X):
    """
    Leave-one-out cross-validation.
    Returns: LOO RMS, LOO gap%, per-galaxy prediction errors
    """
    n = len(Y)
    Xa = np.column_stack([np.ones(n), X])
    errors = np.zeros(n)

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        Xa_train = Xa[mask]
        Y_train = Y[mask]

        beta_loo = np.linalg.lstsq(Xa_train, Y_train, rcond=None)[0]
        pred_i = Xa[i] @ beta_loo
        errors[i] = Y[i] - pred_i

    loo_rms = np.sqrt(np.mean(errors**2))
    sd_y = np.std(Y, ddof=1)
    gap_pct = 100 * (1 - loo_rms**2 / sd_y**2)

    return loo_rms, gap_pct, errors


def main():
    script_dir = Path(__file__).parent
    csv_path = script_dir / 'N45_final_dataset.csv'

    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found. Place it in the same directory.")
        sys.exit(1)

    data = load_data(csv_path)
    n = len(data)

    logA0 = np.array([float(d['logA0']) for d in data])
    logMHI = np.array([float(d['logMHI']) for d in data])
    rcWig = np.array([float(d['rcWiggliness']) for d in data])
    logMhost = np.array([float(d['logMhost']) for d in data])
    logSigma0 = np.array([float(d['logSigma0']) for d in data])
    logMeanRun = np.array([float(d['logMeanRun']) for d in data])
    logUps = np.array([float(d['logUpsilon_disk']) for d in data])
    morphT = np.array([float(d['morphT']) for d in data])
    names = [d['galaxy_name'] for d in data]

    print()
    print("=" * 70)
    print("  REPLICATION FROM SCRATCH")
    print("  B'' Model — N=45 Quality Subsample")
    print("=" * 70)

    print()
    print("-" * 70)
    print("  STEP 1: Construct Upsilon_perp")
    print("-" * 70)

    X_conf = np.column_stack([logMHI, logSigma0, morphT])
    ups_fit = ols_fit(logUps, X_conf)
    ups_perp = ups_fit['residuals']

    print(f"  Auxiliary regression: logUps ~ logMHI + logSigma0 + morphT")
    print(f"  R^2 = {ups_fit['r2']:.4f} (confounders explain {ups_fit['r2']*100:.1f}% of Upsilon)")
    print(f"  Upsilon_perp: mean={np.mean(ups_perp):.6f}, SD={np.std(ups_perp, ddof=1):.4f}")

    print()
    print("-" * 70)
    print("  STEP 2: Fit M0 (null model)")
    print("-" * 70)

    sd_y = np.std(logA0, ddof=1)
    print(f"  mean(logA0) = {np.mean(logA0):.4f}")
    print(f"  SD(logA0) = {sd_y:.4f} dex")
    loo_m0_rms = sd_y * np.sqrt((n - 1) / n)
    print(f"  M0 tau = {sd_y:.4f} dex")

    print()
    print("-" * 70)
    print("  STEP 3: Fit A' (2-variable model)")
    print("-" * 70)

    X_Ap = np.column_stack([logMHI, logMhost])
    fit_Ap = ols_fit(logA0, X_Ap)
    loo_Ap_rms, loo_Ap_gap, _ = loo_cv(logA0, X_Ap)

    print(f"  logA0 ~ logMHI + logMhost")
    print(f"  Coefficients:")
    for j, name in enumerate(['intercept', 'logMHI', 'logMhost']):
        print(f"    {name:>15s} = {fit_Ap['beta'][j]:+.4f}  (t={fit_Ap['t_stats'][j]:+.2f})")
    print(f"  R^2 = {fit_Ap['r2']:.4f}, Adj R^2 = {fit_Ap['adj_r2']:.4f}")
    print(f"  Residual SD = {np.std(fit_Ap['residuals'], ddof=1):.4f} dex")
    print(f"  LOO RMS = {loo_Ap_rms:.4f}, LOO gap% = {loo_Ap_gap:.1f}%")

    print()
    print("-" * 70)
    print("  STEP 4: Fit B' (4-variable model)")
    print("-" * 70)

    X_Bp = np.column_stack([logMHI, rcWig, logMhost, logSigma0])
    fit_Bp = ols_fit(logA0, X_Bp)
    loo_Bp_rms, loo_Bp_gap, _ = loo_cv(logA0, X_Bp)

    print(f"  logA0 ~ logMHI + rcWig + logMhost + logSigma0")
    print(f"  Coefficients:")
    for j, name in enumerate(['intercept', 'logMHI', 'rcWig', 'logMhost', 'logSigma0']):
        print(f"    {name:>15s} = {fit_Bp['beta'][j]:+.4f}  (t={fit_Bp['t_stats'][j]:+.2f})")
    print(f"  R^2 = {fit_Bp['r2']:.4f}, Adj R^2 = {fit_Bp['adj_r2']:.4f}")
    print(f"  Residual SD = {np.std(fit_Bp['residuals'], ddof=1):.4f} dex")
    print(f"  LOO RMS = {loo_Bp_rms:.4f}, LOO gap% = {loo_Bp_gap:.1f}%")

    print()
    print("-" * 70)
    print("  STEP 5: Fit B'' (6-variable FROZEN BASELINE)")
    print("-" * 70)

    X_Bdp = np.column_stack([logMHI, rcWig, logMhost, logSigma0, logMeanRun, ups_perp])
    fit_Bdp = ols_fit(logA0, X_Bdp)
    loo_Bdp_rms, loo_Bdp_gap, loo_errors = loo_cv(logA0, X_Bdp)

    var_names = ['intercept', 'logMHI', 'rcWig', 'logMhost', 'logSigma0', 'logMeanRun', 'Ups_perp']
    expected_signs = [None, -1, +1, -1, +1, +1, +1]

    print(f"  logA0 ~ logMHI + rcWig + logMhost + logSigma0 + logMeanRun + Ups_perp")
    print()
    print(f"  {'Variable':>15s}  {'Coeff':>8s}  {'SE':>8s}  {'t':>8s}  {'Sign OK?':>8s}")
    print(f"  {'-'*15}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    sign_checks = []
    for j, name in enumerate(var_names):
        sign_ok = ''
        if expected_signs[j] is not None:
            ok = np.sign(fit_Bdp['beta'][j]) == expected_signs[j]
            sign_ok = 'YES' if ok else 'NO !!!'
            sign_checks.append(ok)
        print(f"  {name:>15s}  {fit_Bdp['beta'][j]:+8.4f}  {fit_Bdp['se_beta'][j]:8.4f}  {fit_Bdp['t_stats'][j]:+8.3f}  {sign_ok:>8s}")

    print()
    print(f"  In-sample R^2 = {fit_Bdp['r2']:.4f}")
    print(f"  Adj R^2 = {fit_Bdp['adj_r2']:.4f}")
    print(f"  Residual SD = {np.std(fit_Bdp['residuals'], ddof=1):.4f} dex")
    print(f"  LOO RMS = {loo_Bdp_rms:.4f}")
    print(f"  LOO gap% = {loo_Bdp_gap:.1f}%")

    print()
    print("-" * 70)
    print("  STEP 6: Residual diagnostics")
    print("-" * 70)

    resid = fit_Bdp['residuals']
    print(f"  Residual mean = {np.mean(resid):.6f} (should be ~0)")
    print(f"  Residual SD = {np.std(resid, ddof=1):.4f} dex")

    skew = np.mean(((resid - np.mean(resid)) / np.std(resid, ddof=1))**3) * n / ((n-1)*(n-2))
    kurt = np.mean(((resid - np.mean(resid)) / np.std(resid, ddof=1))**4) - 3
    print(f"  Skewness = {skew:.4f}")
    print(f"  Excess kurtosis = {kurt:.4f}")

    signs = resid >= 0
    runs = 1 + np.sum(signs[1:] != signs[:-1])
    n_pos = np.sum(signs)
    n_neg = n - n_pos
    exp_runs = (2 * n_pos * n_neg) / n + 1
    var_runs = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n**2 * (n - 1))
    runs_z = (runs - exp_runs) / np.sqrt(var_runs)
    print(f"  Runs test: runs={runs}, expected={exp_runs:.1f}, z={runs_z:.3f}")

    print()
    print("-" * 70)
    print("  STEP 7: Tau progression")
    print("-" * 70)
    print(f"  M0:  tau = {sd_y:.4f} dex")
    print(f"  A':  tau = {np.std(fit_Ap['residuals'], ddof=1):.4f} dex")
    print(f"  B':  tau = {np.std(fit_Bp['residuals'], ddof=1):.4f} dex")
    print(f"  B'': tau = {np.std(fit_Bdp['residuals'], ddof=1):.4f} dex")

    print()
    print("=" * 70)
    print("  VERIFICATION SUMMARY")
    print("=" * 70)
    print()

    all_signs_ok = all(sign_checks)
    loo_in_range = 45.0 <= loo_Bdp_gap <= 55.0
    r2_in_range = 0.65 <= fit_Bdp['r2'] <= 0.80

    print(f"  All 6 coefficient signs correct?  {'PASS' if all_signs_ok else 'FAIL'}  ({sum(sign_checks)}/6)")
    print(f"  LOO gap% in [45,55]?              {'PASS' if loo_in_range else 'FAIL'}  ({loo_Bdp_gap:.1f}%)")
    print(f"  R^2 in [0.65,0.80]?               {'PASS' if r2_in_range else 'FAIL'}  ({fit_Bdp['r2']:.4f})")
    print()

    if all_signs_ok and loo_in_range and r2_in_range:
        print("  ===================================================")
        print("  REPLICATION CONFIRMED")
        print("  All verification checks passed.")
        print("  ===================================================")
    else:
        print("  ===================================================")
        print("  REPLICATION ISSUES DETECTED")
        print("  Check failed items above.")
        print("  ===================================================")

    print()
    print("-" * 70)
    print("  Per-galaxy predictions and residuals")
    print("-" * 70)
    print(f"  {'Galaxy':>15s}  {'logA0_obs':>10s}  {'logA0_pred':>10s}  {'resid':>8s}  {'LOO_err':>8s}")
    for i in range(n):
        print(f"  {names[i]:>15s}  {logA0[i]:10.4f}  {fit_Bdp['predictions'][i]:10.4f}  {resid[i]:+8.4f}  {loo_errors[i]:+8.4f}")


if __name__ == '__main__':
    main()
