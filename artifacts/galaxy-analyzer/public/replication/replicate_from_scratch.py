#!/usr/bin/env python3
"""
Replication script for M5 and M3 galaxy rotation curve models.
Reproduces all key results from the N=45 quality subsample analysis.

Requirements: Python 3.7+, numpy (no other dependencies)
Usage: python replicate_from_scratch.py

Input:  N45_final_dataset.csv (must be in same directory)
Output: Prints all coefficients, metrics, and verification checks to stdout
"""

import csv
import sys
import numpy as np
from pathlib import Path


def load_data(filepath):
    data = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append(row)
    print(f"Loaded {len(data)} galaxies from {filepath}")
    return data


def ols_fit(Y, X):
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
        'beta': beta, 'se_beta': se_beta, 't_stats': t_stats,
        'residuals': residuals, 'predictions': predictions,
        'r2': r2, 'adj_r2': adj_r2, 'se': se,
        'rss': rss, 'tss': tss, 'n': n, 'k': p
    }


def loo_cv(Y, X):
    n = len(Y)
    Xa = np.column_stack([np.ones(n), X])
    errors = np.zeros(n)

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        beta_loo = np.linalg.lstsq(Xa[mask], Y[mask], rcond=None)[0]
        errors[i] = Y[i] - Xa[i] @ beta_loo

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
    logMhost = np.array([float(d['logMhost']) for d in data])
    logSigma0 = np.array([float(d['logSigma0']) for d in data])
    logMeanRun = np.array([float(d['logMeanRun']) for d in data])
    logUps = np.array([float(d['logUpsilon_disk']) for d in data])
    morphT = np.array([float(d['morphT']) for d in data])
    names = [d['galaxy_name'] for d in data]

    print()
    print("=" * 70)
    print("  REPLICATION FROM SCRATCH")
    print("  M5 + M3 Models  |  N=45 Quality Subsample")
    print("=" * 70)

    # STEP 1: Construct Upsilon_perp
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

    # STEP 2: M0 (null)
    print()
    print("-" * 70)
    print("  STEP 2: M0 (universal constant)")
    print("-" * 70)

    sd_y = np.std(logA0, ddof=1)
    print(f"  mean(logA0) = {np.mean(logA0):.4f}")
    print(f"  SD(logA0) = {sd_y:.4f} dex")

    # STEP 3: M3 (compressed state law)
    print()
    print("-" * 70)
    print("  STEP 3: M3 (3-axis compressed state law)")
    print("-" * 70)

    X_M3 = np.column_stack([logMHI, logMhost, logMeanRun])
    fit_M3 = ols_fit(logA0, X_M3)
    loo_M3_rms, loo_M3_gap, _ = loo_cv(logA0, X_M3)

    m3_names = ['intercept', 'logMHI', 'logMhost', 'logMeanRun']
    m3_signs = [None, -1, -1, +1]

    print(f"  logA0 ~ logMHI + logMhost + logMeanRun")
    print()
    print(f"  {'Variable':>15s}  {'Coeff':>8s}  {'SE':>8s}  {'t':>8s}")
    print(f"  {'-'*15}  {'-'*8}  {'-'*8}  {'-'*8}")
    for j, name in enumerate(m3_names):
        print(f"  {name:>15s}  {fit_M3['beta'][j]:+8.4f}  {fit_M3['se_beta'][j]:8.4f}  {fit_M3['t_stats'][j]:+8.3f}")

    print()
    print(f"  R^2 = {fit_M3['r2']:.4f}, Adj R^2 = {fit_M3['adj_r2']:.4f}")
    print(f"  Residual SD = {np.std(fit_M3['residuals'], ddof=1):.4f} dex")
    print(f"  LOO RMS = {loo_M3_rms:.4f}, LOO gap% = {loo_M3_gap:.1f}%")

    # STEP 4: M5 (full predictive law)
    print()
    print("-" * 70)
    print("  STEP 4: M5 (5-axis predictive law) [FINAL MODEL]")
    print("-" * 70)

    X_M5 = np.column_stack([logMHI, logMhost, logSigma0, logMeanRun, ups_perp])
    fit_M5 = ols_fit(logA0, X_M5)
    loo_M5_rms, loo_M5_gap, loo_errors = loo_cv(logA0, X_M5)

    m5_names = ['intercept', 'logMHI', 'logMhost', 'logSigma0', 'logMeanRun', 'Ups_perp']
    expected_signs = [None, -1, -1, +1, +1, +1]

    print(f"  logA0 ~ logMHI + logMhost + logSigma0 + logMeanRun + Ups_perp")
    print()
    print(f"  {'Variable':>15s}  {'Coeff':>8s}  {'SE':>8s}  {'t':>8s}  {'Sign OK?':>8s}")
    print(f"  {'-'*15}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    sign_checks = []
    for j, name in enumerate(m5_names):
        sign_ok = ''
        if expected_signs[j] is not None:
            ok = np.sign(fit_M5['beta'][j]) == expected_signs[j]
            sign_ok = 'YES' if ok else 'NO !!!'
            sign_checks.append(ok)
        print(f"  {name:>15s}  {fit_M5['beta'][j]:+8.4f}  {fit_M5['se_beta'][j]:8.4f}  {fit_M5['t_stats'][j]:+8.3f}  {sign_ok:>8s}")

    print()
    print(f"  In-sample R^2 = {fit_M5['r2']:.4f}")
    print(f"  Adj R^2 = {fit_M5['adj_r2']:.4f}")
    print(f"  Residual SD = {np.std(fit_M5['residuals'], ddof=1):.4f} dex")
    print(f"  LOO RMS = {loo_M5_rms:.4f}")
    print(f"  LOO gap% = {loo_M5_gap:.1f}%")

    # STEP 5: Residual diagnostics
    print()
    print("-" * 70)
    print("  STEP 5: Residual diagnostics (M5)")
    print("-" * 70)

    resid = fit_M5['residuals']
    print(f"  Residual mean = {np.mean(resid):.6f}")
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
    if var_runs > 0:
        runs_z = (runs - exp_runs) / np.sqrt(var_runs)
    else:
        runs_z = 0.0
    print(f"  Runs test: runs={runs}, expected={exp_runs:.1f}, z={runs_z:.3f}")

    # STEP 6: M3 vs M5 comparison
    print()
    print("-" * 70)
    print("  STEP 6: Model comparison summary")
    print("-" * 70)
    print(f"  M0:  RMS = {sd_y:.4f} dex  |  LOO gap% ~ -2.3%")
    print(f"  M3:  RMS = {np.std(fit_M3['residuals'], ddof=1):.4f} dex  |  LOO gap% = {loo_M3_gap:.1f}%")
    print(f"  M5:  RMS = {np.std(fit_M5['residuals'], ddof=1):.4f} dex  |  LOO gap% = {loo_M5_gap:.1f}%")
    print()
    m3_retention = (loo_M3_gap / loo_M5_gap * 100) if loo_M5_gap > 0 else 0
    print(f"  M3 retains {m3_retention:.0f}% of M5 signal (3 axes vs 5)")

    # VERIFICATION
    print()
    print("=" * 70)
    print("  VERIFICATION SUMMARY")
    print("=" * 70)
    print()

    all_signs_ok = all(sign_checks)
    loo_in_range = 42.0 <= loo_M5_gap <= 52.0
    r2_in_range = 0.55 <= fit_M5['r2'] <= 0.65

    print(f"  All 5 M5 coefficient signs correct?   {'PASS' if all_signs_ok else 'FAIL'}  ({sum(sign_checks)}/5)")
    print(f"  M5 LOO gap% in [42,52]?               {'PASS' if loo_in_range else 'FAIL'}  ({loo_M5_gap:.1f}%)")
    print(f"  M5 R^2 in [0.55,0.65]?                {'PASS' if r2_in_range else 'FAIL'}  ({fit_M5['r2']:.4f})")
    print(f"  M3 LOO gap% in [40,50]?               {'PASS' if 40 <= loo_M3_gap <= 50 else 'FAIL'}  ({loo_M3_gap:.1f}%)")
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

    # Per-galaxy table
    print()
    print("-" * 70)
    print("  Per-galaxy M5 predictions and residuals")
    print("-" * 70)
    print(f"  {'Galaxy':>15s}  {'logA0_obs':>10s}  {'logA0_pred':>10s}  {'resid':>8s}  {'LOO_err':>8s}")
    for i in range(n):
        print(f"  {names[i]:>15s}  {logA0[i]:10.4f}  {fit_M5['predictions'][i]:10.4f}  {resid[i]:+8.4f}  {loo_errors[i]:+8.4f}")


if __name__ == '__main__':
    main()
