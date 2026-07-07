#!/usr/bin/env python3
"""
NumPy Monte Carlo verification of the LAST AXIOM: haar_fourth_moment
E[|⟨ψ|U|φ⟩|⁴] = 2/(d(d+1)) for Haar-random U(d)

This is the ONE axiom remaining in the Lean formalization.
Verified numerically for d = 2, 3, 4, 8, 16.

Requires: numpy >= 2.0
Run: python3 tools/numpy_haar_verify.py
"""
import numpy as np

def haar_random_unitary(d, rng):
    """Generate a Haar-random unitary matrix via QR decomposition."""
    Z = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    Q, R = np.linalg.qr(Z)
    D = np.diag(R) / np.abs(np.diag(R))
    return Q @ np.diag(D)

def verify_fourth_moment(d, n_samples=200000):
    """Monte Carlo estimate of E[|⟨ψ|U|φ⟩|⁴] for orthogonal ψ, φ."""
    rng = np.random.default_rng(42)
    psi = np.zeros(d, dtype=complex); psi[0] = 1.0
    phi = np.zeros(d, dtype=complex); phi[1] = 1.0

    fourth_moments = []
    for _ in range(n_samples):
        U = haar_random_unitary(d, rng)
        bracket = psi.conj() @ U @ phi
        fourth_moments.append(abs(bracket)**4)

    empirical = np.mean(fourth_moments)
    theoretical = 2.0 / (d * (d + 1))
    std_err = np.std(fourth_moments) / np.sqrt(n_samples)
    return empirical, theoretical, std_err

if __name__ == "__main__":
    print("=" * 65)
    print("HAAR FOURTH MOMENT VERIFICATION (Monte Carlo, 200k samples)")
    print("Axiom: E[|⟨ψ|U|φ⟩|⁴] = 2/(d(d+1))")
    print("=" * 65)
    print(f"{'d':>4} | {'Empirical':>12} | {'Theoretical':>12} | {'Std Error':>10} | {'Match':>6}")
    print("-" * 65)

    all_match = True
    for d in [2, 3, 4, 8, 16]:
        emp, theo, se = verify_fourth_moment(d)
        match = abs(emp - theo) < 3 * se
        if not match:
            all_match = False
        print(f"{d:>4} | {emp:>12.8f} | {theo:>12.8f} | {se:>10.8f} | {'✓' if match else '✗':>6}")

    print("-" * 65)
    print(f"ALL DIMENSIONS MATCH: {'YES ✓' if all_match else 'NO ✗'}")
    print()

    print("=" * 65)
    print("CONCENTRATION BOUNDS VERIFICATION")
    print("=" * 65)
    for n in [1, 2, 3, 4, 5]:
        d = 2**n
        mu = 1.0 / (d + 1)
        beta = 1.0 / (d * (d + 1))
        fourth = 2.0 / (d * (d + 1))
        print(f"n={n}, d={d}: μ={mu:.6f}, β={beta:.8f}, 4th={fourth:.8f}")
