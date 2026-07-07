#!/usr/bin/env python3
"""
Purity Regularization: The Quantum SIGReg

Mathematical analysis of purity-constrained quantum states.
For a d-dimensional state with Tr[ρ]=1 and Tr[ρ²]=p:

  a = (1 + √((d²-d)p - d + 1)) / d    (largest eigenvalue)
  b = (1 - a) / (d - 1)                (remaining eigenvalues)

Key results:
  p = 1/d: a = b = 1/d (maximally mixed = barren plateau)
  p = 1: a = 1, b = 0 (pure state)
  p ∈ (1/d, 1): Goldilocks (structured but not collapsed)

SWAP test measures purity in O(1/ε²) shots.
Purity regularizer: L_reg = (Tr[ρ²] - p_target)²

Requires: conda env jepa-verify
Run: conda run -n jepa-verify python3 tools/sage_purity_regularization.py
"""
import math

def optimal_spectrum(d, p):
    """Compute the max-entropy spectrum with trace=1, purity=p."""
    disc = (d*d - d) * p - d + 1
    if disc < 0:
        return None, None
    a = (1 + math.sqrt(disc)) / d
    b = (1 - a) / (d - 1) if d > 1 else 0
    return a, b

if __name__ == "__main__":
    print("Purity Regularization: Quantum SIGReg")
    print("=" * 55)
    print()
    for d in [4, 16, 256]:
        n = int(math.log2(d))
        print(f"d = {d} ({n} qubits):")
        print(f"  {'p':>8} {'a (max λ)':>10} {'b (rest)':>10} {'a/b ratio':>10}")
        print(f"  {'-'*42}")
        for p in [1/d, 0.1, 0.2, 0.5, 0.8]:
            if p >= 1/d:
                a, b = optimal_spectrum(d, p)
                if a is not None and b > 1e-10:
                    print(f"  {p:>8.4f} {a:>10.4f} {b:>10.4f} {a/b:>10.1f}x")
                elif a is not None:
                    print(f"  {p:>8.4f} {a:>10.4f} {'≈0':>10} {'∞':>10}")
        print()

    print("THEOREM: purity > 1/d + ε → spectrum non-degenerate")
    print("COMPUTABLE: SWAP test, O(1/ε²) shots")
    print("OPEN: does purity control → quantum advantage?")
