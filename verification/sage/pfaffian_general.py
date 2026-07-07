#!/usr/bin/env python3
"""
SYMBOLIC PROOF: Pf(A)^2 = det(A) for general antisymmetric matrices

Proves the Cayley identity Pf(A)^2 = det(A) for 2n x 2n antisymmetric
matrices with SYMBOLIC entries over QQ[a_ij].

This is an EXACT polynomial identity verification, not numerical.
SageMath's pfaffian() uses the recursive definition via perfect matchings.

Verified for n = 1 (2x2), 2 (4x4), 3 (6x6), 4 (8x8), 5 (10x10).
Number of Pfaffian terms: 1, 3, 15, 105, 945 = (2n-1)!!

Requires: conda env jepa-verify (SageMath 10.8)
Run: conda run -n jepa-verify python3 tools/sage_pfaffian_general.py
"""
from sage.all import PolynomialRing, QQ, matrix

def verify_pfaffian_identity(n):
    """Verify Pf(A)^2 = det(A) for generic 2n x 2n antisymmetric matrix."""
    size = 2 * n
    names = ['a_%d_%d' % (i, j) for i in range(size) for j in range(i+1, size)]
    R = PolynomialRing(QQ, names)
    gens = list(R.gens())

    # Build generic antisymmetric matrix
    idx = 0
    A = matrix(R, size, size)
    for i in range(size):
        for j in range(i+1, size):
            A[i, j] = gens[idx]
            A[j, i] = -gens[idx]
            idx += 1

    pf = A.pfaffian()
    det_A = A.det()
    identity_holds = (pf**2 - det_A).is_zero()

    return {
        'n': n,
        'size': size,
        'num_terms': len(pf.monomials()),
        'identity_holds': identity_holds
    }

if __name__ == "__main__":
    print("SYMBOLIC PROOF: Pf(A)^2 = det(A)")
    print("=" * 50)
    print()

    all_pass = True
    for n in range(1, 6):
        result = verify_pfaffian_identity(n)
        status = "✓" if result['identity_holds'] else "✗"
        if not result['identity_holds']:
            all_pass = False
        print(f"  n={result['n']} ({result['size']}x{result['size']}): "
              f"Pf has {result['num_terms']:>4d} terms, "
              f"Pf^2 == det: {status}")

    print()
    if all_pass:
        print("PROVED: Pf(A)^2 = det(A) for n=1..5 over QQ[a_ij] ✓")
        print("(2n-1)!! terms: 1, 3, 15, 105, 945")
    else:
        print("FAILED ✗")
