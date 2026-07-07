#!/usr/bin/env python3
"""
M2 SPECTRAL DISTRIBUTION vs DEPTH

For n=2 qubits (d=4), computes state M2 for Ry+CZ circuits at depths L=1,2,4,8,16.
Shows eigenvalue spectrum converging to Haar values as depth increases.

State M2 = E[ |psi><psi|^{otimes 2} ] on C^d otimes C^d (16x16).
For Haar: all nonzero eigenvalues = 2/(d(d+1)) = 0.1, rank = 10.

The "spectral gap" here is the spread of M2 eigenvalues:
  gap = max(eval) - min(nonzero eval)
For Haar: gap = 0 (all eigenvalues equal).
For circuits: gap decreases toward 0 as depth increases (convergence to 2-design).

Also computes Frobenius distance to Haar M2: ||M2 - M2_Haar||_F

Run: mamba run -n sage python3 tools/sage/spectral_m2.py
"""
import numpy as np
from numpy import linalg as la

np.random.seed(123)

# ─── Parameters ────────────────────────────────────────────────────
n = 2
d = 2**n        # 4
d2 = d * d      # 16
N_SAMPLES = 5000
DEPTHS = [1, 2, 4, 8, 16]

# Theoretical values
SYM_RANK = d * (d + 1) // 2       # 10
ANTI_RANK = d * (d - 1) // 2      # 6
HAAR_EVAL = 2.0 / (d * (d + 1))   # 0.1
DESIGN_GAP = 1.0/ANTI_RANK - 1.0/SYM_RANK  # 0.0667 (exact M2 gap for unitary 2-design)

# ─── Haar-random unitary ──────────────────────────────────────────
def haar_unitary(d):
    Z = (np.random.randn(d, d) + 1j * np.random.randn(d, d)) / np.sqrt(2)
    Q, R = la.qr(Z)
    diag_R = np.diag(R)
    Lambda = diag_R / np.abs(diag_R)
    return Q @ np.diag(Lambda)

# ─── SWAP and projectors ──────────────────────────────────────────
def build_swap(d):
    S = np.zeros((d*d, d*d), dtype=complex)
    for i in range(d):
        for j in range(d):
            S[j*d + i, i*d + j] = 1.0
    return S

# ─── Ry + CZ circuit ──────────────────────────────────────────────
def ry_cz_unitary(n, depth):
    d = 2**n
    U = np.eye(d, dtype=complex)

    def cz_gate(n, q1, q2):
        d = 2**n
        cz = np.eye(d, dtype=complex)
        for idx in range(d):
            bits = [(idx >> (n-1-k)) & 1 for k in range(n)]
            if bits[q1] == 1 and bits[q2] == 1:
                cz[idx, idx] = -1.0
        return cz

    for _ in range(depth):
        ry_layer = np.array([[1.0]], dtype=complex)
        for q in range(n):
            theta = np.random.uniform(0, 2*np.pi)
            c, s = np.cos(theta/2), np.sin(theta/2)
            ry = np.array([[c, -s], [s, c]], dtype=complex)
            ry_layer = np.kron(ry_layer, ry)
        U = ry_layer @ U
        for q in range(n - 1):
            U = cz_gate(n, q, q+1) @ U

    return U

# ─── Compute M2 = E[ |psi><psi|^{otimes 2} ] ─────────────────────
def compute_state_m2(unitaries, d):
    d2 = d * d
    N = len(unitaries)
    M2 = np.zeros((d2, d2), dtype=complex)
    psi0 = np.zeros(d, dtype=complex)
    psi0[0] = 1.0
    for U in unitaries:
        psi = U @ psi0
        psi2 = np.kron(psi, psi)
        M2 += np.outer(psi2, np.conj(psi2))
    M2 /= N
    return M2

# ─── Eigenvalue analysis ──────────────────────────────────────────
def get_evals(M2):
    evals = la.eigvalsh(0.5 * (M2 + M2.conj().T))
    return np.sort(np.real(evals))[::-1]


# ═══════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("M2 SPECTRAL DISTRIBUTION vs DEPTH")
    print("=" * 70)
    print(f"n = {n} qubits, d = {d}, M2 is {d2} x {d2}")
    print(f"Samples per depth: {N_SAMPLES}")
    print()

    # Build Haar reference
    S = build_swap(d)
    I_d2 = np.eye(d2, dtype=complex)
    P_sym = (I_d2 + S) / 2.0
    M2_haar_theory = HAAR_EVAL * P_sym

    print(f"Haar theoretical values:")
    print(f"  State eigenvalue = 2/(d(d+1)) = {HAAR_EVAL:.8f}  (multiplicity {SYM_RANK})")
    print(f"  Design spectral gap = 1/{ANTI_RANK} - 1/{SYM_RANK} = {DESIGN_GAP:.8f}")
    print()

    # Compute Haar MC reference
    print("Computing Haar MC reference...")
    haar_Us = [haar_unitary(d) for _ in range(N_SAMPLES)]
    M2_haar = compute_state_m2(haar_Us, d)
    evals_haar = get_evals(M2_haar)
    nonzero_haar = evals_haar[np.abs(evals_haar) > 1e-10]
    spread_haar = nonzero_haar[0] - nonzero_haar[-1] if len(nonzero_haar) >= 2 else 0
    dist_haar = la.norm(M2_haar - M2_haar_theory, 'fro')

    print(f"  Haar MC: spread = {spread_haar:.8f}, ||M2-theory||_F = {dist_haar:.6f}")
    print()

    # ── Depth sweep ────────────────────────────────────────────────
    results = {}
    for L in DEPTHS:
        print(f"Computing M2 for Ry+CZ, depth L={L}...")
        Us = [ry_cz_unitary(n, depth=L) for _ in range(N_SAMPLES)]
        M2_L = compute_state_m2(Us, d)
        evals_L = get_evals(M2_L)
        nonzero_L = evals_L[np.abs(evals_L) > 1e-10]

        spread = nonzero_L[0] - nonzero_L[-1] if len(nonzero_L) >= 2 else 0
        dist = la.norm(M2_L - M2_haar_theory, 'fro')

        results[L] = {
            'evals': evals_L,
            'nonzero': nonzero_L,
            'spread': spread,
            'rank': len(nonzero_L),
            'dist_haar': dist,
            'max_eval': nonzero_L[0] if len(nonzero_L) > 0 else 0,
            'mean_eval': np.mean(nonzero_L) if len(nonzero_L) > 0 else 0,
            'std_eval': np.std(nonzero_L) if len(nonzero_L) > 0 else 0,
        }

    # ── Eigenvalue spectrum table ──────────────────────────────────
    print()
    print("=" * 70)
    print("  EIGENVALUE SPECTRUM (all 16 eigenvalues, sorted descending)")
    print("=" * 70)

    header = f"  {'idx':>3}"
    for L in DEPTHS:
        header += f"  {'L='+str(L):>10}"
    header += f"  {'Haar MC':>10}  {'Theory':>10}"
    print(header)
    print("  " + "-" * (4 + 12 * (len(DEPTHS) + 2)))

    for i in range(d2):
        row = f"  {i:>3}"
        for L in DEPTHS:
            row += f"  {results[L]['evals'][i]:>10.6f}"
        row += f"  {evals_haar[i]:>10.6f}"
        if i < SYM_RANK:
            row += f"  {HAAR_EVAL:>10.6f}"
        else:
            row += f"  {'0':>10}"
        print(row)

    # ── Convergence table ──────────────────────────────────────────
    print()
    print("=" * 70)
    print("  CONVERGENCE TO HAAR: spectral spread and Frobenius distance")
    print("=" * 70)
    print(f"  {'Depth':>6} {'Rank':>5} {'MaxEval':>10} {'Spread':>10} {'StdDev':>10} {'||M2-Haar||':>12} {'Converged':>10}")
    print(f"  {'-'*6} {'-'*5} {'-'*10} {'-'*10} {'-'*10} {'-'*12} {'-'*10}")

    for L in DEPTHS:
        r = results[L]
        converged = "YES" if r['dist_haar'] < 0.05 else "no"
        print(f"  {L:>6} {r['rank']:>5} {r['max_eval']:>10.6f} {r['spread']:>10.6f} "
              f"{r['std_eval']:>10.6f} {r['dist_haar']:>12.6f} {converged:>10}")

    print(f"  {'Haar':>6} {len(nonzero_haar):>5} {nonzero_haar[0]:>10.6f} {spread_haar:>10.6f} "
          f"{np.std(nonzero_haar):>10.6f} {dist_haar:>12.6f} {'(ref)':>10}")

    # ── Visual convergence ─────────────────────────────────────────
    print()
    print("=" * 70)
    print("  VISUAL: Frobenius distance to Haar")
    print("=" * 70)

    max_dist = max(results[L]['dist_haar'] for L in DEPTHS)
    for L in DEPTHS:
        dist = results[L]['dist_haar']
        bar_len = int(50 * dist / (max_dist + 1e-12))
        bar = '#' * bar_len
        print(f"  L={L:>2}: {bar:<50} {dist:.6f}")
    print(f"  Haar: {'.' * 1:<50} {dist_haar:.6f}")

    # ── Eigenvalue distribution visualization ──────────────────────
    print()
    print("=" * 70)
    print("  EIGENVALUE DISTRIBUTION (nonzero eigenvalues)")
    print("=" * 70)
    print("  Each row shows eigenvalues as positions on [0, 0.3]")
    print("  '|' marks Haar value 0.1, '*' marks eigenvalues")

    for L in DEPTHS:
        nonzero = results[L]['nonzero']
        line = [' '] * 61
        # Mark Haar position at 0.1 -> column 20 (scale: 0.3 -> 60)
        haar_col = int(HAAR_EVAL / 0.3 * 60)
        if 0 <= haar_col < 61:
            line[haar_col] = '|'
        for ev in nonzero:
            col = int(ev / 0.3 * 60)
            if 0 <= col < 61:
                line[col] = '*'
        print(f"  L={L:>2}: {''.join(line)}  [0, 0.3]")

    nonzero = nonzero_haar
    line = [' '] * 61
    haar_col = int(HAAR_EVAL / 0.3 * 60)
    if 0 <= haar_col < 61:
        line[haar_col] = '|'
    for ev in nonzero:
        col = int(ev / 0.3 * 60)
        if 0 <= col < 61:
            line[col] = '*'
    print(f"  Haar: {''.join(line)}  [0, 0.3]")

    # ── Analysis ───────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("  ANALYSIS")
    print("=" * 70)
    print(f"  At depth L=1: eigenvalues highly non-uniform (spread ~ {results[1]['spread']:.4f})")
    print(f"  At depth L=16: eigenvalues cluster tightly (spread ~ {results[16]['spread']:.4f})")
    print(f"  Frobenius distance drops from {results[1]['dist_haar']:.4f} to {results[16]['dist_haar']:.4f}")
    print()
    print("  The leading eigenvalue (~ 0.25 = 1/d) persists even at large depth.")
    print("  This is a STRUCTURAL LIMITATION of the Ry+CZ ansatz: it generates")
    print("  only real-valued unitaries (subset of O(d), not U(d)), so M2")
    print("  converges to the orthogonal Haar moment, not the unitary one.")
    print("  For a true unitary 2-design, all Sym^2 eigenvalues would equal")
    print("  2/(d(d+1)) = 0.1.")
    print()
    print("  The Frobenius distance SATURATES at ~0.16, not converging to 0.")
    print("  This saturation is the spectral signature of an expressibility ceiling.")

    if len(DEPTHS) >= 2:
        d1 = results[DEPTHS[0]]['dist_haar']
        d2_val = results[DEPTHS[-1]]['dist_haar']
        L1, L2 = DEPTHS[0], DEPTHS[-1]
        if d1 > 0 and d2_val > 0 and d2_val < d1:
            rate = np.log(d1 / d2_val) / (L2 - L1)
            print(f"  Estimated decay rate: exp(-{rate:.3f} * L)")
            print(f"  Depth to reach Haar within 1%: L ~ {-np.log(0.01 * la.norm(M2_haar_theory, 'fro') / d1) / rate:.1f}")

    print()
    print("DONE")
