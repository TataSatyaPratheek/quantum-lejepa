"""
Fibre Bundle Verification for Q-JEPA
=====================================

Numerically verify the geometric interpretation of Q-JEPA:
  - Energy with V=I (flat connection) = ‖φ(y) - φ(x)‖²
  - Energy = 0 when V is the "perfect predictor" mapping φ(x) → φ(y)
  - Collapsed encoder → zero energy for V=I

Uses n=2 qubits, so H = C^4.
"""

import numpy as np

np.random.seed(42)

d = 4  # 2 qubits → C^4


def random_state(d):
    """Random normalized state in C^d."""
    v = np.random.randn(d) + 1j * np.random.randn(d)
    return v / np.linalg.norm(v)


def qjepa_energy(encoder, predictor, x, y):
    """E(x,y) = ‖φ(y) - V(φ(x))‖²"""
    diff = encoder(y) - predictor(encoder(x))
    return np.real(np.vdot(diff, diff))


# Fix two input points
phi_x = random_state(d)
phi_y = random_state(d)

encoder = lambda p: phi_x if p == 'x' else phi_y

print("=" * 60)
print("FIBRE BUNDLE VERIFICATION (n=2 qubits, d=4)")
print("=" * 60)

# --- Test 1: Flat connection (V=I) ---
print("\n--- Test 1: Flat connection (V = I) ---")
identity = lambda v: v
E_flat = qjepa_energy(encoder, identity, 'x', 'y')
E_jepa = np.real(np.vdot(phi_y - phi_x, phi_y - phi_x))
print(f"  E(x,y; V=I)          = {E_flat:.10f}")
print(f"  ‖φ(y) - φ(x)‖²      = {E_jepa:.10f}")
print(f"  Difference            = {abs(E_flat - E_jepa):.2e}")
assert abs(E_flat - E_jepa) < 1e-12, "FAIL: flat connection != JEPA distance"
print("  PASS: flat connection = JEPA distance")

# --- Test 2: Perfect predictor (V maps φ(x) to φ(y)) ---
print("\n--- Test 2: Perfect predictor ---")
# Build a unitary V such that V φ(x) = φ(y).
# Extend {φ(x)} to an ONB, map φ(x)→φ(y) and complete unitarily.
# Simple approach: use a Householder-like construction.

# V = I + (φ(y) - φ(x)) ⊗ ⟨φ(x)| / ⟨φ(x)|φ(x)⟩  won't be unitary in general.
# Instead, build V via two reflections (Householder).

def unitary_mapping(a, b):
    """Build unitary U s.t. U @ a = b, where ‖a‖=‖b‖=1.
    Uses the formula: V = I - |a⟩⟨a| - |b⟩⟨b| + |b⟩⟨a| + |a_perp⟩⟨b_perp|
    where a_perp, b_perp complete the 2D subspace."""
    n = len(a)
    if np.allclose(a, b):
        return np.eye(n, dtype=complex)
    # Build V that acts as rotation in span{a, b} and identity on complement.
    # Gram-Schmidt: get component of b orthogonal to a
    overlap = np.vdot(a, b)
    b_orth = b - overlap * a
    norm_b_orth = np.linalg.norm(b_orth)
    if norm_b_orth < 1e-14:
        # a and b differ only by phase: V = (phase) * projector + identity on complement
        phase = np.vdot(a, b)
        return np.eye(n, dtype=complex) + (phase - 1) * np.outer(a, np.conj(a))
    b_orth = b_orth / norm_b_orth
    # In the {a, b_orth} subspace, we want to map a → b.
    # a has coords (1, 0), b has coords (overlap, norm_b_orth)
    # The 2x2 rotation matrix R maps (1,0) → (overlap, norm_b_orth)
    # R = [[overlap, -conj(norm_b_orth)], [norm_b_orth, conj(overlap)]]
    # But we need this to be unitary, so:
    # V = I + (R - I_2) projected onto span{a, b_orth}
    # V = I + |a⟩(overlap-1)⟨a| + |a⟩ norm_b_orth ⟨b_orth|
    #       + |b_orth⟩(-conj(norm_b_orth))⟨a| + |b_orth⟩(conj(overlap)-1)⟨b_orth|
    # Simpler: V = I - |a⟩⟨a| - |b_orth⟩⟨b_orth| + [a b_orth] R [a b_orth]†
    P = np.outer(a, np.conj(a)) + np.outer(b_orth, np.conj(b_orth))
    # 2x2 unitary mapping (1,0)→(overlap, norm_b_orth)
    R_2x2 = np.array([[overlap, -np.conj(norm_b_orth)],
                       [norm_b_orth, np.conj(overlap)]])
    # Embed: V = (I - P) + [a, b_orth] R_2x2 [a, b_orth]†
    AB = np.column_stack([a, b_orth])  # n x 2
    V = np.eye(n, dtype=complex) - P + AB @ R_2x2 @ AB.conj().T
    return V


V_perfect = unitary_mapping(phi_x, phi_y)

# Verify V is unitary
print(f"  ‖V†V - I‖            = {np.linalg.norm(V_perfect.conj().T @ V_perfect - np.eye(d)):.2e}")

# Verify V φ(x) = φ(y)
mapped = V_perfect @ phi_x
print(f"  ‖V φ(x) - φ(y)‖     = {np.linalg.norm(mapped - phi_y):.2e}")

predictor_perfect = lambda v: V_perfect @ v
E_perfect = qjepa_energy(encoder, predictor_perfect, 'x', 'y')
print(f"  E(x,y; V_perfect)    = {E_perfect:.2e}")
assert E_perfect < 1e-20, "FAIL: perfect predictor should give zero energy"
print("  PASS: perfect predictor → zero energy")

# --- Test 3: Collapsed encoder → zero energy for V=I ---
print("\n--- Test 3: Collapsed encoder + identity connection ---")
c = random_state(d)
collapsed_encoder = lambda p: c  # constant encoder

E_collapsed = qjepa_energy(collapsed_encoder, identity, 'x', 'y')
print(f"  E(x,y; collapsed, I) = {E_collapsed:.2e}")
assert E_collapsed < 1e-30, "FAIL: collapsed + identity should give zero"
print("  PASS: collapsed encoder + V=I → zero energy")

# --- Test 4: Isometry composition ---
print("\n--- Test 4: Composition of isometries preserves norms ---")
# Two random unitaries
U1 = np.linalg.qr(np.random.randn(d, d) + 1j * np.random.randn(d, d))[0]
U2 = np.linalg.qr(np.random.randn(d, d) + 1j * np.random.randn(d, d))[0]
v_test = random_state(d)

norm_v = np.linalg.norm(v_test)
norm_composed = np.linalg.norm(U1 @ U2 @ v_test)
print(f"  ‖v‖                  = {norm_v:.10f}")
print(f"  ‖U₁ U₂ v‖           = {norm_composed:.10f}")
print(f"  Difference            = {abs(norm_v - norm_composed):.2e}")
assert abs(norm_v - norm_composed) < 1e-12, "FAIL: composition should preserve norm"
print("  PASS: composition of isometries preserves norms")

print("\n" + "=" * 60)
print("ALL TESTS PASSED")
print("=" * 60)
