"""NumPy: Gradient variance scaling n=4,6,8. Run: python3 tools/numpy/gradient_scaling.py"""
import numpy as np
np.random.seed(42)

def partial_trace_single(rho, n, keep):
    rho_t = rho.reshape([2]*2*n)
    traced = sorted(set(range(n)) - {keep}, reverse=True)
    nr, k = n, keep
    for q in traced:
        rho_t = np.trace(rho_t, axis1=q, axis2=q+nr); nr -= 1
        if q < k: k -= 1
    return rho_t.reshape(2, 2)

def local_purity_sum(state, n):
    rho = np.outer(state, state.conj())
    return sum(np.real(np.trace(partial_trace_single(rho, n, k) @ partial_trace_single(rho, n, k))) for k in range(n))

def circuit(theta, n, layers=2):
    d = 2**n; state = np.zeros(d, dtype=complex); state[0] = 1.0
    Ry = lambda a: np.array([[np.cos(a/2),-np.sin(a/2)],[np.sin(a/2),np.cos(a/2)]], dtype=complex)
    idx = 0
    for _ in range(layers):
        U = np.eye(1, dtype=complex)
        for k in range(n): U = np.kron(U, Ry(theta[idx])); idx += 1
        state = U @ state
        for k in range(n-1):
            for i in range(d):
                bits = [(i>>(n-1-q))&1 for q in range(n)]
                if bits[k]==1 and bits[k+1]==1: state[i] *= -1
    return state

print("="*60)
print("GRADIENT VARIANCE SCALING: THE MONEY PLOT")
print("="*60)
print(f"{'n':>3} {'d':>5} {'Var_global':>12} {'Var_local':>12} {'Ratio':>7} {'1/d^2':>12}")
print("-"*55)

results = {}
for n in [4, 6, 8, 10]:
    d = 2**n; layers = 2; np_ = n*layers; shift = np.pi/2
    ref = np.zeros(d, dtype=complex); ref[0] = 1.0
    N = 200 if n <= 6 else (80 if n <= 8 else 50)
    gg, gl = [], []
    for s in range(N):
        if (s+1) % 50 == 0: print(f"  n={n}: {s+1}/{N}")
        theta = np.random.uniform(0, 2*np.pi, np_)
        tp = theta.copy(); tp[0] += shift
        tm = theta.copy(); tm[0] -= shift
        sp, sm = circuit(tp, n, layers), circuit(tm, n, layers)
        gg.append((abs(np.vdot(ref,sp))**2 - abs(np.vdot(ref,sm))**2)/(2*np.sin(shift)))
        gl.append((local_purity_sum(sp,n) - local_purity_sum(sm,n))/(2*np.sin(shift)))
    vg, vl = np.var(gg), np.var(gl)
    results[n] = (vg, vl)
    print(f"{n:>3d} {d:>5d} {vg:>12.4e} {vl:>12.4e} {vl/vg if vg>0 else 0:>7.1f} {1/d**2:>12.4e}")

ns = sorted(results.keys())
cg = np.polyfit(ns, [np.log(results[n][0]) for n in ns], 1)
cl = np.polyfit(ns, [np.log(results[n][1]) for n in ns], 1)
print(f"\nGlobal: log(Var) ~ {cg[0]:.2f}*n  (1/4^n = {-2*np.log(2):.2f})")
print(f"Local:  log(Var) ~ {cl[0]:.2f}*n  (poly if near 0)")
print(f"Global concentrates: {'YES' if cg[0] < -0.8 else 'NO'}")
print(f"Local trainable:     {'YES' if abs(cl[0]) < 0.5 else 'MAYBE'}")
print("DONE")
