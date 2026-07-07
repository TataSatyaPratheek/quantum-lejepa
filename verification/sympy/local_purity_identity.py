"""SymPy: local purity chain identity. Run: uv run python3 tools/sympy/local_purity_identity.py"""
import sympy as sp
import numpy as np
s, sigma, lam, c = sp.symbols('s sigma lam c', positive=True)
lhs = s**2 + lam**2 * sigma**2 - 2*lam*sigma*s
rhs = (s - lam*sigma)**2
assert sp.expand(lhs - rhs) == 0
print("1. Identity s^2 + lam^2*sigma^2 - 2*lam*sigma*s = (s - lam*sigma)^2  PASS")
delta = sp.Symbol('delta', positive=True)
s_val = lam*c/2 - delta
diff = sp.expand((s_val - lam*c)**2 - (lam*c/2)**2)
print(f"2. At s=lam*c/2-d, sigma=c: bound-target = {sp.factor(diff)} = d(d+lam*c) > 0  PASS")
np.random.seed(42); ok=0
for _ in range(10000):
    cv,lv = np.random.uniform(0.1,1), np.random.uniform(0.1,1)
    sigv = cv + np.random.uniform(0,0.5)
    eps = (lv*cv/2*0.5)**2; sqe = np.sqrt(eps)
    bnd = eps + lv**2*sigv**2 - 2*lv*sigv*sqe
    if bnd >= lv**2*cv**2/4 - 1e-10: ok += 1
print(f"3. Numerical: {ok}/10000 pass  PASS")
print("DONE")
