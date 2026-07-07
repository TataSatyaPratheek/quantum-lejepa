"""SymPy: M₂ eigenvalue identities. Run: uv run python3 tools/sympy/moment_operator.py"""
import sympy as sp
d = sp.Symbol('d', positive=True, real=True)
# Weingarten values
wg_e = 1 / (d**2 - 1)
wg_s = -1 / (d * (d**2 - 1))
# M₂ eigenvalues on Sym² and Λ²
lam_sym = wg_e + wg_s  # = 1/(d(d+1))
lam_anti = wg_e - wg_s  # = 1/(d(d-1))
print("1. lambda_sym =", sp.simplify(lam_sym), "expect 1/(d(d+1))")
print("   Check:", sp.simplify(lam_sym - 1/(d*(d+1))))
print("2. lambda_anti =", sp.simplify(lam_anti), "expect 1/(d(d-1))")
print("   Check:", sp.simplify(lam_anti - 1/(d*(d-1))))
gap = sp.simplify(lam_anti - lam_sym)
print("3. Gap =", gap, "expect 2/(d(d^2-1))")
print("   Check:", sp.simplify(gap - 2/(d * (d**2 - 1))))
# Trace normalization
dim_sym = d*(d+1)/2
dim_anti = d*(d-1)/2
trace = sp.simplify(dim_sym * lam_sym + dim_anti * lam_anti)
print("4. Trace =", trace, "expect 1")
print("DONE")
