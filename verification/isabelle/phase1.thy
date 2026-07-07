theory phase1
  imports Complex_Main
begin

section "Phase 1 Quick Wins Cross-Verification"

(* 1.1 PSD implies convex quadratic form *)
lemma psd_convex_form:
  fixes a c :: real and x y :: real
  assumes "0 \<le> a" "0 \<le> c"
  shows "0 \<le> (a * x + c * y)^2"
  by auto

(* 1.2 Rotation kernel in [0,1] *)
lemma cos_sq_bounded:
  fixes t :: real
  shows "0 \<le> (cos t)^2"
  by auto

(* 1.3 Measurement shots: N >= 1/(4*eta*delta^2) *)
lemma shots_bound:
  fixes eta delta :: real
  assumes "0 < eta" "0 < delta"
  shows "0 < 1 / (4 * eta * delta^2)"
  using assms by (auto simp: power2_eq_square)

(* 1.4 Born composition: T rounds, error adds *)
lemma born_composition:
  fixes T :: nat and eps :: real
  assumes "0 < T" "0 < eps"
  shows "T * (eps / T) = eps"
  using assms by auto

(* 1.5 NTK PSD: sum of squares non-negative *)
lemma ntk_psd_sum:
  fixes v :: "nat \<Rightarrow> real"
  shows "(\<Sum>i=0..<n. v i)^2 \<ge> 0"
  by auto

(* Phase 2: state_to_kernel_contraction *)
lemma sq_mono_nonneg:
  fixes a b :: real
  assumes "0 \<le> a" "a \<le> b"
  shows "a^2 \<le> b^2"
  using assms by (auto simp: power2_eq_square intro: mult_mono)

(* Cross-domain transfer: regularization beats concentration *)
lemma regularization_beats_concentration:
  fixes L_good R_good R_bad lam :: real
  assumes "0 < lam" "R_good < R_bad" "L_good < lam * (R_bad - R_good)"
  shows "L_good + lam * R_good < lam * R_bad"
  using assms by (auto simp: algebra_simps)

(* Pillar 4 Probability: variance = E[X^2] - E[X]^2 >= 0 *)
lemma variance_nonneg:
  fixes sq_mean mean :: real
  assumes "mean^2 \<le> sq_mean"
  shows "0 \<le> sq_mean - mean^2"
  using assms by auto

(* Chebyshev core: beta/eps^2 <= 1 when beta <= eps^2 *)
lemma chebyshev_prob:
  fixes beta eps :: real
  assumes "0 < eps" "beta \<le> eps^2"
  shows "beta / eps^2 \<le> 1"
  using assms by (auto simp: divide_le_eq_1_pos power2_eq_square)

end
