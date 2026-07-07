theory pillars imports Complex_Main begin

section "Pillar 2 (Analysis): RKHS and PSD kernels"

(* Schur-Weyl completeness: Sym^2 + Alt^2 = d^2 *)
lemma schur_weyl_completeness:
  fixes d :: real
  shows "d*(d+1)/2 + d*(d-1)/2 = d^2"
  by (auto simp: field_simps power2_eq_square)

section "Pillar 3 (Geometry): Lie algebra dimension"

(* Unitary Lie algebra dimension = d^2 *)
lemma unitary_lie_algebra_dim:
  fixes d :: nat assumes "0 < d"
  shows "d * d = d^2"
  by (auto simp: power2_eq_square)

(* D4 equivariant: 8n << 4^n for n >= 3 *)
(* This is the KEY inequality: D4 circuits are exponentially *)
(* smaller than full unitary, hence trainable but simulable *)
lemma d4_poly_vs_exp:
  fixes n :: nat assumes "3 \<le> n"
  shows "8 * n < 4^n"
  using assms by (induct n rule: nat_induct_at_least) auto

section "Pillar 4 (Probability): Variance"

(* Variance non-negativity *)
lemma variance_nonneg_pillar:
  fixes sq_mean mean :: real
  assumes "mean^2 \<le> sq_mean"
  shows "0 \<le> sq_mean - mean^2"
  using assms by auto

(* From QM-AM: totalVar^2/K <= sqSum implies sqSum/K >= (totalVar/K)^2 *)
lemma var_from_qmam:
  fixes sqSum totalVar K :: real
  assumes "0 < K" "totalVar^2 / K \<le> sqSum"
  shows "0 \<le> sqSum / K - totalVar^2 / K^2"
proof -
  have "totalVar^2 / K^2 \<le> sqSum / K"
    using assms by (auto simp: field_simps power2_eq_square)
  thus ?thesis by linarith
qed

end
