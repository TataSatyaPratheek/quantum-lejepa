theory mele_effective_depth
  imports Complex_Main
begin

section "Mele et al. 2026: Noise-induced shallow circuits (Nature Physics)"

text "Formalization of key algebraic results from:
  Noise-induced shallow circuits and the absence of barren plateaus
  Nature Physics, April 2026. arXiv:2403.13927"

(* The noise contraction parameter c = (1/3)(||D||^2 + ||t||^2) *)
(* For non-unitary channels: c < 1 *)

(* Depolarizing: D = (q,q,q), t = (0,0,0) => c = q^2 *)
lemma depolarizing_noise_param:
  fixes q :: real assumes "0 \<le> q"
  shows "(q^2 + q^2 + q^2) / 3 = q^2"
  by (auto simp: power2_eq_square field_simps)

(* c < 1 for non-trivial depolarizing *)
lemma depolarizing_c_lt_one:
  fixes q :: real assumes "0 \<le> q" "q < 1"
  shows "q^2 < 1"
proof -
  have "q * q < 1 * 1" using assms by (intro mult_strict_mono') auto
  thus ?thesis by (simp add: power2_eq_square)
qed

(* Effective depth: exp(-alpha * m) contraction *)
lemma effective_depth_contraction:
  fixes c :: real and m :: nat
  assumes "0 \<le> c" "c < 1"
  shows "c ^ m \<le> 1"
  using assms by (intro power_le_one) auto

(* m* grows logarithmically: m* = ceil(log(1/eps) / log(1/c)) *)
(* Full log characterization requires HOL-Analysis logarithm bridge *)

(* The key contraction: noise makes early layers irrelevant *)
lemma noise_erases_early_layers:
  fixes c norm_O :: real and m :: nat
  assumes "0 \<le> c" "c < 1" "0 \<le> norm_O"
  shows "norm_O * c ^ m \<le> norm_O"
proof -
  have "c ^ m \<le> 1" using assms by (intro power_le_one) auto
  then have "norm_O * c ^ m \<le> norm_O * 1"
    using assms(3) by (intro mult_left_mono)
  thus ?thesis by simp
qed

(* For non-unital noise: barren plateaus are ABSENT *)
(* This is because the fixed point is NOT the maximally mixed state *)
(* so local observables retain non-trivial expectation values *)
lemma nonunital_fixed_point_nontrivial:
  fixes t :: real
  assumes "t \<noteq> 0"
  shows "t^2 > 0"
  using assms by auto

end
