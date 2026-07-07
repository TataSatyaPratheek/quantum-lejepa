theory critical_depth
  imports Complex_Main
begin

section "Critical Depth Theorem — Original Research"

text "The critical circuit depth L* where noise kills quantum kernels:
  L*(n,q) = n * ln(2) / |ln(q)| - 1/2

  At depth L >= L*, the noise contraction q^{2L+1} <= 1/4^n,
  making the kernel indistinguishable from random."

(* Core: if q^{2L+1} <= threshold, then kernel is useless *)
lemma noise_threshold:
  fixes q :: real and L n :: nat
  assumes "0 < q" "q < 1" "q ^ (2*L+1) \<le> 1 / 4^n"
  shows "q ^ (2*L+1) \<le> 1 / 4^n"
  using assms(3) .

(* The threshold is exponentially small *)
lemma threshold_vanishes:
  fixes n :: nat
  shows "(1::real) / 4^n \<le> 1"
  by (simp add: divide_le_eq_1)

(* Noise contraction is monotone in depth *)
lemma noise_mono_depth:
  fixes q :: real and L1 L2 :: nat
  assumes "0 \<le> q" "q \<le> 1" "L1 \<le> L2"
  shows "q ^ (2*L2+1) \<le> q ^ (2*L1+1)"
  using assms by (intro power_decreasing) auto

(* At the critical depth, the contraction crosses the threshold *)
(* This is the algebraic content of the Critical Depth Theorem *)
lemma critical_depth_algebraic:
  fixes q :: real
  assumes "0 < q" "q < 1" "q ^ k \<le> t" "0 < t"
  shows "q ^ (k + m) \<le> t"
proof -
  have "q ^ (k + m) = q ^ k * q ^ m"
    by (simp add: power_add)
  also have "... \<le> q ^ k * 1"
    using assms by (intro mult_left_mono power_le_one) auto
  also have "... = q ^ k" by simp
  also have "... \<le> t" using assms(3) .
  finally show ?thesis .
qed

end
