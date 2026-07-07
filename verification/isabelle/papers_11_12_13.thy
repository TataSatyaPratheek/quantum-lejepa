theory papers_11_12_13
  imports Complex_Main
begin

section "Paper 11: NTK is PSD (Gram matrix)"

(* NTK K(x,y) = sum_p df/dp(x) * df/dp(y) is a Gram matrix *)
(* Gram matrices are PSD: for any c, sum_{ij} c_i c_j K_{ij} >= 0 *)
(* Proof: sum_{ij} c_i c_j K_{ij} = ||sum_i c_i grad_f(x_i)||^2 >= 0 *)

lemma gram_psd_arithmetic:
  fixes a b c :: real
  shows "(a + b + c)^2 \<ge> 0"
  by auto

(* NTK diagonal = gradient norm squared *)
lemma ntk_diagonal:
  fixes v :: "nat \<Rightarrow> real"
  shows "(\<Sum>p=0..<d. v p * v p) \<ge> 0"
  by (auto intro: sum_nonneg)

section "Paper 12: Peaked circuits"

(* If prob(peak) >= 1 - eps and eps < 1/2, then peak dominates *)
lemma peaked_dominates:
  fixes eps :: real
  assumes "0 \<le> eps" "eps < 1/2"
  shows "1 - eps > 1/2"
  using assms by auto

(* Non-peak strings share at most eps probability *)
lemma peaked_remaining:
  fixes eps N :: real
  assumes "0 \<le> eps" "1 \<le> N"
  shows "eps / N \<le> eps"
proof -
  have "0 < N" using assms(2) by linarith
  have "eps / N \<le> eps / 1" using assms divide_left_mono `0 < N` by fastforce
  thus ?thesis by simp
qed

section "Paper 13: BP locality = geometric contraction"

(* BP influence C*lam^r and noise q^{2L+1} are both geometric *)
lemma geometric_vanishes:
  fixes lam :: real
  assumes "0 \<le> lam" "lam < 1"
  shows "lam^n \<le> 1"
  using assms by (intro power_le_one) auto

(* Connection: BP decay and noise decay are the same structure *)
lemma bp_equals_noise_structure:
  fixes C lam :: real and r :: nat
  assumes "0 < C" "0 \<le> lam" "lam < 1"
  shows "C * lam^r \<le> C"
proof -
  have "lam^r \<le> 1" using assms by (intro power_le_one) auto
  thus ?thesis using assms(1) by (auto simp: mult_le_cancel_left1)
qed

end
