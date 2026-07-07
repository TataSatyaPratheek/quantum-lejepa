theory chebyshev
  imports Complex_Main
begin

section "Chebyshev Concentration Bridge"

lemma chebyshev_arithmetic:
  fixes beta eps :: real
  assumes "0 < eps" "0 \<le> beta" "beta \<le> eps^2"
  shows "0 \<le> 1 - beta / eps^2"
proof -
  have "eps^2 > 0" using assms(1) by (auto simp: power2_eq_square)
  then have "beta / eps^2 \<le> eps^2 / eps^2"
    using assms(3) by (intro divide_right_mono) auto
  then have "beta / eps^2 \<le> 1" using `eps^2 > 0` by auto
  thus ?thesis by linarith
qed

lemma concentration_vanishes:
  fixes d :: real
  assumes "d \<ge> 2"
  shows "1 / (d * (d + 1)) \<le> 1 / d^2"
proof -
  have "d^2 \<le> d * (d + 1)" using assms
    by (auto simp: power2_eq_square algebra_simps)
  moreover have "0 < d^2" using assms by (auto simp: power2_eq_square)
  moreover have "0 < d * (d + 1)" using assms by (auto intro: mult_pos_pos)
  ultimately show ?thesis by (auto intro: divide_left_mono)
qed

lemma kernel_dev_nonneg:
  fixes variance eps :: real
  assumes "0 < eps" "0 \<le> variance"
  shows "0 \<le> variance / eps^2"
  using assms by (auto simp: power2_eq_square)

end
