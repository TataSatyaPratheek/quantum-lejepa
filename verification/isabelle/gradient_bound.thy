theory gradient_bound imports Complex_Main begin

section "Theorem 1: QM-AM (Cauchy-Schwarz) on finite sums"

lemma qm_am_finite:
  fixes x :: "nat \<Rightarrow> real" and N :: nat
  assumes "N > 0"
  shows "(\<Sum>i<N. x i)\<^sup>2 \<le> real N * (\<Sum>i<N. (x i)\<^sup>2)"
proof -
  define mu where "mu = (\<Sum>i<N. x i) / real N"
  have var_nonneg: "0 \<le> (\<Sum>i<N. (x i - mu)\<^sup>2)"
    by (intro sum_nonneg) auto
  have expand: "(\<Sum>i<N. (x i - mu)\<^sup>2)
      = (\<Sum>i<N. (x i)\<^sup>2) - 2 * mu * (\<Sum>i<N. x i) + real N * mu\<^sup>2"
  proof -
    have "(\<Sum>i<N. (x i - mu)\<^sup>2) = (\<Sum>i<N. (x i)\<^sup>2 - 2 * x i * mu + mu\<^sup>2)"
      by (intro sum.cong) (auto simp: power2_eq_square ring_distribs)
    also have "\<dots> = (\<Sum>i<N. (x i)\<^sup>2) - (\<Sum>i<N. 2 * x i * mu) + (\<Sum>i<N. mu\<^sup>2)"
      by (simp add: sum.distrib sum_subtractf)
    also have "\<dots> = (\<Sum>i<N. (x i)\<^sup>2) - 2 * mu * (\<Sum>i<N. x i) + real N * mu\<^sup>2"
      by (simp add: sum_distrib_left mult.commute mult.left_commute)
    finally show ?thesis .
  qed
  have simplify: "(\<Sum>i<N. (x i)\<^sup>2) - 2 * mu * (\<Sum>i<N. x i) + real N * mu\<^sup>2
      = (\<Sum>i<N. (x i)\<^sup>2) - (\<Sum>i<N. x i)\<^sup>2 / real N"
  proof -
    have Npos: "real N > 0" using assms by simp
    have "mu = (\<Sum>i<N. x i) / real N" by (simp add: mu_def)
    hence "2 * mu * (\<Sum>i<N. x i) = 2 * (\<Sum>i<N. x i)\<^sup>2 / real N"
      by (simp add: power2_eq_square field_simps)
    moreover have "real N * mu\<^sup>2 = (\<Sum>i<N. x i)\<^sup>2 / real N"
      using Npos by (simp add: mu_def power2_eq_square field_simps)
    ultimately show ?thesis by linarith
  qed
  from var_nonneg expand simplify
  have "0 \<le> (\<Sum>i<N. (x i)\<^sup>2) - (\<Sum>i<N. x i)\<^sup>2 / real N" by linarith
  hence "(\<Sum>i<N. x i)\<^sup>2 / real N \<le> (\<Sum>i<N. (x i)\<^sup>2)" by linarith
  thus ?thesis using assms by (simp add: field_simps)
qed


section "Theorem 2: Gradient variance lower bound"

text \<open>THE KEY IDENTITY:
  eps + lam^2 * sigma^2 - 2*lam*sigma*sqrt(eps) = (sqrt(eps) - lam*sigma)^2

  Proof strategy: expand the RHS and use sqrt(eps)^2 = eps.\<close>

lemma gradient_variance_identity:
  fixes eps sigma lam :: real
  assumes "eps \<ge> 0" "sigma \<ge> 0" "lam > 0"
  shows "eps + lam\<^sup>2 * sigma\<^sup>2 - 2 * lam * sigma * sqrt eps
       = (sqrt eps - lam * sigma)\<^sup>2"
proof -
  have se: "(sqrt eps)\<^sup>2 = eps" using assms(1) by simp
  show ?thesis using se by (simp add: power2_eq_square ring_distribs)
qed

corollary gradient_variance_nonneg:
  fixes eps sigma lam :: real
  assumes "eps \<ge> 0" "sigma \<ge> 0" "lam > 0"
  shows "eps + lam\<^sup>2 * sigma\<^sup>2 - 2 * lam * sigma * sqrt eps \<ge> 0"
  using gradient_variance_identity[OF assms] by simp

corollary gradient_variance_positive:
  fixes eps sigma lam :: real
  assumes "eps \<ge> 0" "sigma > 0" "lam > 0" "sqrt eps < lam * sigma"
  shows "eps + lam\<^sup>2 * sigma\<^sup>2 - 2 * lam * sigma * sqrt eps > 0"
proof -
  have eq: "eps + lam\<^sup>2 * sigma\<^sup>2 - 2 * lam * sigma * sqrt eps
          = (sqrt eps - lam * sigma)\<^sup>2"
    using gradient_variance_identity assms by simp
  have ne: "sqrt eps - lam * sigma \<noteq> 0" using assms(4) by linarith
  have "sqrt eps - lam * sigma < 0 \<or> sqrt eps - lam * sigma > 0"
    using ne by linarith
  thus ?thesis unfolding eq
    by (auto simp: power2_eq_square zero_less_mult_iff)
qed


section "Theorem 3: Local cost algebraic structure"

text \<open>Variance additivity for averaged independent terms:
  If R = (1/n) * sum_k R_k with variances v_k >= 0,
  total variance = (sum v_k) / n^2 >= 0.\<close>

lemma local_cost_variance_additive:
  fixes v :: "nat \<Rightarrow> real" and n :: nat
  assumes "n > 0"
  assumes "\<And>k. k < n \<Longrightarrow> v k \<ge> 0"
  shows "(\<Sum>k<n. v k) / (real n)\<^sup>2 \<ge> 0"
proof -
  have "(\<Sum>k<n. v k) \<ge> 0"
    by (intro sum_nonneg) (simp add: assms(2))
  moreover have "(real n)\<^sup>2 > 0"
    using assms(1) by simp
  ultimately show ?thesis by (simp add: divide_nonneg_pos)
qed

text \<open>Full chain: combining the local cost bound with gradient_variance_positive.
  If sigma^2 >= c > 0 and eps >= 0 and sqrt(eps) < lam * c,
  then eps + lam^2 * sigma^2 - 2*lam*sigma*sqrt(eps) > 0.\<close>

lemma local_gradient_chain:
  fixes eps sigma lam c :: real
  assumes "eps \<ge> 0" "sigma\<^sup>2 \<ge> c" "c > 0" "lam > 0" "sqrt eps < lam * c"
  shows "eps + lam\<^sup>2 * sigma\<^sup>2 - 2 * lam * sigma * sqrt eps > 0"
proof -
  have sigma_pos: "sigma > 0"
  proof -
    have "sigma\<^sup>2 \<ge> c" using assms(2) .
    hence "sigma\<^sup>2 > 0" using assms(3) by linarith
    thus "sigma > 0" by (metis power2_eq_square zero_less_mult_iff)
  qed
  have sigma_ge_c: "lam * sigma \<ge> lam * c"
  proof -
    have "sigma\<^sup>2 \<ge> c" using assms(2) .
    have "c > 0" using assms(3) .
    have "sigma > 0" using sigma_pos .
    have "sigma \<ge> sqrt c"
      using \<open>sigma\<^sup>2 \<ge> c\<close> \<open>sigma > 0\<close>
      by (metis real_le_rsqrt abs_of_pos less_imp_le)
    moreover have "sqrt c \<ge> c"
    proof (cases "c \<ge> 1")
      case True
      hence "c \<le> c\<^sup>2" by (simp add: power2_eq_square le_self_mult)
      hence "sqrt c \<le> sqrt (c\<^sup>2)" by (intro real_sqrt_le_mono) auto
      hence "sqrt c \<le> c" using \<open>c > 0\<close> by simp
      thus ?thesis using True by linarith
    next
      case False
      hence "c < 1" by simp
      hence "c\<^sup>2 < c" using \<open>c > 0\<close>
        by (simp add: power2_eq_square mult_strict_left_mono)
      hence "sqrt (c\<^sup>2) < sqrt c" by (intro real_sqrt_less_mono) auto
      hence "c < sqrt c" using \<open>c > 0\<close> by simp
      thus ?thesis by linarith
    qed
    ultimately have "sigma \<ge> c" by linarith
    thus "lam * sigma \<ge> lam * c" using assms(4) by (simp add: mult_left_mono)
  qed
  hence "sqrt eps < lam * sigma" using assms(5) by linarith
  thus ?thesis using gradient_variance_positive[OF assms(1) sigma_pos assms(4)] by auto
qed


section "Theorem 5: M2 spectral gap algebra"

text \<open>M2 eigenvalues on Sym^2 and Lambda^2 for Haar circuits:
  lambda_sym = 1/(d(d+1)), lambda_anti = 1/(d(d-1)).
  Spectral gap = lambda_anti - lambda_sym = 2/(d(d^2-1)).
  Trace normalization: dim_sym * lambda_sym + dim_anti * lambda_anti = 1.\<close>

lemma m2_spectral_gap:
  fixes d :: real
  assumes "d > 1"
  shows "1/(d*(d-1)) - 1/(d*(d+1)) = 2 / (d * (d\<^sup>2 - 1))"
proof -
  have d0: "d > 0" using assms by simp
  have dm: "d - 1 > 0" using assms by simp
  have dp: "d + 1 > 0" using assms by simp
  \<comment> \<open>Step 1: simplify each fraction\<close>
  have lhs1: "1/(d*(d-1)) = 1/d * (1/(d-1))"
    using d0 dm by (simp add: field_simps)
  have lhs2: "1/(d*(d+1)) = 1/d * (1/(d+1))"
    using d0 dp by (simp add: field_simps)
  \<comment> \<open>Step 2: factor out 1/d\<close>
  have "1/(d*(d-1)) - 1/(d*(d+1)) = 1/d * (1/(d-1) - 1/(d+1))"
    using lhs1 lhs2 by (simp add: algebra_simps)
  \<comment> \<open>Step 3: compute 1/(d-1) - 1/(d+1)\<close>
  also have "1/(d-1) - 1/(d+1) = 2/((d-1)*(d+1))"
    using dm dp by (simp add: field_simps)
  \<comment> \<open>Step 4: multiply\<close>
  also have "1/d * (2/((d-1)*(d+1))) = 2/(d*((d-1)*(d+1)))"
    using d0 dm dp by (simp add: field_simps)
  \<comment> \<open>Step 5: refactor denominator\<close>
  also have "d*((d-1)*(d+1)) = d*(d\<^sup>2 - 1)"
    by (simp add: power2_eq_square algebra_simps)
  hence "2/(d*((d-1)*(d+1))) = 2/(d*(d\<^sup>2 - 1))"
    by simp
  finally show ?thesis .
qed

lemma m2_trace_normalization:
  fixes d :: real
  assumes "d > 1"
  shows "d*(d+1)/2 * (1/(d*(d+1))) + d*(d-1)/2 * (1/(d*(d-1))) = 1"
proof -
  have d0: "d \<noteq> 0" using assms by simp
  have dp: "d + 1 \<noteq> 0" using assms by simp
  have dm: "d - 1 \<noteq> 0" using assms by simp
  show ?thesis using d0 dp dm
    by (simp add: divide_simps)
qed

end
