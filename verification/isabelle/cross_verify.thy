theory cross_verify
  imports Complex_Main
begin

section \<open>Cross-Verification of Lean Theorems in Isabelle/HOL\<close>

text \<open>
  Each lemma independently verifies a key theorem from the Lean 4
  formalization at jepa-mathlib/QML/. This provides a second proof
  assistant's confirmation of mathematical correctness.

  Run: isabelle build -D tools/isabelle
\<close>

subsection \<open>1. QM-AM = K * Variance (DeepConnections.lean)\<close>

lemma qm_am_is_variance:
  fixes S T K :: real assumes "K > 0"
  shows "S - T\<^sup>2/K = K * (S/K - (T/K)\<^sup>2)"
  using assms by (auto simp: field_simps power2_eq_square)

subsection \<open>2. Bridge Inequality (Bridge.lean)\<close>

lemma bridge_inequality:
  fixes k r :: real
  assumes "r\<^sup>2 \<le> k" "r \<le> 1" "0 \<le> r"
  shows "1 - k \<le> 2*(1-r)"
proof -
  have "1 - k \<le> 1 - r\<^sup>2" using assms(1) by auto
  also have "... = (1-r)*(1+r)" by (auto simp: power2_eq_square algebra_simps)
  also have "... \<le> (1-r)*2"
    using assms(2,3) by (intro mult_left_mono) auto
  finally show ?thesis by auto
qed

subsection \<open>3. Noise Contraction (NoiseChain.lean)\<close>

lemma noise_contraction:
  fixes q tvd :: real
  assumes "0 \<le> q" "q \<le> 1" "0 \<le> tvd"
  shows "q^b * tvd \<le> tvd"
proof -
  have "q^b \<le> 1" using assms(1,2) by (intro power_le_one) auto
  then have "q^b * tvd \<le> 1 * tvd"
    using assms(3) by (intro mult_right_mono)
  thus ?thesis by simp
qed

subsection \<open>4. EMA Contraction (JEPA.lean)\<close>

lemma ema_contraction:
  fixes tau :: real assumes "0 \<le> tau" "tau \<le> 1"
  shows "tau^t \<le> 1"
  using assms by (intro power_le_one) auto

subsection \<open>5. Spectral Gap (SpectralBridge.lean)\<close>

lemma spectral_gap:
  fixes n mu :: real
  shows "(1 + (n-1)*mu) - (1-mu) = n*mu"
  by algebra

subsection \<open>6. Pfaffian Expansion (Pfaffian.lean)\<close>

lemma pfaffian_sq_expand:
  fixes a b c d e f :: real
  shows "(a*f - b*e + c*d)^2 =
    a^2*f^2 + b^2*e^2 + c^2*d^2
    - 2*a*f*b*e + 2*a*f*c*d - 2*b*e*c*d"
  by algebra

subsection \<open>7. Dimension Formula (SchurWeyl.lean)\<close>

lemma sym2_plus_alt2:
  fixes d :: real
  shows "d*(d+1)/2 + d*(d-1)/2 = d^2"
  by (auto simp: power2_eq_square field_simps)

end
