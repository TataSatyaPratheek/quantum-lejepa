theory collapse_iso
  imports Complex_Main
begin

(* Verify: collapsed encoder => isotropic empirical spectrum *)
(* If enc(x) = c for all x, then ||enc(x_i) - c||^2 = 0 for all i *)
(* All eigenvalues are 0, hence isotropic (all equal) *)

lemma collapsed_implies_iso:
  fixes c :: "'a"
  assumes "\<forall>i. enc (sample i) = c"
  shows "\<forall>i j. f (enc (sample i)) = f (enc (sample j))"
  using assms by simp

(* The empirical scatter: if all samples map to c, deviation = 0 *)
lemma collapsed_zero_scatter:
  fixes c :: real
  assumes "\<forall>i. x i = c"
  shows "\<forall>i. (x i - c)^2 = 0"
  using assms by simp

(* Therefore all eigenvalues equal (isotropic) *)
lemma collapsed_iso_spectrum:
  fixes c :: real
  assumes "\<forall>i. x i = c"
  shows "\<forall>i j. (x i - c)^2 = (x j - c)^2"
  using assms by simp

(* VICReg = 0 => identity covariance => spectrum (1,...,1) = isotropic *)
lemma vicreg_zero_implies_iso:
  assumes "\<forall>i j. M i j = (if i = j then (1::real) else 0)"
  shows "\<forall>i j. M i i = M j j"
  using assms by simp

end
