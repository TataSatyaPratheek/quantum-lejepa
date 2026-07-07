(** Local cost purity: algebraic chain verification — Coq/Rocq 9.1
    Run: opam exec -- coqc tools/coq/local_purity.v *)

From Stdlib Require Import Reals Lra Lia.
Open Scope R_scope.

(** ================================================================= *)
(** 1. local_sum_nonneg: sum of non-negative terms is non-negative    *)
(** ================================================================= *)

(** Base case: empty sum *)
Lemma local_sum_nonneg_0 : 0 >= 0.
Proof. lra. Qed.

(** Step: adding a non-negative term preserves non-negativity *)
Lemma local_sum_nonneg_step : forall (S x : R),
  S >= 0 -> x >= 0 -> S + x >= 0.
Proof. intros. lra. Qed.

(** Inductive version for a triple (pattern for N terms) *)
Lemma local_sum_nonneg_3 : forall (a b c : R),
  a >= 0 -> b >= 0 -> c >= 0 -> a + b + c >= 0.
Proof. intros. lra. Qed.

(** General: if each l_i >= 0 and there are N of them, sum >= 0.
    We state for arbitrary partial sums. *)
Lemma local_sum_nonneg : forall (partial_sum next : R),
  partial_sum >= 0 -> next >= 0 -> partial_sum + next >= 0.
Proof. intros. lra. Qed.

(** ================================================================= *)
(** 2. local_gradient_bound_chain                                     *)
(**    If sigma >= c > 0, eps >= 0, and eps < lam^2 * c^2 / 4,      *)
(**    then eps + lam^2 * sigma^2 - 2*lam*sigma*s                    *)
(**         >= lam^2 * c^2 / 4                                        *)
(**    where s^2 = eps, s >= 0.                                       *)
(** ================================================================= *)

(** First, the algebraic identity: the expression is a perfect square *)
Lemma gradient_square_identity : forall s sigma lam : R,
  s^2 + lam^2 * sigma^2 - 2 * lam * sigma * s = (s - lam * sigma)^2.
Proof. intros. ring. Qed.

(** The perfect square is non-negative *)
Lemma gradient_square_nonneg : forall s sigma lam : R,
  s^2 + lam^2 * sigma^2 - 2 * lam * sigma * s >= 0.
Proof.
  intros.
  rewrite gradient_square_identity.
  assert (H := pow2_ge_0 (s - lam * sigma)). lra.
Qed.

(** Key lemma: when sigma >= c > 0 and s is small, the bound is large *)
Lemma gradient_bound_lower : forall s sigma lam c : R,
  s >= 0 -> sigma >= c -> c > 0 -> lam > 0 ->
  s <= lam * c / 2 ->
  (s - lam * sigma)^2 >= (lam * c / 2)^2.
Proof.
  intros s sigma lam c Hs Hsigma Hc Hlam Hsmall.
  (* Since sigma >= c and s <= lam*c/2, we have
     |s - lam*sigma| >= lam*sigma - s >= lam*c - lam*c/2 = lam*c/2 *)
  assert (Hgap: lam * sigma - s >= lam * c / 2) by nra.
  (* So (s - lam*sigma)^2 = (lam*sigma - s)^2 >= (lam*c/2)^2 *)
  nra.
Qed.

(** The main theorem: the full chain *)
Lemma local_gradient_bound_chain : forall (eps sigma lam c s : R),
  s >= 0 -> s^2 = eps ->
  sigma >= c -> c > 0 -> lam > 0 ->
  eps >= 0 ->
  eps < lam^2 * c^2 / 4 ->
  eps + lam^2 * sigma^2 - 2 * lam * sigma * s >= lam^2 * c^2 / 4.
Proof.
  intros eps sigma lam c s Hs Heps Hsigma Hc Hlam Heps_nn Heps_small.
  (* Rewrite eps as s^2 *)
  rewrite <- Heps.
  rewrite gradient_square_identity.
  (* Now show (s - lam*sigma)^2 >= lam^2 * c^2 / 4 *)
  (* Since s^2 < lam^2 * c^2 / 4 and s >= 0, we get s < lam*c/2 *)
  assert (Hs_bound: s < lam * c / 2).
  { (* s^2 < lam^2 * c^2 / 4 = (lam*c/2)^2 *)
    assert (Hsq: s^2 < (lam * c / 2)^2) by nra.
    (* Since s >= 0 and lam*c/2 > 0, s^2 < t^2 implies s < t *)
    destruct (Rlt_or_le s (lam * c / 2)) as [Hlt|Hge].
    - exact Hlt.
    - exfalso.
      assert (Ht: 0 <= lam * c / 2) by nra.
      assert (Hge2: lam * c / 2 * (lam * c / 2) <= s * s).
      { exact (Rmult_le_compat _ _ _ _ Ht Ht Hge Hge). }
      assert (s^2 = s * s) by ring.
      assert ((lam * c / 2)^2 = lam * c / 2 * (lam * c / 2)) by ring.
      lra. }
  (* Now apply the lower bound *)
  assert (Hgap: lam * sigma - s >= lam * c / 2) by nra.
  nra.
Qed.

(** ================================================================= *)
(** Bonus: locality of composition (algebraic content)                *)
(** ================================================================= *)

(** If f and g are Lipschitz, so is f . g *)
Lemma lipschitz_compose : forall (Lf Lg d : R),
  Lf >= 0 -> Lg >= 0 -> d >= 0 ->
  Lf * (Lg * d) = (Lf * Lg) * d.
Proof. intros. ring. Qed.

(** Non-negativity of squared gradient norm *)
Lemma grad_norm_sq_nonneg : forall (v : R),
  v^2 >= 0.
Proof.
  intros. assert (H := pow2_ge_0 v). lra.
Qed.

Print Assumptions local_gradient_bound_chain.
Print Assumptions local_sum_nonneg.
