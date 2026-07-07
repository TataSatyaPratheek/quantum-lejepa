(** Gradient variance lower bound — Coq/Rocq 9.1 verification
    Run: opam exec -- coqc tools/coq/gradient_bound.v *)

From Stdlib Require Import Reals Lra Lia.
From Stdlib Require Import Reals.Rsqrt_def.
Open Scope R_scope.

(** Theorem 1: QM-AM core identity
    For two values: (a + b)² ≤ 2(a² + b²)
    This is the N=2 case of QM-AM: (Σx)² ≤ N·Σx². *)
Lemma qm_am_two : forall a b : R,
  (a + b)^2 <= 2 * (a^2 + b^2).
Proof.
  intros a b.
  (* (a+b)² = a² + 2ab + b² ≤ 2a² + 2b² ↔ 0 ≤ (a-b)² *)
  assert (H: 0 <= (a - b)^2) by (apply pow2_ge_0).
  nra.
Qed.

(** Theorem 1b: QM-AM for N values (stated for triples to show pattern) *)
Lemma qm_am_three : forall a b c : R,
  (a + b + c)^2 <= 3 * (a^2 + b^2 + c^2).
Proof.
  intros a b c.
  assert (H1: 0 <= (a - b)^2) by (apply pow2_ge_0).
  assert (H2: 0 <= (b - c)^2) by (apply pow2_ge_0).
  assert (H3: 0 <= (a - c)^2) by (apply pow2_ge_0).
  nra.
Qed.

(** Theorem 2: Gradient variance identity (algebraic core)
    ε + λ²σ² - 2λσ√ε = (√ε - λσ)²

    We prove the algebraic identity without sqrt:
    if s² = ε (i.e., s = √ε), then ε + λ²σ² - 2λσs = (s - λσ)² *)
Lemma gradient_variance_identity : forall s sigma lam : R,
  s >= 0 -> sigma >= 0 -> lam > 0 ->
  s^2 + lam^2 * sigma^2 - 2 * lam * sigma * s = (s - lam * sigma)^2.
Proof.
  intros s sigma lam Hs Hsigma Hlam.
  ring.
Qed.

(** Corollary: the bound is non-negative *)
Lemma gradient_variance_nonneg : forall s sigma lam : R,
  s >= 0 -> sigma >= 0 -> lam > 0 ->
  s^2 + lam^2 * sigma^2 - 2 * lam * sigma * s >= 0.
Proof.
  intros s sigma lam Hs Hsigma Hlam.
  rewrite gradient_variance_identity; try lra.
  assert (H := pow2_ge_0 (s - lam * sigma)). lra.
Qed.

(** Corollary: strict positivity when s < λσ *)
Lemma gradient_variance_positive : forall s sigma lam : R,
  s >= 0 -> sigma > 0 -> lam > 0 -> s < lam * sigma ->
  s^2 + lam^2 * sigma^2 - 2 * lam * sigma * s > 0.
Proof.
  intros s sigma lam Hs Hsigma Hlam Hlt.
  rewrite gradient_variance_identity; try lra.
  assert (Hne: s - lam * sigma < 0) by lra.
  nra.
Qed.

(** Theorem 2b: Quantitative bound
    When s < λσ/2 (i.e., ε < λ²σ²/4):
    bound ≥ (λσ/2)² = λ²σ²/4 *)
Lemma gradient_variance_quantitative : forall s sigma lam : R,
  s >= 0 -> sigma > 0 -> lam > 0 -> s < lam * sigma / 2 ->
  s^2 + lam^2 * sigma^2 - 2 * lam * sigma * s >= lam^2 * sigma^2 / 4.
Proof.
  intros s sigma lam Hs Hsigma Hlam Hlt.
  rewrite gradient_variance_identity; try lra.
  (* (s - λσ)² ≥ (λσ/2)² when s < λσ/2 *)
  assert (Hge: lam * sigma - s >= lam * sigma / 2) by lra.
  assert (Hpos: s - lam * sigma <= -(lam * sigma / 2)) by lra.
  nra.
Qed.
