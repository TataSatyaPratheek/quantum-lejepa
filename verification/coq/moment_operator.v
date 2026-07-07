(** M₂ moment operator: spectral gap and trace — Coq/Rocq 9.1
    Run: opam exec -- coqc tools/coq/moment_operator.v *)
From Stdlib Require Import Reals Lra.
Open Scope R_scope.

(** M₂ spectral gap: 1/(d(d-1)) - 1/(d(d+1)) = 2/(d(d²-1)) *)
Lemma m2_spectral_gap : forall d : R,
  d > 1 -> d <> 0 -> d - 1 <> 0 -> d + 1 <> 0 ->
  1/(d*(d-1)) - 1/(d*(d+1)) = 2 / (d * (d^2 - 1)).
Proof.
  intros d Hd Hd0 Hdm Hdp. field. repeat split; nra.
Qed.

(** M₂ trace normalization *)
Lemma m2_trace_one : forall d : R,
  d > 1 -> d <> 0 -> d - 1 <> 0 -> d + 1 <> 0 ->
  d*(d+1)/2 * (1/(d*(d+1))) + d*(d-1)/2 * (1/(d*(d-1))) = 1.
Proof.
  intros d Hd Hd0 Hdm Hdp. field. repeat split; nra.
Qed.

(** Spectral gap positive when d > 1 *)
Lemma m2_gap_positive : forall d : R, d > 1 -> 2 / (d * (d^2 - 1)) > 0.
Proof. intros. apply Rdiv_lt_0_compat; nra. Qed.
