# Quantum LeJEPA

**A joint-embedding predictive architecture for quantum machine learning — and the operator that
tells you, without training, which quantum circuits can actually learn.**

![Domain](https://img.shields.io/badge/quantum%20ML-representation%20learning-1f6feb)
![Object](https://img.shields.io/badge/M%E2%82%82-second%20moment%20operator-6f42c1)
![Verification](https://img.shields.io/badge/Lean%204-572%20decls%2C%200%20sorry-brightgreen)
![Status](https://img.shields.io/badge/personal%20research-independent-orange)

Personal research applying the **LeJEPA** self-supervised framework (Balestriero & LeCun, 2025) to
quantum feature maps, and identifying the single operator that governs quantum-kernel trainability.
Independent of my work at Accelequant.

## Papers

| | |
|:--|:--|
| 📄 **[Second-Moment-Operator.pdf](./Second-Moment-Operator.pdf)** | The headline article — *the second moment operator M₂ determines quantum kernel trainability, expressivity, and sample efficiency.* |
| 📄 **[Quantum-LeJEPA-Theory.pdf](./Quantum-LeJEPA-Theory.pdf)** | The full theory — *The Compact Topology of Quantum Representation Learning: A Fugue on Six Themes.* |
| 📄 **[Empirical-Results.pdf](./Empirical-Results.pdf)** | Circuit selection, gradient validation, classification, and RL on board games. |
| 📄 **[Benchmarking-Plan.pdf](./Benchmarking-Plan.pdf)** | Baselines and published-result comparisons (Bowles/Schuld protocol). |

## The core idea: M₂, a training-free circuit diagnostic

The **second moment operator** `M₂ = 𝔼_θ[U(θ)^⊗2 ⊗ U(θ)*^⊗2]` is a positive operator on the doubled
Hilbert space whose spectrum encodes three things at once: its **spectral gap** governs
concentration / barren plateaus (trainability), its **rank** equals the dimension of the dynamical
Lie algebra (expressivity), and its **local restriction** bounds the gradient. So one computable
object predicts whether a circuit can be trained — *before* you train it.

## What the experiments show (honest, including where it fails)

- **M₂ ranks circuits without training.** Screening **216 circuit configurations in 75 seconds**,
  M₂'s spectral rank predicts downstream accuracy: full-rank designs (256/256, Haar distance <0.07)
  reach **88–91%** classification, half-rank designs (124/256) only **78%**.
- **A concrete discovery.** M₂ exposed that the ubiquitous **Ry+CZ ansatz generates only the real
  orthogonal group O(d), not the full U(d)** — both gates are real-valued, so it accesses just half
  the doubled Hilbert space (rank 124/256). Two rotation axes (Rx+Ry, Ry+Rz, or U3) are required.
- **A universal gradient bound.** The identity `Var[∂L_reg/∂θ] ≥ (√ε − λσ)²` plus a *local* purity
  regularizer keeps the gradient alive: local variance stays flat (~0.016) while the global one
  decays as 1/4ⁿ — a **local-to-global ratio exceeding 10⁹ at n = 20**.
- **Low-data advantage.** M₂-selected kernels hit **99.3% from N = 50** samples (vs RBF 90.3%,
  XGBoost 86.0%); classical overtakes only past N ≈ 500.
- **Where it does *not* help (stated plainly).** Continuous regression (R² 0.29 vs RBF 0.67),
  temporal prediction (kernel–task correlation r = 0.019), and contact-rich control all fail
  pre-verification. Quantum kernels help for **static, discrete, geometric classification in the
  low-data regime** — and the paper says exactly where the boundary is.

## Formal verification

The mathematics is formalized in **Lean 4 (572 declarations, 0 `sorry`, 4 axioms)** and
cross-checked across six independent tools (SageMath, Isabelle, Coq/Rocq, Agda, SymPy, NumPy). *This
repository contains the papers; the Lean development is kept separate.*

## Author

**Satya Pratheek Tata** — [github.com/TataSatyaPratheek](https://github.com/TataSatyaPratheek)
· [LinkedIn](https://linkedin.com/in/satyapratheek-tata)

Personal research project, independent of my work at Accelequant.

## License

© 2026 Satya Pratheek Tata. The papers are shared for scholarly reading and citation; please cite
rather than redistribute modified copies.
