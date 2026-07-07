#!/bin/bash
# ================================================================
# MULTI-TOOL VERIFICATION SUITE
# Verifies the jepa-mathlib formalization using 6 independent tools:
#   1. Lean 4 + Mathlib (primary formalization)
#   2. SageMath 10.8 (symbolic computation)
#   3. Isabelle/HOL 2025 (algebraic cross-verification)
#   4. Agda 2.7.0 (constructive type-checking)
#   5. Coq/Rocq 9.1 (independent kernel verification)
#   6. NumPy/SciPy (numerical verification)
#
# Run from project root: bash tools/verify_all.sh
# ================================================================

set -e
cd "$(dirname "$0")/.."

echo "================================================================"
echo "MULTI-TOOL VERIFICATION SUITE — jepa-mathlib"
echo "================================================================"
echo ""

# Tool 1: Lean 4
echo "--- [1/6] Lean 4 + Mathlib ---"
if lake build 2>&1 | tail -1 | grep -q "Build completed successfully"; then
    SORRY_COUNT=$(grep -rn "sorry" QML/ --include="*.lean" | grep -v "^.*:.*--" | wc -l | tr -d ' ')
    echo "  BUILD: ✓ (zero errors)"
    echo "  SORRY: $SORRY_COUNT"
else
    echo "  BUILD: ✗ FAILED"
fi
echo ""

# Tool 2: SageMath
echo "--- [2/6] SageMath 10.8 ---"
if command -v mamba &>/dev/null; then
    for script in tools/sage/*.py; do
        name=$(basename "$script" .py)
        OUT=$(mamba run -n sage python3 "$script" 2>&1)
        if echo "$OUT" | grep -qE "(✓|DONE|pass)"; then
            echo "  $name: ✓"
        else
            echo "  $name: ✗ FAILED"
        fi
    done
else
    echo "  SageMath: not found, skipping"
fi
echo ""

# Tool 3: Isabelle/HOL
echo "--- [3/6] Isabelle/HOL 2025 ---"
ISABELLE=/Applications/Isabelle2025.app/bin/isabelle
if [ -f "$ISABELLE" ]; then
    ISA_OUT=$($ISABELLE build -D tools/isabelle 2>&1)
    if echo "$ISA_OUT" | grep -q "Finished\|elapsed time"; then
        THEORY_COUNT=$(ls tools/isabelle/*.thy | wc -l | tr -d ' ')
        echo "  $THEORY_COUNT theories: ✓ VERIFIED"
    else
        echo "  Isabelle: ✗ FAILED"
    fi
else
    echo "  Isabelle: not found at $ISABELLE, skipping"
fi
echo ""

# Tool 4: Agda
echo "--- [4/6] Agda 2.7.0 ---"
if command -v agda &>/dev/null; then
    for script in tools/agda/*.agda; do
        name=$(basename "$script" .agda)
        AGDA_OUT=$(agda --include-path=tools/agda "$script" 2>&1)
        if [ $? -eq 0 ]; then
            echo "  $name: ✓ TYPE-CHECKED"
        else
            echo "  $name: ✗ FAILED"
        fi
    done
else
    echo "  Agda: not found, skipping"
fi
echo ""

# Tool 5: Coq/Rocq
echo "--- [5/6] Coq/Rocq 9.1 ---"
if opam exec -- coqc --version &>/dev/null; then
    for script in tools/coq/*.v; do
        name=$(basename "$script" .v)
        COQ_OUT=$(opam exec -- coqc "$script" 2>&1)
        if [ $? -eq 0 ]; then
            echo "  $name: ✓ COMPILED"
        else
            echo "  $name: ✗ FAILED"
        fi
    done
else
    echo "  Coq/Rocq: not found, skipping"
fi
echo ""

# Tool 6: NumPy/SciPy
echo "--- [6/6] NumPy/SciPy ---"
if python3 -c "import numpy" 2>/dev/null; then
    for script in tools/numpy/*.py; do
        name=$(basename "$script" .py)
        OUT=$(python3 "$script" 2>&1)
        if echo "$OUT" | grep -qE "(MATCH|VERIFIED|pass|✓)"; then
            echo "  $name: ✓"
        else
            echo "  $name: ✗ FAILED"
        fi
    done
fi
if python3 -c "import scipy" 2>/dev/null; then
    for script in tools/scipy/*.py; do
        name=$(basename "$script" .py)
        OUT=$(python3 "$script" 2>&1)
        if echo "$OUT" | grep -qE "(VERIFIED|isotropic|pass|✓)"; then
            echo "  $name: ✓"
        else
            echo "  $name: ✗ FAILED"
        fi
    done
fi
echo ""

echo "================================================================"
echo "VERIFICATION COMPLETE"
echo "================================================================"
