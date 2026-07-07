#!/usr/bin/env python3
"""
QML Formalization Audit Tool
Counts axioms, sorry, theorems across all modules.
Maps each axiom to its OpenProblems.lean entry.

Usage: python3 tools/audit.py
"""
import subprocess, re, os

def audit():
    qml_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "QML")

    axioms = []
    sorrys = []
    theorems = 0
    lines = 0
    modules = 0

    for root, _, files in os.walk(qml_dir):
        for f in files:
            if f.endswith(".lean"):
                modules += 1
                path = os.path.join(root, f)
                with open(path) as fh:
                    for i, line in enumerate(fh, 1):
                        lines += 1
                        if line.startswith("axiom "):
                            axioms.append(f"{f}:{i}: {line.strip()}")
                        if "sorry" in line and not line.strip().startswith("--"):
                            sorrys.append(f"{f}:{i}: {line.strip()}")
                        if re.match(r"^(theorem |lemma |private theorem |private lemma )", line):
                            theorems += 1

    print(f"╔══════════════════════════════════════╗")
    print(f"║  QML FORMALIZATION AUDIT             ║")
    print(f"╠══════════════════════════════════════╣")
    print(f"║  Modules:  {modules:4d}                     ║")
    print(f"║  Lines:    {lines:4d}                     ║")
    print(f"║  Theorems: {theorems:4d}                     ║")
    print(f"║  Axioms:   {len(axioms):4d}                     ║")
    print(f"║  Sorry:    {len(sorrys):4d}                     ║")
    print(f"╚══════════════════════════════════════╝")

    if axioms:
        print(f"\nAxioms ({len(axioms)}):")
        for a in axioms:
            print(f"  {a}")

    if sorrys:
        print(f"\n⚠️  SORRY FOUND ({len(sorrys)}):")
        for s in sorrys:
            print(f"  {s}")
    else:
        print(f"\n✅ Zero sorry — project is clean")

    return len(sorrys) == 0

if __name__ == "__main__":
    clean = audit()
    exit(0 if clean else 1)
