#!/usr/bin/env python3
"""End-to-end smoke test of the protein minimization pipeline with REAL
GROMACS (through biobb). Not part of the unit tests: it needs `gmx` on PATH
and the biobb packages installed.

    python scripts/smoke_minimization.py [--keep DIR]

Uses tests/data/trpcage_protein.pdb (20-residue trp-cage, first NMR model).
Exit code 0 only if the minimization converged according to GROMACS' log.
"""
import argparse
import os
import re
import shutil
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from ligro.core.pipelines.minimization import (  # noqa: E402
    ProteinMinimization,
    ProteinMinimizationParams,
)
from ligro.core.runners import BiobbRunner  # noqa: E402

INPUT = os.path.join(ROOT, "tests", "data", "trpcage_protein.pdb")


def compare_gro_pdb(gro_path, pdb_path, tol_angstrom=0.01):
    """Coordinates in em.pdb must equal em.gro (nm -> A) within tolerance.
    This is what would have caught the editconf recentring shift."""
    from ligro.core.io.gro import read_gro
    _, atoms, _ = read_gro(gro_path)
    pdb = [l for l in open(pdb_path) if l.startswith("ATOM")]
    if len(pdb) != len(atoms):
        return ["atom count differs: gro={} pdb={}".format(len(atoms), len(pdb))]
    worst = 0.0
    for atom, line in zip(atoms, pdb):
        for key, col in (("x", 30), ("y", 38), ("z", 46)):
            diff = abs(atom[key] * 10.0 - float(line[col:col + 8]))
            worst = max(worst, diff)
    print("max coordinate deviation gro vs pdb: {:.4f} A".format(worst))
    if worst > tol_angstrom:
        return ["em.pdb coordinates differ from em.gro by {:.3f} A".format(worst)]
    return []


def check(workdir):
    problems = []
    gro = os.path.join(workdir, "em.gro")
    pdb_path = os.path.join(workdir, "em.pdb")
    if os.path.exists(gro) and os.path.exists(pdb_path):
        problems += compare_gro_pdb(gro, pdb_path)
    for name in ("em.gro", "em.pdb", "em.log", "em.tpr"):
        path = os.path.join(workdir, name)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            problems.append("missing or empty: " + name)
    log = os.path.join(workdir, "em.log")
    if os.path.exists(log):
        text = open(log, encoding="utf-8", errors="replace").read()
        if not re.search(r"converged to Fmax", text) and "Steepest Descents converged" not in text:
            problems.append("em.log does not report convergence")
        m = re.search(r"Potential Energy\s*=\s*([-\d.eE+]+)", text)
        if m:
            print("final potential energy:", m.group(1), "kJ/mol")
    pdb = os.path.join(workdir, "em.pdb")
    if os.path.exists(pdb):
        atoms = [l for l in open(pdb) if l.startswith("ATOM") or l.startswith("HETATM")]
        protein = [l for l in atoms if l[17:20].strip() not in ("SOL", "HOH", "WAT", "NA", "CL", "Na+", "Cl-")]
        residues = {l[22:26] for l in protein}
        print("em.pdb atoms:", len(atoms), "| protein residues:", len(residues))
        if len(residues) != 20:
            problems.append("expected 20 protein residues in em.pdb, found {}".format(len(residues)))
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", help="directory to keep the run in (default: temp dir, removed)")
    ap.add_argument("--steps", type=int, default=5000, help="minimization steps (small for a smoke test)")
    args = ap.parse_args()

    workdir = args.keep or tempfile.mkdtemp(prefix="ligro_smoke_")
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(INPUT, os.path.join(workdir, "protein.pdb"))

    params = ProteinMinimizationParams(
        ignore_h=True,
        ion_mode="Concentration", ion_value="0",   # neutralize only
        min_algorithm="SD Algorithm", min_steps=args.steps,
    )
    print("plan:", " -> ".join(ProteinMinimization(params).plan()))
    print("workdir:", workdir)

    try:
        ProteinMinimization(params).run(workdir, BiobbRunner(binary="gmx"))
    except Exception as error:
        print("PIPELINE FAILED:", error)
        return 2

    problems = check(workdir)
    if problems:
        print("CHECKS FAILED:")
        for p in problems:
            print("  -", p)
        return 1
    print("OK: minimization finished and checks passed")
    if not args.keep:
        shutil.rmtree(workdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
