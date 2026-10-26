"""Protein preparation + energy minimization.

Port of the protein branch of ligro.py (mount_simulation_prot and the
minimization block after it). Same file names, same order, same MDP
values -- but every GROMACS call goes through a Runner, and every file
write goes through Python instead of `cat << EOF`.

Usage:
    p = ProteinMinimization(ProteinMinimizationParams(...))
    p.plan()                              # list of step names, no side effects
    p.run(workdir, BiobbRunner())         # real run
    p.run(workdir, FakeRunner())          # tests
"""
import os
from dataclasses import dataclass, field
from typing import Callable, List, Optional

from ligro.core import config as C
from ligro.core.io.pdb import clean_protein_pdb, ion_option
from ligro.core.mdp.builder import MinimizationMdp, render_minimization_mdp


@dataclass
class ProteinMinimizationParams:
    force_field: str = "amber99sb-ildn"
    water: str = "tip3p"
    ignore_h: bool = False
    box_type: str = "cubic"
    box_distance: float = 1.0
    metal: Optional[str] = None            # e.g. "ZN"; None or "None" = no metal
    ion_mode: str = "Na (Number)"          # "Na (Number)", "Cl (Number)" or concentration
    ion_value: str = "0"
    min_algorithm: str = "SD Algorithm"    # "SD Algorithm" or "CG Algorithm"
    min_steps: int = 50000

    def __post_init__(self):
        if self.min_algorithm not in ("SD Algorithm", "CG Algorithm"):
            raise ValueError("min_algorithm must be 'SD Algorithm' or 'CG Algorithm'")
        if self.min_steps < 0:
            raise ValueError("min_steps must be >= 0")


@dataclass
class Context:
    workdir: str
    params: ProteinMinimizationParams

    def path(self, name):
        return os.path.join(self.workdir, name)


@dataclass
class Step:
    name: str
    action: Callable[["Context", object], None]


def _write(ctx, name, text):
    with open(ctx.path(name), "w", encoding="utf-8") as fh:
        fh.write(text)


def _steps() -> List[Step]:
    """Ordered steps. Each action gets (ctx, runner)."""

    def clean(ctx, r):
        clean_protein_pdb(ctx.path(C.PROTEIN_INPUT), ctx.path(C.PROTEIN_CLEAN),
                          metal=ctx.params.metal)

    def pdb2gmx(ctx, r):
        p = ctx.params
        r.pdb2gmx(ctx.path(C.PROTEIN_CLEAN), ctx.path(C.PROTEIN_GRO), ctx.path(C.TOPOLOGY),
                  force_field=p.force_field, water=p.water, ignore_h=p.ignore_h)

    def editconf(ctx, r):
        r.editconf(ctx.path(C.PROTEIN_GRO), ctx.path(C.BOX_GRO),
                   box_type=ctx.params.box_type, distance=ctx.params.box_distance)

    def solvate(ctx, r):
        r.solvate(ctx.path(C.BOX_GRO), ctx.path(C.TOPOLOGY),
                  ctx.path(C.SOLVATED_GRO), ctx.path(C.TOPOLOGY))

    def ionization_grompp(ctx, r):
        _write(ctx, "em.mdp", render_minimization_mdp(MinimizationMdp(
            integrator="steep", nsteps=50000, energygrps="system",
            cutoff_scheme="Verlet")))
        r.grompp(ctx.path("em.mdp"), ctx.path(C.SOLVATED_GRO), ctx.path(C.TOPOLOGY),
                 ctx.path(C.IONIZED_TPR))

    def genion(ctx, r):
        p = ctx.params
        r.genion(ctx.path(C.IONIZED_TPR), ctx.path(C.TOPOLOGY),
                 ctx.path(C.EM_INPUT_GRO), ctx.path(C.TOPOLOGY),
                 ion_opt=ion_option(p.ion_mode, p.ion_value))

    def steep(ctx, r):
        """SD branch / first pass of the CG branch."""
        _write(ctx, "em_real.mdp", render_minimization_mdp(MinimizationMdp(
            integrator="steep", nsteps=ctx.params.min_steps, energygrps="Protein")))
        r.grompp(ctx.path("em_real.mdp"), ctx.path(C.EM_INPUT_GRO), ctx.path(C.TOPOLOGY),
                 ctx.path(C.EM_DEFFNM + ".tpr"))
        r.mdrun(ctx.path(C.EM_DEFFNM + ".tpr"), ctx.path(C.EM_DEFFNM))

    def cg(ctx, r):
        """CG pass: continues from the steep result (em.gro)."""
        _write(ctx, "em_real.mdp", render_minimization_mdp(MinimizationMdp(
            integrator="cg", nsteps=ctx.params.min_steps, energygrps="Protein")))
        r.grompp(ctx.path("em_real.mdp"), ctx.path(C.EM_FINAL_GRO), ctx.path(C.TOPOLOGY),
                 ctx.path(C.EM_DEFFNM + ".tpr"))
        r.mdrun(ctx.path(C.EM_DEFFNM + ".tpr"), ctx.path(C.EM_DEFFNM))

    def postprocess(ctx, r):
        """Writes em.pdb from em.gro with editconf (no centering).
        Replaces the old trjconv/mv sequence; see biobb_gmx.gro_to_pdb."""
        r.gro_to_pdb(ctx.path(C.EM_FINAL_GRO), ctx.path(C.EM_FINAL_PDB))

    steps = [
        Step("clean_pdb", clean),
        Step("pdb2gmx", pdb2gmx),
        Step("editconf", editconf),
        Step("solvate", solvate),
        Step("ionization_grompp", ionization_grompp),
        Step("genion", genion),
        Step("steep_minimization", steep),
    ]
    return steps, cg, postprocess


class ProteinMinimization:
    def __init__(self, params: ProteinMinimizationParams):
        self.params = params

    def _all_steps(self) -> List[Step]:
        steps, cg, postprocess = _steps()
        if self.params.min_algorithm == "CG Algorithm":
            steps.append(Step("cg_minimization", cg))
        steps.append(Step("postprocess", postprocess))
        return steps

    def plan(self) -> List[str]:
        """Names of the steps, in order. Does not touch the filesystem."""
        return [s.name for s in self._all_steps()]

    def run(self, workdir, runner):
        """Runs every step. Stops at the first exception; the exception is
        re-raised with the step name attached so the GUI can report it."""
        os.makedirs(workdir, exist_ok=True)
        ctx = Context(workdir=workdir, params=self.params)
        for step in self._all_steps():
            try:
                step.action(ctx, runner)
            except Exception as error:
                raise RuntimeError("minimization step '{}' failed: {}".format(step.name, error)) from error
        return ctx.path(C.EM_FINAL_PDB)
