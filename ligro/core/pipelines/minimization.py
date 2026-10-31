"""Protein preparation + energy minimization.

Port of the protein branch of main.py (mount_simulation_prot and the
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
import shutil
import tempfile
import zipfile
from dataclasses import dataclass, field
from typing import Callable, List, Optional

from ligro.core import config as C
from ligro.core.io.pdb import clean_protein_pdb, ion_option
from ligro.core.mdp.builder import MinimizationMdp, render_minimization_mdp


MARTINI_3_FORCE_FIELD = "Martini 3 (CG)"


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
    protein_input: Optional[str] = None

    def __post_init__(self):
        if self.min_algorithm not in ("SD Algorithm", "CG Algorithm"):
            raise ValueError("min_algorithm must be 'SD Algorithm' or 'CG Algorithm'")
        if self.min_steps < 0:
            raise ValueError("min_steps must be >= 0")


@dataclass
class Context:
    workdir: str
    params: ProteinMinimizationParams
    scratch_dir: Optional[str] = None

    def path(self, name):
        return os.path.join(self.workdir, name)

    def scratch_path(self, name):
        if self.scratch_dir is None:
            raise RuntimeError("temporary work directory is unavailable")
        return os.path.join(self.scratch_dir, name)


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


def _martini_steps(min_algorithm) -> List[Step]:
    def clean(ctx, _runner):
        protein_input = ctx.params.protein_input or ctx.path(C.PROTEIN_INPUT)
        if not os.path.isfile(protein_input):
            raise FileNotFoundError(
                "Protein PDB input was not found: {}. Select a PDB file or place "
                "protein.pdb in the project directory.".format(protein_input)
            )
        clean_protein_pdb(protein_input, ctx.path(C.PROTEIN_CLEAN))

    def martinize(ctx, runner):
        coarse_structure = ctx.scratch_path("protein_cg.pdb")
        runner.martinize2(
            ctx.path(C.PROTEIN_CLEAN),
            coarse_structure,
            ctx.scratch_path("protein.top"),
            force_field="martini3001",
        )
        shutil.copyfile(coarse_structure, ctx.path(C.PROTEIN_GRO))

    def box(ctx, runner):
        runner.editconf(
            ctx.scratch_path("protein_cg.pdb"),
            ctx.path(C.BOX_GRO),
            box_type=ctx.params.box_type,
            distance=ctx.params.box_distance,
        )

    def steep_minimize(ctx, runner):
        _write(ctx, "em_real.mdp", _render_martini_minimization_mdp(
            ctx.params.min_steps
        ))
        topology_archive = _archive_martini_topology(
            ctx.scratch_dir,
            ctx.scratch_path("protein.top"),
            ctx.path(C.TOPOLOGY),
        )
        tpr_path = ctx.path(C.EM_DEFFNM + ".tpr")
        runner.grompp(
            ctx.path("em_real.mdp"),
            ctx.path(C.BOX_GRO),
            topology_archive,
            tpr_path,
        )
        runner.mdrun(tpr_path, ctx.path(C.EM_DEFFNM))

    def cg_minimize(ctx, runner):
        _write(ctx, "em_cg.mdp", _render_martini_minimization_mdp(
            ctx.params.min_steps, integrator="cg"
        ))
        runner.grompp(
            ctx.path("em_cg.mdp"),
            ctx.path(C.EM_FINAL_GRO),
            ctx.path(C.TOPOLOGY),
            ctx.path(C.EM_DEFFNM + ".tpr"),
        )
        runner.mdrun(
            ctx.path(C.EM_DEFFNM + ".tpr"),
            ctx.path(C.EM_DEFFNM),
        )

    def postprocess(ctx, runner):
        runner.gro_to_pdb(ctx.path(C.EM_FINAL_GRO), ctx.path(C.EM_FINAL_PDB))

    steps = [
        Step("clean_pdb", clean),
        Step("martinize2", martinize),
        Step("editconf", box),
        Step("steep_minimization", steep_minimize),
    ]
    if min_algorithm == "CG Algorithm":
        steps.append(Step("cg_minimization", cg_minimize))
    steps.extend([
        Step("postprocess", postprocess),
    ])
    return steps


def _render_martini_minimization_mdp(nsteps: int, integrator="steep") -> str:
    return (
        "; Martini 3 coarse-grained energy minimization\n"
        "integrator = {}\n"
        "emtol = 10.0\n"
        "emstep = 0.01\n"
        "nsteps = {}\n"
        "cutoff-scheme = Verlet\n"
        "nstlist = 20\n"
        "coulombtype = Reaction-Field\n"
        "rcoulomb = 1.1\n"
        "rvdw = 1.1\n"
        "epsilon_r = 15\n"
        "pbc = xyz\n"
    ).format(integrator, nsteps)


class ProteinMinimization:
    def __init__(self, params: ProteinMinimizationParams):
        self.params = params

    def _all_steps(self) -> List[Step]:
        if self.params.force_field == MARTINI_3_FORCE_FIELD:
            return _martini_steps(self.params.min_algorithm)
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
        if self.params.force_field == MARTINI_3_FORCE_FIELD:
            with tempfile.TemporaryDirectory(
                prefix=".martini3-", dir=workdir
            ) as scratch_dir:
                ctx = Context(
                    workdir=workdir,
                    params=self.params,
                    scratch_dir=scratch_dir,
                )
                self._run_steps(ctx, runner)
        else:
            ctx = Context(workdir=workdir, params=self.params)
            self._run_steps(ctx, runner)
        return ctx.path(C.EM_FINAL_PDB)

    def _run_steps(self, ctx, runner):
        for step in self._all_steps():
            try:
                step.action(ctx, runner)
            except Exception as error:
                raise RuntimeError("minimization step '{}' failed: {}".format(step.name, error)) from error


def _archive_martini_topology(scratch_dir, top_path, output_path):
    if not os.path.isfile(top_path):
        raise RuntimeError("martinize2 topology was not produced: " + top_path)

    with zipfile.ZipFile(
        output_path, "w", compression=zipfile.ZIP_DEFLATED
    ) as archive:
        archive.write(top_path, "trp.top")
        for current, _dirs, files in os.walk(scratch_dir):
            for name in files:
                if not name.endswith(".itp"):
                    continue
                full_path = os.path.join(current, name)
                relative_path = os.path.relpath(full_path, scratch_dir)
                archive.write(
                    full_path,
                    relative_path.replace(os.sep, "/"),
                )
    return output_path
