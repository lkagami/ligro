"""Builds GROMACS .mdp text for minimization runs.

Replaces the ~10 copies of `cat << EOF >| em*.mdp` template strings in
ligro.py. Defaults are the values the protein branch used, so a builder
call with no arguments reproduces the original steep-descent input.
"""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class MinimizationMdp:
    integrator: str = "steep"          # "steep" (SD) or "cg" (conjugate gradient)
    nsteps: int = 50000
    emtol: float = 1000.0
    emstep: float = 0.01
    energygrps: str = "system"         # "system" for em.mdp, "Protein" for em_real
    nstlist: int = 1
    ns_type: str = "grid"
    rlist: float = 1.0
    coulombtype: str = "PME"
    rcoulomb: float = 1.0
    rvdw: float = 1.0
    pbc: str = "xyz"
    cutoff_scheme: Optional[str] = None  # the ionization mdp sets "Verlet"
    title: str = "Minimization"
    extra: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.integrator not in ("steep", "cg"):
            raise ValueError("integrator must be 'steep' or 'cg', got {!r}".format(self.integrator))
        if self.nsteps < 0:
            raise ValueError("nsteps must be >= 0")


def render_minimization_mdp(p: MinimizationMdp) -> str:
    """Returns the .mdp file contents as a string."""
    lines = [
        "; LINES STARTING WITH ';' ARE COMMENTS",
        "title       = {}".format(p.title),
        "",
        "; Parameters describing what to do, when to stop and what to save",
        "integrator  = {}".format(p.integrator),
        "emtol       = {}".format(p.emtol),
        "emstep      = {}".format(p.emstep),
        "nsteps      = {}".format(p.nsteps),
        "energygrps  = {}".format(p.energygrps),
        "",
        "; Parameters describing how to find the neighbors of each atom and how to calculate the interactions",
        "nstlist     = {}".format(p.nstlist),
    ]
    if p.cutoff_scheme:
        lines.append("cutoff-scheme = {}".format(p.cutoff_scheme))
    lines += [
        "ns_type     = {}".format(p.ns_type),
        "rlist       = {}".format(p.rlist),
        "coulombtype = {}".format(p.coulombtype),
        "rcoulomb    = {}".format(p.rcoulomb),
        "rvdw        = {}".format(p.rvdw),
        "pbc         = {}".format(p.pbc),
    ]
    for key, value in p.extra.items():
        lines.append("{:<12}= {}".format(key, value))
    return "\n".join(lines) + "\n"
