"""The Runner interface: the only thing the pipelines use to talk to GROMACS.

Two implementations:
  BiobbRunner -- production; delegates to ligro.biobb_gmx (biobb_gromacs).
  FakeRunner  -- tests; records calls and writes nothing.

Paths are strings. Topologies are zipped (biobb convention).
"""
from typing import Protocol


class Runner(Protocol):
    def pdb2gmx(self, pdb, gro, top_zip, force_field, water, ignore_h): ...
    def editconf(self, gro_in, gro_out, box_type, distance): ...
    def solvate(self, gro_in, top_zip_in, gro_out, top_zip_out): ...
    def grompp(self, mdp, gro, top_zip, tpr, ref=None): ...
    def genion(self, tpr, top_zip, gro_out, top_zip_out, ion_opt): ...
    def mdrun(self, tpr, deffnm, write_xtc=False, write_cpt=False): ...
    def gro_to_pdb(self, gro_in, pdb_out): ...


class BiobbRunner:
    """Adapter from the Runner interface to ligro.biobb_gmx."""

    def __init__(self, binary="gmx"):
        self.binary = binary

    def _bg(self):
        from ligro import biobb_gmx  # lazy: biobb is only needed here
        return biobb_gmx

    def pdb2gmx(self, pdb, gro, top_zip, force_field, water, ignore_h):
        return self._bg().pdb2gmx(pdb, gro, top_zip, force_field=force_field,
                                  water=water, binary=self.binary, ignore_h=ignore_h)

    def editconf(self, gro_in, gro_out, box_type, distance):
        return self._bg().editconf(gro_in, gro_out, box_type=box_type,
                                   distance=distance, binary=self.binary)

    def solvate(self, gro_in, top_zip_in, gro_out, top_zip_out):
        return self._bg().solvate(gro_in, top_zip_in, gro_out, top_zip_out,
                                  binary=self.binary)

    def grompp(self, mdp, gro, top_zip, tpr, ref=None):
        return self._bg().grompp(mdp, gro, top_zip, tpr, ref=ref, binary=self.binary)

    def genion(self, tpr, top_zip, gro_out, top_zip_out, ion_opt):
        key, value = ion_opt   # 'num_pos' -> np, 'num_neg' -> nn, 'conc' -> concentration
        biobb_key = {"num_pos": "np", "num_neg": "nn", "conc": "concentration"}[key]
        return self._bg().genion(tpr, gro_out, top_zip, top_zip_out,
                                 binary=self.binary, **{biobb_key: value})

    def mdrun(self, tpr, deffnm, write_xtc=False, write_cpt=False):
        return self._bg().mdrun(tpr, deffnm, binary=self.binary,
                                 write_xtc=write_xtc, write_cpt=write_cpt)

    def gro_to_pdb(self, gro_in, pdb_out):
        # Pure Python on purpose: editconf recentres, which shifts coordinates.
        from ligro.core.io.gro import gro_to_pdb as _convert
        return _convert(gro_in, pdb_out)


class FakeRunner:
    """Records every call as (method_name, args, kwargs) in `self.calls`.
    `fail_on` names a method that raises, to test error reporting."""

    def __init__(self, fail_on=None):
        self.calls = []
        self.fail_on = fail_on

    def _record(self, name, args, kwargs):
        self.calls.append((name, args, kwargs))
        if self.fail_on == name:
            raise RuntimeError("FakeRunner: simulated failure in " + name)

    def pdb2gmx(self, *a, **k):  self._record("pdb2gmx", a, k)
    def editconf(self, *a, **k): self._record("editconf", a, k)
    def solvate(self, *a, **k):  self._record("solvate", a, k)
    def grompp(self, *a, **k):   self._record("grompp", a, k)
    def genion(self, *a, **k):   self._record("genion", a, k)
    def mdrun(self, *a, **k):    self._record("mdrun", a, k)
    def gro_to_pdb(self, *a, **k): self._record("gro_to_pdb", a, k)

    def names(self):
        return [name for name, _, _ in self.calls]
