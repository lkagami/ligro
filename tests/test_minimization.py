"""Tests for the core minimization pipeline. No GROMACS, no biobb, no Qt."""
import os
import shutil
import sys
import tempfile
import unittest

from ligro.core.io.pdb import clean_protein_pdb, ion_option
from ligro.core.mdp.builder import MinimizationMdp, render_minimization_mdp
from ligro.core.pipelines.minimization import (
    ProteinMinimization,
    ProteinMinimizationParams,
)
from ligro.core.runners import FakeRunner


def _read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class MdpBuilderTests(unittest.TestCase):
    def test_defaults_match_original_protein_em_mdp(self):
        text = render_minimization_mdp(MinimizationMdp())
        self.assertIn("integrator  = steep", text)
        self.assertIn("nsteps      = 50000", text)
        self.assertIn("energygrps  = system", text)
        self.assertNotIn("cutoff-scheme", text)

    def test_cutoff_scheme_only_when_requested(self):
        text = render_minimization_mdp(MinimizationMdp(cutoff_scheme="Verlet"))
        self.assertIn("cutoff-scheme = Verlet", text)

    def test_invalid_integrator_rejected(self):
        with self.assertRaises(ValueError):
            MinimizationMdp(integrator="leapfrog")


class PdbHelperTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.src = os.path.join(self.dir, "in.pdb")
        with open(self.src, "w") as fh:
            fh.write("ATOM      1  N   ALA A   1\n"
                     "HETATM    2 ZN    ZN A 200\n"
                     "HETATM    3  O   HOH A 300\n"
                     "CONECT    1    2\n")
        self.dst = os.path.join(self.dir, "out.pdb")

    def tearDown(self):
        shutil.rmtree(self.dir)

    def read(self):
        with open(self.dst) as fh:
            return fh.read()

    def test_keeps_only_atoms_without_metal(self):
        n = clean_protein_pdb(self.src, self.dst, metal=None)
        self.assertEqual(n, 1)
        self.assertIn("ATOM", self.read())
        self.assertNotIn("ZN", self.read())

    def test_keeps_metal_line_when_requested(self):
        clean_protein_pdb(self.src, self.dst, metal="ZN")
        self.assertIn("ZN A 200", self.read())

    def test_string_none_means_no_metal(self):
        clean_protein_pdb(self.src, self.dst, metal="None")
        self.assertNotIn("ZN", self.read())

    def test_ion_option_mapping(self):
        self.assertEqual(ion_option("Na (Number)", "5"), ("num_pos", 5))
        self.assertEqual(ion_option("Cl (Number)", "3"), ("num_neg", 3))
        self.assertEqual(ion_option("Concentration", "0.15"), ("conc", 0.15))

    def test_ion_option_rejects_bad_number(self):
        with self.assertRaises(ValueError):
            ion_option("Na (Number)", "abc")


class PlanTests(unittest.TestCase):
    def test_sd_plan_has_no_cg_step(self):
        plan = ProteinMinimization(ProteinMinimizationParams(min_algorithm="SD Algorithm")).plan()
        self.assertNotIn("cg_minimization", plan)
        self.assertEqual(plan[-1], "postprocess")

    def test_cg_plan_runs_cg_before_postprocess(self):
        plan = ProteinMinimization(ProteinMinimizationParams(min_algorithm="CG Algorithm")).plan()
        self.assertEqual(plan[-2:], ["cg_minimization", "postprocess"])

    def test_plan_has_no_side_effects(self):
        p = ProteinMinimization(ProteinMinimizationParams())
        p.plan()  # must not need a workdir or a runner


class RunTests(unittest.TestCase):
    def setUp(self):
        self.wd = tempfile.mkdtemp()
        with open(os.path.join(self.wd, "protein.pdb"), "w") as fh:
            fh.write("ATOM      1  N   ALA A   1\nHETATM    2  O   HOH A 300\n")

    def tearDown(self):
        shutil.rmtree(self.wd)

    def test_sd_run_call_order_and_arguments(self):
        r = FakeRunner()
        params = ProteinMinimizationParams(ion_mode="Na (Number)", ion_value="5",
                                           min_algorithm="SD Algorithm", min_steps=1234)
        ProteinMinimization(params).run(self.wd, r)

        self.assertEqual(r.names(), [
            "pdb2gmx", "editconf", "solvate",
            "grompp",            # ionization tpr from em.mdp
            "genion",
            "grompp", "mdrun",   # steep minimization
            "gro_to_pdb",
        ])
        genion_kwargs = [k for n, _, k in r.calls if n == "genion"][0]
        self.assertEqual(genion_kwargs, {"ion_opt": ("num_pos", 5)})
        self.assertNotIn("trjconv", r.names())

        em_real = _read(os.path.join(self.wd, "em_real.mdp"))
        self.assertIn("integrator  = steep", em_real)
        self.assertIn("nsteps      = 1234", em_real)
        self.assertIn("energygrps  = Protein", em_real)

    def test_cg_run_runs_steep_then_cg_from_em_gro(self):
        r = FakeRunner()
        ProteinMinimization(ProteinMinimizationParams(min_algorithm="CG Algorithm")).run(self.wd, r)
        names = r.names()
        self.assertEqual(names.count("mdrun"), 2)
        grompp_calls = [(a, k) for n, a, k in r.calls if n == "grompp"]
        # second minimization grompp must start from em.gro
        self.assertTrue(grompp_calls[-1][0][1].endswith("em.gro"))
        em_real = _read(os.path.join(self.wd, "em_real.mdp"))
        self.assertIn("integrator  = cg", em_real)

    def test_failure_names_the_step(self):
        r = FakeRunner(fail_on="genion")
        with self.assertRaises(RuntimeError) as ctx:
            ProteinMinimization(ProteinMinimizationParams()).run(self.wd, r)
        self.assertIn("genion", str(ctx.exception))

    def test_core_does_not_import_gui_toolkits(self):
        # Checked in a fresh interpreter: other test modules stub PyQt5 in
        # this process, which would make an in-process check meaningless.
        import subprocess
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        code = ("import sys; import ligro.core.pipelines.minimization, ligro.core.runners;"
                "bad=[m for m in ('PyQt5','gi') if m in sys.modules];"
                "print(bad); sys.exit(1 if bad else 0)")
        result = subprocess.run([sys.executable, "-c", code], cwd=root,
                                env=dict(os.environ, PYTHONPATH=root),
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class ParamValidationTests(unittest.TestCase):
    def test_bad_algorithm_rejected(self):
        with self.assertRaises(ValueError):
            ProteinMinimizationParams(min_algorithm="Newton")


if __name__ == "__main__":
    unittest.main()
