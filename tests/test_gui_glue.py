"""Tests the glue between ligro.py's widgets and the core pipeline.

PyQt5 and the other GUI/analysis imports are replaced by stubs, so the
window class can be imported without a display. The widgets are fakes with
fixed values. The runner is a simulated GROMACS that writes the files the
pipeline expects. This checks the wiring (parameters, working directory,
topology restoration, error handling), not the Qt interface itself.
"""
import importlib
import os
import shutil
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

for name in ("PyQt5", "PyQt5.QtCore", "PyQt5.QtGui", "PyQt5.QtWidgets",
             "matplotlib", "matplotlib.pyplot", "pylab", "biopandas",
             "biopandas.mol2", "plip", "plip.modules", "acpype_lib"):
    sys.modules.setdefault(name, mock.MagicMock())

sys.modules["plip.modules"].report.__version__ = "2.4.0"   # read by ligro.py at import
ui = importlib.import_module("ligro.ligro")

from ligro.core.io.gro import gro_to_pdb  # noqa: E402
from ligro.core.runners import FakeRunner  # noqa: E402

GRO_ONE_ATOM = "t\n    1\n    1ALA     CA    1   1.000   1.000   1.000\n   3.00000   3.00000   3.00000\n"


class Widget:
    def __init__(self, value):
        self.value = value

    def currentText(self):
        return self.value

    def text(self):
        return self.value

    def isChecked(self):
        return True


class SimulatedGromacs(FakeRunner):
    """FakeRunner that writes the files each step is expected to produce."""

    def pdb2gmx(self, pdb, gro, top_zip, *a, **k):
        super().pdb2gmx(pdb, gro, top_zip, *a, **k)
        with zipfile.ZipFile(top_zip, "w") as z:
            z.writestr("p2g.top", "[ system ]\n")
            z.writestr("posre.itp", "[ position_restraints ]\n")

    def mdrun(self, tpr, deffnm):
        super().mdrun(tpr, deffnm)
        with open(deffnm + ".gro", "w") as fh:
            fh.write(GRO_ONE_ATOM)

    def gro_to_pdb(self, gro_in, pdb_out):
        super().gro_to_pdb(gro_in, pdb_out)
        return gro_to_pdb(gro_in, pdb_out)


class GromacsDetectionTests(unittest.TestCase):
    """gromacs_flag must only accept the GROMACS driver, never biobb's
    bare-name command-line tools (pdb2gmx, mdrun, ...)."""

    def test_bare_biobb_names_are_not_gromacs(self):
        for name in ("mdrun", "pdb2gmx", "grompp", "genion", "editconf"):
            self.assertFalse(ui.gromacs_flag(name), name)

    def test_gmx_is_probed_with_gmx_version(self):
        with mock.patch.object(ui.subprocess, "run") as run:
            run.return_value = mock.Mock(returncode=0)
            self.assertTrue(ui.gromacs_flag("gmx"))
            self.assertEqual(run.call_args[0][0], ["gmx", "--version"])

    def test_missing_gmx_returns_false(self):
        with mock.patch.object(ui.subprocess, "run", side_effect=FileNotFoundError):
            self.assertFalse(ui.gromacs_flag("gmx"))


class GlueTests(unittest.TestCase):
    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        shutil.copy(os.path.join(ROOT, "tests", "data", "trpcage_protein.pdb"),
                    os.path.join(self.cwd, "protein.pdb"))
        self.old_cwd = os.getcwd()
        os.chdir(self.cwd)
        self.dialogs = []
        self.runner = SimulatedGromacs()

    def tearDown(self):
        os.chdir(self.old_cwd)
        shutil.rmtree(self.cwd)

    def window(self, min_alg="SD Algorithm", ioniz="Concentration (mol/liter)", conc="0"):
        w = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        w.sel_metal = Widget("None")
        w.sel_ff = Widget("amber99sb-ildn")
        w.sel_water = Widget("tip3p")
        w.ig_hyd = Widget(True)
        w.comboBox_5 = Widget("cubic")
        w.sel_dist = Widget("1.0")
        w.sel_ioniz = Widget(ioniz)
        w.sel_conc = Widget(conc)
        w.sel_min_alg = Widget(min_alg)
        w.sel_min_step = Widget("100")
        return w

    def run_glue(self, window):
        def fake_dialog(title, text):
            self.dialogs.append((title, text))
        with mock.patch.object(ui, "BiobbRunner", lambda binary="gmx": self.runner), \
             mock.patch.object(ui, "showdialog", fake_dialog):
            return window._run_protein_minimization()

    def test_success_writes_topology_and_returns_true(self):
        ok = self.run_glue(self.window())
        self.assertTrue(ok, self.dialogs)
        self.assertEqual(self.dialogs, [])
        self.assertTrue(os.path.exists(os.path.join(self.cwd, "trp.top")))
        self.assertTrue(os.path.exists(os.path.join(self.cwd, "posre.itp")))
        self.assertTrue(os.path.exists(os.path.join(self.cwd, "em.pdb")))

    def test_widget_values_reach_the_pipeline(self):
        self.run_glue(self.window())
        pdb2gmx = [k for n, _, k in self.runner.calls if n == "pdb2gmx"][0]
        self.assertEqual(pdb2gmx, {"force_field": "amber99sb-ildn", "water": "tip3p", "ignore_h": True})
        genion = [k for n, _, k in self.runner.calls if n == "genion"][0]
        self.assertEqual(genion, {"ion_opt": ("conc", 0.0)})

    def test_sd_plus_cg_runs_the_second_minimization(self):
        self.run_glue(self.window(min_alg="SD + CG"))
        self.assertEqual(self.runner.names().count("mdrun"), 2)

    def test_sd_only_runs_one_minimization(self):
        self.run_glue(self.window(min_alg="SD Algorithm"))
        self.assertEqual(self.runner.names().count("mdrun"), 1)

    def test_sodium_count_maps_to_num_pos(self):
        self.run_glue(self.window(ioniz="Na (Number)", conc="3"))
        genion = [k for n, _, k in self.runner.calls if n == "genion"][0]
        self.assertEqual(genion, {"ion_opt": ("num_pos", 3)})

    def test_pipeline_failure_is_shown_and_returns_false(self):
        self.runner.fail_on = "genion"
        ok = self.run_glue(self.window())
        self.assertFalse(ok)
        self.assertEqual(len(self.dialogs), 1)
        self.assertEqual(self.dialogs[0][0], "Minimization failed")
        self.assertIn("genion", self.dialogs[0][1])

    def test_bad_numeric_input_is_reported_not_raised(self):
        w = self.window()
        w.sel_min_step = Widget("many")
        ok = self.run_glue(w)
        self.assertFalse(ok)
        self.assertEqual(self.dialogs[0][0], "Minimization failed")


if __name__ == "__main__":
    unittest.main()
