"""Tests for the BioExcel-coordinated NVT/NPT/production workflow."""
import os
import shutil
import tempfile
import unittest
import zipfile

from ligro.core.pipelines.dynamics import DynamicsParams, MolecularDynamics


class RecordingRunner:
    def __init__(self):
        self.calls = []

    def grompp(self, mdp, gro, top_zip, tpr, ref=None):
        self.calls.append(("grompp", os.path.basename(tpr), ref))
        with open(tpr, "w", encoding="utf-8") as fh:
            fh.write("fake tpr")

    def mdrun(self, tpr, deffnm, write_xtc=False, write_cpt=False):
        self.calls.append(("mdrun", os.path.basename(deffnm), write_xtc, write_cpt))
        for suffix in (".gro", ".edr", ".log", ".trr"):
            with open(deffnm + suffix, "w", encoding="utf-8") as fh:
                fh.write("fake output")
        if write_cpt:
            with open(deffnm + ".cpt", "w", encoding="utf-8") as fh:
                fh.write("fake checkpoint")
        if write_xtc:
            with open(deffnm + ".xtc", "w", encoding="utf-8") as fh:
                fh.write("fake trajectory")


class DynamicsTests(unittest.TestCase):
    def setUp(self):
        self.wd = tempfile.mkdtemp()
        for name in ("trp.top", "nvt.mdp", "npt.mdp", "md.mdp", "em.gro"):
            with open(os.path.join(self.wd, name), "w", encoding="utf-8") as fh:
                fh.write("test\n")
        with zipfile.ZipFile(os.path.join(self.wd, "trp.zip"), "w") as z:
            z.writestr("p2g.top", "[ system ]\n")
            z.writestr("amber99sb-ildn.ff/forcefield.itp", "; force field\n")

    def tearDown(self):
        shutil.rmtree(self.wd)

    def test_runs_nvt_npt_then_production_in_order(self):
        runner = RecordingRunner()
        result = MolecularDynamics().run(self.wd, runner)
        self.assertEqual([call[:2] for call in runner.calls], [
            ("grompp", "nvt.tpr"), ("mdrun", "nvt"),
            ("grompp", "npt.tpr"), ("mdrun", "npt"),
            ("grompp", "md.tpr"), ("mdrun", "md"),
        ])
        self.assertTrue(os.path.isfile(result["xtc"]))
        for key in ("gro", "edr", "trr", "log", "cpt"):
            self.assertTrue(os.path.isfile(result[key]), key)
        with zipfile.ZipFile(os.path.join(self.wd, "trp_runtime.zip")) as z:
            self.assertIn("p2g.top", z.namelist())
            self.assertIn("amber99sb-ildn.ff/forcefield.itp", z.namelist())
            self.assertIn("test", z.read("p2g.top").decode())

    def test_can_prepare_tpr_without_production_md(self):
        runner = RecordingRunner()
        result = MolecularDynamics(DynamicsParams(run_production=False)).run(self.wd, runner)
        self.assertTrue(os.path.isfile(result["tpr"]))
        self.assertIsNone(result["xtc"])
        self.assertEqual([call[0] for call in runner.calls],
                         ["grompp", "mdrun", "grompp", "mdrun", "grompp"])

    def test_stops_if_required_input_is_missing(self):
        os.remove(os.path.join(self.wd, "md.mdp"))
        with self.assertRaisesRegex(RuntimeError, "md.mdp"):
            MolecularDynamics().run(self.wd, RecordingRunner())


if __name__ == "__main__":
    unittest.main()
