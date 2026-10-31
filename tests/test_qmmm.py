import os
import tempfile
import unittest
import zipfile
from unittest import mock

from ligro.core.pipelines.qmmm import (
    QmmmParams,
    QmmmSimulation,
    render_qmmm_awh_mdp,
)
from ligro.core.runners import BiobbRunner, FakeRunner
from ligro import biobb_gmx


class QmmmParamsTests(unittest.TestCase):
    def test_rejects_invalid_qm_group_and_awh_interval(self):
        with self.assertRaisesRegex(ValueError, "QM index group"):
            QmmmParams(
                "system.gro", "topology.zip", "index.ndx", qm_group="QM atoms",
                awh_group1="QM", awh_group2="Ligand",
            )
        with self.assertRaisesRegex(ValueError, "lower than"):
            QmmmParams(
                "system.gro", "topology.zip", "index.ndx",
                awh_group1="QM", awh_group2="Ligand",
                awh_start_nm=1.0, awh_end_nm=0.5,
            )

    def test_renders_cp2k_and_awh_settings(self):
        params = QmmmParams(
            "system.gro",
            "topology.zip",
            "index.ndx",
            qm_group="ReactiveRegion",
            qm_charge=-1,
            qm_multiplicity=2,
            qm_method="BLYP",
            awh_group1="ReactiveRegion",
            awh_group2="Substrate",
            awh_start_nm=0.15,
            awh_end_nm=0.45,
        )

        mdp = render_qmmm_awh_mdp(params)

        self.assertIn("qmmm-cp2k-active = yes", mdp)
        self.assertIn("qmmm-cp2k-qmgroup = ReactiveRegion", mdp)
        self.assertIn("qmmm-cp2k-qmcharge = -1", mdp)
        self.assertIn("qmmm-cp2k-qmmultiplicity = 2", mdp)
        self.assertIn("qmmm-cp2k-qmmethod = BLYP", mdp)
        self.assertIn("pull-group2-name = Substrate", mdp)
        self.assertIn("pull-coord1-type = external-potential", mdp)
        self.assertIn("pull-coord1-potential-provider = awh", mdp)
        self.assertIn("awh1-dim1-coord-provider = pull", mdp)
        self.assertIn("awh1-dim1-diffusion = 1e-05", mdp)
        self.assertIn("nstxout-compressed = 100", mdp)
        self.assertIn("awh = yes", mdp)


class QmmmPipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.workdir = self.temp_dir.name
        self.structure = os.path.join(self.workdir, "system.gro")
        self.topology = os.path.join(self.workdir, "topology.zip")
        self.index = os.path.join(self.workdir, "index.ndx")
        with open(self.structure, "w", encoding="ascii") as output:
            output.write(
                "equilibrated system\n"
                "4\n"
                "    1ALA      CA    1   0.100   0.100   0.100\n"
                "    1ALA      CB    2   0.200   0.100   0.100\n"
                "    2LIG      C1    3   0.300   0.100   0.100\n"
                "    2LIG      C2    4   0.400   0.100   0.100\n"
                "   2.00000   2.00000   2.00000\n"
            )
        with zipfile.ZipFile(self.topology, "w") as archive:
            archive.writestr("topol.top", "[ system ]\nTest\n")
        with open(self.index, "w", encoding="ascii") as output:
            output.write(
                "[ QMatoms ]\n1 2\n[ Protein ]\n1 2 3\n"
                "[ Substrate ]\n3 4\n"
            )

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_runs_grompp_with_index_and_mdrun_and_checks_outputs(self):
        class OutputRunner(FakeRunner):
            def grompp(self, mdp, gro, top_zip, tpr, ref=None, ndx=None):
                super().grompp(mdp, gro, top_zip, tpr, ref=ref, ndx=ndx)
                with open(tpr, "w", encoding="ascii") as output:
                    output.write("run input")

            def mdrun(self, tpr, deffnm, write_xtc=False, write_cpt=False):
                super().mdrun(
                    tpr, deffnm, write_xtc=write_xtc, write_cpt=write_cpt
                )
                for suffix in (".gro", ".edr", ".log", ".trr", ".xtc", ".cpt"):
                    with open(deffnm + suffix, "w", encoding="ascii") as output:
                        output.write("output")

        params = QmmmParams(
            self.structure,
            self.topology,
            self.index,
            awh_group1="QMatoms",
            awh_group2="Substrate",
        )
        runner = OutputRunner()

        outputs = QmmmSimulation(params).run(
            os.path.join(self.workdir, "results"), runner
        )

        self.assertEqual(runner.names(), ["grompp", "mdrun"])
        self.assertEqual(runner.calls[0][2]["ndx"], self.index)
        self.assertTrue(all(os.path.isfile(path) for path in outputs.values()))
        self.assertTrue(runner.calls[1][2]["write_xtc"])
        self.assertTrue(runner.calls[1][2]["write_cpt"])

    def test_requires_index_groups_to_exist(self):
        params = QmmmParams(
            self.structure,
            self.topology,
            self.index,
            qm_group="MissingQM",
            awh_group1="MissingQM",
            awh_group2="Substrate",
        )

        with self.assertRaisesRegex(ValueError, "MissingQM"):
            QmmmSimulation(params).run(os.path.join(self.workdir, "results"), FakeRunner())

    def test_rejects_index_atom_ids_outside_the_selected_gro(self):
        with open(self.index, "w", encoding="ascii") as output:
            output.write("[ QMatoms ]\n1 5\n[ Substrate ]\n3 4\n")
        params = QmmmParams(
            self.structure,
            self.topology,
            self.index,
            awh_group1="QMatoms",
            awh_group2="Substrate",
        )

        with self.assertRaisesRegex(ValueError, "between 1 and 4"):
            QmmmSimulation(params).run(os.path.join(self.workdir, "results"), FakeRunner())

    def test_reports_missing_simulation_outputs(self):
        class NoOutputRunner(FakeRunner):
            def grompp(self, *args, **kwargs):
                super().grompp(*args, **kwargs)
                with open(args[3], "w", encoding="ascii") as output:
                    output.write("run input")

            def mdrun(self, *args, **kwargs):
                super().mdrun(*args, **kwargs)

        params = QmmmParams(
            self.structure,
            self.topology,
            self.index,
            awh_group1="QMatoms",
            awh_group2="Substrate",
        )

        with self.assertRaisesRegex(RuntimeError, "did not produce"):
            QmmmSimulation(params).run(
                os.path.join(self.workdir, "results"), NoOutputRunner()
            )

    def test_refuses_to_overwrite_previous_run_files(self):
        params = QmmmParams(
            self.structure,
            self.topology,
            self.index,
            awh_group1="QMatoms",
            awh_group2="Substrate",
        )
        output_dir = os.path.join(self.workdir, "results")
        os.makedirs(output_dir)
        with open(os.path.join(output_dir, "qmmm_awh.tpr"), "w") as tpr:
            tpr.write("previous run")

        with self.assertRaisesRegex(FileExistsError, "already exist"):
            QmmmSimulation(params).run(output_dir, FakeRunner())

    def test_biobb_runner_passes_index_file_to_grompp(self):
        runner = BiobbRunner(binary="/opt/gromacs-cp2k/bin/gmx")
        with mock.patch.object(biobb_gmx, "grompp") as grompp:
            runner.grompp(
                "qmmm.mdp", self.structure, self.topology, "qmmm.tpr",
                ndx=self.index,
            )

        self.assertEqual(grompp.call_args.kwargs["ndx"], self.index)
        self.assertEqual(
            grompp.call_args.kwargs["binary"], "/opt/gromacs-cp2k/bin/gmx"
        )
