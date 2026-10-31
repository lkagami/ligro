"""Tests for the core minimization pipeline. No GROMACS, no biobb, no Qt."""
import os
import shutil
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

from ligro import biobb_gmx
from ligro.core.io.pdb import clean_protein_pdb, ion_option
from ligro.core.mdp.builder import MinimizationMdp, render_minimization_mdp
from ligro.core.pipelines.minimization import (
    MARTINI_3_FORCE_FIELD,
    ProteinMinimization,
    ProteinMinimizationParams,
)
from ligro.core.runners import FakeRunner


def _read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class FakeMartiniRunner(FakeRunner):
    def martinize2(self, pdb, gro, top, **kwargs):
        super().martinize2(pdb, gro, top, **kwargs)
        with open(gro, "w", encoding="ascii") as structure:
            structure.write("CG structure\n")
        with open(top, "w", encoding="ascii") as topology:
            topology.write(
                '#include "martini_v3.0.0.itp"\n'
                '#include "Protein_A.itp"\n'
                "[ system ]\nProtein\n[ molecules ]\nProtein_A 1\n"
            )
        with open(os.path.join(os.path.dirname(top), "Protein_A.itp"),
                  "w", encoding="ascii") as include:
            include.write("[ moleculetype ]\nProtein_A 1\n")

    def editconf(self, gro_in, gro_out, *args, **kwargs):
        super().editconf(gro_in, gro_out, *args, **kwargs)
        shutil.copyfile(gro_in, gro_out)

    def mdrun(self, tpr, deffnm, *args, **kwargs):
        super().mdrun(tpr, deffnm, *args, **kwargs)
        with open(deffnm + ".gro", "w", encoding="ascii") as structure:
            structure.write("minimized CG structure\n")

    def gro_to_pdb(self, gro_in, pdb_out):
        super().gro_to_pdb(gro_in, pdb_out)
        with open(pdb_out, "w", encoding="ascii") as structure:
            structure.write("minimized CG PDB\n")


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


class Martinize2AdapterTests(unittest.TestCase):
    def test_runs_martinize2_and_checks_generated_files(self):
        with tempfile.TemporaryDirectory() as workdir:
            pdb_path = os.path.join(workdir, "protein.pdb")
            gro_path = os.path.join(workdir, "protein_cg.pdb")
            top_path = os.path.join(workdir, "protein.top")
            with open(pdb_path, "w", encoding="ascii") as pdb:
                pdb.write("PDB\n")

            def run(command, **kwargs):
                with open(command[command.index("-x") + 1], "w") as output:
                    output.write("CG structure\n")
                with open(command[command.index("-o") + 1], "w") as output:
                    output.write("CG topology\n")
                return mock.Mock(returncode=0, stdout="", stderr="")

            with mock.patch.object(biobb_gmx.subprocess, "run", side_effect=run) as call:
                result = biobb_gmx.martinize2(
                    pdb_path, gro_path, top_path, force_field="martini3001"
                )

        self.assertEqual(result, gro_path)
        command = call.call_args.args[0]
        self.assertEqual(command[0], "martinize2")
        self.assertEqual(command[command.index("-ff") + 1], "martini3001")
        self.assertTrue(os.path.isabs(command[command.index("-o") + 1]))

    def test_reports_martinize2_failure(self):
        with mock.patch.object(
            biobb_gmx.subprocess,
            "run",
            return_value=mock.Mock(
                returncode=2, stdout="", stderr="invalid protein structure"
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "invalid protein structure"):
                biobb_gmx.martinize2(
                    "protein.pdb", "protein_cg.pdb", "protein.top"
                )


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

    def test_martini_plan_converts_boxes_and_minimizes_without_atomistic_steps(self):
        plan = ProteinMinimization(
            ProteinMinimizationParams(force_field=MARTINI_3_FORCE_FIELD)
        ).plan()
        self.assertEqual(
            plan,
            [
                "clean_pdb",
                "martinize2",
                "editconf",
                "steep_minimization",
                "postprocess",
            ],
        )

    def test_martini_sd_plus_cg_plan_includes_both_minimizers(self):
        plan = ProteinMinimization(
            ProteinMinimizationParams(
                force_field=MARTINI_3_FORCE_FIELD,
                min_algorithm="CG Algorithm",
            )
        ).plan()
        self.assertEqual(
            plan[-3:],
            ["steep_minimization", "cg_minimization", "postprocess"],
        )


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

    def test_martini_run_uses_biobb_gromacs_and_archives_generated_topology(self):
        runner = FakeMartiniRunner()
        params = ProteinMinimizationParams(
            force_field=MARTINI_3_FORCE_FIELD,
            min_steps=1500,
        )

        result = ProteinMinimization(params).run(self.wd, runner)

        self.assertEqual(result, os.path.join(self.wd, "em.pdb"))
        self.assertEqual(
            runner.names(),
            ["martinize2", "editconf", "grompp", "mdrun", "gro_to_pdb"],
        )
        martinize_call = next(
            call for call in runner.calls if call[0] == "martinize2"
        )
        self.assertEqual(martinize_call[2], {"force_field": "martini3001"})
        with open(os.path.join(self.wd, "em_real.mdp"), encoding="utf-8") as mdp:
            contents = mdp.read()
        self.assertIn("integrator = steep", contents)
        self.assertIn("coulombtype = Reaction-Field", contents)
        self.assertIn("epsilon_r = 15", contents)
        self.assertNotIn("PME", contents)
        with zipfile.ZipFile(os.path.join(self.wd, "trp.zip")) as topology:
            self.assertIn("trp.top", topology.namelist())
            self.assertIn("Protein_A.itp", topology.namelist())
        self.assertFalse(os.path.exists(os.path.join(self.wd, "Protein_A.itp")))

    def test_martini_uses_selected_pdb_outside_project_directory(self):
        source_dir = tempfile.mkdtemp()
        try:
            selected_pdb = os.path.join(source_dir, "selected.pdb")
            with open(selected_pdb, "w", encoding="ascii") as pdb:
                pdb.write("ATOM      1  CA  ALA A   1\n")

            runner = FakeMartiniRunner()
            ProteinMinimization(
                ProteinMinimizationParams(
                    force_field=MARTINI_3_FORCE_FIELD,
                    protein_input=selected_pdb,
                )
            ).run(self.wd, runner)

            martinize_call = next(
                call for call in runner.calls if call[0] == "martinize2"
            )
            self.assertEqual(
                martinize_call[1][0],
                os.path.join(self.wd, "protein_clean.pdb"),
            )
            with open(martinize_call[1][0], encoding="ascii") as cleaned:
                self.assertIn("ALA", cleaned.read())
        finally:
            shutil.rmtree(source_dir)

    def test_martini_missing_input_reports_the_expected_file(self):
        missing_pdb = os.path.join(self.wd, "not-selected.pdb")
        with self.assertRaisesRegex(
            RuntimeError, "Protein PDB input was not found"
        ):
            ProteinMinimization(
                ProteinMinimizationParams(
                    force_field=MARTINI_3_FORCE_FIELD,
                    protein_input=missing_pdb,
                )
            ).run(self.wd, FakeMartiniRunner())

    def test_martini_ignores_selected_metal_and_keeps_only_protein_records(self):
        source = os.path.join(self.wd, "with_ligand_and_zinc.pdb")
        with open(source, "w", encoding="ascii") as pdb:
            pdb.write(
                "ATOM      1  CA  ALA A   1\n"
                "HETATM    2 ZN    ZN A 200\n"
                "HETATM    3  C1  LIG A 501\n"
            )
        runner = FakeMartiniRunner()

        ProteinMinimization(
            ProteinMinimizationParams(
                force_field=MARTINI_3_FORCE_FIELD,
                protein_input=source,
                metal="ZN",
            )
        ).run(self.wd, runner)

        with open(os.path.join(self.wd, "protein_clean.pdb"), encoding="ascii") as clean:
            cleaned = clean.read()
        self.assertIn("ALA", cleaned)
        self.assertNotIn("ZN", cleaned)
        self.assertNotIn("LIG", cleaned)

    def test_martini_cg_minimization_uses_steep_then_cg(self):
        runner = FakeMartiniRunner()
        ProteinMinimization(
            ProteinMinimizationParams(
                force_field=MARTINI_3_FORCE_FIELD,
                min_algorithm="CG Algorithm",
            )
        ).run(self.wd, runner)

        self.assertEqual(runner.names().count("mdrun"), 2)
        with open(os.path.join(self.wd, "em_cg.mdp"), encoding="utf-8") as mdp:
            self.assertIn("integrator = cg", mdp.read())

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
