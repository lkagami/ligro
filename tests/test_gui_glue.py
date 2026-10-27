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


class LigandTopologyTests(unittest.TestCase):
    def test_ligand_include_and_molecule_are_added_to_valid_sections(self):
        with tempfile.TemporaryDirectory() as workdir:
            protein_top = os.path.join(workdir, "trp.top")
            ligand_itp = os.path.join(workdir, "Ligand.itp")
            with open(protein_top, "w", encoding="utf-8") as topology:
                topology.write(
                    '#include "amber99sb.ff/forcefield.itp"\n'
                    '[ moleculetype ]\nProtein_chain_A 3\n'
                    '[ molecules ]\nProtein_chain_A 1\n'
                )
            with open(ligand_itp, "w", encoding="utf-8") as ligand:
                ligand.write("[ moleculetype ]\nMOL 3\n")

            ui._write_complex_topology(protein_top, ligand_itp, protein_top)

            with open(protein_top, encoding="utf-8") as topology:
                contents = topology.read()
            self.assertIn(
                '#include "amber99sb.ff/forcefield.itp"\n#include "Ligand.itp"\n',
                contents,
            )
            self.assertRegex(contents, r"(?m)^\[ molecules \]\nProtein_chain_A 1\nMOL 1$")

    def test_existing_ligand_molecule_entry_is_not_duplicated(self):
        with tempfile.TemporaryDirectory() as workdir:
            protein_top = os.path.join(workdir, "trp.top")
            ligand_itp = os.path.join(workdir, "Ligand.itp")
            with open(protein_top, "w", encoding="utf-8") as topology:
                topology.write(
                    '#include "amber99sb.ff/forcefield.itp"\n'
                    '[ molecules ]\nProtein_chain_A 1\nMOL 1\n'
                )
            with open(ligand_itp, "w", encoding="utf-8") as ligand:
                ligand.write("[ moleculetype ]\nMOL 3\n")

            ui._write_complex_topology(protein_top, ligand_itp, protein_top)

            with open(protein_top, encoding="utf-8") as topology:
                contents = topology.read()
            self.assertEqual(contents.count("MOL 1"), 1)

    def test_missing_molecules_section_is_reported(self):
        with tempfile.TemporaryDirectory() as workdir:
            protein_top = os.path.join(workdir, "trp.top")
            ligand_itp = os.path.join(workdir, "Ligand.itp")
            with open(protein_top, "w", encoding="utf-8") as topology:
                topology.write('#include "amber99sb.ff/forcefield.itp"\n')
            with open(ligand_itp, "w", encoding="utf-8") as ligand:
                ligand.write("[ moleculetype ]\nLIG 3\n")

            with self.assertRaisesRegex(ValueError, r"\[ molecules \]"):
                ui._write_complex_topology(protein_top, ligand_itp, protein_top)


class HeredocWriterTests(unittest.TestCase):
    def test_indented_delimiter_is_not_written_into_mdp_file(self):
        with tempfile.TemporaryDirectory() as workdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(workdir)
                status = ui._run_heredoc(
                    "cat << EOF > em_real.mdp\n"
                    "integrator = steep\n"
                    "                EOF\n"
                )
                self.assertEqual(status, 0)
                with open("em_real.mdp", encoding="utf-8") as mdp:
                    contents = mdp.read()
                self.assertEqual(contents, "integrator = steep\n")
            finally:
                os.chdir(old_cwd)


class ViewerTests(unittest.TestCase):
    def test_dots_representation_sets_vertex_shader_point_size(self):
        from types import SimpleNamespace
        from vismol.libgl import representations

        dot_representation = representations.DotsRepresentation.__new__(
            representations.DotsRepresentation
        )
        dot_representation._check_vao_and_vbos = mock.Mock()
        dot_representation._enable_anti_alias_to_lines = mock.Mock()
        dot_representation._disable_anti_alias_to_lines = mock.Mock()
        dot_representation.vm_glcore = SimpleNamespace(
            vm_config=SimpleNamespace(gl_parameters={"dots_size": 2}),
            height=420,
            dist_cam_zrp=10,
            _get_uniform_location=mock.Mock(return_value=7),
            load_matrices=mock.Mock(),
            load_fog=mock.Mock(),
            load_lights=mock.Mock(),
        )
        dot_representation.vm_object = SimpleNamespace(model_mat=object())
        dot_representation.shader_program = 3
        dot_representation.vao = 4
        dot_representation.elements = 1
        dot_representation.was_rep_coord_modified = False
        dot_representation.was_rep_ind_modified = False
        dot_representation.was_col_modified = False

        uniform_upload = mock.Mock()
        with mock.patch.multiple(
            representations.GL,
            glEnable=mock.Mock(),
            glUseProgram=mock.Mock(),
            glPointSize=mock.Mock(),
            glUniform1f=uniform_upload,
            glBindVertexArray=mock.Mock(),
            glDrawElements=mock.Mock(),
            glDisable=mock.Mock(),
        ):
            dot_representation.draw_representation()

        uniform_upload.assert_called_once_with(7, 42.0)

    def test_optional_builder_pick_fallback_returns_a_pair(self):
        from vismol.libgl.vismol_glcore import _OptionalBuilderFallback

        fallback = _OptionalBuilderFallback()
        self.assertEqual(fallback._read_depth_and_atom_at_pixel(None, 0, 0), (None, None))

    def test_structure_review_loads_the_file_in_vismol(self):
        class Signal:
            def connect(self, callback):
                self.callback = callback

        class Dialog:
            Accepted = 1

            def __init__(self, _parent):
                self.buttons = []

            def setWindowTitle(self, _title):
                pass

            def resize(self, _width, _height):
                pass

            def setModal(self, _modal):
                pass

            def show(self):
                pass

            def accept(self):
                pass

            def reject(self):
                pass

            def exec_(self):
                return self.Accepted

        class Layout:
            def __init__(self, _dialog):
                pass

            def addWidget(self, _widget):
                pass

        class ButtonBox:
            Ok = 1
            Cancel = 2

            def __init__(self, _flags, parent):
                self.accepted = Signal()
                self.rejected = Signal()

        initialization_order = []
        qt_widgets = mock.Mock(
            QDialog=Dialog,
            QVBoxLayout=Layout,
            QDialogButtonBox=ButtonBox,
            QApplication=mock.Mock(
                processEvents=mock.Mock(
                    side_effect=lambda: initialization_order.append("events")
                )
            ),
        )
        viewer = mock.Mock()
        viewer.load_structure.side_effect = lambda _path: initialization_order.append("load")
        pyvista = mock.Mock(create_pyvista_viewer=mock.Mock(return_value=viewer))
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        structure = os.path.abspath("em.pdb")

        with mock.patch.object(ui, "QtWidgets", qt_widgets), \
             mock.patch.dict(sys.modules, {"ligro.viewer": pyvista}):
            self.assertTrue(window._confirm_structure(structure))

        qt_widgets.QApplication.processEvents.assert_called_once_with()
        viewer.load_structure.assert_called_once_with(os.path.abspath(structure))
        pyvista.create_pyvista_viewer.assert_called_once_with(mock.ANY)
        self.assertEqual(initialization_order, ["events", "load"])


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

    def test_run_simulation_copies_all_production_outputs(self):
        w = self.window()
        w.sel_proj_name = Widget("trpcage")
        w.sel_sav_direc = Widget(self.cwd)
        outputs = {}
        for key in ("xtc", "edr", "gro", "tpr", "trr", "log", "cpt"):
            output_path = os.path.join(self.cwd, "md." + key)
            with open(output_path, "w", encoding="utf-8") as output_file:
                output_file.write("test")
            outputs[key] = output_path

        with mock.patch.object(w, "_run_dynamics_pipeline", return_value=outputs), \
             mock.patch.object(ui, "showdialog", lambda *_args: None):
            results_path = w.run_simulation()

        self.assertEqual(results_path, os.path.join(self.cwd, "trpcage_results"))
        for key in outputs:
            self.assertTrue(
                os.path.isfile(os.path.join(results_path, "trpcage." + key)),
                key,
            )

    def test_run_simulation_reports_missing_production_output(self):
        w = self.window()
        w.sel_proj_name = Widget("trpcage")
        w.sel_sav_direc = Widget(self.cwd)
        outputs = {key: None for key in ("xtc", "edr", "gro", "tpr", "trr", "log", "cpt")}
        outputs["xtc"] = os.path.join(self.cwd, "missing.xtc")
        with mock.patch.object(w, "_run_dynamics_pipeline", return_value=outputs), \
             mock.patch.object(ui, "showdialog", side_effect=lambda *args: self.dialogs.append(args)):
            self.assertIsNone(w.run_simulation())

        self.assertEqual(self.dialogs[0][0], "Error")
        self.assertIn("missing.xtc", self.dialogs[0][1])


if __name__ == "__main__":
    unittest.main()
