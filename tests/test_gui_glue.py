"""Tests the glue between main.py's widgets and the core pipeline.

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

sys.modules["plip.modules"].report.__version__ = "2.4.0"   # read by main.py at import
ui = importlib.import_module("ligro.main")

from ligro.core.io.gro import gro_to_pdb  # noqa: E402
from ligro.core.pipelines.minimization import MARTINI_3_FORCE_FIELD  # noqa: E402
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
    def test_qt_default_context_requests_desktop_opengl_32(self):
        surface_format = mock.Mock()
        surface_format_type = mock.Mock(return_value=surface_format)
        surface_format_type.OpenGL = 1
        surface_format_type.CoreProfile = 2
        qt_gui = mock.Mock(QSurfaceFormat=surface_format_type)

        with mock.patch.object(ui, "QtGui", qt_gui):
            ui._configure_qt_opengl()

        surface_format.setRenderableType.assert_called_once_with(1)
        surface_format.setProfile.assert_called_once_with(2)
        surface_format.setVersion.assert_called_once_with(3, 2)
        surface_format.setDepthBufferSize.assert_called_once_with(8)
        surface_format_type.setDefaultFormat.assert_called_once_with(surface_format)

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

    def test_structure_review_loads_the_file_in_the_main_window(self):
        initialization_order = []
        qt_widgets = mock.Mock(QApplication=mock.Mock(
            processEvents=mock.Mock(
                side_effect=lambda *_args: initialization_order.append("events")
            )
        ))
        event_loop = mock.Mock()
        event_loop.exec_.side_effect = lambda: window._continue_structure_review()
        viewer = mock.Mock()
        viewer.load_structure.side_effect = lambda _path: initialization_order.append("load")
        pyvista = mock.Mock(create_pyvista_viewer=mock.Mock(return_value=viewer))
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        structure = os.path.abspath("em.pdb")
        window.viewer_layout = mock.Mock()
        window.viewer_container = mock.Mock()
        window.viewer_placeholder = mock.Mock()
        window.structure_viewer = None
        window.monitor_tab = mock.Mock()
        window.tabWidget = mock.Mock()
        window.tabWidget.indexOf.return_value = 5
        window.centralwidget = mock.Mock()
        window.process_log = mock.Mock()
        window.structure_review_controls = mock.Mock()
        window._structure_review_loop = None
        window._structure_review_result = False
        window._structure_review_cancelled = False

        with mock.patch.object(ui, "QtWidgets", qt_widgets), \
             mock.patch.object(ui.QtCore, "QEventLoop", return_value=event_loop), \
             mock.patch.dict(sys.modules, {"ligro.viewer": pyvista}):
            self.assertTrue(window._confirm_structure(structure, "complex.pdb"))

        self.assertEqual(qt_widgets.QApplication.processEvents.call_count, 3)
        viewer.load_structure.assert_called_once_with(os.path.abspath(structure))
        pyvista.create_pyvista_viewer.assert_called_once_with(window.viewer_container)
        window.tabWidget.setCurrentIndex.assert_called_once_with(5)
        window.viewer_layout.addWidget.assert_called_once_with(viewer)
        self.assertEqual(
            initialization_order, ["events", "load", "events", "events"]
        )
        event_loop.exec_.assert_called_once_with()
        event_loop.quit.assert_called_once_with()
        self.assertEqual(window.structure_review_controls.hide.call_count, 2)

    def test_restore_reloads_original_structure_and_stops_pipeline(self):
        qt_widgets = mock.Mock(QApplication=mock.Mock(processEvents=mock.Mock()))
        event_loop = mock.Mock()
        viewer = mock.Mock()
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        window.viewer_layout = mock.Mock()
        window.viewer_container = mock.Mock()
        window.viewer_placeholder = mock.Mock()
        window.structure_viewer = None
        window.monitor_tab = mock.Mock()
        window.tabWidget = mock.Mock()
        window.tabWidget.indexOf.return_value = 5
        window.centralwidget = mock.Mock()
        window.process_log = mock.Mock()
        window.structure_review_controls = mock.Mock()
        window._structure_review_loop = None
        window._structure_review_result = False
        window._structure_review_cancelled = False

        with mock.patch.object(ui, "QtWidgets", qt_widgets), \
             mock.patch.object(ui.QtCore, "QEventLoop", return_value=event_loop), \
             mock.patch.dict(
                 sys.modules,
                 {"ligro.viewer": mock.Mock(
                     create_pyvista_viewer=mock.Mock(return_value=viewer)
                 )},
             ):
            event_loop.exec_.side_effect = lambda: window._restore_review_structure()
            accepted = window._confirm_structure("em.pdb", "complex.pdb")

        self.assertFalse(accepted)
        self.assertTrue(window._structure_review_cancelled)
        viewer.load_structure.assert_has_calls([
            mock.call(os.path.abspath("em.pdb")),
            mock.call(os.path.abspath("complex.pdb")),
        ])
        event_loop.quit.assert_called_once_with()
        self.assertEqual(window.structure_review_controls.hide.call_count, 2)

    def test_opening_a_local_pdb_displays_it_in_the_viewer(self):
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        window.sel_pdb = mock.Mock()
        window._display_loaded_pdb = mock.Mock()
        dialog = mock.Mock()

        with tempfile.TemporaryDirectory() as workdir, \
             mock.patch.object(ui, "path", workdir), \
             mock.patch.object(ui.QtWidgets, "QFileDialog", return_value=dialog):
            pdb_path = os.path.join(workdir, "input.pdb")
            dialog.getOpenFileName.return_value = (pdb_path, "")
            with open(pdb_path, "w", encoding="ascii") as pdb_file:
                pdb_file.write("PDB")

            window.openpdbfile()

            window.sel_pdb.setText.assert_called_once_with(pdb_path)
            window._display_loaded_pdb.assert_called_once_with(pdb_path)
            self.assertTrue(os.path.isfile(os.path.join(workdir, "protein.pdb")))

    def test_downloading_a_pdb_displays_it_in_the_viewer(self):
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        window.pdb_id_input = Widget("1abc")
        window.pdb_source_path = None
        window.pdb_chain_combo = mock.Mock()
        window.pdb_ligand_combo = mock.Mock()
        window.prepare_pdb_button = mock.Mock()
        window.pdb_id_status = mock.Mock()
        window._display_loaded_pdb = mock.Mock()
        window._append_process_log = mock.Mock()
        window.update_pdb_ligand_choices = mock.Mock()

        with tempfile.TemporaryDirectory() as workdir, \
             mock.patch.object(ui, "path", workdir), \
             mock.patch.object(
                 ui, "download_pdb", return_value=os.path.join(workdir, "1ABC.pdb")
             ), \
             mock.patch.object(
                 ui, "inspect_pdb", return_value=(["A"], [object()])
             ):
            window.download_pdb_by_id()

        window._display_loaded_pdb.assert_called_once_with(
            os.path.join(workdir, "1ABC.pdb")
        )

    def test_gromacs_log_tail_is_added_to_the_main_window_log(self):
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        window.process_log = mock.Mock()
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as log_file:
            log_file.write("GROMACS completed successfully\n")
            log_file.flush()
            window._append_process_log_file(log_file.name)

        message = window.process_log.appendPlainText.call_args[0][0]
        self.assertIn("GROMACS completed successfully", message)
        self.assertIn("GROMACS log contents", message)


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
        w.tabWidget = mock.Mock()
        w.monitor_tab = mock.Mock()
        w.process_log = mock.Mock()
        w._structure_review_cancelled = False
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

    def test_restore_cancellation_does_not_start_md(self):
        window = self.window()
        window.sel_proj_name = Widget("project")
        window.sel_sav_direc = Widget(self.cwd)
        window.sel_exp_fol = Widget(False)
        window.mount_simulation_prot = mock.Mock(
            side_effect=lambda: setattr(window, "_structure_review_cancelled", True)
        )
        window.run_simulation = mock.Mock()

        with mock.patch.object(ui.os, "chdir"), \
             mock.patch.object(
                 ui.os.path,
                 "exists",
                 side_effect=lambda file_path: file_path == "protein.pdb",
             ), \
             mock.patch.object(ui, "showdialog"):
            window.choose_run_simulation()

        window.mount_simulation_prot.assert_called_once_with()
        window.run_simulation.assert_not_called()

    def test_external_command_output_is_streamed_into_process_log(self):
        window = self.window()

        with mock.patch.object(ui.QtWidgets.QApplication, "processEvents"):
            status = window._run_logged_command(
                "printf 'ACPYPE output\\n'; printf 'GROMACS warning\\n' >&2"
            )

        self.assertEqual(status, 0)
        logged_text = "\n".join(
            call.args[0]
            for call in window.process_log.appendPlainText.call_args_list
        )
        self.assertIn("ACPYPE output", logged_text)
        self.assertIn("GROMACS warning", logged_text)
        self.assertIn("Command finished with exit code 0", logged_text)

    def test_external_command_failure_is_logged_and_preserves_nonzero_status(self):
        window = self.window()

        with mock.patch.object(ui.QtWidgets.QApplication, "processEvents"):
            status = window._run_logged_command(
                "printf 'GROMACS failed\\n'; exit 3"
            )

        self.assertEqual(status, 3 << 8)
        logged_text = "\n".join(
            call.args[0]
            for call in window.process_log.appendPlainText.call_args_list
        )
        self.assertIn("GROMACS failed", logged_text)
        self.assertIn("exit code 3", logged_text)

    def run_glue(self, window):
        def fake_dialog(title, text):
            self.dialogs.append((title, text))
        with mock.patch.object(ui, "BiobbRunner", lambda binary="gmx": self.runner), \
             mock.patch.object(ui, "showdialog", fake_dialog):
            return window._run_protein_minimization()

    def test_martini_selection_uses_coarse_grained_minimization(self):
        window = self.window()
        window.sel_ff = Widget(MARTINI_3_FORCE_FIELD)
        window.sel_metal = Widget("ZN")
        window.sel_pdb = Widget(os.path.join(self.cwd, "protein.pdb"))
        runner = FakeRunner()

        with mock.patch.object(ui, "BiobbRunner", return_value=runner), \
             mock.patch.object(ui, "ProteinMinimization") as pipeline, \
             mock.patch.object(ui, "extract_topology"):
            self.assertTrue(window._run_protein_minimization())

        params = pipeline.call_args.args[0]
        self.assertEqual(params.force_field, MARTINI_3_FORCE_FIELD)
        self.assertEqual(
            params.protein_input,
            os.path.join(self.cwd, "protein.pdb"),
        )
        self.assertIsNone(params.metal)
        pipeline.return_value.run.assert_called_once_with(os.getcwd(), runner)

    def test_martini_protein_route_stops_before_atomistic_dynamics(self):
        window = self.window()
        window.sel_ff = Widget(MARTINI_3_FORCE_FIELD)
        window.sel_metal = Widget("None")
        window.sel_pdb = Widget(os.path.join(self.cwd, "protein.pdb"))
        window._run_protein_minimization = mock.Mock(return_value=True)
        window._display_loaded_pdb = mock.Mock()
        window._append_process_log = mock.Mock()
        window._run_logged_command = mock.Mock()

        with mock.patch.object(
            ui, "inspect_pdb", return_value=(["A"], [object()])
        ), mock.patch.object(ui, "showdialog") as showdialog:
            window.mount_simulation_prot()

        window._display_loaded_pdb.assert_called_once_with(
            os.path.join(os.getcwd(), "em.pdb")
        )
        window._run_logged_command.assert_not_called()
        showdialog.assert_called_once()
        self.assertIn("only the protein will be used", showdialog.call_args.args[1])

    def test_martini_with_ligand_file_still_runs_protein_only(self):
        window = self.window()
        window.sel_ff = Widget(MARTINI_3_FORCE_FIELD)
        window.sel_proj_name = Widget("project")
        window.sel_sav_direc = Widget(self.cwd)
        window.sel_exp_fol = mock.Mock(isChecked=mock.Mock(return_value=False))
        window._structure_review_cancelled = False
        window.mount_simulation_prot = mock.Mock()
        window.run_simulation = mock.Mock()
        dialogs = []

        with mock.patch.object(ui.os, "chdir"), \
             mock.patch.object(
                 ui.os.path, "exists",
                 side_effect=lambda file_path: file_path in {
                     "protein.pdb", "Ligand.mol2"
                 },
             ), \
             mock.patch.object(ui, "showdialog", side_effect=lambda *args: dialogs.append(args)):
            window.choose_run_simulation()

        window.mount_simulation_prot.assert_called_once_with()
        window.run_simulation.assert_not_called()

    def test_direct_ligand_route_uses_protein_only_martini_route(self):
        window = ui.Ui_MainWindo.__new__(ui.Ui_MainWindo)
        window.sel_ff = Widget(MARTINI_3_FORCE_FIELD)
        window.mount_simulation_prot = mock.Mock()
        window._run_logged_command = mock.Mock()

        window.mount_simulation_lig()

        window.mount_simulation_prot.assert_called_once_with()
        window._run_logged_command.assert_not_called()

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


class QmmmUiTests(unittest.TestCase):
    def test_setup_ui_creates_the_qmmm_tab(self):
        window = ui.Ui_MainWindo()
        with mock.patch.object(ui, "gromacs_flag", return_value=True):
            window.setupUi(mock.Mock())

        self.assertIsNotNone(window.qmmm_tab)
        self.assertIsNotNone(window.qmmm_run_button)
        self.assertIsNotNone(window.qmmm_awh_group1)


if __name__ == "__main__":
    unittest.main()
