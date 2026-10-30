# PyVista Qt molecular viewer

LiGRO embeds a PyVista/VTK molecular viewer using `pyvistaqt.QtInteractor` in
the main window's **Complex viewer** tab. A process log remains visible below
the tabs while navigating the application, showing preparation and dynamics
progress plus the final GROMACS production log. Structure review no longer opens
a separate viewer window.
Output from ACPYPE and GROMACS shell commands (both standard output and errors)
is streamed to this panel while each command runs.
LiGRO requests a desktop OpenGL 3.2 core context before creating the Qt
application; this is the minimum needed for VTK's GLSL 1.50 shaders. If the
driver cannot provide that context, update the graphics driver or use a
graphics device that supports OpenGL 3.2 or newer.
The viewer reads the prepared PDB directly, displays protein and non-solvent
ligand atoms as element-colored spheres, and infers nearby covalent bonds for
the stick-like connections. Bulk water and common ions are hidden to keep the
complex legible.

## Install

From the project root and in the same Python environment used to launch LiGRO:

```bash
python -m pip install -r requirements.txt
```

The viewer is implemented in `ligro/viewer/pyvista_qt.py`; it does not depend
on the vendored VisMol renderer or an external molecular-viewer application.
