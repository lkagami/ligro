# PyVista Qt molecular viewer

LiGRO embeds a PyVista/VTK molecular viewer using `pyvistaqt.QtInteractor`.
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
