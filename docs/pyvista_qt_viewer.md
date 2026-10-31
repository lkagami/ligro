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
The viewer reads the PDB directly and draws the protein backbone as a smooth
cartoon tube. Ligand atoms and all atoms of protein residues within 4.5
angstroms of a non-protein heavy atom are displayed as element-colored spheres
with inferred covalent bonds. Labels identify each displayed ligand component
and pocket residue by residue name, chain, and residue number. Bulk water and
common ions are hidden to keep the complex legible. Coarse-grained protein
PDBs without alpha-carbon atoms are displayed as bead spheres.

## Install

From the project root and in the same Python environment used to launch LiGRO:

```bash
python -m pip install -r requirements.txt
```

The viewer is implemented in `ligro/viewer/pyvista_qt.py`; it does not depend
on the vendored VisMol renderer or an external molecular-viewer application.
