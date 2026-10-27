"""Embedded PyVista molecular viewer for LiGRO."""
from __future__ import annotations

import numpy as np

from .structure import (
    PROTEIN_RESIDUES,
    infer_bonds,
    read_pdb_atoms,
    select_complex_atoms,
)

try:
    import pyvista as pv
    from pyvistaqt import QtInteractor
except ImportError as exc:
    missing_module = getattr(exc, "name", None) or str(exc)
    raise ImportError(
        "LiGRO's molecular viewer requires PyVista and pyvistaqt. "
        "Install them with: python -m pip install -r requirements.txt "
        "(missing: {})".format(missing_module)
    ) from exc


_ELEMENT_STYLE = {
    "H": ((235, 239, 244), 0.20),
    "C": ((170, 183, 197), 0.34),
    "N": ((76, 137, 220), 0.33),
    "O": ((229, 83, 87), 0.32),
    "F": ((105, 190, 119), 0.31),
    "P": ((240, 151, 62), 0.38),
    "S": ((235, 204, 76), 0.38),
    "CL": ((92, 190, 112), 0.40),
    "BR": ((165, 91, 67), 0.43),
    "I": ((135, 92, 170), 0.45),
}
_DEFAULT_ELEMENT_STYLE = ((183, 157, 208), 0.36)


class PyVistaStructureViewer(QtInteractor):
    """Qt widget that displays protein/ligand atoms and inferred bonds."""

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.set_background("#111827")
        self.add_axes()

    def load_structure(self, path: str) -> None:
        atoms = select_complex_atoms(read_pdb_atoms(path))
        coordinates = np.asarray([atom.coordinates for atom in atoms], dtype=np.float32)
        atom_style = []
        for atom in atoms:
            if atom.element == "C" and atom.residue not in PROTEIN_RESIDUES:
                atom_style.append(((239, 180, 76), 0.36))
            else:
                atom_style.append(
                    _ELEMENT_STYLE.get(atom.element, _DEFAULT_ELEMENT_STYLE)
                )
        colors = np.asarray([style[0] for style in atom_style], dtype=np.uint8)
        radii = np.asarray([style[1] for style in atom_style], dtype=np.float32)

        self.clear()
        self.set_background("#111827")
        self.add_axes()

        bonds = infer_bonds(atoms)
        if len(bonds):
            bond_mesh = pv.PolyData(coordinates)
            bond_cells = np.column_stack(
                (np.full(len(bonds), 2, dtype=np.int64), bonds)
            )
            bond_mesh.lines = bond_cells.ravel()
            self.add_mesh(
                bond_mesh,
                color="#93A4B8",
                line_width=2.0,
                render_lines_as_tubes=True,
            )

        atom_cloud = pv.PolyData(coordinates)
        atom_cloud["radius"] = radii
        atom_cloud["rgb"] = colors
        sphere = pv.Sphere(
            radius=1.0,
            theta_resolution=12,
            phi_resolution=12,
        )
        atom_mesh = atom_cloud.glyph(
            geom=sphere,
            scale="radius",
            orient=False,
        )
        self.add_mesh(
            atom_mesh,
            scalars="rgb",
            rgb=True,
            smooth_shading=True,
        )
        self.add_text(
            "{} atoms | {} bonds".format(len(atoms), len(bonds)),
            position="upper_left",
            font_size=11,
            color="white",
        )
        self.reset_camera()
        self.view_isometric()
        self.render()


def create_pyvista_viewer(parent=None) -> PyVistaStructureViewer:
    """Create an embedded Qt viewer for a prepared PDB structure."""
    return PyVistaStructureViewer(parent=parent)
