"""Embedded PyVista molecular viewer for LiGRO."""
from __future__ import annotations

import numpy as np

from .structure import (
    POCKET_CUTOFF,
    PROTEIN_RESIDUES,
    coarse_backbone_segments,
    infer_bonds,
    is_coarse_grained_protein,
    protein_backbone_segments,
    read_pdb_atoms,
    residue_label_positions,
    residue_key,
    select_complex_atoms,
    select_pocket_residues,
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
    """Qt widget that displays a protein cartoon and atomic binding pocket."""

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.set_background("#111827")
        self.add_axes()

    def load_structure(self, path: str) -> None:
        atoms = select_complex_atoms(read_pdb_atoms(path))
        coarse_grained = is_coarse_grained_protein(atoms)
        pocket_residues = select_pocket_residues(atoms)
        atomic_atoms = [
            atom for atom in atoms
            if atom.residue not in PROTEIN_RESIDUES
            or residue_key(atom) in pocket_residues
            or coarse_grained
        ]
        coordinates = np.asarray(
            [atom.coordinates for atom in atomic_atoms], dtype=np.float32
        )
        atom_style = []
        for atom in atomic_atoms:
            if coarse_grained:
                atom_style.append(((89, 190, 205), 1.8))
            elif atom.element == "C" and atom.residue not in PROTEIN_RESIDUES:
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

        backbone_segments = (
            coarse_backbone_segments(atoms)
            if coarse_grained
            else protein_backbone_segments(atoms)
        )
        for backbone in backbone_segments:
            cartoon = pv.Spline(
                backbone,
                n_points=max(20, len(backbone) * 12),
            ).tube(radius=0.3 if coarse_grained else 0.55)
            self.add_mesh(
                cartoon,
                color="#397B8C" if coarse_grained else "#7F94B0",
                smooth_shading=True,
            )

        bonds = (
            infer_bonds(atomic_atoms)
            if atomic_atoms and not coarse_grained
            else np.empty((0, 2))
        )
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

        if atomic_atoms:
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
        label_points, labels = residue_label_positions(atoms, pocket_residues)
        if labels:
            self.add_point_labels(
                label_points,
                labels,
                font_size=12,
                text_color="#F9FAFB",
                shape="rounded_rect",
                shape_color="#111827",
                shape_opacity=0.8,
                margin=3,
                always_visible=True,
                bold=True,
                shadow=True,
                pickable=False,
            )
        structure_summary = (
            "{} CG beads | {} backbone links".format(
                len(atomic_atoms),
                sum(len(segment) - 1 for segment in backbone_segments),
            )
            if coarse_grained
            else "{} atoms (ligand/pocket)".format(len(atomic_atoms))
        )
        self.add_text(
            "{} | {} pocket residues | {} bonds".format(
                structure_summary, len(pocket_residues), len(bonds)
            ),
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
