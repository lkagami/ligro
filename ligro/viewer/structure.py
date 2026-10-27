"""Small, dependency-light helpers for molecular PDB structures."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


SOLVENT_AND_ION_RESIDUES = {
    "HOH", "WAT", "H2O", "SOL", "TIP3", "TIP3P", "TIP4", "TIP4P",
    "SPC", "SPCE", "NA", "NA+", "CL", "CL-", "K", "K+", "MG", "MG2",
    "CA", "CA2", "ZN",
}
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "ASH", "CYM", "CYX", "GLH", "HID", "HIE", "HIP",
    "HSD", "HSE", "HSP",
}
COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "BR": 1.20,
    "I": 1.39,
}


@dataclass(frozen=True)
class PDBAtom:
    name: str
    element: str
    residue: str
    chain: str
    residue_id: str
    record_type: str
    coordinates: tuple[float, float, float]
    occupancy: float
    alternate_location: str


def _element_from_line(line: str, atom_name: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    stripped = atom_name.strip().lstrip("0123456789")
    if not stripped:
        return "C"
    if len(atom_name) >= 2 and atom_name[0].isspace():
        return stripped[0].upper()
    candidate = stripped[:2].upper()
    if candidate in COVALENT_RADII:
        return candidate
    return candidate[0]


def read_pdb_atoms(path: str | Path) -> list[PDBAtom]:
    """Read the first model of a fixed-column PDB file."""
    atoms_by_key: dict[tuple[str, ...], PDBAtom] = {}
    with open(path, encoding="ascii", errors="replace") as pdb_file:
        for line_number, line in enumerate(pdb_file, start=1):
            record_type = line[:6].strip().upper()
            if record_type == "ENDMDL":
                break
            if record_type not in {"ATOM", "HETATM"}:
                continue
            if len(line) < 54:
                raise ValueError(
                    "Malformed PDB atom record on line {}".format(line_number)
                )
            atom_name = line[12:16].strip()
            alternate_location = line[16:17].strip()
            residue = line[17:20].strip().upper()
            occupancy_text = line[54:60].strip()
            try:
                coordinates = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
                occupancy = float(occupancy_text) if occupancy_text else 1.0
            except ValueError as exc:
                raise ValueError(
                    "Invalid atom coordinates or occupancy on PDB line {}".format(
                        line_number
                    )
                ) from exc
            if not np.isfinite(coordinates).all():
                raise ValueError(
                    "Non-finite atom coordinates on PDB line {}".format(line_number)
                )

            atom = PDBAtom(
                name=atom_name,
                element=_element_from_line(line, atom_name),
                residue=residue,
                chain=line[21:22].strip(),
                residue_id=line[22:27].strip(),
                record_type=record_type,
                coordinates=coordinates,
                occupancy=occupancy,
                alternate_location=alternate_location,
            )
            key = (atom.chain, atom.residue_id, atom.residue, atom.name)
            previous = atoms_by_key.get(key)
            if previous is None or (
                previous.alternate_location and not atom.alternate_location
            ) or (
                previous.alternate_location == "B"
                and atom.alternate_location == "A"
            ):
                atoms_by_key[key] = atom

    atoms = list(atoms_by_key.values())
    if not atoms:
        raise ValueError("No atoms were found in PDB file: {}".format(path))
    return atoms


def select_complex_atoms(atoms: list[PDBAtom]) -> list[PDBAtom]:
    """Keep protein and non-solvent ligand records; omit bulk solvent and ions."""
    selected = [
        atom for atom in atoms
        if atom.residue not in SOLVENT_AND_ION_RESIDUES
    ]
    if not selected:
        raise ValueError(
            "The PDB contains no protein or ligand atoms after solvent and ion "
            "removal."
        )
    return selected


def infer_bonds(atoms: list[PDBAtom]) -> np.ndarray:
    """Infer close covalent contacts using a spatial grid rather than all pairs."""
    coordinates = np.asarray([atom.coordinates for atom in atoms], dtype=np.float32)
    cell_width = 2.0
    cells: dict[tuple[int, int, int], list[int]] = {}
    bonds: list[tuple[int, int]] = []
    neighbor_offsets = [
        (x, y, z)
        for x in (-1, 0, 1)
        for y in (-1, 0, 1)
        for z in (-1, 0, 1)
    ]

    for index, (atom, coordinate) in enumerate(zip(atoms, coordinates)):
        cell = tuple(np.floor(coordinate / cell_width).astype(int))
        radius = COVALENT_RADII.get(atom.element, 0.77)
        for offset in neighbor_offsets:
            neighbor_cell = tuple(cell[axis] + offset[axis] for axis in range(3))
            for other_index in cells.get(neighbor_cell, ()):
                other = atoms[other_index]
                if atom.element == other.element == "H":
                    continue
                other_radius = COVALENT_RADII.get(other.element, 0.77)
                delta = coordinate - coordinates[other_index]
                distance_squared = float(np.dot(delta, delta))
                max_distance = min(1.95, 1.25 * (radius + other_radius))
                if 0.16 < distance_squared <= max_distance * max_distance:
                    bonds.append((other_index, index))
        cells.setdefault(cell, []).append(index)

    return np.asarray(bonds, dtype=np.int64).reshape((-1, 2))
