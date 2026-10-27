import numpy as np
import pytest

from ligro.viewer.structure import (
    PDBAtom,
    infer_bonds,
    read_pdb_atoms,
    select_complex_atoms,
)


def _atom_line(
    serial, name, residue, residue_id, xyz, element, record="ATOM", altloc=" "
):
    x, y, z = xyz
    return (
        f"{record:<6}{serial:5d} {name:>4s}{altloc}{residue:>3s} A"
        f"{residue_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


def _atom(name, element, residue, xyz, record="ATOM"):
    return PDBAtom(
        name=name,
        element=element,
        residue=residue,
        chain="A",
        residue_id="1",
        record_type=record,
        coordinates=xyz,
        occupancy=1.0,
        alternate_location="",
    )


def test_read_pdb_atoms_parses_coordinates_and_uses_first_model(tmp_path):
    pdb = tmp_path / "structure.pdb"
    pdb.write_text(
        "MODEL        1\n"
        + _atom_line(1, "CA", "ALA", 1, (1.0, 2.0, 3.0), "C")
        + _atom_line(2, "O", "ALA", 1, (2.0, 2.0, 3.0), "O")
        + "ENDMDL\nMODEL        2\n"
        + _atom_line(1, "CA", "ALA", 1, (9.0, 9.0, 9.0), "C")
        + "ENDMDL\n",
        encoding="ascii",
    )

    atoms = read_pdb_atoms(pdb)

    assert len(atoms) == 2
    assert atoms[0].coordinates == (1.0, 2.0, 3.0)
    assert atoms[0].element == "C"


def test_read_pdb_atoms_prefers_primary_alternate_location(tmp_path):
    pdb = tmp_path / "alternate.pdb"
    pdb.write_text(
        _atom_line(1, "CA", "ALA", 1, (4.0, 0.0, 0.0), "C", altloc="B")
        + _atom_line(2, "CA", "ALA", 1, (1.0, 0.0, 0.0), "C", altloc="A")
        + "END\n",
        encoding="ascii",
    )

    atoms = read_pdb_atoms(pdb)

    assert len(atoms) == 1
    assert atoms[0].alternate_location == "A"
    assert atoms[0].coordinates == (1.0, 0.0, 0.0)


def test_select_complex_atoms_excludes_water_and_ions_but_keeps_ligand():
    atoms = [
        _atom("CA", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("O", "O", "HOH", (8.0, 0.0, 0.0), "HETATM"),
        _atom("NA", "NA", "NA", (9.0, 0.0, 0.0), "HETATM"),
        _atom("C1", "C", "LIG", (1.5, 0.0, 0.0), "HETATM"),
    ]

    selected = select_complex_atoms(atoms)

    assert [atom.residue for atom in selected] == ["ALA", "LIG"]


def test_infer_bonds_finds_nearby_covalent_atoms_only():
    atoms = [
        _atom("C", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("O", "O", "ALA", (1.3, 0.0, 0.0)),
        _atom("C", "C", "ALA", (6.0, 0.0, 0.0)),
    ]

    np.testing.assert_array_equal(infer_bonds(atoms), np.array([[0, 1]]))


def test_select_complex_atoms_rejects_solvent_only_structure():
    with pytest.raises(ValueError, match="no protein or ligand atoms"):
        select_complex_atoms([_atom("O", "O", "SOL", (0.0, 0.0, 0.0))])
