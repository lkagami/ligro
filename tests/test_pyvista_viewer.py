import numpy as np
import pytest

from ligro.viewer.structure import (
    PDBAtom,
    coarse_backbone_segments,
    infer_bonds,
    is_coarse_grained_protein,
    protein_backbone_segments,
    read_pdb_atoms,
    residue_label_positions,
    select_pocket_residues,
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


def test_select_pocket_residues_uses_distance_to_non_protein_heavy_atoms():
    atoms = [
        _atom("CA", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("CA", "C", "GLY", (9.0, 0.0, 0.0)),
        _atom("C1", "C", "LIG", (4.4, 0.0, 0.0), "HETATM"),
        _atom("H1", "H", "LIG", (9.0, 0.0, 0.0), "HETATM"),
    ]

    assert select_pocket_residues(atoms) == {("A", "1", "ALA")}


def test_residue_labels_cover_ligand_and_pocket_not_distant_protein():
    atoms = [
        _atom("CA", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("CB", "C", "ALA", (2.0, 0.0, 0.0)),
        _atom("CA", "C", "GLY", (20.0, 0.0, 0.0)),
        _atom("C1", "C", "LIG", (4.0, 0.0, 0.0), "HETATM"),
        _atom("H1", "H", "LIG", (100.0, 0.0, 0.0), "HETATM"),
    ]

    positions, labels = residue_label_positions(
        atoms, {("A", "1", "ALA")}
    )

    assert labels == ["ALA A:1", "LIG A:1"]
    np.testing.assert_allclose(positions, ((1.0, 0.0, 0.0), (4.0, 0.0, 0.0)))


def test_protein_backbone_segments_split_chains_and_large_gaps():
    atoms = [
        _atom("CA", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("CA", "C", "GLY", (3.8, 0.0, 0.0)),
        _atom("CA", "C", "SER", (12.0, 0.0, 0.0)),
        PDBAtom(
            name="CA",
            element="C",
            residue="ALA",
            chain="B",
            residue_id="1",
            record_type="ATOM",
            coordinates=(20.0, 0.0, 0.0),
            occupancy=1.0,
            alternate_location="",
        ),
        PDBAtom(
            name="CA",
            element="C",
            residue="GLY",
            chain="B",
            residue_id="2",
            record_type="ATOM",
            coordinates=(23.8, 0.0, 0.0),
            occupancy=1.0,
            alternate_location="",
        ),
    ]

    segments = protein_backbone_segments(atoms)

    assert [len(segment) for segment in segments] == [2, 2]
    np.testing.assert_allclose(segments[0][1], (3.8, 0.0, 0.0))
    np.testing.assert_allclose(segments[1][0], (20.0, 0.0, 0.0))


def test_coarse_grained_protein_is_detected_from_missing_alpha_carbon_atoms():
    coarse_atoms = [
        _atom("BB", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("SC1", "C", "ALA", (1.0, 0.0, 0.0)),
    ]
    atomistic_atoms = [_atom("CA", "C", "ALA", (0.0, 0.0, 0.0))]

    assert is_coarse_grained_protein(coarse_atoms)
    assert not is_coarse_grained_protein(atomistic_atoms)


def test_coarse_backbone_segments_follow_bb_beads():
    atoms = [
        _atom("BB", "C", "ALA", (0.0, 0.0, 0.0)),
        _atom("SC1", "C", "ALA", (0.0, 2.0, 0.0)),
        _atom("BB", "C", "GLY", (3.8, 0.0, 0.0)),
        _atom("BB", "C", "SER", (12.0, 0.0, 0.0)),
    ]

    segments = coarse_backbone_segments(atoms)

    assert [len(segment) for segment in segments] == [2]
    np.testing.assert_allclose(segments[0], ((0.0, 0.0, 0.0), (3.8, 0.0, 0.0)))
