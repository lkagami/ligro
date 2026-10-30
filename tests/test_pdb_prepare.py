import os
from unittest import mock

import pytest

from ligro.core.io import pdb_prepare


def _atom(serial, record, name, residue, chain, residue_id, xyz, element):
    x, y, z = xyz
    return (
        "{:<6}{:5d} {:>4s} {:>3s} {:1s}{:4d}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
    ).format(
        record, serial, name, residue, chain, residue_id,
        x, y, z, 1.0, 0.0, element,
    )


@pytest.fixture
def pdb_text():
    return (
        "MODEL        1\n"
        + _atom(1, "ATOM", "CA", "ALA", "A", 1, (0, 0, 0), "C")
        + _atom(2, "ATOM", "N", "ALA", "A", 1, (1, 0, 0), "N")
        + _atom(3, "ATOM", "CA", "GLY", "B", 2, (2, 0, 0), "C")
        + _atom(4, "HETATM", "C1", "LIG", "A", 501, (3, 0, 0), "C")
        + _atom(5, "HETATM", "C1", "LIG", "B", 502, (4, 0, 0), "C")
        + _atom(6, "HETATM", "O", "HOH", "A", 900, (5, 0, 0), "O")
        + "ENDMDL\nMODEL        2\n"
        + _atom(7, "ATOM", "CA", "ALA", "C", 3, (6, 0, 0), "C")
        + "ENDMDL\nEND\n"
    )


def test_inspect_pdb_lists_protein_chains_and_ligands_without_solvent(
    tmp_path, pdb_text
):
    structure = tmp_path / "input.pdb"
    structure.write_text(pdb_text, encoding="ascii")

    chains, ligands = pdb_prepare.inspect_pdb(str(structure))

    assert chains == ["A", "B"]
    assert [(item.chain, item.residue_id, item.residue_name) for item in ligands] == [
        ("A", "501", "LIG"),
        ("B", "502", "LIG"),
    ]


def test_prepare_complex_writes_only_selected_chain_and_ligand(
    tmp_path, pdb_text
):
    structure = tmp_path / "input.pdb"
    structure.write_text(pdb_text, encoding="ascii")
    ligand = pdb_prepare.PDBLigand("A", "501", "LIG")

    converted_input = {}

    def convert(command, **_kwargs):
        output = command[command.index("-O") + 1]
        with open(command[2], encoding="ascii") as ligand_pdb:
            converted_input["ligand_pdb"] = ligand_pdb.read()
        with open(output, "w", encoding="ascii") as mol2:
            mol2.write("@<TRIPOS>MOLECULE\nLIG\n@<TRIPOS>ATOM\n")
        return mock.Mock(returncode=0, stdout="", stderr="")

    with mock.patch.object(pdb_prepare.shutil, "which", return_value="/usr/bin/obabel"), \
         mock.patch.object(pdb_prepare.subprocess, "run", side_effect=convert) as run:
        protein_path, mol2_path = pdb_prepare.prepare_complex_inputs(
            str(structure), "A", ligand, str(tmp_path)
        )

    protein = open(protein_path, encoding="ascii").read()
    assert "ALA" in protein
    assert "GLY" not in protein
    assert "LIG" not in protein
    assert os.path.isfile(mol2_path)
    command = run.call_args.args[0]
    ligand_pdb = converted_input["ligand_pdb"]
    assert "LIG" in ligand_pdb
    assert "HOH" not in ligand_pdb
    assert "GLY" not in ligand_pdb
    hetatm_records = [
        line for line in ligand_pdb.splitlines() if line.startswith("HETATM")
    ]
    assert hetatm_records
    assert all(line[21:22].strip() == "A" for line in hetatm_records)
    assert "-h" in command


def test_prepare_complex_requires_openbabel(tmp_path, pdb_text):
    structure = tmp_path / "input.pdb"
    structure.write_text(pdb_text, encoding="ascii")
    ligand = pdb_prepare.PDBLigand("A", "501", "LIG")

    with mock.patch.object(pdb_prepare.shutil, "which", return_value=None):
        with pytest.raises(RuntimeError, match="Open Babel"):
            pdb_prepare.prepare_complex_inputs(
                str(structure), "A", ligand, str(tmp_path)
            )

    assert not (tmp_path / "protein.pdb").exists()
    assert not (tmp_path / "Ligand.mol2").exists()


def test_prepare_complex_rejects_ligand_from_another_chain(tmp_path):
    with pytest.raises(ValueError, match="does not belong"):
        pdb_prepare.prepare_complex_inputs(
            "unused.pdb", "A", pdb_prepare.PDBLigand("B", "1", "LIG"),
            str(tmp_path),
        )


def test_download_pdb_validates_id_before_network_access(tmp_path):
    with mock.patch.object(pdb_prepare, "urlopen") as urlopen:
        with pytest.raises(ValueError, match="exactly four"):
            pdb_prepare.download_pdb("../x", str(tmp_path / "structure.pdb"))

    urlopen.assert_not_called()


def test_download_pdb_writes_rcsb_response(tmp_path):
    response = mock.MagicMock()
    response.__enter__.return_value.read.return_value = b"ATOM  pdb data\n"
    with mock.patch.object(pdb_prepare, "urlopen", return_value=response) as urlopen:
        output = pdb_prepare.download_pdb("1abc", str(tmp_path / "1abc.pdb"))

    assert output == str(tmp_path / "1abc.pdb")
    assert open(output, "rb").read() == b"ATOM  pdb data\n"
    urlopen.assert_called_once_with(
        "https://files.rcsb.org/download/1ABC.pdb", timeout=30
    )
