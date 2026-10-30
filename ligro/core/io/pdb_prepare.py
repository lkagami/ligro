"""Download and clean PDB structures for protein-ligand MD inputs."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "ASH", "CYM", "CYX", "GLH", "HID", "HIE", "HIP",
    "HSD", "HSE", "HSP",
}
SOLVENT_AND_ION_RESIDUES = {
    "HOH", "WAT", "H2O", "SOL", "TIP3", "TIP3P", "TIP4", "TIP4P",
    "SPC", "SPCE", "NA", "NA+", "CL", "CL-", "K", "K+", "MG", "MG2",
    "CA", "CA2", "ZN",
}


@dataclass(frozen=True)
class PDBLigand:
    chain: str
    residue_id: str
    residue_name: str

    @property
    def label(self) -> str:
        chain_label = self.chain or "(blank)"
        return "{} — chain {}, residue {}".format(
            self.residue_name, chain_label, self.residue_id
        )


def _first_model_lines(pdb_path: str) -> list[str]:
    lines = []
    model_seen = False
    in_first_model = True
    with open(pdb_path, encoding="ascii", errors="replace") as pdb_file:
        for line in pdb_file:
            record = line[:6].strip().upper()
            if record == "MODEL":
                if model_seen:
                    in_first_model = False
                else:
                    model_seen = True
                continue
            if record == "ENDMDL":
                if model_seen:
                    break
                continue
            if in_first_model:
                lines.append(line)
    return lines


def inspect_pdb(pdb_path: str) -> tuple[list[str], list[PDBLigand]]:
    """Return protein chains and distinct non-solvent HETATM residues."""
    chains = set()
    ligands = {}
    for line in _first_model_lines(pdb_path):
        record = line[:6].strip().upper()
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 27:
            continue
        chain = line[21:22].strip()
        residue_name = line[17:20].strip().upper()
        residue_id = line[22:27].strip()
        if record == "ATOM" and residue_name in PROTEIN_RESIDUES:
            chains.add(chain)
        elif (
            record == "HETATM"
            and residue_name not in PROTEIN_RESIDUES
            and residue_name not in SOLVENT_AND_ION_RESIDUES
        ):
            key = (chain, residue_id, residue_name)
            ligands[key] = PDBLigand(chain, residue_id, residue_name)
    return sorted(chains, key=lambda value: (value != "", value)), sorted(
        ligands.values(),
        key=lambda ligand: (ligand.chain, ligand.residue_id, ligand.residue_name),
    )


def download_pdb(pdb_id: str, output_path: str) -> str:
    """Download a four-character PDB identifier from RCSB."""
    normalized_id = str(pdb_id).strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", normalized_id):
        raise ValueError("PDB ID must contain exactly four letters or numbers.")

    url = "https://files.rcsb.org/download/{}.pdb".format(normalized_id)
    try:
        with urlopen(url, timeout=30) as response:
            contents = response.read()
    except HTTPError as error:
        raise RuntimeError(
            "RCSB could not find PDB ID {} (HTTP {}).".format(
                normalized_id, error.code
            )
        ) from error
    except (OSError, URLError) as error:
        raise RuntimeError("Could not download PDB ID {}: {}".format(
            normalized_id, error
        )) from error

    if b"ATOM  " not in contents and b"HETATM" not in contents:
        raise RuntimeError(
            "The RCSB response for {} does not contain PDB atom records.".format(
                normalized_id
            )
        )
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "wb") as pdb_file:
        pdb_file.write(contents)
    return os.path.abspath(output_path)


def prepare_complex_inputs(
    pdb_path: str,
    chain: str,
    ligand: PDBLigand,
    output_dir: str,
) -> tuple[str, str]:
    """Write a single-chain receptor PDB and selected ligand MOL2 for LiGRO."""
    if ligand.chain != chain:
        raise ValueError("The selected ligand does not belong to the selected chain.")

    protein_records = []
    ligand_records = []
    target_residue = (ligand.chain, ligand.residue_id, ligand.residue_name)
    for line in _first_model_lines(pdb_path):
        record = line[:6].strip().upper()
        if record not in {"ATOM", "HETATM"} or len(line) < 27:
            continue
        alt_location = line[16:17]
        if alt_location not in {" ", "A"}:
            continue
        atom_chain = line[21:22].strip()
        residue_name = line[17:20].strip().upper()
        residue_id = line[22:27].strip()
        if (
            record == "ATOM"
            and atom_chain == chain
            and residue_name in PROTEIN_RESIDUES
        ):
            protein_records.append(line[:16] + " " + line[17:])
        elif (
            record == "HETATM"
            and (atom_chain, residue_id, residue_name) == target_residue
        ):
            ligand_records.append(line[:16] + " " + line[17:])

    if not protein_records:
        raise ValueError("No protein atoms were found for chain {!r}.".format(chain))
    if not ligand_records:
        raise ValueError(
            "No atoms were found for ligand {}.".format(ligand.label)
        )

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    protein_path = os.path.join(output_dir, "protein.pdb")
    mol2_path = os.path.join(output_dir, "Ligand.mol2")
    with tempfile.TemporaryDirectory(prefix=".ligro-pdb-", dir=output_dir) as temp_dir:
        ligand_pdb = os.path.join(temp_dir, "selected_ligand.pdb")
        temp_mol2 = os.path.join(temp_dir, "Ligand.mol2")
        temp_protein = os.path.join(temp_dir, "protein.pdb")
        with open(ligand_pdb, "w", encoding="ascii") as ligand_file:
            ligand_file.writelines(ligand_records)
            ligand_file.write("END\n")
        with open(temp_protein, "w", encoding="ascii") as protein_file:
            protein_file.writelines(protein_records)
            protein_file.write("TER\nEND\n")

        executable = shutil.which("obabel")
        if executable is None:
            raise RuntimeError(
                "Open Babel (obabel) is required to convert the selected PDB "
                "ligand to MOL2. Install Open Babel and try again."
            )
        try:
            result = subprocess.run(
                [
                    executable, "-ipdb", ligand_pdb, "-omol2",
                    "-O", temp_mol2, "-h",
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError("Open Babel timed out while converting the ligand.") from error
        except OSError as error:
            raise RuntimeError("Could not run Open Babel: {}".format(error)) from error

        if result.returncode != 0:
            details = (result.stderr or result.stdout or "").strip()
            raise RuntimeError(
                "Open Babel could not convert the selected ligand to MOL2"
                + (": " + details[-1000:] if details else ".")
            )
        if not os.path.isfile(temp_mol2):
            raise RuntimeError("Open Babel did not create the ligand MOL2 file.")
        with open(temp_mol2, encoding="utf-8", errors="replace") as mol2_file:
            if "@<TRIPOS>ATOM" not in mol2_file.read():
                raise RuntimeError(
                    "Open Babel produced an invalid MOL2 file without atoms."
                )

        os.replace(temp_protein, protein_path)
        os.replace(temp_mol2, mol2_path)

    return protein_path, mol2_path
