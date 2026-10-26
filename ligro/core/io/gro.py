"""Pure-Python .gro -> .pdb conversion that keeps coordinates exactly.

Why this exists: `gmx editconf` recentres the molecule when it sets a new
box (the smoke test measured a 2.57 nm shift even with center_molecule=False),
so it cannot be used for a coordinate-preserving conversion without a
confirmed option to disable centring. This module does the conversion
directly: no shift, original box kept in CRYST1.

Limitations (raised, not silently wrong):
  - orthorhombic boxes only (triclinic raises NotImplementedError)
  - PDB serial numbers are 5 digits: more than 99999 atoms raises ValueError
"""
import math

ION_RESNAMES = {"NA", "CL", "NA+", "CL-", "K", "MG", "ZN"}


def read_gro(path):
    """Returns (title, atoms, box). atoms: list of dicts with keys
    resnr, resname, name, serial, x, y, z (all in nm, as in the file).
    box: (a, b, c) in nm."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    if len(lines) < 3:
        raise ValueError("gro file too short: " + path)
    title = lines[0]
    natoms = int(lines[1].strip())
    if len(lines) < natoms + 3:
        raise ValueError("gro file declares {} atoms but has fewer lines".format(natoms))
    atoms = []
    for line in lines[2:2 + natoms]:
        atoms.append({
            "resnr": int(line[0:5]),
            "resname": line[5:10].strip(),
            "name": line[10:15].strip(),
            "serial": int(line[15:20]),
            "x": float(line[20:28]),
            "y": float(line[28:36]),
            "z": float(line[36:44]),
        })
    box_values = [float(v) for v in lines[2 + natoms].split()]
    if len(box_values) < 3:
        raise ValueError("gro box line has fewer than 3 values")
    if any(abs(v) > 1e-12 for v in box_values[3:]):
        raise NotImplementedError("triclinic boxes are not supported by gro_to_pdb")
    return title, atoms, tuple(box_values[:3])


def _element(atom):
    if atom["resname"] in ION_RESNAMES and atom["name"].upper().rstrip("+-") in ION_RESNAMES:
        return atom["name"].upper().rstrip("+-")[:2].capitalize()
    letters = "".join(ch for ch in atom["name"] if ch.isalpha())
    return letters[:1].upper() if letters else ""


def write_pdb(atoms, box, path, title=""):
    if len(atoms) > 99999:
        raise ValueError("PDB format supports at most 99999 atoms; got {}".format(len(atoms)))
    a, b, c = (round(v * 10.0, 3) for v in box)   # nm -> Angstrom
    with open(path, "w", encoding="utf-8") as fh:
        if title:
            fh.write("TITLE     {}\n".format(title[:70]))
        fh.write("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n".format(
            a, b, c, 90.0, 90.0, 90.0))
        for i, atom in enumerate(atoms, start=1):
            x, y, z = atom["x"] * 10.0, atom["y"] * 10.0, atom["z"] * 10.0
            if any(math.isnan(v) for v in (x, y, z)):
                raise ValueError("NaN coordinate in atom {}".format(atom["serial"]))
            name = atom["name"]
            # PDB convention: 4-char atom name; 1-letter elements start at column 14
            name_field = name.rjust(3) if len(name) < 4 else name[:4]
            fh.write("ATOM  {:>5d} {:<4s} {:>3s} A{:>4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n".format(
                i % 100000, name_field, atom["resname"][:3], atom["resnr"] % 10000,
                x, y, z, 1.00, 0.00, _element(atom)))
        fh.write("END\n")


def gro_to_pdb(gro_path, pdb_path):
    title, atoms, box = read_gro(gro_path)
    write_pdb(atoms, box, pdb_path, title=title)
    return pdb_path
