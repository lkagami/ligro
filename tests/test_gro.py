"""Tests for the pure-Python gro -> pdb conversion (no GROMACS needed)."""
import os
import shutil
import tempfile
import unittest

from ligro.core.io.gro import gro_to_pdb, read_gro

GRO = (
    "test system\n"
    "    3\n"
    "    1ASN     N    1   1.000   2.000   3.000\n"
    "    1ASN    CA    2   1.500   2.500   3.500\n"
    "    2NA      NA    3   4.250   5.125   6.000\n"
    "   10.00000  11.00000  12.00000\n"
)


class GroConversionTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.gro = os.path.join(self.dir, "in.gro")
        self.pdb = os.path.join(self.dir, "out.pdb")
        with open(self.gro, "w") as fh:
            fh.write(GRO)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def atoms_from_pdb(self):
        out = []
        with open(self.pdb) as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        return out

    def test_coordinates_are_exact_nm_to_angstrom(self):
        gro_to_pdb(self.gro, self.pdb)
        coords = self.atoms_from_pdb()
        self.assertEqual(len(coords), 3)
        self.assertAlmostEqual(coords[0][0], 10.0, places=3)
        self.assertAlmostEqual(coords[0][1], 20.0, places=3)
        self.assertAlmostEqual(coords[2][0], 42.5, places=3)
        self.assertAlmostEqual(coords[2][1], 51.25, places=3)

    def test_no_recentring_shift(self):
        # The regression that motivated this module: editconf moved atoms.
        gro_to_pdb(self.gro, self.pdb)
        _, atoms, _ = read_gro(self.gro)
        coords = self.atoms_from_pdb()
        for atom, (x, y, z) in zip(atoms, coords):
            self.assertAlmostEqual(x, atom["x"] * 10.0, places=3)
            self.assertAlmostEqual(y, atom["y"] * 10.0, places=3)
            self.assertAlmostEqual(z, atom["z"] * 10.0, places=3)

    def test_box_is_kept_in_cryst1(self):
        gro_to_pdb(self.gro, self.pdb)
        with open(self.pdb) as fh:
            cryst = [l for l in fh if l.startswith("CRYST1")][0]
        self.assertAlmostEqual(float(cryst[6:15]), 100.0, places=3)
        self.assertAlmostEqual(float(cryst[15:24]), 110.0, places=3)
        self.assertAlmostEqual(float(cryst[24:33]), 120.0, places=3)

    def test_ion_element_symbol(self):
        gro_to_pdb(self.gro, self.pdb)
        with open(self.pdb) as fh:
            last = [l for l in fh if l.startswith("ATOM")][-1]
        self.assertEqual(last[76:78].strip(), "Na")

    def test_triclinic_box_rejected(self):
        bad = GRO.replace("  10.00000  11.00000  12.00000",
                          "  10.00000  11.00000  12.00000   0.00000   0.00000   0.50000")
        with open(self.gro, "w") as fh:
            fh.write(bad)
        with self.assertRaises(NotImplementedError):
            gro_to_pdb(self.gro, self.pdb)

    def test_truncated_file_rejected(self):
        with open(self.gro, "w") as fh:
            fh.write("title\n    5\n    1ASN     N    1   1.000   2.000   3.000\n")
        with self.assertRaises(ValueError):
            gro_to_pdb(self.gro, self.pdb)


if __name__ == "__main__":
    unittest.main()
