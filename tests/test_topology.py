"""Tests for restoring trp.top and the .itp includes from the biobb zip."""
import io
import os
import shutil
import tempfile
import unittest
import zipfile

from ligro.core.io.topology import extract_topology


def make_zip(path, members):
    with zipfile.ZipFile(path, "w") as z:
        for name, text in members.items():
            z.writestr(name, text)


class TopologyTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.zip = os.path.join(self.dir, "trp.zip")

    def tearDown(self):
        shutil.rmtree(self.dir)

    def read(self, name):
        with open(os.path.join(self.dir, name), encoding="utf-8") as fh:
            return fh.read()

    def test_renames_top_and_keeps_includes(self):
        make_zip(self.zip, {"p2g.top": "[ system ]\n", "posre.itp": "[ position_restraints ]\n"})
        result = extract_topology(self.zip, self.dir)
        self.assertEqual(os.path.basename(result), "trp.top")
        self.assertIn("[ system ]", self.read("trp.top"))
        self.assertIn("[ position_restraints ]", self.read("posre.itp"))

    def test_zip_slip_member_is_flattened(self):
        make_zip(self.zip, {"p2g.top": "x\n", "../../evil.itp": "bad\n"})
        extract_topology(self.zip, self.dir)
        self.assertTrue(os.path.exists(os.path.join(self.dir, "evil.itp")))
        self.assertFalse(os.path.exists(os.path.join(os.path.dirname(self.dir), "evil.itp")))

    def test_missing_top_rejected(self):
        make_zip(self.zip, {"posre.itp": "x\n"})
        with self.assertRaises(RuntimeError):
            extract_topology(self.zip, self.dir)

    def test_two_tops_rejected(self):
        make_zip(self.zip, {"a.top": "x\n", "b.top": "y\n"})
        with self.assertRaises(RuntimeError):
            extract_topology(self.zip, self.dir)

    def test_missing_archive_rejected(self):
        with self.assertRaises(RuntimeError):
            extract_topology(os.path.join(self.dir, "nope.zip"), self.dir)


if __name__ == "__main__":
    unittest.main()
