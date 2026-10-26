#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  residue_bond_orders.py
#
#  Bond orders of residues WITHOUT hydrogens (raw PDB files), 2026-10-15.
#
#  The valence heuristic of vismol_object.perceive_bond_order_for_pairs
#  assumes every hydrogen is present: in a raw PDB a carbon with three heavy
#  neighbours looks unsaturated and gets a spurious double bond. For every
#  residue that has no hydrogen atom this module replaces those orders:
#
#    - standard residues (amino acids and their Amber/CHARMM variants,
#      nucleotides, water, ions): Kekule template by atom name (amino acids
#      from the pDynamo3 PDB component templates, nucleotides from the PDB
#      Chemical Component Dictionary);
#    - any other residue (ligands, modified residues): EasyHybrid's native
#      perception, pdynamo/ligand_perception.py (hydrogen counts from the
#      heavy-atom geometry + Wang & Case orders), when it can be imported;
#      otherwise the valence orders are kept;
#    - bonds between two residues (peptide, phosphodiester, S-S, metal
#      contacts) are single when either residue has no hydrogen.
#
#  Residues that do have hydrogens keep the valence orders (they are right
#  there).
#

import logging

import numpy as np

logger = logging.getLogger(__name__)

_BACKBONE_CO = (("C", "O"), ("C", "OC1"), ("C", "OT1"), ("C", "O1"))

# . non-single bonds of each template (all the other intra-residue bonds are single)
_AMINO = {
    "ALA": [], "CYS": [], "GLY": [], "ILE": [], "LEU": [], "LYS": [], "MET": [],
    "PRO": [], "SER": [], "THR": [], "VAL": [], "UNK": [], "MSE": [],
    "ARG": [("CZ", "NH2", 2)],
    "ASN": [("CG", "OD1", 2)],
    "ASP": [("CG", "OD1", 2)],
    "GLN": [("CD", "OE1", 2)],
    "GLU": [("CD", "OE1", 2)],
    "HIS": [("CG", "CD2", 2), ("ND1", "CE1", 2)],     # . HIE-like (no H to tell the tautomer)
    "HID": [("CG", "CD2", 2), ("CE1", "NE2", 2)],
    "HIE": [("CG", "CD2", 2), ("ND1", "CE1", 2)],
    "HIP": [("CG", "CD2", 2), ("ND1", "CE1", 2)],
    "PHE": [("CG", "CD1", 2), ("CD2", "CE2", 2), ("CE1", "CZ", 2)],
    "TYR": [("CG", "CD1", 2), ("CD2", "CE2", 2), ("CE1", "CZ", 2)],
    "TRP": [("CG", "CD1", 2), ("CD2", "CE2", 2), ("CE3", "CZ3", 2), ("CZ2", "CH2", 2)],
}
_AMINO_ALIAS = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP", "HIN": "HIE", "ASH": "ASP", "GLH": "GLU",
                "CYX": "CYS", "CYM": "CYS", "CYZ": "CYS", "LYN": "LYS", "ARN": "ARG", "TYM": "TYR"}

_PHOSPHATE = [("P", "OP1", 2), ("P", "O1P", 2)]
_NUCLEIC = {
    "A": [("C8", "N7", 2), ("C4", "C5", 2), ("C6", "N1", 2), ("C2", "N3", 2)],
    "G": [("C8", "N7", 2), ("C4", "C5", 2), ("C6", "O6", 2), ("C2", "N3", 2)],
    "C": [("C2", "O2", 2), ("N3", "C4", 2), ("C5", "C6", 2)],
    "U": [("C2", "O2", 2), ("C4", "O4", 2), ("C5", "C6", 2)],
    "T": [("C2", "O2", 2), ("C4", "O4", 2), ("C5", "C6", 2)],
}
_NO_DOUBLE = {"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3", "SPC", "NA", "K", "CL", "MG", "CA", "ZN",
              "MN", "FE", "CU", "CO", "NI", "CD", "BR", "IOD", "LI", "CS", "RB", "SR", "BA", "NA+", "CL-",
              "SOD", "POT", "CLA", "CAL"}


def _nucleic_base(resn):
    """ 'DA', 'DA5', 'RA3', 'A', 'ADE'... -> 'A' (None when not a nucleotide name). """
    r = resn.upper()
    full = {"ADE": "A", "GUA": "G", "CYT": "C", "URA": "U", "THY": "T"}
    if r in full:
        return full[r]
    if len(r) > 1 and r[-1] in "35N":
        r = r[:-1]
    if len(r) == 2 and r[0] in "DR":
        r = r[1]
    return r if r in _NUCLEIC else None


def template(resn, names):
    """ {(name_i, name_j): order} of the non-single bonds of a standard residue,
        or None when resn is not a standard residue. names: atom names present. """
    r = resn.upper().strip()
    if r in _NO_DOUBLE:
        return {}
    if len(r) == 4 and r[0] in "NC" and (r[1:] in _AMINO or r[1:] in _AMINO_ALIAS):
        r = r[1:]                                    # . Amber terminal residues (NALA, CGLY...)
    r = _AMINO_ALIAS.get(r, r)
    out = {}
    if r in _AMINO:
        for a, b in _BACKBONE_CO:
            if a in names and b in names:
                out[(a, b)] = 2
                break
        for a, b, o in _AMINO[r]:
            out[(a, b)] = o
        return out
    base = _nucleic_base(r)
    if base is not None:
        for a, b, o in _NUCLEIC[base]:
            out[(a, b)] = o
        for a, b, o in _PHOSPHATE:
            if a in names and b in names:
                out[(a, b)] = o
                break
        return out
    return None


def _native_perceiver():
    try:
        from pdynamo.ligand_perception import perceive_ligand
        return perceive_ligand
    except Exception:
        return None


def residue_aware_orders(atoms, flat_pairs, orders, coords):
    """ Corrects the valence orders for residues without hydrogens.

        atoms      dict index -> vismol Atom (.symbol, .name, .residue)
        flat_pairs [i0, j0, i1, j1, ...]
        orders     orders of the pairs (valence heuristic), same length
        coords     (N, 3) coordinates (for the native perception)

        Returns a new uint32 array (orders unchanged when every residue has
        hydrogens). """
    ib = np.asarray(flat_pairs, dtype=np.int64).ravel()
    orders = np.array(orders, dtype=np.uint32, copy=True)
    n_bonds = ib.shape[0] // 2
    if n_bonds == 0:
        return orders
    # . residues and whether they have hydrogens
    res_of = {}
    members = {}
    has_h = {}
    for idx, atom in atoms.items():
        res = getattr(atom, "residue", None)
        if res is None:                       # . no residue information: keep the valence orders
            continue
        key = id(res)
        res_of[idx] = key
        members.setdefault(key, []).append(idx)
        if str(getattr(atom, "symbol", "")).upper() == "H":
            has_h[key] = True
    raw = {key for key in members if not has_h.get(key, False)}
    if not raw:
        return orders
    # . bonds of each raw residue; inter-residue bonds touching a raw residue -> single
    intra = {}
    for k in range(n_bonds):
        i, j = int(ib[2 * k]), int(ib[2 * k + 1])
        ri, rj = res_of.get(i), res_of.get(j)
        if ri != rj:
            if ri in raw or rj in raw:
                orders[k] = 1
            continue
        if ri in raw:
            intra.setdefault(ri, []).append(k)
    perceive = None
    n_template = n_native = 0
    for key, bond_ks in intra.items():
        idxs = members[key]
        atom0 = atoms[idxs[0]]
        resn = str(getattr(getattr(atom0, "residue", None), "name", "") or "")
        names = {str(atoms[i].name).strip().upper(): i for i in idxs}
        tpl = template(resn, names)
        if tpl is not None:
            by_pair = {}
            for (a, b), o in tpl.items():
                if a in names and b in names:
                    ia, ibb = names[a], names[b]
                    by_pair[(min(ia, ibb), max(ia, ibb))] = o
            for k in bond_ks:
                i, j = int(ib[2 * k]), int(ib[2 * k + 1])
                orders[k] = by_pair.get((min(i, j), max(i, j)), 1)
            n_template += 1
            continue
        if perceive is None:
            perceive = _native_perceiver() or False
        if not perceive or len(idxs) < 2:
            continue
        local = {a: n for n, a in enumerate(idxs)}
        pairs = [(local[int(ib[2 * k])], local[int(ib[2 * k + 1])]) for k in bond_ks]
        try:
            res = perceive([atoms[a].symbol for a in idxs], np.asarray(coords)[idxs], bonds=pairs,
                           protonate_groups=False)
        except Exception as e:                          # . keep the valence orders
            logger.debug("native perception failed for residue %s: %s", resn, e)
            continue
        got = {(min(i, j), max(i, j)): o for i, j, o in res["bonds"] if i < len(idxs) and j < len(idxs)}
        for k, (i, j) in zip(bond_ks, pairs):
            orders[k] = got.get((min(i, j), max(i, j)), 1)
        n_native += 1
    logger.debug("residue_aware_orders: %d residue(s) without H: %d by template, %d by native perception",
                 len(raw), n_template, n_native)
    return orders
