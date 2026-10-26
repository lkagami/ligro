#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ribbon_backbone.py
#
#  Independent Cartoon representation -- DATA layer (chain-aware backbone
#  extraction + secondary-structure classification). See this project's
#  own plan file for the full design rationale (user's own request: a
#  from-scratch cartoon representation, not influenced by the existing
#  vismol.utils.cartoon_BCK.py it replaces).
#
#  [EN] Deliberately does NOT reuse VismolObject.define_Calpha_backbone()
#  (vismol_object.py) -- that method flattens EVERY chain's CA atoms into
#  ONE global list and only detects a backbone break via consecutive PDB
#  residue numbering (residue.index == previous.index + 1) plus a < 4.0 A
#  distance sanity check, NEVER comparing chain identity -- two different
#  chains whose PDB numbering happens to run consecutively (common) can
#  get bridged by a bogus CA-CA "bond". extract_backbone_chains() below
#  fixes that class of bug by construction: chain identity is the PRIMARY
#  segmentation key, residue-numbering continuity is only ever checked
#  WITHIN an already-chain-scoped run.
#
#  classify_secondary_structure() implements the same published P-SEA
#  (Labesse et al. 1997) Ca-only H/S/C classification thresholds the
#  (being-replaced) cartoon_BCK.py's own calculate_secondary_structure()
#  used -- the THRESHOLDS are cited science, not borrowed code; this is a
#  fresh, independent implementation, operating on one already-chain-
#  scoped Ca array at a time (never a cross-chain-merged list), with its
#  own dihedral/angle helpers and a simple minimum-run-length smoothing
#  pass instead of whatever block-merging the old module did internally.
# ============================================================================
import numpy as np


# [EN] 2026-09-28 user request ("nossa representacao de cartoon esta
# considerando que nao possam existir quebras na cadeia principal"):
# backbone continuity is decided by GEOMETRY, not by residue numbering --
# see _residues_are_bonded() below. Thresholds:
#   PEPTIDE_BOND_MAX_DISTANCE: C(i)-N(i+1); a real peptide bond is ~1.33 A,
#       2.0 A leaves room for poor/strained geometry while any genuinely
#       missing residue puts it several Angstroms away.
#   CA_CA_MAX_DISTANCE: fallback for Ca-only models (no C/N atoms): trans
#       peptides give ~3.80 A, cis ~2.9 A; a missing residue gives >= ~4.5 A.
PEPTIDE_BOND_MAX_DISTANCE = 2.0
CA_CA_MAX_DISTANCE        = 4.3


def _residues_are_bonded ( prev, curr ):
    """ True if residue `curr` continues the backbone of residue `prev`:
    a peptide bond C(prev)-N(curr) when both atoms exist, else the Ca-Ca
    distance. `prev`/`curr` are the per-residue coordinate dicts built in
    extract_backbone_chains() ("ca", "c", "n": np.ndarray or None). """
    # plain squared-distance arithmetic, not np.linalg.norm: this runs once
    # per residue on every rebuild (every trajectory frame), where numpy's
    # per-call overhead on 3-element vectors dominates.
    if prev["c"] is not None and curr["n"] is not None:
        a, b, limit = prev["c"], curr["n"], PEPTIDE_BOND_MAX_DISTANCE
    else:
        a, b, limit = prev["ca"], curr["ca"], CA_CA_MAX_DISTANCE
    dx = float ( b[0] ) - float ( a[0] )
    dy = float ( b[1] ) - float ( a[1] )
    dz = float ( b[2] ) - float ( a[2] )
    return dx * dx + dy * dy + dz * dz <= limit * limit


def extract_backbone_chains ( vismol_object, frame = None ):
    """ Returns a list of per-chain backbone records, one per CONTIGUOUS
    protein run:
        { "chain": Chain, "residues": [Residue, ...],
          "ca": np.ndarray (N,3), "o": [np.ndarray(3) or None, ...],
          "ca_atoms": [Atom, ...] }

    A residue starts a NEW record whenever (a) it belongs to a different
    Chain object than the one currently running, or (b) it is not
    covalently continuous with the previous residue -- no peptide bond
    C(i)-N(i+1), or, for Ca-only models, a Ca-Ca distance too long for
    adjacent residues (see _residues_are_bonded()) -- so a chain with an
    internal gap correctly produces separate ribbon segments, not one
    falsely-continuous one bridging the gap.

    [EN] 2026-09-28: (b) used to be "residue.index != previous + 1", i.e.
    it trusted the NUMBERING. That fails for prepared/renumbered systems
    (AMBER prmtop, PSF, tleap/psfgen output...), which number residues
    sequentially straight across missing loops and often drop chain IDs:
    6VXX renumbered that way became ONE ribbon with 35 bridged gaps (up to
    28 A, across all 3 chains). It also broke the ribbon wherever the
    numbering jumps without any missing residue (insertion codes, antibody
    numbering schemes). Geometry is the ground truth either way.

    Residues that aren't recognised as protein (residue.is_protein False --
    ligands, waters, ions; see Residue._is_protein(), reused as-is here, a
    boolean CLASSIFIER not a rendering concern) or that are missing a CA
    atom simply close out whatever record is running and are skipped, same
    as a chain change -- a ribbon cannot pass through a point that doesn't
    exist. Records shorter than 2 residues (nothing to draw a spline
    through) are dropped entirely.

    `frame` is passed straight to Atom.coords(); None (the default) uses
    whatever frame the session currently has active -- matching how every
    other "just show me the CURRENT structure" representation build
    already behaves in this codebase, NOT a live trajectory subscription
    (this feature's own plan documents that live per-frame rebuild is out
    of scope for v1). """
    records = [ ]
    current = None

    for chain in vismol_object.chains.values ( ):
        # residues is a dict keyed by PDB residue index -- iterate in that
        # NUMERIC order, not dict/insertion order (not guaranteed already
        # sorted for every loader).
        for residue_index in sorted ( chain.residues.keys ( ) ):
            residue = chain.residues[residue_index]

            ca_atom = None
            o_atom  = None
            c_atom  = None
            n_atom  = None
            if residue.is_protein:
                for atom in residue.atoms.values ( ):
                    if atom.name == "CA":
                        ca_atom = atom
                    elif atom.name == "O":
                        o_atom = atom
                    elif atom.name == "C":
                        c_atom = atom
                    elif atom.name == "N":
                        n_atom = atom

            if ca_atom is None:
                current = None
                continue

            # [EN] Continuity is judged on FRAME 0 (the structure's own
            # topology), never on the frame being drawn: where the chain
            # breaks is a property of the molecule, not of one trajectory
            # frame -- judging it per frame would let a single distorted
            # frame split/merge ribbons (and invalidate the static-SS
            # cache, which is keyed on this segmentation).
            backbone = { "ca": ca_atom.coords ( 0 ),
                         "c":  c_atom.coords ( 0 ) if c_atom is not None else None,
                         "n":  n_atom.coords ( 0 ) if n_atom is not None else None }

            starts_new = (
                current is None
                or current["chain"] is not chain
                or not _residues_are_bonded ( current["_last_backbone"], backbone )
            )
            if starts_new:
                current = { "chain": chain, "residues": [ ], "_ca": [ ], "o": [ ], "ca_atoms": [ ] }
                records.append ( current )
            current["_last_backbone"] = backbone

            current["residues"].append ( residue )
            current["_ca"].append ( ca_atom.coords ( frame ) )
            current["o"].append ( o_atom.coords ( frame ) if o_atom is not None else None )
            # [EN] 2026-09-26 (user request: click-to-select on the ribbon,
            # via a background Ca-trace) -- the ACTUAL Atom object, not just
            # its coords, so a picking-geometry builder can read its
            # precomputed .color_id (see Atom._generate_atom_unique_color_id()
            # and CartoonRepresentation._build_pick_trace_geometry()).
            current["ca_atoms"].append ( ca_atom )

    result = [ ]
    for record in records:
        if len ( record["residues"] ) < 2:
            continue
        record["ca"] = np.asarray ( record["_ca"], dtype = np.float64 )
        del record["_ca"]
        del record["_last_backbone"]
        result.append ( record )
    return result


def _dihedral ( p0, p1, p2, p3 ):
    """ Signed dihedral angle (radians, [-pi, pi]) of the p0-p1-p2-p3
    quadruplet -- standard atan2-of-cross/dot formulation (numerically
    stable at angles near 0/180 degrees, unlike an arccos-of-dot
    formulation would be). """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = np.linalg.norm ( b1 )
    if b1_norm < 1e-9:
        return 0.0
    b1 = b1 / b1_norm
    v = b0 - np.dot ( b0, b1 ) * b1
    w = b2 - np.dot ( b2, b1 ) * b1
    x = np.dot ( v, w )
    y = np.dot ( np.cross ( b1, v ), w )
    return float ( np.arctan2 ( y, x ) )


def _angle_between ( v0, v1 ):
    """ Unsigned angle (radians, [0, pi]) between two vectors -- returns
    None if either is (near) zero-length (degenerate, e.g. two
    coincident Ca positions). """
    n0 = np.linalg.norm ( v0 )
    n1 = np.linalg.norm ( v1 )
    if n0 < 1e-9 or n1 < 1e-9:
        return None
    cos_t = np.clip ( np.dot ( v0, v1 ) / ( n0 * n1 ), -1.0, 1.0 )
    return float ( np.arccos ( cos_t ) )


def _signed_angle_delta_deg ( angle_deg, target_deg ):
    """ Shortest signed difference (degrees) between angle_deg and
    target_deg, correctly wrapping across the +-180 degree boundary (a
    plain subtraction would report ~360 degrees "off" for e.g. 179 vs.
    -179, when the true difference is 2 degrees). """
    return ( angle_deg - target_deg + 180.0 ) % 360.0 - 180.0


def _smooth_short_runs ( codes, min_run_helix = 3, min_run_strand = 2 ):
    """ Demotes any H run shorter than `min_run_helix` (or S run shorter
    than `min_run_strand`) residues back to 'C' -- at these Ca-only
    thresholds, a run that short is almost always noise, not a real
    secondary-structure element. Simple run-length pass, independent of
    whatever block-acceptance logic the module being replaced used
    internally.

    [EN] 2026-09-26 user request ("faz sentido um unico residuo estar em
    conformacao de fita?... o mesmo deve valer para helices"): split what
    used to be ONE shared `min_run` into per-SS-type minimums, since real
    secondary-structure convention treats them differently -- an alpha
    helix needs at least ~3-4 residues to complete even one full turn
    (the i, i+4 H-bond pattern that actually DEFINES a helix has nowhere
    to form on fewer residues than that), while a single beta BRIDGE
    (as short as 2 residues) is already a recognised, real feature in
    DSSP/STRIDE-style classification -- so a lower default minimum for
    strand than helix is a deliberate, literature-grounded choice, not a
    leftover asymmetry. """
    codes = codes.copy ( )
    n = len ( codes )
    i = 0
    while i < n:
        j = i
        while j < n and codes[j] == codes[i]:
            j += 1
        run_len = j - i
        if codes[i] == "H" and run_len < min_run_helix:
            codes[i:j] = "C"
        elif codes[i] == "S" and run_len < min_run_strand:
            codes[i:j] = "C"
        i = j
    return codes


def classify_secondary_structure ( ca_coords, min_run_helix = 3, min_run_strand = 2, helix_tolerance_scale = 1.0, strand_tolerance_scale = 1.5 ):
    """ P-SEA-style (Labesse et al. 1997) Ca-only H(elix)/S(trand)/C(oil)
    classification. `ca_coords` must already be a SINGLE chain's own
    contiguous Ca array (see extract_backbone_chains() above) -- this
    function itself has no chain-awareness of its own, by design, so it
    can never accidentally classify across a chain boundary.

    For each residue i (needs 2 residues of context on either side --
    i-1 and i+1..i+3 -- so only i in [1, n-4] are ever classified; the
    first/last few residues of any chain always come back 'C'), computes:
      - d2/d3/d4: Ca(i-1)-Ca(i+1)/(i+2)/(i+3) distances
      - t : bond angle at Ca(i), between the (i-1)->i and (i+1)->i
            directions (i.e. the Ca(i-1)-Ca(i)-Ca(i+1) angle)
      - a : dihedral of the Ca(i-1),Ca(i),Ca(i+1),Ca(i+2) quadruplet
    against the published threshold ranges:

                           Helix (published)     Strand (published)
      Angle t (deg)        89 +/- 12             124 +/- 14
      Dihedral a (deg)     50 +/- 20            -170 +/- 45
      Distance d2 (A)      5.5 +/- 0.5            6.7 +/- 0.6
      Distance d3 (A)      5.3 +/- 0.5            9.9 +/- 0.9
      Distance d4 (A)      6.4 +/- 0.6           12.4 +/- 1.1

    Every "+/-" tolerance above is multiplied by `helix_tolerance_scale` /
    `strand_tolerance_scale` respectively (uniformly across all 5 criteria
    for that SS type) BEFORE comparing -- 1.0 = the published values
    as-is, >1.0 = more permissive (more residues classified as that SS
    type), <1.0 = stricter.

    [EN] TOLERANCE WIDENED (user report: "a tolerancia para considerar ou
    nao uma estrutura em fita esta muito rigorosa"): strand_tolerance_
    scale's own DEFAULT is 1.5 (i.e. the published values +50%), not 1.0,
    since real beta strands are frequently a little distorted/twisted
    relative to the idealized geometry the paper's own thresholds were
    fit to, and the published tolerances were visibly rejecting real
    strand residues as coil. helix_tolerance_scale defaults to 1.0 (only
    strand detection was reported as too strict). Both are now also
    user-adjustable live from Preferences > Cartoon (see DEFAULT_CARTOON_
    SETTINGS in ribbon_geometry.py and setup_interface.py's
    __apply_cartoon_parameters()) -- the user's own follow-up request for
    "um criterio variavel para a atribuicao de helices e fitas". min_run's
    own default (see signature) was also lowered from 3 to 2 for the same
    original reason.

    then runs _smooth_short_runs() to drop isolated single-residue H or S
    calls back to coil. Returns a numpy array of single-character '<U1'
    codes, same length as ca_coords. """
    n = len ( ca_coords )
    codes = np.full ( n, "C", dtype = "<U1" )
    if n < 5:
        return codes

    # [EN] 2026-09-25 PERFORMANCE FIX (user's own report: "para o modelo
    # de pdb 6vxx, a representacao fica muito lenta" -- 6VXX is the
    # SARS-CoV-2 spike trimer, ~1150 residues PER CHAIN x 3 chains, far
    # beyond anything this function had been profiled against before --
    # crambin.pdb is 46 residues, the largest synthetic test in round 7/8
    # of this module's own development history was 500). This is EXACTLY
    # the scenario that history already flagged as the next thing to fix
    # "if a complaint ever comes in for a much larger structure":
    # classify_secondary_structure() was still a plain per-residue Python
    # `for` loop, calling np.linalg.norm()/np.dot()/np.cross()/np.clip()
    # ONE RESIDUE AT A TIME -- the exact same "numpy dispatch overhead
    # dwarfs tiny-array arithmetic" anti-pattern build_ribbon_mesh() was
    # vectorized away from in round 7 (see that round's own ~20-37x
    # numbers), just never applied here since caching (round 8) already
    # covered the common "one build, then reuse across trajectory frames"
    # case -- but a structure this size is slow on its VERY FIRST build,
    # before there's anything to cache yet. Profiled first, not guessed:
    # ~66 ms per 1150-residue chain (cProfile: _dihedral()'s own np.cross/
    # np.linalg.norm calls alone were >50% of that). Fixed by batching
    # EVERY residue's d2/d3/d4/angle/dihedral into ONE set of whole-array
    # numpy calls (index arrays instead of a Python loop over `range`),
    # same technique as round 7 -- algorithm and thresholds are otherwise
    # completely unchanged.
    idx = np.arange ( 1, n - 3 )
    if idx.size == 0:
        return codes

    p_im1 = ca_coords[idx - 1]
    p_i   = ca_coords[idx]
    p_ip1 = ca_coords[idx + 1]
    p_ip2 = ca_coords[idx + 2]
    p_ip3 = ca_coords[idx + 3]

    d2 = np.linalg.norm ( p_im1 - p_ip1, axis = 1 )
    d3 = np.linalg.norm ( p_im1 - p_ip2, axis = 1 )
    d4 = np.linalg.norm ( p_im1 - p_ip3, axis = 1 )

    # Batched _angle_between ( p_i - p_im1, p_i - p_ip1 ) -- a degenerate
    # (near-zero-length) vector means "skip this residue" in the original
    # scalar version (`continue`, leaving it 'C'); `angle_valid` tracks
    # that exactly, masked into is_helix/is_strand below instead of an
    # early `continue`.
    v0 = p_i - p_im1
    v1 = p_i - p_ip1
    n0 = np.linalg.norm ( v0, axis = 1 )
    n1 = np.linalg.norm ( v1, axis = 1 )
    angle_valid = ( n0 >= 1e-9 ) & ( n1 >= 1e-9 )
    denom = np.where ( angle_valid, n0 * n1, 1.0 )
    cos_t = np.clip ( np.sum ( v0 * v1, axis = 1 ) / denom, -1.0, 1.0 )
    t = np.degrees ( np.arccos ( cos_t ) )

    # Batched _dihedral ( p_im1, p_i, p_ip1, p_ip2 ) -- same degenerate-
    # b1 guard as the scalar version (which returns dihedral 0.0 exactly,
    # not "skip"), reproduced via an explicit np.where() override so the
    # batched and scalar versions agree bit-for-bit on every input.
    b0 = p_im1 - p_i
    b1 = p_ip1 - p_i
    b2 = p_ip2 - p_ip1
    b1_norm = np.linalg.norm ( b1, axis = 1 )
    b1_safe = b1_norm >= 1e-9
    b1_unit = b1 / np.where ( b1_safe, b1_norm, 1.0 )[ :, None ]
    v = b0 - np.sum ( b0 * b1_unit, axis = 1 )[ :, None ] * b1_unit
    w = b2 - np.sum ( b2 * b1_unit, axis = 1 )[ :, None ] * b1_unit
    x = np.sum ( v * w, axis = 1 )
    y = np.sum ( np.cross ( b1_unit, v ) * w, axis = 1 )
    dihedral_rad = np.where ( b1_safe, np.arctan2 ( y, x ), 0.0 )
    a = np.degrees ( dihedral_rad )

    is_helix = (
        angle_valid
        & ( np.abs ( t - 89.0 ) <= 12.0 * helix_tolerance_scale )
        & ( np.abs ( _signed_angle_delta_deg ( a, 50.0 ) ) <= 20.0 * helix_tolerance_scale )
        & ( np.abs ( d2 - 5.5 ) <= 0.5 * helix_tolerance_scale )
        & ( np.abs ( d3 - 5.3 ) <= 0.5 * helix_tolerance_scale )
        & ( np.abs ( d4 - 6.4 ) <= 0.6 * helix_tolerance_scale )
    )
    is_strand = (
        angle_valid
        & ( np.abs ( t - 124.0 ) <= 14.0 * strand_tolerance_scale )
        & ( np.abs ( _signed_angle_delta_deg ( a, -170.0 ) ) <= 45.0 * strand_tolerance_scale )
        & ( np.abs ( d2 - 6.7 ) <= 0.6 * strand_tolerance_scale )
        & ( np.abs ( d3 - 9.9 ) <= 0.9 * strand_tolerance_scale )
        & ( np.abs ( d4 - 12.4 ) <= 1.1 * strand_tolerance_scale )
    )

    codes[ idx[ is_helix & ~is_strand ] ] = "H"
    codes[ idx[ is_strand & ~is_helix ] ] = "S"

    return _smooth_short_runs ( codes, min_run_helix = min_run_helix, min_run_strand = min_run_strand )
