#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  bond_perception.py
#
#  Description:
#      Bond ORDER perception (single/double/triple) from geometry alone --
#      the algorithm documented in `roteiro_percepcao_ligacoes.md` (the
#      user's own write-up of how Gabedit's src/Geometry/GeomXYZ.c does
#      this: connected(), connecteds(), reset_connections_XYZ(),
#      set_multiple_bonds()), ported as closely as possible:
#
#          Section 3 (Etapa A)              -> connectivity_from_distance()
#          Section 4 (Etapa B, steps B.1-B.3) -> promote_multiple_bonds()
#          Section 6 (pseudocódigo consolidado) -> perceive_bond_orders()
#
#      Deliberately a PURELY GEOMETRIC/HEURISTIC method (explicitly NOT a
#      quantum bond-order calculation, e.g. Wiberg/Mayer) -- good enough
#      to DRAW a molecule correctly (single/double/triple sticks), not
#      for rigorous chemical analysis. Same documented limitation as the
#      roadmap's own section 5: this does NOT detect rings, does not
#      count pi electrons, does not apply Hückel's rule -- the
#      single/double alternation that emerges in aromatic rings (benzene,
#      etc.) is a SIDE EFFECT of "one unit of valence slack per atom"
#      combined with the sequential, index-order-dependent scan, not real
#      aromaticity detection. Verified against exactly this benzene case
#      with a standalone numpy simulation before wiring this in (see the
#      conversation this was built in) -- confirmed perfect alternation
#      (3 double + 3 single around the ring), matching the roadmap's own
#      predicted "emergent" behaviour. Also verified against ethylene
#      (C=C), acetylene (C#C), and methane (no promotion at all).
#
#      DELIBERATELY PLACED HERE (graphics_engine, not EasyHybrid's own
#      gui/windows/builder/atom_ops.py, where a first version of this
#      lived briefly): this algorithm has ZERO EasyHybrid-specific
#      dependencies -- it only ever needs plain symbol/position/covalent-
#      radius lists -- so keeping it in the engine avoids the exact
#      layering mistake this project already hit once before (vismol_
#      glcore.py importing gui.windows.builder at module level caused a
#      real, reproducible circular-import crash on app startup -- see
#      that fix's own history). find_bonded_and_nonbonded_atoms() below
#      in vismol_object.py calls straight into this module (a same-
#      submodule import, no risk of that repeating), and EasyHybrid's own
#      Builder can still use these same functions for its own "Perceive
#      Bonds"-style features later, by importing THIS module instead of
#      duplicating the algorithm a second time.
#
#      TOGGLE: see AUTO_PERCEIVE_MULTIPLE_BONDS below -- flip that one
#      constant to turn this off globally (every automatic file load
#      falls back to the previous behaviour: every detected bond is a
#      plain single line) without needing to touch any call site. Every
#      call site also accepts its own perceive_multiple_bonds=True/False
#      override, for turning it off (or on) for one specific call only.
#

# [EN] Master on/off switch for AUTOMATIC bond-order perception (used by
# VismolObject.find_bonded_and_nonbonded_atoms(), the path every
# normally-loaded file with no explicit connectivity/bond-order
# information takes -- e.g. a plain .xyz file). Flip this to False to
# go back to the previous behaviour everywhere (every automatically-
# detected bond stays a plain single line) without editing anything
# else. Does NOT affect systems that already HAVE real bond-order data
# from elsewhere (an MM force field's topology, a MOL2/SDF file with
# explicit bond orders, or the Builder's own manually-drawn bonds) --
# those never go through this module at all.
AUTO_PERCEIVE_MULTIPLE_BONDS = True

# [EN] Tolerance added to EACH atom's own covalent radius before summing
# a pair's bonding threshold (Section 2.1 of the roadmap) -- Gabedit's
# own value, replicated exactly (adding a single +0.2 to the SUM instead
# would be a different, smaller total slack).
DEFAULT_TOLERANCE = 0.2

# [EN] Typical valence (NOT the full ~112-element table from Gabedit's
# own AtomsProp.c -- this covers the elements most common in organic
# chemistry/QM-MM; easy to extend later if more are needed). Keys
# uppercase, matched via symbol.upper().
GABEDIT_MAX_VALENCE = {
    'H' : 1, 'HE': 0,
    'LI': 1, 'BE': 2, 'B' : 3, 'C' : 4, 'N' : 3, 'O' : 2, 'F' : 1, 'NE': 0,
    'NA': 1, 'MG': 2, 'AL': 3, 'SI': 4, 'P' : 3, 'S' : 2, 'CL': 1, 'AR': 0,
    'K' : 1, 'CA': 2, 'BR': 1, 'I' : 1,
    'FE': 2, 'ZN': 2, 'CU': 2, 'MN': 2, 'NI': 2, 'CO': 2,
}


def connectivity_from_distance ( symbols, positions, cov_radii, tolerance = DEFAULT_TOLERANCE ):
    """ [EN] Section 3 ("Etapa A") of the roadmap -- purely distance-based
    connectivity, no angles, no hybridisation, no ring logic at all.

    d(i,j) < (cov_radii[i] + tolerance) + (cov_radii[j] + tolerance)

    symbols/positions/cov_radii must all be the same length N (index i of
    each refers to the same atom). Returns an NxN list-of-lists,
    bond_order[i][j] == bond_order[j][i], 0 (no bond) or 1 (bonded) --
    multiple-bond promotion happens separately, in
    promote_multiple_bonds() below. O(N²), matching the roadmap's own
    note in section 3 (a grid/cell-list speed-up is suggested there for
    very large molecules -- VismolObject.find_bonded_and_nonbonded_atoms()
    already HAS exactly that, via c_distances.get_atomic_bonds_from_grid(),
    which is why THAT function calls promote_multiple_bonds() directly on
    its own, already fast, connectivity result instead of using this
    slower function -- this one is here for standalone/testing use, and
    for any future caller that doesn't already have a connectivity result
    of its own to build on). """
    n = len ( symbols )
    bond_order = [ [ 0 ] * n for _ in range ( n ) ]

    for i in range ( n ):
        for j in range ( i + 1, n ):
            dx = positions[i][0] - positions[j][0]
            dy = positions[i][1] - positions[j][1]
            dz = positions[i][2] - positions[j][2]
            distance = ( dx * dx + dy * dy + dz * dz ) ** 0.5

            threshold = ( cov_radii[i] + tolerance ) + ( cov_radii[j] + tolerance )
            if distance < threshold:
                bond_order[i][j] = 1
                bond_order[j][i] = 1

    return bond_order


def promote_multiple_bonds ( bond_order, max_valences ):
    """ [EN] Section 4 ("Etapa B") of the roadmap, steps B.1/B.2/B.3 --
    promotes existing single bonds (bond_order[i][j] == 1) to double,
    then existing double bonds to triple, based ONLY on how much
    "valence slack" (max_valence - current bond count) each atom still
    has -- distance is NEVER looked at again here, only counts.

    Mutates `bond_order` IN PLACE (matching the roadmap's own pseudocode,
    section 6, where n_bonds is updated immediately DURING the scan, not
    from a frozen snapshot -- this is what makes the algorithm greedy/
    sequential and dependent on atom index order, a property explicitly
    documented as a KNOWN, expected limitation in the roadmap's own
    section 5, not a bug to "fix" by e.g. sorting by distance first (see
    that section's own "melhorias possíveis" list for that idea,
    deliberately left as a future, separate option).

    Also returns `bond_order`, for convenience when called as part of an
    expression rather than as a standalone mutating statement. """
    n = len ( bond_order )

    n_bonds = [ sum ( 1 for j in range ( n ) if bond_order[i][j] != 0 ) for i in range ( n ) ]

    # --- B.2: promove para DUPLA ---
    for i in range ( n ):
        for j in range ( i + 1, n ):
            if bond_order[i][j] != 1:
                continue
            if n_bonds[i] < max_valences[i] and n_bonds[j] < max_valences[j]:
                bond_order[i][j] = 2
                bond_order[j][i] = 2
                n_bonds[i] += 1
                n_bonds[j] += 1

    # --- B.3: promove para TRIPLA ---
    for i in range ( n ):
        for j in range ( i + 1, n ):
            if bond_order[i][j] != 2:
                continue
            if n_bonds[i] < max_valences[i] and n_bonds[j] < max_valences[j]:
                bond_order[i][j] = 3
                bond_order[j][i] = 3
                n_bonds[i] += 1
                n_bonds[j] += 1

    return bond_order


def perceive_bond_orders ( symbols, positions, cov_radii, tolerance = DEFAULT_TOLERANCE, max_valence_table = None ):
    """ [EN] Section 6 of the roadmap -- the single consolidated entry
    point: runs connectivity_from_distance() (Etapa A) then
    promote_multiple_bonds() (Etapa B) in sequence, from scratch. See
    promote_multiple_bonds_for_existing_connectivity() below instead when
    connectivity has ALREADY been determined some other (e.g. faster,
    grid-accelerated) way and only the bond-ORDER promotion step is
    needed.

    symbols    : list of N element symbols (e.g. "C", "H", "O", ...).
    positions  : Nx3 array-like of coordinates (any consistent unit, same
                 one cov_radii is expressed in -- Angstrom throughout the
                 rest of this app).
    cov_radii  : list of N covalent radii, ALREADY resolved per atom
                 (e.g. from each Atom's own .cov_rad) -- kept as a
                 separate, explicit argument rather than looking symbols
                 up in a table internally, so this stays decoupled from
                 any one specific data source (this app's periodic table,
                 Gabedit's own ~112-element AtomsProp.c table, or
                 anything else).
    max_valence_table : optional dict SYMBOL (uppercase) -> int, falls
                 back to GABEDIT_MAX_VALENCE above; unknown symbols
                 default to 4 (a permissive middle-ground -- better to
                 risk one extra bond-order promotion than to silently
                 forbid all multiple bonds for an element this table
                 simply doesn't happen to list yet).

    Returns an NxN list-of-lists bond_order matrix, symmetric, values
    0 (not bonded) / 1 / 2 / 3. """
    max_valences = _resolve_max_valences ( symbols, max_valence_table )
    bond_order = connectivity_from_distance ( symbols, positions, cov_radii, tolerance = tolerance )
    promote_multiple_bonds ( bond_order, max_valences )
    return bond_order


def promote_multiple_bonds_for_existing_connectivity ( symbols, index_bonds, max_valence_table = None ):
    """ [EN] Same Etapa B (steps B.1-B.3) as promote_multiple_bonds()
    above, but starting from a CONNECTIVITY ALREADY DETERMINED some other
    way (e.g. VismolObject.find_bonded_and_nonbonded_atoms()'s own fast,
    grid-accelerated c_distances.get_atomic_bonds_from_grid() result) --
    avoids recomputing connectivity a second, slower (O(N²)) time with
    connectivity_from_distance() when a faster method has already found
    it.

    symbols     : list of N element symbols (index i of THIS list is
                  what index_bonds' pair values refer to).
    index_bonds : FLAT list of atom-index pairs, [i0, j0, i1, j1, ...]
                  (same format VismolObject.index_bonds/
                  define_bonds_from_external() already use throughout
                  this app).

    Returns a dict {(i,j): order, ...} (i <= j) -- the SAME "Convention
    C" format VismolObject._bonds_from_pair_of_indexes_list() already
    accepts for its external_orders argument, ready to hand straight to
    define_bonds_from_external(). Deliberately a dict, not a bond_order
    matrix or a positional list: a matrix would need O(N²) memory for
    what's normally a sparse set of bonds on a real molecule, and a
    positional list would reintroduce exactly the "does this line up
    with some OTHER list's order/length" risk this project already found
    and fixed once (see pdynamo/pDynamo2EasyHybrid/session.py's own bond-
    order-extraction fix, same underlying concern). """
    max_valences_by_index = _resolve_max_valences ( symbols, max_valence_table )

    bond_order_by_pair = { }
    for k in range ( 0, len ( index_bonds ) - 1, 2 ):
        i, j = index_bonds[k], index_bonds[k + 1]
        pair = ( i, j ) if i <= j else ( j, i )
        bond_order_by_pair[pair] = 1

    n_bonds = { }
    for pair in bond_order_by_pair:
        n_bonds[pair[0]] = n_bonds.get ( pair[0], 0 ) + 1
        n_bonds[pair[1]] = n_bonds.get ( pair[1], 0 ) + 1

    # --- B.2: promove para DUPLA ---
    for pair in sorted ( bond_order_by_pair.keys ( ) ):
        if bond_order_by_pair[pair] != 1:
            continue
        i, j = pair
        if n_bonds.get ( i, 0 ) < max_valences_by_index[i] and n_bonds.get ( j, 0 ) < max_valences_by_index[j]:
            bond_order_by_pair[pair] = 2
            n_bonds[i] += 1
            n_bonds[j] += 1

    # --- B.3: promove para TRIPLA ---
    for pair in sorted ( bond_order_by_pair.keys ( ) ):
        if bond_order_by_pair[pair] != 2:
            continue
        i, j = pair
        if n_bonds.get ( i, 0 ) < max_valences_by_index[i] and n_bonds.get ( j, 0 ) < max_valences_by_index[j]:
            bond_order_by_pair[pair] = 3
            n_bonds[i] += 1
            n_bonds[j] += 1

    return bond_order_by_pair


def _resolve_max_valences ( symbols, max_valence_table ):
    """ [EN] Shared helper -- builds the per-atom max-valence list (see
    perceive_bond_orders()'s own docstring for the fallback-to-4
    reasoning), indexed the same way as `symbols` itself (a plain list,
    position i == atom_id i -- the same convention index_bonds/
    vismol_object.atoms already use everywhere else in this app). """
    if max_valence_table is None:
        max_valence_table = GABEDIT_MAX_VALENCE
    return [ max_valence_table.get ( symbol.upper ( ), 4 ) for symbol in symbols ]
