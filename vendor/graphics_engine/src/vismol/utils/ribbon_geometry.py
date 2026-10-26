#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ribbon_geometry.py
#
#  Independent Cartoon representation -- GEOMETRY layer (spline, per-
#  residue orientation frames, per-secondary-structure cross-section
#  extrusion into a triangle mesh). See this project's own plan file for
#  the full design rationale. Pure numpy, no GTK/GL dependency -- the
#  Representation class (representations.py) is the only consumer of
#  build_ribbon_mesh_for_object() below.
#
#  [EN] JMOL-like aesthetic (the user's own explicit request): helix and
#  strand are both rendered as WIDE FLAT RIBBONS (not cylinders/tubes --
#  that's the PyMOL/Jmol default "cartoon" look, as opposed to the
#  smoother rounded "putty" tube style some other viewers default to),
#  strand tapers to a triangular arrowhead over its last residue, coil/
#  loop is a slim tube. Colored by secondary-structure type by default
#  (helix=red/magenta, strand=yellow, coil=light grey) -- the classic
#  Richardson convention most viewers, Jmol's own "Structure" colour
#  scheme included, default to.
# ============================================================================
import numpy as np

from vismol.utils.ribbon_backbone import extract_backbone_chains, classify_secondary_structure


# [EN] BUG FIX (user report: "as pontas de helice parecem nao estar se
# conectando corretamente com as estruturas de coil"): coil and helix/
# strand used to be DIFFERENT polygon shapes with DIFFERENT vertex counts
# (an 8-sided octagon for coil, a 4-sided rectangle for helix/strand) --
# build_ribbon_mesh()'s own wall-stitching loop explicitly SKIPPED
# connecting two rings whenever their vertex counts differed (documented
# at the time as a known v1 gap, not silently missed), leaving a real
# hole/seam at every single coil<->helix or coil<->strand boundary in the
# whole protein. Fixed by using this ONE profile shape (an octagon) for
# EVERY secondary-structure type -- only its WIDTH/HEIGHT vary (an
# octagon with width==height reads as a round tube; width >> height reads
# as a flat ribbon, since 6 of its 8 points end up hugging the two long
# flat edges) -- so ring vertex counts always match and every consecutive
# pair of rings now always gets connected, with no shape-topology special
# case left anywhere in the mesh-building loop below.
_RIBBON_PROFILE_SIDES = 8

# [EN] User request (round 5 -- "Preferences > Cartoon" tab, plus a
# per-object "variable criterion" for helix/strand assignment and
# thickness): every number that used to be a hardcoded module constant
# (cross-section half-extents, SS colors, P-SEA tolerance scaling,
# min_run) now lives in this ONE defaults dict instead, and every
# function below takes an optional `settings` dict that overrides
# individual keys on top of these defaults (_resolve_settings() does the
# merge) -- so CartoonRepresentation.__init__ (representations.py) can
# pass in whatever the user set in Preferences (read from
# vm_config.gl_parameters['cartoon_*'], see setup_interface.py's
# __apply_cartoon_parameters()) while every existing caller that passes
# no `settings` at all (the pure-geometry tests, in particular) keeps
# behaving exactly as before, unchanged. Values here are the same
# defaults this file already had hardcoded (half_height 0.18 = the
# "mais planas" flattening request, strand tolerance_scale 1.5 = the
# "tolerancia... muito rigorosa" loosening request, both from earlier
# rounds this session).
DEFAULT_CARTOON_SETTINGS = {
    "coil_half_width":       0.32,
    "coil_half_height":      0.32,
    "helix_half_width":      1.60,
    "helix_half_height":     0.18,
    "strand_half_width":     1.40,
    "strand_half_height":    0.18,
    "arrow_half_width":      2.20,   # widest point of the arrowhead base
    "arrow_tip_half_width":  0.04,   # near-zero, but never exactly 0 (degenerate)
    # [EN] 2026-09-26 aesthetic fix (user report: strand arrowheads looked
    # like sharp "shark fins" compared to PyMOL/VMD's longer, graceful
    # triangular taper) -- see _per_residue_half_extents() below. This is
    # the number of intermediate taper STEPS from the wide arrow base down
    # to the tip (1 = the old single-segment taper, exactly the previous
    # behaviour; 2 spreads it over 2 segments instead, roughly doubling the
    # taper's arc length for a softer point), capped per-run so a short
    # strand never has its ENTIRE length consumed by the arrowhead.
    "arrow_taper_residues":  2,
    # [EN] 2026-09-26 aesthetic fix (user report: strand ribbons looked too
    # "sinuosa" (wavy/sinuous) -- splining directly through the raw Ca
    # trace faithfully shows a real beta strand's own physical PLEAT
    # (Ca positions genuinely alternate slightly above/below the mean
    # strand plane), but PyMOL/VMD-style cartoons deliberately SMOOTH a
    # strand's guide points before splining to turn that into a clean,
    # mostly-straight "arrow" shaft instead -- see
    # _smooth_strand_guide_points() below. 0 = off (the old, unsmoothed
    # behaviour); each iteration is one Laplacian-smoothing pass over the
    # run's own INTERIOR points only (its first/last residues -- the
    # user's own "onde a fita comeca e termina" -- are always held fixed).
    "strand_smoothing_iterations": 4,
    "color_helix":  ( 0.90, 0.20, 0.55 ),   # magenta-ish red
    "color_strand": ( 0.95, 0.85, 0.15 ),   # yellow
    "color_coil":   ( 0.85, 0.85, 0.85 ),   # light grey
    # multiplies ALL 5 P-SEA criteria (angle/dihedral/d2/d3/d4) for that
    # SS type uniformly -- see classify_secondary_structure() in
    # ribbon_backbone.py; 1.0 = the published Labesse et al. values.
    "helix_tolerance_scale":  1.0,
    "strand_tolerance_scale": 1.5,
    # [EN] 2026-09-26 user request ("faz sentido um unico residuo estar em
    # conformacao de fita? o mesmo deve valer para helices") -- split what
    # used to be one shared `min_run` into per-SS-type minimums; see
    # ribbon_backbone.py's own _smooth_short_runs() docstring for why 3
    # and 2 (not equal) are the deliberate, literature-grounded defaults.
    "min_run_helix":  3,
    "min_run_strand": 2,
}


def _resolve_settings ( settings ):
    """ Merges a partial user-supplied settings dict on top of
    DEFAULT_CARTOON_SETTINGS -- every key is optional, missing ones fall
    back to the default. `settings=None` (the default for every function
    below) returns DEFAULT_CARTOON_SETTINGS itself unchanged. """
    if not settings:
        return DEFAULT_CARTOON_SETTINGS
    merged = dict ( DEFAULT_CARTOON_SETTINGS )
    merged.update ( settings )
    return merged


def _color_for_code ( code, settings ):
    return { "H": settings["color_helix"], "S": settings["color_strand"], "C": settings["color_coil"] }[ code ]


# [EN] PERFORMANCE (user report, round 7: "ficou muito lento... como esta
# sendo feito o calculo do cartoon?"): these used to be recomputed from
# scratch (np.linspace + cos/sin) on every single _octagon_profile() call
# -- cheap in isolation, but that function used to run once per SPLINE
# POINT (hundreds to thousands of times per mesh, and, since round 6,
# once per TRAJECTORY FRAME too). Hoisted to module level, computed once
# at import time -- see build_ribbon_mesh()'s own perf note below for
# the much bigger win this was bundled with.
_PROFILE_ANGLES = np.linspace ( 0.0, 2.0 * np.pi, _RIBBON_PROFILE_SIDES, endpoint = False )
_PROFILE_COS = np.cos ( _PROFILE_ANGLES )
_PROFILE_SIN = np.sin ( _PROFILE_ANGLES )


def _octagon_profile ( half_width, half_height ):
    return np.stack ( [ half_width * _PROFILE_COS, half_height * _PROFILE_SIN ], axis = 1 )


def catmull_rom ( points, subdivisions = 5 ):
    """ Catmull-Rom spline through `points` (N, 3) -- the classic closed-
    form 4-point matrix basis (tension fixed at the standard 0.5), NOT a
    port of the Hermite-tangent parameterisation the module being
    replaced used internally; independently written from the well-known
    formula every ribbon/curve renderer in this space shares.

    Endpoints are clamped by duplicating the first/last control point as
    a phantom neighbour (`p[0] + (p[0]-p[1])` / mirror at the far end) --
    the standard, simplest way to define tangents at an open curve's own
    ends without requiring a circular/periodic point set.

    Returns (spline_points, t_values): spline_points is (M, 3); t_values
    is (M,) giving each output point's FRACTIONAL original-point index
    (e.g. 2.4 = 40% of the way from points[2] to points[3]) -- callers
    use this to map a spline point back to "which residue" for secondary-
    structure colouring/cross-section-shape lookups and arrowhead
    tapering, without re-deriving it from scratch at every call site. """
    points = np.asarray ( points, dtype = np.float64 )
    n = len ( points )
    if n < 2:
        return points.copy ( ), np.zeros ( len ( points ), dtype = np.float64 )

    p_ext = np.vstack ( [
        points[0] + ( points[0] - points[1] ),
        points,
        points[-1] + ( points[-1] - points[-2] ),
    ] )

    # [EN] PERFORMANCE (see build_ribbon_mesh()'s own perf note for the
    # full story -- this is one of the two remaining hotspots left after
    # vectorizing that function, per a follow-up cProfile run): the inner
    # per-SUBDIVISION-STEP loop below used to build one 3-vector at a
    # time via plain Python arithmetic. The outer per-SEGMENT loop stays
    # (only n-1 iterations, i.e. one per residue -- already cheap), but
    # each segment's `steps` output points are now computed in a single
    # batched expression (np.outer over the step-parameter array) instead
    # of one Python-level point per step. Same closed-form matrix basis,
    # bit-for-bit identical output.
    out_points = [ ]
    out_t      = [ ]
    for i in range ( 1, n ):
        p0, p1, p2, p3 = p_ext[i - 1], p_ext[i], p_ext[i + 1], p_ext[i + 2]
        # [EN] Include the segment's own END point only on the LAST
        # segment -- otherwise it would be emitted twice (once as this
        # segment's t=1, once as the next segment's t=0).
        steps = subdivisions if i < n - 1 else subdivisions + 1
        t = np.arange ( steps, dtype = np.float64 ) / float ( subdivisions )
        t2 = t * t
        t3 = t2 * t
        seg_points = 0.5 * (
            ( 2.0 * p1 )[None, :]
            + np.outer ( t,  ( -p0 + p2 ) )
            + np.outer ( t2, ( 2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 ) )
            + np.outer ( t3, ( -p0 + 3.0 * p1 - 3.0 * p2 + p3 ) )
        )
        out_points.append ( seg_points )
        out_t.append ( ( i - 1 ) + t )

    return np.concatenate ( out_points, axis = 0 ), np.concatenate ( out_t, axis = 0 )


def _arbitrary_perpendicular ( v ):
    """ Any unit vector perpendicular to `v` -- used only as a last-resort
    fallback when BOTH the O-atom offset and the Ca(i-1)-Ca(i)-Ca(i+1)
    bisector normal degenerate (near-zero length, e.g. 3 collinear Ca
    atoms with no O recorded at all). """
    v = np.asarray ( v, dtype = np.float64 )
    reference = np.array ( [ 1.0, 0.0, 0.0 ] ) if abs ( v[0] ) < 0.9 else np.array ( [ 0.0, 1.0, 0.0 ] )
    perpendicular = np.cross ( v, reference )
    norm = np.linalg.norm ( perpendicular )
    return perpendicular / norm if norm > 1e-6 else np.array ( [ 0.0, 0.0, 1.0 ] )


def _compute_side_vectors ( ca_coords, o_coords ):
    """ One unit "side" (orientation) vector per residue -- the local
    frame's cross-sections are built from this + the spline tangent, so
    an inconsistent/twisting side-vector sequence is exactly what causes
    a ribbon to visibly flip/twist along its length (the classic ribbon-
    rendering failure mode). Primary source: the carbonyl O -> Ca offset
    (the standard Carson & Bugg / Richardson reference vector). Fallback,
    for residues missing an O atom (common in real, incomplete PDBs): the
    Ca(i-1)-Ca(i)-Ca(i+1) bisector-plane normal, via a cross product of
    the two neighbour-facing bond vectors.

    Consistency pass: any side vector whose dot product with the
    PREVIOUS one is negative gets flipped -- standard practice, prevents
    the frame from flickering 180 degrees between adjacent residues (an
    O-atom offset vector, in particular, can point to alternating sides
    of the backbone in a way that's chemically normal but would look like
    constant twisting if used un-corrected). """
    n = len ( ca_coords )
    sides = np.zeros ( ( n, 3 ), dtype = np.float64 )

    for i in range ( n ):
        vector = None
        if o_coords[i] is not None:
            candidate = np.asarray ( o_coords[i], dtype = np.float64 ) - ca_coords[i]
            if np.linalg.norm ( candidate ) > 1e-6:
                vector = candidate

        if vector is None:
            prev_i = max ( i - 1, 0 )
            next_i = min ( i + 1, n - 1 )
            b1 = ca_coords[prev_i] - ca_coords[i]
            b2 = ca_coords[next_i] - ca_coords[i]
            candidate = np.cross ( b1, b2 )
            vector = candidate if np.linalg.norm ( candidate ) > 1e-6 else _arbitrary_perpendicular ( b2 - b1 )

        norm = np.linalg.norm ( vector )
        sides[i] = vector / norm if norm > 1e-6 else np.array ( [ 1.0, 0.0, 0.0 ] )

    for i in range ( 1, n ):
        if np.dot ( sides[i], sides[i - 1] ) < 0.0:
            sides[i] = -sides[i]

    return sides


def _interpolate_side_vectors ( sides, t_values ):
    """ Normalized-linear-interpolation ("nlerp") of the per-residue side
    vectors at each spline t_value -- simpler than a full spherical
    interpolation (slerp), and sufficient here since the flip-consistency
    pass above already keeps consecutive residues' side vectors within a
    small angle of each other in any real structure.

    [EN] PERFORMANCE: vectorized over the WHOLE t_values array at once
    (one of the two remaining hotspots after build_ribbon_mesh()'s own
    vectorization -- see its perf note) instead of a Python-level call
    per spline point; same floor+clip+nlerp+fallback semantics as
    before, bit-for-bit. """
    n = len ( sides )
    i0 = np.clip ( np.floor ( t_values ).astype ( np.int64 ), 0, n - 1 )
    i1 = np.clip ( i0 + 1, 0, n - 1 )
    frac = ( t_values - i0 )[:, None]
    vector = sides[i0] * ( 1.0 - frac ) + sides[i1] * frac
    norms = np.linalg.norm ( vector, axis = 1 )
    degenerate = norms <= 1e-6
    out = vector / np.where ( degenerate, 1.0, norms )[:, None]
    out[degenerate] = sides[i0][degenerate]
    return out


def _find_ss_runs ( ss_codes ):
    """ Returns [(code, start_residue_index, end_residue_index), ...] --
    end is INCLUSIVE, matching how the arrowhead-taper lookup below wants
    to ask "is residue i within the last residue of ITS OWN strand run". """
    runs = [ ]
    n = len ( ss_codes )
    i = 0
    while i < n:
        j = i
        while j < n and ss_codes[j] == ss_codes[i]:
            j += 1
        runs.append ( ( ss_codes[i], i, j - 1 ) )
        i = j
    return runs


def _smooth_strand_guide_points ( ca_coords, ss_codes, iterations ):
    """ [EN] 2026-09-26 (user report): "acho que seria mais interessante
    que ela tivesse um formato menos sinuosa, em formato de flecha...
    os carbonos alfa onde a fita comeca e termina sao muito importantes
    para gerar as fitas em forma de flechas, que as vezes sao torcidas,
    mas nunca sinuosas." Exactly the distinction every mainstream
    "PyMOL-style" cartoon renderer (Carson & Bugg 1986 and every
    descendant) makes: a beta strand's raw Ca trace genuinely zigzags
    (the real physical pleat of a beta sheet), and splining straight
    through it -- which is what this module did before this fix, and
    still does for helix/coil -- faithfully reproduces that zigzag as a
    visibly wavy ribbon. PyMOL/VMD instead SMOOTH the guide points used
    for a strand's own centerline before splining, producing a clean,
    mostly-straight "arrow" shaft. Deliberately position-only: a strand's
    TWIST comes from _compute_side_vectors(), called on the ORIGINAL
    (unsmoothed) ca_coords/o_coords -- so an arrow can still end up
    "torcida" (twisted, real chemistry, left alone) while never being
    "sinuosa" (a wavy centerline, smoothed away here), matching the
    user's own distinction exactly.

    Iterative Laplacian smoothing of a strand run's INTERIOR guide points
    ONLY -- run[0] and run[-1] (the residues "onde a fita comeca e
    termina") are held perfectly fixed at every iteration, so the ribbon
    still threads through the true start/end Ca position and connects
    seamlessly to whatever coil/helix precedes or follows it. Each
    iteration nudges every interior point toward the average of its two
    immediate neighbours (a classic, minimal Laplacian smoothing pass) --
    a few iterations noticeably straighten the strand while a fixed
    start/end anchor keeps it from ever collapsing onto a dead-straight
    line (real sheets do have their own gentle large-scale curvature/
    twist, which this preserves; only the small, residue-to-residue pleat
    wiggle is removed, since Laplacian smoothing damps high-frequency
    variation much faster than the strand's own overall trend).
    `iterations=0` reproduces the OLD unsmoothed behaviour exactly (a
    plain copy, no-op). """
    smoothed = np.array ( ca_coords, dtype = np.float64, copy = True )
    if iterations <= 0:
        return smoothed
    runs = _find_ss_runs ( ss_codes )
    for code, start, end in runs:
        if code != "S" or ( end - start ) < 2:
            continue   # < 3 residues -- no interior point to smooth
        for _ in range ( iterations ):
            prev_pts = smoothed[start:end - 1]
            next_pts = smoothed[start + 2:end + 1]
            interior = smoothed[start + 1:end]
            smoothed[start + 1:end] = 0.5 * interior + 0.25 * prev_pts + 0.25 * next_pts
    return smoothed


def _per_residue_half_extents ( ss_codes, runs_by_residue, settings ):
    """ Precomputes ONE (half_width, half_height) TARGET per residue --
    the arrowhead is expressed here too: a strand run's own last
    `n_taper` residues (see below) ramp smoothly from the arrow's widest
    BASE extent down to a near-zero TIP extent at the very last residue.

    build_ribbon_mesh() below then LINEARLY INTERPOLATES between two
    adjacent residues' targets at every spline point in between (see
    _interpolate_half_extents(), same nlerp-style scheme
    _interpolate_side_vectors() already uses for orientation) instead of
    switching abruptly at each residue boundary the way the previous
    per-spline-point lookup did.

    [EN] BUG FIX (user report: helix/coil connection gap, and "a
    tolerancia para fita esta muito rigorosa" prompted re-checking this
    whole taper mechanism too): this ALSO smooths the strand's own
    arrowhead taper and, more importantly, every coil<->helix/strand
    boundary -- previously an instant width jump (e.g. coil's 0.32 straight
    to helix's 1.60 on the very next spline sub-point); now spread evenly
    across the one spline segment between those two residues, matching
    the topology fix above (build_ribbon_mesh()'s wall loop no longer
    needs to skip anything, since a constant octagon vertex COUNT was
    already guaranteed by using _octagon_profile() everywhere -- this
    additional interpolation is what makes the WIDTH itself flare/taper
    smoothly too, rather than merely being topologically connected but
    still visually abrupt).

    [EN] 2026-09-26 aesthetic fix (user report: "nao estou gostando do
    resultado estetico das fitas beta" -- compared against PyMOL/VMD).
    Grounded in a real screenshot of ubiquitin (examples/pdb1ubq.pdb,
    strand run lengths measured as 5/7/3/3/8 residues): the OLD taper
    above (widest BASE at the run's second-to-last residue, near-zero TIP
    at its last) crammed the entire base->tip transition into a SINGLE
    spline segment (~3.8 A of arc length) -- visually a sharp "shark fin"
    kink rather than PyMOL/VMD's longer, graceful triangular arrowhead,
    made worse on short runs (a 3-residue run had 2 of its 3 residues
    consumed by the arrowhead, leaving almost no flat "body" at all).
    Fixed by spreading the base->tip ramp over `arrow_taper_residues`
    SEGMENTS instead of always exactly 1 -- `n_taper` below is that value
    capped at `max(1, run_length // 2)` so a short run never has more than
    half its own length consumed by the arrowhead (this cap also means
    `arrow_taper_residues=2`'s new default reproduces the OLD 1-segment
    behaviour bit-for-bit on any 2- or 3-residue run, where
    `run_length // 2` is already 1 -- only strands long enough to afford a
    longer taper actually get one). The widen-then-narrow "shoulder" where
    the flat body meets the arrow's wide base is still a single abrupt
    step, deliberately -- that hard base is a real, visible feature of
    PyMOL/VMD's own classic arrowhead too, not something to smooth away. """
    n = len ( ss_codes )
    half_widths  = np.zeros ( n, dtype = np.float64 )
    half_heights = np.zeros ( n, dtype = np.float64 )
    for i in range ( n ):
        code = ss_codes[i]
        if code == "H":
            half_widths[i]  = settings["helix_half_width"]
            half_heights[i] = settings["helix_half_height"]
        elif code == "S":
            run_start, run_end = runs_by_residue[i]
            run_length = run_end - run_start + 1
            n_taper = min ( settings["arrow_taper_residues"], max ( 1, run_length // 2 ) )
            dist_from_end = run_end - i   # 0 at the run's own last residue
            if run_length > 1 and dist_from_end <= n_taper:
                t = dist_from_end / float ( n_taper )   # 0 at the tip, 1 at the base
                base_w = settings["arrow_half_width"]
                tip_w  = settings["arrow_tip_half_width"]
                half_widths[i]  = tip_w + t * ( base_w - tip_w )
                half_heights[i] = settings["strand_half_height"]
            else:
                half_widths[i]  = settings["strand_half_width"]
                half_heights[i] = settings["strand_half_height"]
        else:
            half_widths[i]  = settings["coil_half_width"]
            half_heights[i] = settings["coil_half_height"]
    return half_widths, half_heights


def _interpolate_half_extents ( half_widths, half_heights, t ):
    n = len ( half_widths )
    i0 = int ( np.floor ( t ) )
    i0 = max ( 0, min ( i0, n - 1 ) )
    i1 = min ( i0 + 1, n - 1 )
    frac = t - i0
    half_width  = half_widths[i0]  * ( 1.0 - frac ) + half_widths[i1]  * frac
    half_height = half_heights[i0] * ( 1.0 - frac ) + half_heights[i1] * frac
    return half_width, half_height


def _interpolate_half_extents_batch ( half_widths, half_heights, t_values ):
    """ Vectorized form of _interpolate_half_extents() above -- ONE call
    over the WHOLE t_values array (every spline point at once) instead of
    one Python-level call per point. Bit-for-bit the same floor+clip+
    nlerp semantics, just batched -- see build_ribbon_mesh()'s own perf
    note for why this matters. The scalar version above is kept as-is
    (still used by the end-cap code, which only ever needs 2 points, and
    by this project's own earlier verification scripts). """
    n = len ( half_widths )
    i0 = np.clip ( np.floor ( t_values ).astype ( np.int64 ), 0, n - 1 )
    i1 = np.clip ( i0 + 1, 0, n - 1 )
    frac = t_values - i0
    half_width  = half_widths[i0]  * ( 1.0 - frac ) + half_widths[i1]  * frac
    half_height = half_heights[i0] * ( 1.0 - frac ) + half_heights[i1] * frac
    return half_width, half_height


def _colors_for_codes_batch ( ss_codes, settings ):
    """ Vectorized _color_for_code() -- builds the (n_residues, 3) color
    array with 3 boolean-mask assignments total, instead of one Python
    dict lookup per residue. """
    ss_codes = np.asarray ( ss_codes )
    out = np.empty ( ( len ( ss_codes ), 3 ), dtype = np.float64 )
    out[ ss_codes == "H" ] = settings["color_helix"]
    out[ ss_codes == "S" ] = settings["color_strand"]
    out[ ss_codes == "C" ] = settings["color_coil"]
    return out


def build_ribbon_mesh ( chain_record, ss_codes, spline_detail = 5, settings = None, residue_colors = None ):
    """ Builds the triangle mesh for ONE chain record (see extract_
    backbone_chains()) -- coords/normals/colors are per-vertex, `indices`
    are triangle-list indices (3 per triangle) INTO those per-vertex
    arrays, all already offset-adjusted for this record alone (callers
    combining several records must offset `indices` by their own running
    vertex count -- see build_ribbon_mesh_for_object() below, which does
    exactly that).

    Vertex-normal history (4 rounds of real, user-reported/self-found
    lighting bugs, each verified numerically, not just by eye):
      1) originally one shared, radially-AVERAGED normal per profile
         POINT -- visibly wrong on a flat ribbon face (a diagonal light/
         dark gradient across what should be flat), fixed by duplicating
         vertices per EDGE with an analytic per-edge outward normal;
      2) that analytic normal turned out to be anti-parallel to the
         actual triangle winding for the wall geometry (and, separately,
         for one of the two end caps) -- fixed by aligning signs to the
         real winding order;
      3) even sign-aligned, a single analytic normal per ring is only
         correct for a near-planar quad -- a sharp backbone bend could
         still twist a quad enough that its real geometric normal
         diverged from (occasionally even inverted against) the ring's
         own analytic normal (found via this fix's own whole-mesh
         winding-vs-normal sweep, not a user report: 15/3616 triangles on
         a real 46-residue protein, all at sharp coil turns). Fixed by
         giving every WALL triangle its own 3 fresh, never-shared
         vertices with a normal computed directly from ITS OWN 3D
         corners -- winding and stored normal identical by construction;
      4) that FIXED the sign, but "never shared" also means every
         triangle facet renders hard-edged against its neighbours -- an
         8-sided octagon then visibly reads as faceted rather than a
         smooth tube/ribbon (user report: "as faces de cada poligono
         estao evidentes... as normais nao estao combinadas para
         amenizar as faces"). Fixed by going back to ONE shared vertex
         per (ring, profile-point) -- see the wall-mesh block below --
         with that shared vertex's normal set to the normalized SUM of
         the true geometric face-normals of every triangle touching it
         (standard "smooth shading group" technique). This still uses
         round 3's own correctly-signed per-triangle face normals as its
         only input (summing same-signed vectors can't flip a sign), so
         it keeps round 3's fix while additionally blending across
         shared corners for a smooth Gouraud look; a genuinely flat run
         of faces (the ribbon's own big top/bottom) still reads flat,
         since its own incident face normals are already near-identical
         and averaging them barely moves anything -- only real corners
         get rounded off, which is the desired look. End caps were never
         affected by any of this (their own fan winding was verified
         consistent back in round 2) and remain unchanged.

    `residue_colors` (2026-09-28, per-representation colour schemes):
    optional (n_residues, 3) RGB array, one row per residue of THIS
    record, used instead of the settings' SS-type colours -- see
    residue_colors_for_records() below. None keeps the SS colouring. """
    settings = _resolve_settings ( settings )
    ca_coords = chain_record["ca"]
    o_coords  = chain_record["o"]
    n_residues = len ( ca_coords )

    # [EN] 2026-09-26: the SPLINE CENTERLINE is built from smoothed strand
    # guide points (see _smooth_strand_guide_points()'s own docstring for
    # the full "sinuosa vs flecha" story) -- but orientation/twist below
    # still reads the ORIGINAL, unsmoothed ca_coords/o_coords, since twist
    # is real per-residue chemistry the user explicitly wants preserved
    # ("as vezes torcidas"), not something to smooth away.
    smoothed_ca_coords = _smooth_strand_guide_points (
        ca_coords, ss_codes, settings["strand_smoothing_iterations"] )
    spline_points, t_values = catmull_rom ( smoothed_ca_coords, subdivisions = spline_detail )
    side_per_residue = _compute_side_vectors ( ca_coords, o_coords )
    side_per_point    = _interpolate_side_vectors ( side_per_residue, t_values )

    runs = _find_ss_runs ( ss_codes )
    runs_by_residue = [ None ] * n_residues
    for code, start, end in runs:
        for residue_index in range ( start, end + 1 ):
            runs_by_residue[residue_index] = ( start, end )
    half_widths, half_heights = _per_residue_half_extents ( ss_codes, runs_by_residue, settings )

    m = len ( spline_points )

    # [EN] PERFORMANCE (user report, round 7: "ficou muito lento... como
    # esta sendo feito o calculo do cartoon?"): round 6 made this whole
    # function re-run on EVERY trajectory/MD frame, which turned a
    # pre-existing inefficiency into an actually noticeable stall.
    # Profiled (cProfile) on a real protein: ~100ms for a tiny 46-residue
    # chain, over 1 SECOND for a 500-residue one, PER rebuild -- almost
    # entirely numpy per-call dispatch overhead, from calling np.cross()/
    # np.linalg.norm() ONE SPLINE POINT OR ONE WALL TRIANGLE AT A TIME
    # inside a Python for loop (np.cross(), called once per wall
    # triangle in the old normal-accumulation loop, alone accounted for
    # >60% of total time -- tens of thousands of individual 3-element
    # cross products, each paying numpy's full array-function dispatch
    # overhead for work that's microseconds' worth of actual math).
    #
    # Below is the SAME algorithm (same maths, same winding, same vertex
    # layout/order as before -- re-verified byte-identical against the
    # pre-vectorization version on a real protein), just computed with
    # numpy operating on WHOLE ARRAYS (all m spline points, or all wall
    # triangles) at once instead of one Python-level call per point/
    # triangle. This is the single highest-value fix -- see this
    # project's own memory notes for the full before/after benchmark.
    tangents = np.zeros ( ( m, 3 ), dtype = np.float64 )
    tangents[:-1] = spline_points[1:] - spline_points[:-1]
    tangents[-1]  = tangents[-2] if m > 1 else np.array ( [ 1.0, 0.0, 0.0 ] )
    tangent_norms = np.linalg.norm ( tangents, axis = 1 )
    degenerate_tangent = tangent_norms <= 1e-6
    tangents = tangents / np.where ( degenerate_tangent, 1.0, tangent_norms )[:, None]
    tangents[degenerate_tangent] = np.array ( [ 1.0, 0.0, 0.0 ] )

    # re-orthogonalise side against tangent (nlerp doesn't guarantee
    # perpendicularity) and derive the local "up" axis from both -- for
    # every spline point at once.
    side = side_per_point - np.einsum ( 'ij,ij->i', side_per_point, tangents )[:, None] * tangents
    side_norms = np.linalg.norm ( side, axis = 1 )
    degenerate_side = side_norms <= 1e-6
    side = side / np.where ( degenerate_side, 1.0, side_norms )[:, None]
    if np.any ( degenerate_side ):
        # rare (side exactly parallel to tangent) -- fall back per-point,
        # same as the pre-vectorization code; not worth batching since it
        # essentially never triggers on a real backbone.
        for k in np.nonzero ( degenerate_side )[0]:
            side[k] = _arbitrary_perpendicular ( tangents[k] )
    normal_axis = np.cross ( tangents, side )

    half_widths_pt, half_heights_pt = _interpolate_half_extents_batch ( half_widths, half_heights, t_values )

    residue_idx_for_color = np.clip ( np.round ( t_values ).astype ( np.int64 ), 0, n_residues - 1 )
    if residue_colors is None:
        residue_colors = _colors_for_codes_batch ( ss_codes, settings )
    residue_colors = np.asarray ( residue_colors, dtype = np.float64 )
    colors_per_point = residue_colors[ residue_idx_for_color ]

    # [EN] side-wall triangles between consecutive rings. Every ring uses
    # the SAME octagon topology (see _RIBBON_PROFILE_SIDES's own
    # docstring above -- the fix for the user-reported "pontas de helice
    # nao conectando com coil" gap), so every consecutive pair of rings
    # always connects, with no shape-mismatch case to special-case here.
    #
    # [EN] SMOOTH-SHADING PASS (user report, round 4: "as faces de cada
    # poligono estao evidentes... as normais nao estao combinadas para
    # amenizar as faces"): round 3 gave every wall TRIANGLE its own 3
    # fresh, never-shared vertices with a normal taken straight from that
    # triangle's own winding -- this made the winding-vs-normal alignment
    # perfect (see the bug fix 2 paragraphs up), but it also means EVERY
    # triangle facet is hard-edged against its neighbours, so a
    # coarse-in-cross-section shape (an 8-sided octagon) reads as visibly
    # faceted rather than a smoothly rounded tube/ribbon -- exactly the
    # complaint. Standard fix (the same technique every 3D modelling tool
    # calls "smooth shading group"): go back to ONE shared vertex per
    # (ring, profile-point) -- so a corner is the SAME vertex for both
    # triangles that meet there -- then give that vertex a SINGLE normal
    # equal to the normalized SUM of the true geometric face-normals of
    # every triangle touching it. This keeps the exact per-triangle
    # geometric normals round 3 fixed (no analytic formula, no possible
    # winding-sign disagreement -- a sum of already-consistent triangle
    # normals can't flip sign), just BLENDS them at shared corners so
    # Gouraud interpolation produces a smooth gradient across each facet
    # boundary instead of a hard edge -- a genuinely flat run of faces
    # (the ribbon's own big top/bottom faces) still reads as flat, since
    # their own face normals are already nearly identical and averaging
    # them changes almost nothing; only the polygon's actual corners
    # (where adjacent faces really do meet at an angle) get rounded off,
    # which is exactly the desired look.
    n_edges = _RIBBON_PROFILE_SIDES

    # (m, n_edges, 3) -- every ring's own n_edges profile points, ALL at
    # once via broadcasting (replaces the old per-(k,i) nested loop).
    # Flattened to (m*n_edges, 3) in ring-major/edge-minor order, exactly
    # matching the old nested loop's append order (k outer, i inner).
    offsets = (
        ( half_widths_pt[:, None]  * _PROFILE_COS[None, :] )[:, :, None] * side[:, None, :]
        + ( half_heights_pt[:, None] * _PROFILE_SIN[None, :] )[:, :, None] * normal_axis[:, None, :]
    )
    wall_coords_arr = ( spline_points[:, None, :] + offsets ).reshape ( m * n_edges, 3 )
    wall_colors_arr = np.repeat ( colors_per_point, n_edges, axis = 0 )

    n_pairs = m - 1
    if n_pairs > 0:
        # wall triangle indices -- (m-1) ring pairs * n_edges edges * 2
        # triangles, built via array arithmetic instead of a nested
        # Python loop (interleaved tri1/tri2 to match the original
        # append order exactly, though the renderer itself doesn't care
        # about triangle order).
        k_idx  = np.repeat ( np.arange ( n_pairs ), n_edges )
        e_idx  = np.tile ( np.arange ( n_edges ), n_pairs )
        e2_idx = ( e_idx + 1 ) % n_edges
        start_a = k_idx * n_edges
        start_b = ( k_idx + 1 ) * n_edges
        v0 = start_a + e_idx
        v1 = start_a + e2_idx
        v2 = start_b + e_idx
        v3 = start_b + e2_idx
        wall_tris_arr = np.empty ( ( 2 * n_pairs * n_edges, 3 ), dtype = np.int64 )
        wall_tris_arr[0::2] = np.stack ( [ v0, v2, v1 ], axis = 1 )
        wall_tris_arr[1::2] = np.stack ( [ v1, v2, v3 ], axis = 1 )

        # ONE batched cross product over every wall triangle at once --
        # this single call replaces what used to be tens of thousands of
        # individual np.cross() calls, the actual bottleneck found above.
        pa = wall_coords_arr[ wall_tris_arr[:, 0] ]
        pb = wall_coords_arr[ wall_tris_arr[:, 1] ]
        pc = wall_coords_arr[ wall_tris_arr[:, 2] ]
        face_normals = np.cross ( pb - pa, pc - pa )
        face_norms = np.linalg.norm ( face_normals, axis = 1 )
        degenerate_face = face_norms <= 1e-9
        face_normals = face_normals / np.where ( degenerate_face, 1.0, face_norms )[:, None]
        face_normals[degenerate_face] = 0.0

        # scatter-add each triangle's face normal into its 3 vertices --
        # np.add.at handles the repeated-index accumulation correctly
        # (unlike plain fancy-index assignment), same result as the old
        # per-triangle "normal_accum[i] += face_normal" loop.
        normal_accum = np.zeros_like ( wall_coords_arr )
        np.add.at ( normal_accum, wall_tris_arr[:, 0], face_normals )
        np.add.at ( normal_accum, wall_tris_arr[:, 1], face_normals )
        np.add.at ( normal_accum, wall_tris_arr[:, 2], face_normals )
    else:
        wall_tris_arr = np.zeros ( ( 0, 3 ), dtype = np.int64 )
        normal_accum = np.zeros_like ( wall_coords_arr )

    accum_norms = np.linalg.norm ( normal_accum, axis = 1, keepdims = True )
    accum_norms[accum_norms < 1e-9] = 1.0
    wall_normals_arr = normal_accum / accum_norms

    # end caps (real chain termini only -- ring 0 and ring m-1) -- built
    # from the profile's own plain (non-duplicated) points, one fresh
    # vertex per corner, all sharing the SAME flat cap normal (+-tangent)
    # -- entirely separate from the per-edge wall vertices above. Only 2
    # rings' worth of work (never more), so left as plain Python lists --
    # not a meaningful cost next to the wall pass above, not worth the
    # extra complexity of vectorizing.
    wall_offset = len ( wall_coords_arr )
    cap_coords  = [ ]
    cap_normals = [ ]
    cap_colors  = [ ]
    cap_indices = [ ]

    for k in ( 0, m - 1 ):
        side_k, normal_axis_k = side[k], normal_axis[k]
        half_width, half_height = _interpolate_half_extents ( half_widths, half_heights, t_values[k] )
        profile = _octagon_profile ( half_width, half_height )
        count = len ( profile )
        if count < 3:
            continue
        color = residue_colors[ int ( round ( np.clip ( t_values[k], 0, n_residues - 1 ) ) ) ]
        # [EN] BUG FIX (found alongside the wall-normal sign fix above, by
        # the same winding-vs-normal dot-product check): the two caps use
        # OPPOSITE fan winding (`winding` below, reversed for k == m-1, so
        # both caps' triangles face outward despite sharing the same index
        # pattern) -- the normal must follow that SAME per-cap sign, not a
        # single formula that (as originally written here) happened to
        # reduce to "-tangents[k]" for BOTH branches. Confirmed
        # numerically: with the winding as-is, the k == 0 cap needs
        # +tangents[0], the k == m-1 cap needs -tangents[m-1].
        cap_normal = tangents[k] if k == 0 else -tangents[k]

        start = wall_offset + len ( cap_coords )
        for ( a, b ) in profile:
            world = spline_points[k] + a * side_k + b * normal_axis_k
            cap_coords.append ( world )
            cap_normals.append ( cap_normal )
            cap_colors.append ( color )
        local_start = start - wall_offset
        center = np.mean ( np.array ( cap_coords[local_start:local_start + count] ), axis = 0 )
        center_index = wall_offset + len ( cap_coords )
        cap_coords.append ( center )
        cap_normals.append ( cap_normal )
        cap_colors.append ( color )

        winding = range ( count ) if k == 0 else range ( count - 1, -1, -1 )
        winding = list ( winding )
        for idx in range ( count ):
            e0 = winding[idx]
            e1 = winding[ ( idx + 1 ) % count ]
            cap_indices.extend ( [ center_index, start + e0, start + e1 ] )

    if cap_coords:
        coords  = np.concatenate ( [ wall_coords_arr,  np.asarray ( cap_coords,  dtype = np.float64 ) ], axis = 0 )
        normals = np.concatenate ( [ wall_normals_arr, np.asarray ( cap_normals, dtype = np.float64 ) ], axis = 0 )
        colors  = np.concatenate ( [ wall_colors_arr,  np.asarray ( cap_colors,  dtype = np.float64 ) ], axis = 0 )
    else:
        coords, normals, colors = wall_coords_arr, wall_normals_arr, wall_colors_arr

    if n_pairs > 0 or cap_indices:
        indices = np.concatenate ( [ wall_tris_arr.reshape ( -1 ), np.asarray ( cap_indices, dtype = np.uint32 ) ] )
    else:
        indices = np.zeros ( ( 0, ), dtype = np.uint32 )

    return coords.astype ( np.float32 ), normals.astype ( np.float32 ), colors.astype ( np.float32 ), indices.astype ( np.uint32 )


# [EN] 2026-09-28 user request ("o usuario pode escolher o padrao de cores
# para o cartoon (cada representacao pode ter um padrao diferente)... a
# possibilidade de colorir tipo rainbow"). The colour scheme is a property
# of each CartoonRepresentation (rep.color_mode), passed down here --
# NOT a global gl_parameters key -- so two objects can show different
# schemes at the same time. Same short, hand-picked chain cycle the
# molecular-surface tool already uses for its own "Chain" mode (util/
# molecular_surface.py's _CHAIN_COLOR_CYCLE -- duplicated, not imported:
# that module lives in the EasyHybrid tree, outside the vismol package).
CARTOON_COLOR_MODES = ( "ss", "rainbow", "chain", "atom", "uniform" )

_CHAIN_COLOR_CYCLE = (
    ( 0.90, 0.30, 0.30 ), ( 0.30, 0.55, 0.90 ), ( 0.35, 0.75, 0.35 ), ( 0.95, 0.70, 0.20 ),
    ( 0.65, 0.35, 0.85 ), ( 0.30, 0.80, 0.80 ), ( 0.90, 0.45, 0.70 ), ( 0.55, 0.55, 0.55 ),
)


def rainbow_colors ( fractions ):
    """ Blue (N-terminus, fraction 0) -> cyan -> green -> yellow -> red
    (C-terminus, fraction 1) -- the usual "rainbow"/"spectrum" cartoon
    convention (PyMOL's spectrum, Chimera's rainbow): a plain HSV hue
    sweep from 240 deg down to 0 deg at full saturation/value. """
    f = np.clip ( np.asarray ( fractions, dtype = np.float64 ), 0.0, 1.0 )
    h = ( 1.0 - f ) * 4.0            # hue in sextants: 4 = blue ... 0 = red
    x = 1.0 - np.abs ( np.mod ( h, 2.0 ) - 1.0 )
    out = np.zeros ( ( len ( f ), 3 ), dtype = np.float64 )
    s0 = h < 1.0                     # red -> yellow
    s1 = ( h >= 1.0 ) & ( h < 2.0 )  # yellow -> green
    s2 = ( h >= 2.0 ) & ( h < 3.0 )  # green -> cyan
    s3 = h >= 3.0                    # cyan -> blue
    out[s0] = np.stack ( [ np.ones ( s0.sum ( ) ), x[s0], np.zeros ( s0.sum ( ) ) ], axis = 1 )
    out[s1] = np.stack ( [ x[s1], np.ones ( s1.sum ( ) ), np.zeros ( s1.sum ( ) ) ], axis = 1 )
    out[s2] = np.stack ( [ np.zeros ( s2.sum ( ) ), np.ones ( s2.sum ( ) ), x[s2] ], axis = 1 )
    out[s3] = np.stack ( [ np.zeros ( s3.sum ( ) ), x[s3], np.ones ( s3.sum ( ) ) ], axis = 1 )
    return out


def residue_colors_for_records ( chain_records, ss_codes_per_record, settings = None,
                                 color_mode = "ss", uniform_color = None ):
    """ One (n_residues, 3) RGB array per chain record, for `color_mode`:
        "ss"      -- by secondary structure (the Preferences > Cartoon
                     helix/strand/coil colours; the original behaviour);
        "rainbow" -- blue->red from N- to C-terminus, PER CHAIN: the
                     position is counted over ALL of that chain's records
                     (a chain split by a sequence gap still reads as one
                     continuous spectrum), computed BEFORE any partial-
                     display filtering, so hiding part of a chain never
                     shifts the colours of what stays visible;
        "chain"   -- one colour per chain (_CHAIN_COLOR_CYCLE);
        "atom"    -- each residue's own Ca atom colour (atom.color), so a
                     system colour / "Change Color" on the selection is
                     reflected on the ribbon;
        "uniform" -- one flat colour (`uniform_color`, RGB 0-1).
    Unknown modes fall back to "ss". """
    settings = _resolve_settings ( settings )
    if color_mode == "rainbow":
        totals = { }
        for record in chain_records:
            key = id ( record["chain"] )
            totals[key] = totals.get ( key, 0 ) + len ( record["ca"] )
        seen = { }
        out = [ ]
        for record in chain_records:
            key = id ( record["chain"] )
            n = len ( record["ca"] )
            start = seen.get ( key, 0 )
            seen[key] = start + n
            denominator = max ( totals[key] - 1, 1 )
            out.append ( rainbow_colors ( ( start + np.arange ( n ) ) / denominator ) )
        return out
    if color_mode == "chain":
        chain_slot = { }
        out = [ ]
        for record in chain_records:
            key = id ( record["chain"] )
            if key not in chain_slot:
                chain_slot[key] = len ( chain_slot )
            color = _CHAIN_COLOR_CYCLE[ chain_slot[key] % len ( _CHAIN_COLOR_CYCLE ) ]
            out.append ( np.tile ( np.asarray ( color, dtype = np.float64 ), ( len ( record["ca"] ), 1 ) ) )
        return out
    if color_mode == "atom":
        return [ np.asarray ( [ np.asarray ( atom.color, dtype = np.float64 ).ravel ( )[:3]
                                for atom in record["ca_atoms"] ], dtype = np.float64 )
                 for record in chain_records ]
    if color_mode == "uniform":
        color = np.asarray ( uniform_color if uniform_color is not None else settings["color_coil"],
                             dtype = np.float64 ).ravel ( )[:3]
        return [ np.tile ( color, ( len ( record["ca"] ), 1 ) ) for record in chain_records ]
    return [ _colors_for_codes_batch ( codes, settings ) for codes in ss_codes_per_record ]


def split_record_by_mask ( record, mask ):
    """ [EN] 2026-09-28 user request ("representar apenas parte da
    estrutura na forma de cartoon, uma vez selecionada e atribuida pelo
    usuario"): splits one chain record into the contiguous runs where
    `mask` (one bool per residue) is True. Returns a list of
    (sub_record, start, end) with end exclusive, so the caller can slice
    that record's own per-residue arrays (SS codes, colours) with the
    same bounds. Runs shorter than 2 residues are dropped -- same rule as
    extract_backbone_chains() (nothing to spline through). Deliberately
    applied AFTER secondary-structure classification of the FULL chain:
    P-SEA is a windowed Ca criterion, so classifying only the visible
    fragment would mislabel residues near every cut (half a helix shown
    must still be drawn as helix). """
    mask = np.asarray ( mask, dtype = bool )
    pieces = [ ]
    n = len ( mask )
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < n and mask[j]:
            j += 1
        if j - i >= 2:
            sub = dict ( record )
            sub["residues"] = record["residues"][i:j]
            sub["ca"]       = record["ca"][i:j]
            sub["o"]        = record["o"][i:j]
            sub["ca_atoms"] = record["ca_atoms"][i:j]
            pieces.append ( ( sub, i, j ) )
        i = j
    return pieces


def build_ribbon_mesh_for_object ( vismol_object, frame = None, spline_detail = 5, settings = None,
                                    cached_ss_codes = None, return_ss_codes = False,
                                    color_mode = "ss", uniform_color = None, ca_filter = None ):
    """ Single entry point CartoonRepresentation calls -- analogous in
    ROLE to the old cartoon.cartoon(vismol_object, ...) one-liner it
    replaces, but combining extract_backbone_chains() + classify_
    secondary_structure() + build_ribbon_mesh() over EVERY chain record
    in the object, concatenating each record's mesh into one flat set of
    buffers (index values offset by the running vertex count so far).
    Returns (coords, normals, colors, indices), all numpy arrays, ready
    to flatten and upload exactly like the representation being
    replaced's own coords2/colors2/normals2/indexes2 did.

    `settings` (optional, see DEFAULT_CARTOON_SETTINGS above): resolved
    ONCE here and passed down to both classify_secondary_structure() (its
    tolerance-scale/min_run keys) and build_ribbon_mesh() (its cross-
    section/color keys), so a caller building several chains gets
    consistent settings across all of them from a single resolve.

    [EN] User request (round 8 -- "Preferences: cartoon dinamico ou nao"):
    `cached_ss_codes`, if given, is a list of per-chain H/S/C code arrays
    (in the SAME order extract_backbone_chains() returns chain records --
    the same order every call gets for a given object/topology) to REUSE
    instead of re-running classify_secondary_structure() for that chain --
    the point being that P-SEA classification is itself now a comparable
    cost to the whole rest of the mesh build (~45% of total time on a
    real protein, once build_ribbon_mesh() itself got vectorized -- see
    this project's own perf notes), and re-running it on EVERY trajectory
    frame is both wasted work AND a source of visible flicker for a
    single, structurally-rigid protein whose SS genuinely doesn't change
    frame to frame in the overwhelming majority of real use. A cached
    entry is used ONLY if it's non-None AND its own length matches that
    chain's CURRENT residue count (a structural/topology change --
    unusual mid-trajectory, but not impossible -- falls back to a fresh
    classification rather than silently misapplying a stale, wrong-length
    array). `return_ss_codes=True` additionally returns the per-chain
    codes actually used (freshly computed OR reused) as a 5th value, so
    the caller can cache them for the NEXT call -- default False keeps
    the original 4-tuple return for every existing caller unchanged.

    [EN] 2026-09-28: `color_mode`/`uniform_color` -- see
    residue_colors_for_records(). `ca_filter` -- optional callable
    (Ca Atom -> bool); only residues whose Ca it accepts are meshed
    (partial cartoon, see split_record_by_mask()). Classification, the SS cache
    and rainbow positions are all still computed over the FULL chain, so
    changing what is shown never changes how the rest is drawn. """
    settings = _resolve_settings ( settings )
    chain_records = extract_backbone_chains ( vismol_object, frame = frame )

    all_coords  = [ ]
    all_normals = [ ]
    all_colors  = [ ]
    all_indices = [ ]
    ss_codes_per_chain = [ ]
    vertex_offset = 0

    for chain_index, record in enumerate ( chain_records ):
        reused = None
        if cached_ss_codes is not None and chain_index < len ( cached_ss_codes ):
            candidate = cached_ss_codes[chain_index]
            if candidate is not None and len ( candidate ) == len ( record["ca"] ):
                reused = candidate

        if reused is not None:
            ss_codes = reused
        else:
            ss_codes = classify_secondary_structure (
                record["ca"],
                min_run_helix          = settings["min_run_helix"],
                min_run_strand         = settings["min_run_strand"],
                helix_tolerance_scale  = settings["helix_tolerance_scale"],
                strand_tolerance_scale = settings["strand_tolerance_scale"],
            )
        ss_codes_per_chain.append ( ss_codes )

    residue_colors_per_chain = residue_colors_for_records (
        chain_records, ss_codes_per_chain, settings = settings,
        color_mode = color_mode, uniform_color = uniform_color )

    for record, ss_codes, residue_colors in zip ( chain_records, ss_codes_per_chain, residue_colors_per_chain ):
        if ca_filter is None:
            pieces = [ ( record, 0, len ( record["ca"] ) ) ]
        else:
            mask = [ bool ( ca_filter ( ca_atom ) ) for ca_atom in record["ca_atoms"] ]
            pieces = split_record_by_mask ( record, mask )

        for piece, start, end in pieces:
            coords, normals, colors, indices = build_ribbon_mesh (
                piece, np.asarray ( ss_codes )[start:end], spline_detail = spline_detail,
                settings = settings, residue_colors = residue_colors[start:end] )
            if len ( coords ) == 0:
                continue
            all_coords.append ( coords )
            all_normals.append ( normals )
            all_colors.append ( colors )
            all_indices.append ( indices + vertex_offset )
            vertex_offset += len ( coords )

    if not all_coords:
        empty3 = np.zeros ( ( 0, 3 ), dtype = np.float32 )
        result = ( empty3, empty3.copy ( ), empty3.copy ( ), np.zeros ( ( 0, ), dtype = np.uint32 ) )
    else:
        result = (
            np.concatenate ( all_coords,  axis = 0 ),
            np.concatenate ( all_normals, axis = 0 ),
            np.concatenate ( all_colors,  axis = 0 ),
            np.concatenate ( all_indices, axis = 0 ),
        )

    if return_ss_codes:
        return ( *result, ss_codes_per_chain )
    return result
