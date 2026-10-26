# -*- coding: utf-8 -*-
#
#  LiGRO: thin wrappers over BioExcel Building Blocks (biobb_gromacs).
#
#  Every GROMACS step that used to be launched with os.system() goes
#  through one of these functions. Each one:
#    - takes plain file paths and a dict of options,
#    - runs the biobb block (which calls gmx internally),
#    - raises RuntimeError on failure instead of silently continuing.
#
#  biobb packages are imported lazily, so this module can be imported
#  (and py_compile'd) on machines where biobb is not installed yet.
#
#  Install: pip install biobb-common biobb-gromacs
#
#  NOTE: the keyword names below follow the biobb_gromacs documentation
#  as known when this file was written. Confirm them against the
#  installed version with `help(biobb_gromacs.gromacs.grompp.grompp)`
#  before relying on them (see the TODO markers).

import os


def _run(block_fn, label, **kwargs):
    """Calls a biobb block and turns a non-zero return code into an
    exception that carries the step name."""
    rc = block_fn(**kwargs)
    if rc not in (None, 0):
        raise RuntimeError("biobb step '{}' failed (return code {})".format(label, rc))
    return rc


def _gmx_props(extra=None, binary="gmx"):
    props = {"binary_path": binary}
    if extra:
        props.update(extra)
    return props


def grompp(mdp, gro, top_zip, tpr, ndx=None, cpt=None, ref=None,
           maxwarn=1000, binary="gmx"):
    """gmx grompp: builds a .tpr run input.

    mdp      -- path to the .mdp file
    gro      -- input coordinates (.gro)
    top_zip  -- zipped topology directory (biobb convention)
    tpr      -- output .tpr path
    ndx, cpt, ref -- optional index / checkpoint / restraint reference
    """
    from biobb_gromacs.gromacs.grompp import grompp as _grompp  # confirmed signature

    io = {
        "input_gro_path": gro,
        "input_top_zip_path": top_zip,
        "input_mdp_path": mdp,          # confirmed: mdp path is a separate argument
        "output_tpr_path": tpr,
    }
    if ndx:
        io["input_ndx_path"] = ndx
    if cpt:
        io["input_cpt_path"] = cpt
    if ref:
        io["input_ref_path"] = ref
    props = _gmx_props({"maxwarn": maxwarn}, binary)
    _run(_grompp, "grompp", properties=props, **io)
    return tpr


def mdrun(tpr, deffnm, binary="gmx", cpt_in=None, write_xtc=False, write_cpt=False):
    """gmx mdrun for a .tpr. Writes <deffnm>.{trr,gro,edr,log} and, only
    when asked, <deffnm>.xtc and <deffnm>.cpt. Minimization produces neither
    trajectory, so passing them made biobb warn about missing files.
    Replaces `gmx mdrun -v -deffnm <deffnm>`.
    """
    from biobb_gromacs.gromacs.mdrun import mdrun as _mdrun  # confirmed signature

    io = {
        "input_tpr_path": tpr,
        "output_trr_path": deffnm + ".trr",
        "output_gro_path": deffnm + ".gro",
        "output_edr_path": deffnm + ".edr",
        "output_log_path": deffnm + ".log",
    }
    if write_xtc:
        io["output_xtc_path"] = deffnm + ".xtc"
    if write_cpt:
        io["output_cpt_path"] = deffnm + ".cpt"
    if cpt_in:
        io["input_cpt_path"] = cpt_in
    _run(_mdrun, "mdrun", properties=_gmx_props(None, binary), **io)
    return deffnm + ".gro"


def pdb2gmx(pdb, gro, top_zip, force_field="amber99sb-ildn",
            water="tip3p", binary="gmx", ignore_h=False):
    """gmx pdb2gmx: builds the topology for a protein PDB.
    Replaces `gmx pdb2gmx -ff ... -water ...` calls.
    ignore_h=True corresponds to the GUI's '-ignh' option."""
    from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx as _pdb2gmx  # confirmed signature

    io = {
        "input_pdb_path": pdb,
        "output_gro_path": gro,
        "output_top_zip_path": top_zip,
    }
    extra = {"force_field": force_field, "water_type": water}
    if ignore_h:
        extra["ignh"] = True          # confirmed property name
    props = _gmx_props(extra, binary)
    _run(_pdb2gmx, "pdb2gmx", properties=props, **io)
    return gro


def editconf(gro_in, gro_out, box_type="cubic", distance=1.0, binary="gmx"):
    """gmx editconf -bt <box_type> -d <distance>. Replaces `gmx editconf`."""
    from biobb_gromacs.gromacs.editconf import editconf as _editconf  # confirmed signature

    props = _gmx_props({"box_type": box_type,
                        "distance_to_molecule": distance,
                        "center_molecule": False}, binary)  # old command had no -c
    _run(_editconf, "editconf", properties=props,
         input_gro_path=gro_in, output_gro_path=gro_out)
    return gro_out


def solvate(gro_in, top_zip_in, gro_out, top_zip_out, binary="gmx"):
    """gmx solvate with the standard spc216 water box."""
    from biobb_gromacs.gromacs.solvate import solvate as _solvate  # confirmed signature

    _run(_solvate, "solvate", properties=_gmx_props(None, binary),
         input_solute_gro_path=gro_in,
         input_top_zip_path=top_zip_in,
         output_gro_path=gro_out,
         output_top_zip_path=top_zip_out)
    return gro_out


def genion(tpr, gro_out, top_zip, top_zip_out, neutral=True,
           binary="gmx", np=None, nn=None, concentration=None):
    """gmx genion. Confirmed biobb_gromacs signature:
    genion(input_tpr_path, output_gro_path, input_top_zip_path,
           output_top_zip_path, properties={...}).

    The system is read from the .tpr, so there is no separate input .gro.
    Exactly one of np (positive ions), nn (negative ions) or
    concentration (mol/L) is expected.
    replaced_group defaults to 'SOL' in biobb, matching the old pipeline.
    """
    from biobb_gromacs.gromacs.genion import genion as _genion  # confirmed signature

    chosen = [v for v in (np, nn, concentration) if v is not None]
    if len(chosen) != 1:
        raise ValueError("genion needs exactly one of np, nn, concentration")
    extra = {"neutral": neutral, "replaced_group": "SOL"}
    if np is not None:
        extra["np"] = np
    if nn is not None:
        extra["nn"] = nn
    if concentration is not None:
        extra["concentration"] = concentration
    _run(_genion, "genion", properties=_gmx_props(extra, binary),
         input_tpr_path=tpr,
         output_gro_path=gro_out,
         input_top_zip_path=top_zip,
         output_top_zip_path=top_zip_out)
    return gro_out


def make_ndx(structure, ndx_out, selection, ndx_in=None, binary="gmx"):
    """gmx make_ndx for ONE selection (confirmed biobb_gromacs signature:
    make_ndx(input_structure_path, output_ndx_path, input_ndx_path=None,
    properties={'selection': ...})).

    To build several groups, call this repeatedly, feeding each output
    back in as ndx_in:
        make_ndx(s, "a.ndx", "Protein")
        make_ndx(s, "b.ndx", "Protein | Other", ndx_in="a.ndx")
    Group names come from the selection text, as in GROMACS itself.
    """
    from biobb_gromacs.gromacs.make_ndx import make_ndx as _make_ndx  # confirmed signature

    io = {
        "input_structure_path": structure,
        "output_ndx_path": ndx_out,
    }
    if ndx_in:
        io["input_ndx_path"] = ndx_in
    _run(_make_ndx, "make_ndx", properties=_gmx_props({"selection": selection}, binary), **io)
    return ndx_out


def ensure_dir(path):
    """Helper used by callers that create run folders before biobb writes."""
    os.makedirs(path, exist_ok=True)
    return path
