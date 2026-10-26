"""Turns the biobb topology zip back into the plain files the rest of the
LiGRO workflow expects (trp.top plus the .itp includes such as posre.itp).

The old pipeline wrote trp.top directly. biobb writes trp.zip, whose
topology file is named p2g.top. This module restores the old layout.
"""
import os
import zipfile


def extract_topology(zip_path, workdir, top_name="trp.top"):
    """Extracts every member of the topology zip into workdir.

    The main .top member (exactly one is expected) is written as top_name.
    All other members (.itp files, posre.itp, ...) keep their own names.
    Only base names are used, so a zip member cannot write outside workdir.
    Returns the path of the written topology file.
    """
    if not os.path.exists(zip_path):
        raise RuntimeError("topology archive not found: " + zip_path)
    top_path = os.path.join(workdir, top_name)
    with zipfile.ZipFile(zip_path) as archive:
        tops = [m for m in archive.namelist() if m.endswith(".top")]
        if len(tops) != 1:
            raise RuntimeError("expected one .top file in {}, found {}".format(zip_path, tops))
        for member in archive.namelist():
            if member.endswith("/"):
                continue
            base = os.path.basename(member)
            if not base:
                continue
            target = top_path if member == tops[0] else os.path.join(workdir, base)
            with archive.open(member) as src, open(target, "wb") as dst:
                dst.write(src.read())
    return top_path



def pack_topology(workdir, top_name="trp.top", zip_name="trp_runtime.zip",
                  base_zip="trp.zip"):
    """Create a BioBB topology archive using the current working topology.

    Starts from the BioBB-generated archive (to preserve force-field files),
    replaces its main .top with the edited working copy, and adds local .itp
    files such as ligand/cofactor includes. Archive member paths are preserved.
    """
    import tempfile

    top_path = os.path.join(workdir, top_name)
    out_path = os.path.join(workdir, zip_name)
    base_path = os.path.join(workdir, base_zip)
    if not os.path.isfile(top_path):
        raise RuntimeError("topology file not found: " + top_path)

    # Preserve BioBB's force-field directory tree where available.
    members = {}
    if os.path.isfile(base_path):
        with zipfile.ZipFile(base_path, "r") as source:
            for name in source.namelist():
                if name.endswith("/"):
                    continue
                members[name] = source.read(name)

    # The edited topology must override the original topology member.
    top_members = [name for name in members if name.endswith(".top")]
    top_member = top_members[0] if len(top_members) == 1 else top_name
    with open(top_path, "rb") as fh:
        members[top_member] = fh.read()

    # Add current-directory include files; don't flatten any archived force
    # field members. The local top file's includes normally refer to basenames.
    for current, dirs, files in os.walk(workdir):
        # Do not archive VCS/cache directories or the output archive itself.
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__"}]
        for name in files:
            if not name.endswith(".itp"):
                continue
            full = os.path.join(current, name)
            rel = os.path.relpath(full, workdir).replace(os.sep, "/")
            with open(full, "rb") as fh:
                members[rel] = fh.read()

    # BioBB expects a zip archive; write atomically to avoid half-written input.
    fd, tmp_path = tempfile.mkstemp(prefix=".topology-", suffix=".zip", dir=workdir)
    os.close(fd)
    try:
        with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_DEFLATED) as target:
            for name, data in members.items():
                target.writestr(name, data)
        os.replace(tmp_path, out_path)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
    return out_path
