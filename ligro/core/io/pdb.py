"""Pure-Python replacements for the grep/shell helpers used by ligro.py."""
import os


def clean_protein_pdb(src, dst, metal=None):
    """Keeps ATOM records from `src` and, if given, every line containing
    the metal code (e.g. 'ZN', 'FE'). Replaces:
        grep 'ATOM ' protein.pdb > protein_clean.pdb
        grep <metal> protein.pdb >> protein_clean.pdb
    Returns the number of lines written.
    """
    if metal is not None and metal.strip() in ("", "None"):
        metal = None
    written = 0
    with open(src, "r", encoding="utf-8", errors="replace") as fin, \
         open(dst, "w", encoding="utf-8") as fout:
        for line in fin:
            if "ATOM " in line or (metal is not None and metal in line):
                fout.write(line)
                written += 1
    return written


def ion_option(mode, value):
    """Maps the GUI's ionization choice to the gmx genion option.

    'Na (Number)' -> ('num_pos', n)   (old: -np n)
    'Cl (Number)' -> ('num_neg', n)   (old: -nn n)
    anything else -> ('conc', value)  (old: -conc value)
    """
    value = str(value).strip()
    if mode == "Na (Number)":
        return ("num_pos", int(value))
    if mode == "Cl (Number)":
        return ("num_neg", int(value))
    return ("conc", float(value))


def join(workdir, name):
    return os.path.join(workdir, name)
