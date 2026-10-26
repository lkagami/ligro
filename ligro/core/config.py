"""File names shared by the pipelines.

These match the names the original ligro.py wrote, so the rest of the
application (analysis, trajectory tools) keeps finding the same files.
"""

PROTEIN_INPUT = "protein.pdb"
PROTEIN_CLEAN = "protein_clean.pdb"
PROTEIN_GRO = "trp.pdb"
TOPOLOGY = "trp.zip"            # biobb convention: zipped topology directory
BOX_GRO = "trpb4solv.pdb"
SOLVATED_GRO = "trpb4ion.pdb"
IONIZED_TPR = "ion.tpr"
EM_INPUT_GRO = "trpb4em.pdb"
EM_DEFFNM = "em"
EM_FINAL_GRO = "em.gro"
EM_FINAL_PDB = "em.pdb"

SPC_WATER_BOX = "spc216.gro"
