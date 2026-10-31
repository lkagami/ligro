# LiGRO - Version 1.1

This software is available to you under the terms of the GPL-3. See ~/ligro/LICENCE for more informations.
Software is created and maintained by Laboratorio de Sintese Organica Medicinal-LaSOM at
Universidade Federal do Rio Grande do Sul.

## Contributors:

MSc. Luciano Porto Kagami

[luciano_dot_kagami_at_ufrgs_dot_br](mailto:luciano.kagami@ufrgs.br)

Universidade Federal do Rio Grande do Sul
Laboratorio de Sintese Organica Medicinal -LaSOM

Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

## Prepare a protein-ligand input from a PDB ID

In LiGRO's **Input Files** tab, enter a four-character PDB ID and choose
**Download / Inspect**. Select the protein chain and ligand residue, then
choose **Prepare MD inputs**. LiGRO removes non-selected chains and solvent,
writes `protein.pdb` and converts the selected ligand to `Ligand.mol2`; these
files are then used by the existing MD workflow. This conversion requires
Open Babel (`obabel`) on `PATH`.

## Martini 3 coarse-grained protein preparation

The **System Parametrization** force-field selector also offers **Martini 3
(CG)** for protein-only preparation and energy minimization. This option uses
`martinize2` (provided by the `vermouth` Python package) to create the coarse-
grained protein and topology, then uses the BioExcel GROMACS blocks for box
setup and steepest-descent minimization. The minimized coarse-grained PDB is
written as `em.pdb`; its topology is restored as `trp.top`.

Install `martinize2` in the same environment as LiGRO with
`python -m pip install vermouth`, and install the Martini 3 GROMACS topology
files (including `martini_v3.0.0.itp`) where GROMACS can find them, for example
through its topology directory or `GMXLIB`. Martini does not parameterize
ligands, cofactors, or metals: if these are present, LiGRO warns that they will
be omitted and processes only the protein. This initial
integration stops after protein minimization; it does not add coarse-grained
solvent or run molecular dynamics. The existing atomistic force-field options
and workflows remain available.

## GROMACS–CP2K QM/MM with AWH

The **QM/MM (CP2K)** tab runs a one-dimensional AWH QM/MM simulation from an
already equilibrated atomistic system. Select its `.gro` coordinates, the
BioBB-compatible topology `.zip`, and an `.ndx` file containing the QM atom
group and the two groups that define the reaction-coordinate distance. Set
the QM region's formal charge and spin multiplicity, choose PBE or BLYP, and
provide the AWH distance interval in nm.

The selected GROMACS executable must be a compatible build with the
GROMACS–CP2K interface, such as the build described by the
[BioExcel GROMACS–CP2K tutorial](https://github.com/bioexcel/gromacs-2022-cp2k-tutorial);
an ordinary GROMACS installation is not sufficient. CP2K basis/potential data
must also be available to that build. BioBB supplies generic
GROMACS blocks used for preprocessing and running, but does not provide a
dedicated integrated QM/MM workflow. This tab does not prepare or equilibrate
the MM system, create the QM index group, or run umbrella windows. AWH requires
a constant-temperature ensemble. The AWH data are stored in the energy output
and can be inspected with `gmx awh`. See the [GROMACS AWH documentation](https://manual.gromacs.org/2025.4/reference-manual/special/awh.html)
for sampling details.

#### For Antechamber, please cite:

>1.  Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and bond type perception in molecular mechanical calculations". Journal of Molecular Graphics and Modelling , 25, 2006, 247260.
>2.  Wang, J., Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A. "Development and testing of a general AMBER force field". Journal of Computational Chemistry, 25, 2004, 1157-1174.

#### For ACPYPE, please cite:

>SOUSA DA SILVA, A. W. & VRANKEN, W. F. ACPYPE - AnteChamber PYthon Parser interfacE. BMC Research Notes 2012, 5:367 doi:10.1186/1756-0500-5-367 http://www.biomedcentral.com/1756-0500/5/367

#### For GROMACS please cite:

>Pronk, Sander, et al. "GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation toolkit." Bioinformatics 29.7 (2013): 845-854.
        
#### For PLIP please cite:

>Salentin, Sebastian, et al. "PLIP: fully automated protein–ligand interaction profiler."
>Nucleic acids research 43.W1 (2015): W443-W447.

#### For this code plese cite:

>Kagami, L.P., das Neves, G.M., da Silva, A.W.S. et al. J Mol Model (2017) 23: 304. https://doi.org/10.1007/s00894-017-3475-9

For install use Anaconda package Python 3.6:
```bash
conda install -c acpype ligro
```
    
