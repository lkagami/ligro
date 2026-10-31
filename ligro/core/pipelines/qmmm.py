"""GROMACS QM/MM simulation setup using the CP2K interface."""
import os
import re
import zipfile
from dataclasses import dataclass


@dataclass
class QmmmParams:
    structure: str
    topology_archive: str
    index_file: str
    qm_group: str = "QMatoms"
    qm_charge: int = 0
    qm_multiplicity: int = 1
    qm_method: str = "PBE"
    timestep_ps: float = 0.001
    steps: int = 1000
    temperature_k: float = 300.0
    awh_group1: str = ""
    awh_group2: str = ""
    awh_start_nm: float = 0.2
    awh_end_nm: float = 0.5
    awh_force_constant: float = 1000.0
    awh_diffusion: float = 1e-5

    def __post_init__(self):
        if not self.qm_group or any(char.isspace() for char in self.qm_group):
            raise ValueError("The QM index group name must be non-empty and contain no spaces.")
        if self.qm_method not in ("PBE", "BLYP"):
            raise ValueError("QM method must be PBE or BLYP.")
        if self.qm_multiplicity < 1:
            raise ValueError("QM multiplicity must be at least 1.")
        if self.steps < 1:
            raise ValueError("Simulation steps must be positive.")
        if self.timestep_ps <= 0 or self.temperature_k <= 0:
            raise ValueError("Time step and temperature must be positive.")
        if any(
            not group or any(char.isspace() for char in group)
            for group in (self.awh_group1, self.awh_group2)
        ):
            raise ValueError(
                "AWH reaction-coordinate group names must be non-empty and contain no spaces."
            )
        if self.awh_group1 == self.awh_group2:
            raise ValueError("AWH reaction-coordinate groups must be different.")
        if self.awh_start_nm < 0 or self.awh_start_nm >= self.awh_end_nm:
            raise ValueError(
                "AWH distance interval must be non-negative with start lower than end."
            )
        if self.awh_force_constant <= 0 or self.awh_diffusion <= 0:
            raise ValueError("AWH force constant and diffusion must be positive.")


def _index_groups(index_path):
    groups = {}
    current_group = None
    with open(index_path, "r", encoding="utf-8") as index_file:
        for line in index_file:
            match = re.match(r"^\s*\[\s*(.*?)\s*\]\s*$", line)
            if match:
                current_group = match.group(1).strip()
                if not current_group or current_group in groups:
                    raise ValueError("Index file contains an empty or duplicate group name.")
                groups[current_group] = []
                continue
            if current_group is None:
                if line.strip():
                    raise ValueError("Index file data appears before the first group header.")
                continue
            try:
                groups[current_group].extend(int(atom) for atom in line.split())
            except ValueError as error:
                raise ValueError(
                    "Index file contains a non-integer atom index in group {}."
                    .format(current_group)
                ) from error
    return groups


def _gro_atom_count(structure_path):
    with open(structure_path, "r", encoding="utf-8") as structure_file:
        structure_file.readline()
        try:
            return int(structure_file.readline().strip())
        except ValueError as error:
            raise ValueError(
                "The equilibrated structure is not a valid GRO file."
            ) from error


def render_qmmm_awh_mdp(params):
    """Render the GROMACS 2022/CP2K QM/MM AWH input."""
    return (
        "; QM/MM MD using GROMACS with the CP2K interface and AWH\n"
        "integrator = md\n"
        "dt = {timestep}\n"
        "nsteps = {steps}\n"
        "nstxout = 100\n"
        "nstvout = 100\n"
        "nstxout-compressed = 100\n"
        "compressed-x-grps = System\n"
        "nstlog = 100\n"
        "nstcalcenergy = 100\n"
        "nstenergy = 100\n"
        "cutoff-scheme = Verlet\n"
        "nstlist = 10\n"
        "rlist = 1.2\n"
        "coulombtype = PME\n"
        "rcoulomb = 1.2\n"
        "vdwtype = Cut-off\n"
        "rvdw = 1.2\n"
        "tcoupl = V-rescale\n"
        "tc-grps = System\n"
        "tau-t = 0.1\n"
        "ref-t = {temperature}\n"
        "continuation = yes\n"
        "gen-vel = no\n"
        "\n"
        "; CP2K QM/MM settings from the BioExcel GROMACS-CP2K tutorial\n"
        "qmmm-cp2k-active = yes\n"
        "qmmm-cp2k-qmgroup = {qm_group}\n"
        "qmmm-cp2k-qmmethod = {qm_method}\n"
        "qmmm-cp2k-qmcharge = {qm_charge}\n"
        "qmmm-cp2k-qmmultiplicity = {qm_multiplicity}\n"
        "\n"
        "; Reaction-coordinate groups used by AWH\n"
        "pull = yes\n"
        "pull-ngroups = 2\n"
        "pull-ncoords = 1\n"
        "pull-group1-name = {awh_group1}\n"
        "pull-group2-name = {awh_group2}\n"
        "pull-coord1-groups = 1 2\n"
        "pull-coord1-geometry = distance\n"
        "pull-coord1-type = external-potential\n"
        "pull-coord1-potential-provider = awh\n"
        "\n"
        "; AWH free-energy sampling\n"
        "awh = yes\n"
        "awh-nbias = 1\n"
        "awh-potential = convolved\n"
        "awh1-ndim = 1\n"
        "awh1-dim1-coord-provider = pull\n"
        "awh1-dim1-coord-index = 1\n"
        "awh1-dim1-start = {awh_start}\n"
        "awh1-dim1-end = {awh_end}\n"
        "awh1-dim1-diffusion = {diffusion}\n"
        "awh1-dim1-force-constant = {force_constant}\n"
    ).format(
        timestep=params.timestep_ps,
        steps=params.steps,
        temperature=params.temperature_k,
        qm_group=params.qm_group,
        qm_method=params.qm_method,
        qm_charge=params.qm_charge,
        qm_multiplicity=params.qm_multiplicity,
        awh_group1=params.awh_group1,
        awh_group2=params.awh_group2,
        awh_start=params.awh_start_nm,
        awh_end=params.awh_end_nm,
        diffusion=params.awh_diffusion,
        force_constant=params.awh_force_constant,
    )


class QmmmSimulation:
    """Generate an AWH QM/MM MDP and execute it through a GROMACS runner."""

    def __init__(self, params):
        self.params = params

    def run(self, workdir, runner):
        if not str(workdir).strip():
            raise ValueError("Select an output directory for the QM/MM run.")
        workdir = os.path.abspath(workdir)
        p = self.params
        for path, description in (
            (p.structure, "equilibrated structure"),
            (p.topology_archive, "topology archive"),
            (p.index_file, "index file"),
        ):
            if not os.path.isfile(path):
                raise FileNotFoundError("{} was not found: {}".format(description, path))

        groups = _index_groups(p.index_file)
        required_groups = {p.qm_group, p.awh_group1, p.awh_group2}
        missing_groups = sorted(required_groups - set(groups))
        if missing_groups:
            raise ValueError(
                "Index file is missing required group(s): {}".format(
                    ", ".join(missing_groups)
                )
            )
        atom_count = _gro_atom_count(p.structure)
        for name in required_groups:
            atom_indices = groups[name]
            if not atom_indices:
                raise ValueError("Index group {} is empty.".format(name))
            if min(atom_indices) < 1 or max(atom_indices) > atom_count:
                raise ValueError(
                    "Atom indices in group {} must be between 1 and {}."
                    .format(name, atom_count)
                )
        if not zipfile.is_zipfile(p.topology_archive):
            raise ValueError("The GROMACS topology archive is not a valid ZIP file.")
        try:
            with zipfile.ZipFile(p.topology_archive) as topology:
                if len([name for name in topology.namelist() if name.endswith(".top")]) != 1:
                    raise ValueError(
                        "The topology archive must contain exactly one .top file."
                    )
        except zipfile.BadZipFile as error:
            raise ValueError("The GROMACS topology archive is corrupt.") from error

        os.makedirs(workdir, exist_ok=True)
        mdp_path = os.path.join(workdir, "qmmm_awh.mdp")
        tpr_path = os.path.join(workdir, "qmmm_awh.tpr")
        deffnm = os.path.join(workdir, "qmmm_awh")
        output_paths = (
            tpr_path,
            *(
                deffnm + suffix
                for suffix in (".gro", ".edr", ".log", ".trr", ".xtc", ".cpt")
            ),
        )
        existing_outputs = [path for path in output_paths if os.path.exists(path)]
        if existing_outputs:
            raise FileExistsError(
                "QM/MM outputs already exist. Choose an empty output directory or "
                "move the previous results: {}".format(", ".join(existing_outputs))
            )
        with open(mdp_path, "w", encoding="utf-8") as mdp_file:
            mdp_file.write(render_qmmm_awh_mdp(p))

        runner.grompp(
            mdp_path,
            p.structure,
            p.topology_archive,
            tpr_path,
            ndx=p.index_file,
        )
        runner.mdrun(tpr_path, deffnm, write_xtc=True, write_cpt=True)

        outputs = {
            "mdp": mdp_path,
            "tpr": tpr_path,
            "gro": deffnm + ".gro",
            "edr": deffnm + ".edr",
            "log": deffnm + ".log",
            "trr": deffnm + ".trr",
            "xtc": deffnm + ".xtc",
            "cpt": deffnm + ".cpt",
        }
        missing = [path for path in outputs.values() if not os.path.isfile(path)]
        if missing:
            raise RuntimeError(
                "QM/MM simulation did not produce expected output(s): {}".format(
                    ", ".join(missing)
                )
            )
        return outputs
