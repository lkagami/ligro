#!/usr/bin/python
"""
---LiGRO - Version 0.2 ---

This software is available to you under the terms of the GPL-3. See ~/ligro/LICENCE for more informations.
Software is created and maintained by Laboratorio de Sintese Organica Medicinal-LaSOM at
Universidade Federal do Rio Grande do Sul.
Contributors:

Luciano Porto Kagami
luciano_dot_kagami_at_ufrgs_dot_br
Universidade Federal do Rio Grande do Sul
Laboratorio de Sintese Organica Medicinal -LaSOM
Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

               For Antechamber, please cite:

               1.  Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and
               bond type perception in molecular mechanical calculations". Journal of
               Molecular Graphics and Modelling , 25, 2006, 247260.
               2.  Wang, J., Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A.
               "Development and testing of a general AMBER force field". Journal of
               Computational Chemistry, 25, 2004, 1157-1174.

               For ACPYPE, please cite:

               SOUSA DA SILVA, A. W. & VRANKEN, W. F.
               ACPYPE - AnteChamber PYthon Parser interfacE.
               BMC Research Notes 2012, 5:367 doi:10.1186/1756-0500-5-367
               http://www.biomedcentral.com/1756-0500/5/367

               Alan Wilter Sousa da Silva, D.Sc.
               Bioinformatician, UniProt, EMBL-EBI
               Hinxton, Cambridge CB10 1SD, UK.
               >>http://www.ebi.ac.uk/~awilter<<

               alanwilter _at_ gmail _dot_ com

               For GROMACS please cite:

               For PLIP please cite:

               1. Salentin,S. et al. PLIP: fully automated protein-
               ligand interaction profiler. Nucl. Acids Res. (1 July 2015) 43 (W1):
               W443-W447. doi: 10.1093/nar/gkv315

               For GROMACS please cite:


               For this code plese cite:
               Manuscript to be submitted.

-------------------
"""
title = 'LiGRO: Version 0.2 - For GROMACS 5.x'

# Import Pmw from this directory tree.
import sys

sys.path[:0] = ['../../..']
import Tkinter, Tkconstants, tkFileDialog
import Tkinter
import Pmw
import os
import shutil
import webbrowser
import tkMessageBox as mbox
import tempfile
import matplotlib.pyplot as plt
import numpy
import subprocess
from datetime import date
import re
# add to installer
import pandas
import pylab



path2 = os.environ.get('HOME')

try:
    path = tempfile.mkdtemp()

except:
    pass


class GUI:
    def __init__(self, parent):
        #
        notebook = Pmw.NoteBook(parent)
        notebook.pack(fill='both', expand=2 ,padx=100, pady=10)

        #
        tab1 = notebook.add('System Preparation')
        notebook.tab('System Preparation').focus_set()

        #
        group_pdb = Pmw.Group(tab1, tag_text='Select PDB file:')
        group_pdb.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_pdb.interior(), text='                     ').grid(row=0, column=0)
        self.pdb1 = Tkinter.Label(group_pdb.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.pdb1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_pdb.interior(), text='').grid(row=0, column=2)
        pdb1_bt = Tkinter.Button(group_pdb.interior(), text='Browse...', command=self.openpdbfile)
        pdb1_bt.grid(row=0, column=3)

        group_conf = Pmw.Group(tab1, tag_text='Select Cofactor MOL2 file:')
        group_conf.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_conf.interior(), text='                     ').grid(row=0, column=0)
        self.conf1 = Tkinter.Label(group_conf.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.conf1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_conf.interior(), text='').grid(row=0, column=2)
        conf1_bt = Tkinter.Button(group_conf.interior(), text='Browse...', command=self.openconffile)
        conf1_bt.grid(row=0, column=3)

        self.b1 = Tkinter.Checkbutton(group_pdb.interior(), text='Ignore hydrogen atoms', variable='var1',
                                      onvalue="-ignh", offvalue="")
        self.b1.select()
        self.b1.grid(row=2, column=1, padx=10, pady=10)
        self.metal = Pmw.OptionMenu(group_pdb.interior(),
                                    labelpos='w',
                                    label_text='Metal:',
                                    menubutton_textvariable=None,
                                    items=['None', 'ZN', 'MG', 'CA'],
                                    menubutton_width=10,
                                    )
        self.metal.grid(row=2, column=2, padx=10, pady=10)
        #
        group_mol2 = Pmw.Group(tab1, tag_text='Select MOL2 file')
        group_mol2.pack(fill='both', expand=1, padx=100, pady=10)
        self.mol21 = Tkinter.Label(group_mol2.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.mol21.grid(row=0, column=0, padx=100, pady=10)
        self.b = Tkinter.Checkbutton(group_mol2.interior(), text='     ACPYPE user mode', variable='var5',
                                     onvalue="-c user", offvalue="")
        self.b.grid(row=1, column=0, padx=10, pady=10)
        Tkinter.Label(group_mol2.interior(), text='                        ').grid(row=0, column=1)
        mol21_bt = Tkinter.Button(group_mol2.interior(), text='Browse...', command=self.openmol2file)
        mol21_bt.grid(row=0, column=2)
        #
        self.ff_menu = Pmw.OptionMenu(tab1,
                                      labelpos='w',
                                      label_text='Select Force Field:',
                                      menubutton_textvariable=None,
                                      items=['amber03', 'amber94', 'amber99', 'amber99sb', 'amber99sb-ildn', 'amberGS',
                                             'oplsaa'],
                                      menubutton_width=10,
                                      )
        self.ff_menu.pack(anchor='w', padx=10, pady=10)
        #
        tab2 = notebook.add('Solvation')
        group_solv = Pmw.Group(tab2, tag_text='Solvate')
        group_solv.pack(fill='both', expand=1, padx=100, pady=10)
        self.wt_menu = Pmw.OptionMenu(group_solv.interior(),
                                      labelpos='w',
                                      label_text='Select Water Model:',
                                      menubutton_textvariable=None,
                                      items=['spc', 'spce', 'tip3p'],
                                      menubutton_width=10,
                                      )
        self.wt_menu.grid(row=1, column=0, padx=10, pady=10)
        self.bx_menu = Pmw.OptionMenu(group_solv.interior(),
                                      labelpos='w',
                                      label_text='Select Box Type:     ',
                                      menubutton_textvariable=None,
                                      items=['triclinic', 'cubic', 'dodecahedron', 'octahedron'],
                                      menubutton_width=10,
                                      )
        self.bx_menu.grid(row=2, column=0, padx=10, pady=10)
        self.dist = Pmw.EntryField(group_solv.interior(), labelpos='w', label_text='Distance:     ', value='1')
        self.dist.grid(row=3, column=0, padx=10, pady=10)
        #
        group_ion = Pmw.Group(tab2, tag_text='Neutralization / Ionize')
        group_ion.pack(fill='both', expand=1, padx=100, pady=10)

        self.ion_menu = Pmw.OptionMenu(group_ion.interior(),
                                       labelpos='w',
                                       menubutton_textvariable=None,
                                       items=['CONC', 'Na', 'Cl'],
                                       menubutton_width=10,
                                       )
        self.ion_menu.grid(row=6, column=0, padx=10, pady=10)

        self.ion_cc = Pmw.EntryField(group_ion.interior(), labelpos='w', value='0.15')
        self.ion_cc.grid(row=6, column=1, padx=10, pady=10)
        #
        tab3 = notebook.add('MDP parameters')
        group_others = Pmw.Group(tab3, tag_text='General options')
        group_others.grid(row=1, column=1, padx=10, pady=10)
        self.time_st = Pmw.EntryField(group_others.interior(), labelpos='w', value='0.002',
                                      label_text='Time step for integration (ps):')
        self.time_st.grid(row=1, column=1, padx=10, pady=10)
        self.temperature = Pmw.EntryField(group_others.interior(), labelpos='w', value='310.15',
                                          label_text='Temperature (K):                   ')
        self.temperature.grid(row=2, column=1, padx=10, pady=10)

        group_min = Pmw.Group(tab3, tag_text='Minimization')
        group_min.grid(row=1, column=2, padx=10, pady=10)
        self.step = Pmw.EntryField(group_min.interior(), labelpos='w', value='1000', label_text='Steps:      ')
        self.step.grid(row=2, column=0, padx=10, pady=10)
        self.min_menu = Pmw.OptionMenu(group_min.interior(),
                                       labelpos='w',
                                       items=['SD Algorithm', 'SD + CG'],
                                       menubutton_width=10,
                                       )
        self.min_menu.grid(row=4, column=0, padx=10, pady=10)
        group_NVT = Pmw.Group(tab3, tag_text='NVT')
        group_NVT.grid(row=3, column=1, padx=10, pady=10)
        self.time_nvt = Pmw.EntryField(group_NVT.interior(), labelpos='w', value='500000', label_text='Steps:                               ')
        self.time_nvt.grid(row=4, column=0, padx=10, pady=10)
        group_NPT = Pmw.Group(tab3, tag_text='NPT')
        group_NPT.grid(row=3, column=2, padx=10, pady=10)
        self.time_npt = Pmw.EntryField(group_NPT.interior(), labelpos='w', value='500000',
                                       label_text='Steps:           ')
        self.time_npt.grid(row=4, column=1, padx=10, pady=10)
        group_md = Pmw.Group(tab3, tag_text='Molecular Dynamics')
        group_md.grid(row=4, column=1, padx=10, pady=10)
        self.time_md = Pmw.EntryField(group_md.interior(), labelpos='w', value='500000', label_text='Steps:                                ')
        self.time_md.grid(row=6, column=0, padx=10, pady=10)

        #
        tab4 = notebook.add('Run MD')
        group_options = Pmw.Group(tab4, tag_text='Options')
        group_options.pack(fill='both', expand=1, padx=100, pady=10)
        hj = str(date.today())
        self.project = Pmw.EntryField(group_options.interior(), labelpos='w', value='MD_'+hj, label_text='Project Name')
        self.project.grid(row=1, column=3, padx=10, pady=100)
        save_bt = Tkinter.Button(group_options.interior(), text='Save TPR file', command=self.choose_simulation_save_tpr_file)
        save_bt.grid(row=3, column=0)
        run_bt = Tkinter.Button(group_options.interior(), text='Run Dynamics', command=self.choose_run_simulation)
        run_bt.grid(row=3, column=4)
        self.save = Pmw.EntryField(group_options.interior(), labelpos='w', value=path2, label_text='Save directory:')
        self.save.grid(row=2, column=3, padx=10, pady=5)
        #
        tab5 = notebook.add('LIE Binding Enegy')
        lie_group_min = Pmw.Group(tab5, tag_text='LIE Binding Energy Calculation')
        lie_group_min.grid(row=4, column=4, padx=10, pady=10)
        Tkinter.Label(lie_group_min.interior(), text='DG (kJ/mol)').grid(row=0, column=0)
        self.lie_ener= Tkinter.Label(lie_group_min.interior(), text='Empty', bg='white', padx=100, pady=5)
        self.lie_ener.grid(row=0, column=1, padx=10, pady=10)
        lie_run_bt = Tkinter.Button(tab5, text='Run LIE calculation', command=self.mount_simulation_LIE_complex)
        lie_run_bt.grid(row=7, column=4)


        #
        tab6 = notebook.add('Analysis')
        group_tpr = Pmw.Group(tab6, tag_text='Select TPR file:')
        group_tpr.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_tpr.interior(), text='                      ').grid(row=0, column=0)
        self.tpr1 = Tkinter.Label(group_tpr.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.tpr1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_tpr.interior(), text='                     ').grid(row=0, column=2)
        tpr1_bt = Tkinter.Button(group_tpr.interior(), text='Browse...', command=self.opentprfile)
        tpr1_bt.grid(row=0, column=3)
        #
        group_xtc = Pmw.Group(tab6, tag_text='Select XTC file')
        group_xtc.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_xtc.interior(), text='').grid(row=0, column=0)
        self.xtc = Tkinter.Label(group_xtc.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.xtc.grid(row=0, column=1, padx=100, pady=10)
        xtc1_bt = Tkinter.Button(group_xtc.interior(), text='Browse...', command=self.openxtcfile)
        xtc1_bt.grid(row=0, column=2)
        #
        group_edr = Pmw.Group(tab6, tag_text='Select EDR file:')
        group_edr.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_edr.interior(), text='                      ').grid(row=0, column=0)
        self.edr1 = Tkinter.Label(group_edr.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.edr1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_edr.interior(), text='                     ').grid(row=0, column=2)
        edr1_bt = Tkinter.Button(group_edr.interior(), text='Browse...', command=self.openedrfile)
        edr1_bt.grid(row=0, column=3)
        #
        group_analysis = Pmw.Group(tab6, tag_text='Analysis')
        group_analysis.pack(fill='both', expand=1, padx=100, pady=10)

        self.structure_menu = Pmw.OptionMenu(group_analysis.interior(),
                                      labelpos='w',
                                      label_text='Select Structure:',
                                      menubutton_textvariable=None,
                                      items=['Protein','Protein-H','C-alpha','Backbone','LIG'],
                                      menubutton_width=10,
                                      )
        self.structure_menu.grid(row=0, column=0)
        self.analysis_menu = Pmw.OptionMenu(group_analysis.interior(),
                                             labelpos='w',
                                             label_text='Select Analysis:  ',
                                             menubutton_textvariable=None,
                                             items=['RMSD', 'RMSF', 'RG','MSD', 'H_bond','LJSR-CoulSR IE'],
                                             menubutton_width=10,
                                             )
        self.analysis_menu.grid(row=0, column=2)
        run_analysis_bt = Tkinter.Button(group_analysis.interior(), text='Run', command=self.RMSD)
        run_analysis_bt.grid(row=0, column=4)
        Tkinter.Label(group_analysis.interior(),text='           ').grid(row=0, column=3)
        group_plip = Pmw.Group(tab6, tag_text='Protein-Ligand Interaction Profiler (PLIP) v1.3.3')
        group_plip.pack(fill='both', expand=1, padx=100, pady=10)
        self.ft = Pmw.EntryField(group_plip.interior(), labelpos='w', value='100', label_text='Frame time (ps):     ')
        self.ft.grid(row=0, column=0, padx=10, pady=10)
        run_analysis_bt = Tkinter.Button(group_plip.interior(), text='Run', command=self.plip)
        run_analysis_bt.grid(row=0, column=4)

    def openpdbfile(self):
        try:
            self.pdbfile = tkFileDialog.askopenfilename(initialdir=path2, title="Select PDB file",
                                                        filetypes=(("PDB files", "*.pdb"), ("all files", "*.*")))
            self.pdb1['text'] = self.pdbfile
            shutil.copy(self.pdbfile, path+'/protein.pdb')
        except:
            self.pdb1['text'] = 'No select file'

            #

    def openmol2file(self):
        try:
            self.mol2file = tkFileDialog.askopenfilename(initialdir=path2, title="Select MOL2 file",
                                                         filetypes=(("MOL2 files", "*.mol2"), ("all files", "*.*")))
            self.mol21['text'] = self.mol2file
            shutil.copy(self.mol2file, path+'/Ligand.mol2')
        except:
            self.mol21['text'] = 'No select file'

    def openconffile(self):
        try:
            self.conffile = tkFileDialog.askopenfilename(initialdir=path2, title="Select MOL2 file",
                                                         filetypes=(("MOL2 files", "*.mol2"), ("all files", "*.*")))
            self.conf1['text'] = self.conffile
            shutil.copy(self.conffile, path+'/Cofactor.mol2')
        except:
            self.conf1['text'] = 'No select file'        #

    def choose_simulation_save_tpr_file(self):
        try:
            os.chdir(path)
        except:
            pass
        find1=os.path.exists('Ligand.mol2')
        find2=os.path.exists('Cofactor.mol2')
        find3=os.path.exists('protein.pdb')
        if find1 == False:
            mbox.showinfo("INFO", "Run Protein simulation")
            self.mount_simulation_prot()
            self.save_tprfile()
        elif find2 == False:
            mbox.showinfo("INFO","Run Protein-Ligand simulation")
            self.mount_simulation_lig()
            self.save_tprfile()
            quit()
        elif find3 == False:
            mbox.showerror("Error", "PDB file not found.")
            quit()

        elif find2 == True and find1 == False:
            mbox.showerror("Error", "It is not possible run Protein-Cofactor")
            quit()

        else:
            mbox.showinfo("INFO", "Run Protein-Ligand simulation with Cofactor")
            self.mount_simulation_cof()
            self.save_tprfile()
            quit()

    def choose_run_simulation(self):
        try:
            os.chdir(path)
        except:
            pass
        find1=os.path.exists('Ligand.mol2')
        find2=os.path.exists('Cofactor.mol2')
        find3=os.path.exists('protein.pdb')
        if find1 == False:
            mbox.showinfo("INFO", "Run Protein simulation")
            self.mount_simulation_prot()
            self.run_simulation()
        elif find2 == False:
            mbox.showinfo("INFO","Run Protein-Ligand simulation")
            self.mount_simulation_lig()
            self.run_simulation()
            quit()
        elif find3 == False:
            mbox.showerror("Error", "PDB file not found.")
            quit()

        elif find2 == True and find1 == False:
            mbox.showerror("Error", "It is not possible run Protein-Cofactor")
            quit()

        else:
            mbox.showinfo("INFO", "Run Protein-Ligand simulation with Cofactor")
            self.mount_simulation_cof()
            self.run_simulation()
            quit()

    def mount_simulation_lig(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.metal.getcurselection())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.ff_menu.getcurselection())
        wt = str(self.wt_menu.getcurselection())
        ig = self.b1.getvar('var1')
        cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        os.system(cmd)
        os.system('babel Ligand.mol2 Ligand.sdf')
        os.system('babel Ligand.sdf Ligand.mol2')
        os.system("sed -i 's/ARG/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PHE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLN/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TYR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HIS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/CYS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/MET/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TRP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ASX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLU/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PCA/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HYP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/A1/LIG/g' Ligand.mol2")
        os.system("sed -i 's/UNK/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ACE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/FOR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HOH/LIG/g' Ligand.mol2")
        os.system("sed -i 's/DOD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/SO4/LIG/g' Ligand.mol2")
        os.system("sed -i 's/P04/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NAD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NDP/LIG/g' Ligand.mol2")
        usrm = self.b.getvar('var5')
        cmdd = 'acpype -i Ligand.mol2 {0}'.format(usrm)
        os.system(cmdd)
        os.system('grep -h ATOM trp.pdb Ligand.acpype/Ligand_NEW.pdb > complex.pdb')

        if ff != 'oplsaa':
            os.system('cp Ligand.acpype/Ligand_GMX.itp Ligand.itp')
        else:
            os.system('cp Ligand.acpype/Ligand_GMX_OPLS.itp Ligand.itp')

        os.system('cp trp.top Complex.top')
        os.system("cat Complex.top | sed '/forcefield\.itp\"/a\
                   #include \"Ligand.itp\" ' >| Complex2.top")

        os.system('echo "Ligand              1" >> Complex2.top')
        os.system('mv Complex2.top trp.top')
        bx = str(self.bx_menu.getcurselection())
        dst = str(self.dist.getvalue())

        cmd1 = 'gmx editconf -bt {0} -f complex.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
        os.system(cmd1)
        os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
        os.system('''
                    cat << EOF >| em.mdp
                    ; LINES STARTING WITH ';' ARE COMMENTS
                    title		= Minimization	; Title of run

                    ; Parameters describing what to do, when to stop and what to save
                    integrator	= steep		; Algorithm (steep = steepest descent minimization)
                    emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
                    emstep      = 0.01      ; Energy step size
                    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
                    energygrps	= system	; Which energy group(s) to write to disk

                    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
                    cutoff-scheme   = Verlet
                    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
                    rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
                    coulombtype	    = PME		; Treatment of long range electrostatic interactions
                    rcoulomb	    = 1.0		; long range electrostatic cut-off
                    rvdw		    = 1.0		; long range Van der Waals cut-off
                    pbc             = xyz 		; Periodic Boundary Conditions
                    ''')
        os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr')

        if self.ion_menu.getcurselection() == 'Na':
            c_ion = '-np ' + self.ion_cc.getvalue()
        elif self.ion_menu.getcurselection() == 'Cl':
            c_ion = '-nn ' + self.ion_cc.getvalue()
        else:
            c_ion = '-conc ' + self.ion_cc.getvalue()
        cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
        os.system(cmd2)

        inte = self.min_menu.getcurselection()

        nst = self.step.getvalue()

        cmd3 = '''
            cat << EOF >| em_real.mdp
            ; LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein LIG	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
            '''.format(nst)
        cmd4 = '''
            cat << EOF >| em_real.mdp
            ; LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= cg		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein LIG	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
            '''.format(nst)

        if inte == 'SD Algorithm':
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
            gmx make_ndx -f em.gro -o index2.ndx << EOF
            "Protein" | "Other"
            q
            EOF ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
            Protein_Other
            System
            EOF ''')
            os.system('''
            mv em_2.gro em.gro ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
            Protein_Other
            System
            EOF
            ''')
            os.system('''
                    mv em_2.pdb em.pdb ''')

        else:
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system(cmd4)
            os.system('gmx grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
                        gmx make_ndx -f em.gro -o index2.ndx << EOF
                        "Protein" | "Other"
                        q
                        EOF ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                        Protein_Other
                        System
                        EOF ''')
            os.system('''
                        mv em_2.gro em.gro ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                        Protein_Other
                        System
                        EOF
                        ''')
            os.system('mv em_2.pdb em.pdb')


        stnvt = str(self.time_nvt.getvalue())
        stnpt = str(self.time_npt.getvalue())
        stmd = str(self.time_md.getvalue())
        temp = str(self.temperature.getvalue())
        t_st = str(self.time_st.getvalue())

        cmd5 = '''cat << EOF >| nvt.mdp
        title       = Protein-ligand complex NVT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein LIG
        ; Bond parameters
        continuation    = no            ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = no        ; no pressure coupling in NVT
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = yes       ; assign velocities from Maxwell distribution
        gen_temp    = {2}       ; temperature for Maxwell distribution
        gen_seed    = -1        ; generate a random seed
        EOF
        '''.format(stnvt, t_st, temp)

        os.system(cmd5)

        cmd6 = '''
        cat << EOF >| npt.mdp
        title       = Protein-ligand complex NPT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein LIG
        ; Bond parameters
        continuation    = yes           ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                     ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
        pcoupltype  = isotropic                     ; uniform scaling of box vectors
        tau_p       = 2.0                           ; time constant, in ps
        ref_p       = 1.0                           ; reference pressure, in bar
        compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
        refcoord_scaling    = com
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = no        ; velocity generation off after NVT
        EOF
        '''.format(stnpt, t_st, temp)

        os.system(cmd6)

        cmd7 = '''
        cat << EOF >| md.mdp
        title       = Protein-ligand complex MD simulation
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout             = 0         ; suppress .trr output
        nstvout             = 0         ; suppress .trr output
        nstenergy           = 5000      ; save energies every 10.0 ps
        nstlog              = 5000      ; update log file every 10.0 ps
        nstxout-compressed  = 5000      ; write .xtc trajectory every 10.0 ps
        compressed-x-grps   = System
        energygrps          = Protein LIG
        ; Bond parameters
        continuation    = yes           ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                     ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
        pcoupltype  = isotropic                     ; uniform scaling of box vectors
        tau_p       = 2.0                           ; time constant, in ps
        ref_p       = 1.0                           ; reference pressure, in bar
        compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = no        ; assign velocities from Maxwell distribution
        EOF
        '''.format(stmd, t_st, temp)

        os.system(cmd7)

        cmd8='pymol '+path+'/em.pdb'

        os.system(cmd8)

        if mbox.askyesno('View Complex', 'Is complex OK??'):
            pass
        else:
            mbox.showinfo('No', 'Process has been cancelled')
            sys.exit(0)

        os.system('echo 2 | gmx genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
        os.system(r'''
        sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
        ''')
        os.system('mv trp2.top trp.top')

        os.system('cat << EOF > | queue.sh')

        os.system("echo 'gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr' >> queue.sh")
        os.system("echo 'gmx mdrun -v -deffnm nvt' >> queue.sh")
        os.system("echo 'gmx grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr' >> queue.sh")
        os.system("echo 'gmx mdrun -v -deffnm npt' >> queue.sh")
        os.system("echo 'gmx grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr' >> queue.sh")

    def mount_simulation_cof(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.metal.getcurselection())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.ff_menu.getcurselection())
        wt = str(self.wt_menu.getcurselection())
        ig = self.b1.getvar('var1')
        cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        os.system(cmd)
        os.system('babel Ligand.mol2 Ligand.sdf')
        os.system('babel Ligand.sdf Ligand.mol2')
        os.system("sed -i 's/ARG/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PHE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLN/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TYR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HIS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/CYS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/MET/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TRP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ASX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLU/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PCA/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HYP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/A1/LIG/g' Ligand.mol2")
        os.system("sed -i 's/UNK/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ACE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/FOR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HOH/LIG/g' Ligand.mol2")
        os.system("sed -i 's/DOD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/SO4/LIG/g' Ligand.mol2")
        os.system("sed -i 's/P04/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NAD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NDP/LIG/g' Ligand.mol2")
        os.system('babel Cofactor.mol2 Cofactor.sdf')
        os.system('babel Cofactor.sdf Cofactor.mol2')
        os.system("sed -i 's/ARG/COF/g' Cofactor.mol2")
        os.system("sed -i 's/PHE/COF/g' Cofactor.mol2")
        os.system("sed -i 's/GLN/COF/g' Cofactor.mol2")
        os.system("sed -i 's/TYR/COF/g' Cofactor.mol2")
        os.system("sed -i 's/HIS/COF/g' Cofactor.mol2")
        os.system("sed -i 's/CYS/COF/g' Cofactor.mol2")
        os.system("sed -i 's/MET/COF/g' Cofactor.mol2")
        os.system("sed -i 's/TRP/COF/g' Cofactor.mol2")
        os.system("sed -i 's/ASX/COF/g' Cofactor.mol2")
        os.system("sed -i 's/GLX/COF/g' Cofactor.mol2")
        os.system("sed -i 's/GLU/COF/g' Cofactor.mol2")
        os.system("sed -i 's/PCA/COF/g' Cofactor.mol2")
        os.system("sed -i 's/HYP/COF/g' Cofactor.mol2")
        os.system("sed -i 's/A1/COF/g' Cofactor.mol2")
        os.system("sed -i 's/UNK/COF/g' Cofactor.mol2")
        os.system("sed -i 's/ACE/COF/g' Cofactor.mol2")
        os.system("sed -i 's/FOR/COF/g' Cofactor.mol2")
        os.system("sed -i 's/HOH/COF/g' Cofactor.mol2")
        os.system("sed -i 's/DOD/COF/g' Cofactor.mol2")
        os.system("sed -i 's/SO4/COF/g' Cofactor.mol2")
        os.system("sed -i 's/P04/COF/g' Cofactor.mol2")
        os.system("sed -i 's/NAD/COF/g' Cofactor.mol2")
        os.system("sed -i 's/NDP/COF/g' Cofactor.mol2")
        usrm = self.b.getvar('var5')
        cmdd0 = 'acpype -i Ligand.mol2 {0}'.format(usrm)
        cmdd = 'acpype -i Cofactor.mol2 {0}'.format(usrm)
        os.system(cmdd0)
        os.system(cmdd)
        os.system("grep -h 'ATOM ' trp.pdb Ligand.acpype/Ligand_NEW.pdb > complex.pdb")
        os.system("grep -h 'ATOM ' complex.pdb Cofactor.acpype/Cofactor_NEW.pdb > complex2.pdb")

        if ff != 'oplsaa':
            os.system('cp Ligand.acpype/Ligand_GMX.itp Ligand.itp')
            os.system('cp Cofactor.acpype/Cofactor_GMX.itp Cofactor.itp')
        else:
            os.system('cp Ligand.acpype/Ligand_GMX_OPLS.itp Ligand.itp')
            os.system('cp Cofactor.acpype/Cofactor_GMX_OPLS.itp Cofactor.itp')



        with open('Ligand.itp') as infile, open('atomtype1.txt', 'w') as outfile:
            outfile.write('[ atomtypes ]\n')
            copy = False
            for line in infile:
                if line.strip() == "[ atomtypes ]":
                    copy = True
                elif line.strip() == "[ moleculetype ]":
                    copy = False
                elif copy:
                    outfile.write(line)

        with open('Cofactor.itp') as infile, open('atomtype2.txt', 'w') as outfile:

            copy = False
            for line in infile:
                if line.strip() == ";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb":
                    copy = True
                elif line.strip() == "[ moleculetype ]":
                    copy = False
                elif copy:
                     outfile.write(line)

        filenames = ['atomtype1.txt', 'atomtype2.txt']
        with open('atomtype.itp', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())

        with open('Ligand.itp') as infile, open('Ligand_clean.itp', 'w') as outfile:
            outfile.write('[ moleculetype ]\n')
            copy = False
            for line in infile:
                if line.strip() == "[ moleculetype ]":
                    copy = True
                elif copy:
                    outfile.write(line)

        with open('Cofactor.itp') as infile, open('Cofactor_clean.itp', 'w') as outfile:
            outfile.write('[ moleculetype ]\n')
            copy = False
            for line in infile:
                if line.strip() == "[ moleculetype ]":
                    copy = True
                elif copy:
                    outfile.write(line)
        os.system('cp Ligand_clean.itp Ligand.itp')
        os.system('cp Cofactor_clean.itp Cofactor.itp')
        os.system('cp trp.top Complex.top')
        os.system("cat Complex.top | sed '/forcefield\.itp\"/a\
                   #include \"atomtype.itp\" ' >| Complex2.top")
        os.system("cat Complex2.top | sed '/atomtype\.itp\"/a\
                   #include \"Ligand.itp\" ' >| Complex3.top")
        os.system("cat Complex3.top | sed '/Ligand\.itp\"/a\
                   #include \"Cofactor.itp\" ' >| Complex4.top")
        os.system('echo "Ligand              1" >> Complex4.top')
        os.system('echo "Cofactor            1" >> Complex4.top')
        os.system('mv Complex4.top trp.top')

        bx = str(self.bx_menu.getcurselection())
        dst = str(self.dist.getvalue())

        cmd1 = 'gmx editconf -bt {0} -f complex2.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
        os.system(cmd1)
        os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
        os.system('''
                    cat << EOF >| em.mdp
                    ; LINES STARTING WITH ';' ARE COMMENTS
                    title		= Minimization	; Title of run

                    ; Parameters describing what to do, when to stop and what to save
                    integrator	= steep		; Algorithm (steep = steepest descent minimization)
                    emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
                    emstep      = 0.01      ; Energy step size
                    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
                    energygrps	= system	; Which energy group(s) to write to disk

                    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
                    cutoff-scheme   = Verlet
                    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
                    rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
                    coulombtype	    = PME		; Treatment of long range electrostatic interactions
                    rcoulomb	    = 1.0		; long range electrostatic cut-off
                    rvdw		    = 1.0		; long range Van der Waals cut-off
                    pbc             = xyz 		; Periodic Boundary Conditions
                    ''')
        os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')

        if self.ion_menu.getcurselection() == 'Na':
            c_ion = '-np ' + self.ion_cc.getvalue()
        elif self.ion_menu.getcurselection() == 'Cl':
            c_ion = '-nn ' + self.ion_cc.getvalue()
        else:
            c_ion = '-conc ' + self.ion_cc.getvalue()
        cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
        os.system(cmd2)

        inte = self.min_menu.getcurselection()

        nst = self.step.getvalue()

        cmd3 = '''
            cat << EOF >| em_real.mdp
            ; LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein non-Protein	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
            '''.format(nst)
        cmd4 = '''cat << EOF >| em_real.mdp
            ; LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= cg		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein non-Protein	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF    = xyz 		; Periodic Boundary Conditions (yes/no)
            '''.format(nst)

        if inte == 'SD Algorithm':
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
            gmx make_ndx -f em.gro -o index2.ndx << EOF
            "Protein" | "Other"
            q
            EOF ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
            Protein_Other
            System
            EOF ''')
            os.system('''
            mv em_2.gro em.gro ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
            Protein_Other
            System
            EOF
            ''')
            os.system('''
                    mv em_2.pdb em.pdb ''')

        else:
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
            os.system('gmx mdrun -v -deffnm em')
            os.system(cmd4)
            os.system('gmx grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr -maxwarn 1000')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
                        gmx make_ndx -f em.gro -o index2.ndx << EOF
                        "Protein" | "Other"
                        q
                        EOF ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                        Protein_Other
                        System
                        EOF ''')
            os.system('''
                        mv em_2.gro em.gro ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                        Protein_Other
                        System
                        EOF
                        ''')
            os.system('mv em_2.pdb em.pdb')


        stnvt = str(self.time_nvt.getvalue())
        stnpt = str(self.time_npt.getvalue())
        stmd = str(self.time_md.getvalue())
        temp = str(self.temperature.getvalue())
        t_st = str(self.time_st.getvalue())

        cmd5 = '''cat << EOF >| nvt.mdp
        title       = Protein-ligand complex NVT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein LIG
        ; Bond parameters
        continuation    = no            ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = no        ; no pressure coupling in NVT
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = yes       ; assign velocities from Maxwell distribution
        gen_temp    = {2}       ; temperature for Maxwell distribution
        gen_seed    = -1        ; generate a random seed
        EOF
        '''.format(stnvt, t_st, temp)

        os.system(cmd5)

        cmd6 = '''cat << EOF >| npt.mdp
        title       = Protein-ligand complex NPT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein LIG
        ; Bond parameters
        continuation    = yes           ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                     ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
        pcoupltype  = isotropic                     ; uniform scaling of box vectors
        tau_p       = 2.0                           ; time constant, in ps
        ref_p       = 1.0                           ; reference pressure, in bar
        compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
        refcoord_scaling    = com
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = no        ; velocity generation off after NVT
        EOF
        '''.format(stnpt, t_st, temp)

        os.system(cmd6)

        cmd7 = '''cat << EOF >| md.mdp
        title       = Protein-ligand complex MD simulation
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout             = 0         ; suppress .trr output
        nstvout             = 0         ; suppress .trr output
        nstenergy           = 5000      ; save energies every 10.0 ps
        nstlog              = 5000      ; update log file every 10.0 ps
        nstxout-compressed  = 5000      ; write .xtc trajectory every 10.0 ps
        compressed-x-grps   = System
        energygrps          = Protein LIG
        ; Bond parameters
        continuation    = yes           ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                     ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
        pcoupltype  = isotropic                     ; uniform scaling of box vectors
        tau_p       = 2.0                           ; time constant, in ps
        ref_p       = 1.0                           ; reference pressure, in bar
        compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = no        ; assign velocities from Maxwell distribution
        EOF
        '''.format(stmd, t_st, temp)

        os.system(cmd7)

        cmd8='pymol '+path+'/em.pdb'

        os.system(cmd8)

        if mbox.askyesno('View Complex', 'Is complex OK??'):
            pass
        else:
            mbox.showinfo('No', 'Process has been cancelled')
            sys.exit(0)

        os.system('echo 2 | gmx genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
        os.system(r'''
        sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
        ''')
        os.system('mv trp2.top trp.top')

        os.system('cat << EOF > | queue.sh')

        os.system("echo 'gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -maxwarn 1000' >> queue.sh")
        os.system("echo 'gmx mdrun -v -deffnm nvt' >> queue.sh")
        os.system("echo 'gmx grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -maxwarn 1000' >> queue.sh")
        os.system("echo 'gmx mdrun -v -deffnm npt' >> queue.sh")
        os.system("echo 'gmx grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -maxwarn 1000' >> queue.sh")


    def mount_simulation_prot(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.metal.getcurselection())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.ff_menu.getcurselection())
        wt = str(self.wt_menu.getcurselection())
        ig = self.b1.getvar('var1')
        cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        os.system(cmd)
        bx = str(self.bx_menu.getcurselection())
        dst = str(self.dist.getvalue())
        cmd1 = 'gmx editconf -bt {0} -f trp.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
        os.system(cmd1)
        os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
        os.system('''
                    cat << EOF >| em.mdp
                    ; LINES STARTING WITH ';' ARE COMMENTS
                    title		= Minimization	; Title of run

                    ; Parameters describing what to do, when to stop and what to save
                    integrator	= steep		; Algorithm (steep = steepest descent minimization)
                    emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
                    emstep      = 0.01      ; Energy step size
                    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
                    energygrps	= system	; Which energy group(s) to write to disk

                    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
                    cutoff-scheme   = Verlet
                    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
                    rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
                    coulombtype	    = PME		; Treatment of long range electrostatic interactions
                    rcoulomb	    = 1.0		; long range electrostatic cut-off
                    rvdw		    = 1.0		; long range Van der Waals cut-off
                    pbc             = xyz 		; Periodic Boundary Conditions
                    ''')
        os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr')

        if self.ion_menu.getcurselection() == 'Na':
            c_ion = '-np ' + self.ion_cc.getvalue()
        elif self.ion_menu.getcurselection() == 'Cl':
            c_ion = '-nn ' + self.ion_cc.getvalue()
        else:
            c_ion = '-conc ' + self.ion_cc.getvalue()
        cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
        os.system(cmd2)

        inte = self.min_menu.getcurselection()

        nst = self.step.getvalue()

        cmd3 = '''cat << EOF >| em_real.mdp
            ; LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein 	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
            '''.format(nst)
        cmd4 = '''cat << EOF >| em_real.mdp
            ; LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= cg		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF		    = xyz 		; Periodic Boundary Conditions (yes/no)
                    '''.format(nst)

        if inte == 'SD Algorithm':
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro << EOF
            Protein
            System
            EOF ''')
            os.system('''
            mv em_2.gro em.gro ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb << EOF
            Protein
            System
            EOF
            ''')
            os.system('''
                    mv em_2.pdb em.pdb ''')

        else:
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system(cmd4)
            os.system('gmx grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro << EOF
                        Protein
                        System
                        EOF ''')
            os.system('''
                        mv em_2.gro em.gro ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb << EOF
                        Protein
                        System
                        EOF
                        ''')
            os.system('mv em_2.pdb em.pdb')


        stnvt = str(self.time_nvt.getvalue())
        stnpt = str(self.time_npt.getvalue())
        stmd = str(self.time_md.getvalue())
        temp = str(self.temperature.getvalue())
        t_st = str(self.time_st.getvalue())

        cmd5 = '''
        cat << EOF >| nvt.mdp
        title       = Protein NVT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein
        ; Bond parameters
        continuation    = no            ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein Water_and_ions    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = no        ; no pressure coupling in NVT
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = yes       ; assign velocities from Maxwell distribution
        gen_temp    = {2}       ; temperature for Maxwell distribution
        gen_seed    = -1        ; generate a random seed
        EOF
        '''.format(stnvt, t_st, temp)

        os.system(cmd5)

        cmd6 = '''
        cat << EOF >| npt.mdp
        title       = Protein NPT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein
        ; Bond parameters
        continuation    = yes           ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein Water_and_ions    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                     ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
        pcoupltype  = isotropic                     ; uniform scaling of box vectors
        tau_p       = 2.0                           ; time constant, in ps
        ref_p       = 1.0                           ; reference pressure, in bar
        compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
        refcoord_scaling    = com
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = no        ; velocity generation off after NVT
        EOF
        '''.format(stnpt, t_st, temp)

        os.system(cmd6)

        cmd7 = '''
        cat << EOF >| md.mdp
        title       = Protein MD simulation
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout             = 0         ; suppress .trr output
        nstvout             = 0         ; suppress .trr output
        nstenergy           = 5000      ; save energies every 10.0 ps
        nstlog              = 5000      ; update log file every 10.0 ps
        nstxout-compressed  = 5000      ; write .xtc trajectory every 10.0 ps
        compressed-x-grps   = System
        energygrps          = Protein
        ; Bond parameters
        continuation    = yes           ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein Water_and_ions    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                     ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
        pcoupltype  = isotropic                     ; uniform scaling of box vectors
        tau_p       = 2.0                           ; time constant, in ps
        ref_p       = 1.0                           ; reference pressure, in bar
        compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = no        ; assign velocities from Maxwell distribution
        EOF
        '''.format(stmd, t_st, temp)

        os.system(cmd7)

        cmd8='pymol '+path+'/em.pdb'

        os.system(cmd8)

        if mbox.askyesno('View Complex', 'Is protein OK??'):
            pass
        else:
            mbox.showinfo('No', 'Process has been cancelled')
            sys.exit(0)

        os.system('cat << EOF > | queue.sh')

        os.system("echo 'gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr' >> queue.sh")
        os.system("echo 'gmx mdrun -v -deffnm nvt' >> queue.sh")
        os.system("echo 'gmx grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr' >> queue.sh")
        os.system("echo 'gmx mdrun -v -deffnm npt' >> queue.sh")
        os.system("echo 'gmx grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr' >> queue.sh")


    def mount_simulation_LIE_complex(self):
        try:
            os.chdir(path)
        except:
            pass
        find1=os.path.exists('Ligand.mol2')
        find2=os.path.exists('Cofactor.mol2')
        find3=os.path.exists('protein.pdb')
        if find1 == False:
            mbox.showerror("INFO","Ligand not found")
            quit()
        elif find2 == False:
            mbox.showinfo("INFO","Run LIE calculation")
            self.simulation_LIE_complex()
            self.simulation_LIE_lig()
            self.lie_calculation()
            mbox.showinfo("INFO","LIE calculation is finished")
            return self.c
        elif find3 == False:
            mbox.showerror("Error", "PDB file not found.")
            quit()

        else:
            mbox.showerror("Error", "LIE calculation can not be performed with Cofactor.")
            quit()

    def simulation_LIE_complex(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.metal.getcurselection())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.ff_menu.getcurselection())
        wt = str(self.wt_menu.getcurselection())
        ig = self.b1.getvar('var1')
        cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        os.system(cmd)
        os.system('babel Ligand.mol2 Ligand.sdf')
        os.system('babel Ligand.sdf Ligand.mol2')
        os.system("sed -i 's/ARG/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PHE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLN/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TYR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HIS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/CYS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/MET/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TRP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ASX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLU/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PCA/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HYP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/A1/LIG/g' Ligand.mol2")
        os.system("sed -i 's/UNK/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ACE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/FOR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HOH/LIG/g' Ligand.mol2")
        os.system("sed -i 's/DOD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/SO4/LIG/g' Ligand.mol2")
        os.system("sed -i 's/P04/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NAD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NDP/LIG/g' Ligand.mol2")
        usrm = self.b.getvar('var5')
        cmdd = 'acpype -i Ligand.mol2 {0}'.format(usrm)
        os.system(cmdd)
        os.system('grep -h ATOM trp.pdb Ligand.acpype/Ligand_NEW.pdb > complex.pdb')

        if ff != 'oplsaa':
            os.system('cp Ligand.acpype/Ligand_GMX.itp Ligand.itp')
        else:
            os.system('cp Ligand.acpype/Ligand_GMX_OPLS.itp Ligand.itp')

        os.system('cp trp.top Complex.top')
        os.system("cat Complex.top | sed '/forcefield\.itp\"/a\
                   #include \"Ligand.itp\" ' >| Complex2.top")

        os.system('echo "Ligand              1" >> Complex2.top')
        os.system('mv Complex2.top trp.top')
        bx = str(self.bx_menu.getcurselection())
        dst = str(self.dist.getvalue())

        cmd1 = 'gmx editconf -bt {0} -f complex.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
        os.system(cmd1)
        os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
        os.system('''
                    cat << EOF >| em.mdp
                    ; LINES STARTING WITH ';' ARE COMMENTS
                    title		= Minimization	; Title of run

                    ; Parameters describing what to do, when to stop and what to save
                    integrator	= steep		; Algorithm (steep = steepest descent minimization)
                    emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
                    emstep      = 0.01      ; Energy step size
                    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
                    energygrps	= system	; Which energy group(s) to write to disk

                    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
                    cutoff-scheme   = Verlet
                    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
                    rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
                    coulombtype	    = PME		; Treatment of long range electrostatic interactions
                    rcoulomb	    = 1.0		; long range electrostatic cut-off
                    rvdw		    = 1.0		; long range Van der Waals cut-off
                    pbc             = xyz 		; Periodic Boundary Conditions
                    ''')
        os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr')

        if self.ion_menu.getcurselection() == 'Na':
            c_ion = '-np ' + self.ion_cc.getvalue()
        elif self.ion_menu.getcurselection() == 'Cl':
            c_ion = '-nn ' + self.ion_cc.getvalue()
        else:
            c_ion = '-conc ' + self.ion_cc.getvalue()
        cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
        os.system(cmd2)

        inte = self.min_menu.getcurselection()

        nst = self.step.getvalue()

        cmd3 = '''
            cat << EOF >| em_real.mdp
            ;LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein 	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
            '''.format(nst)
        cmd4 = '''cat << EOF >| em_real.mdp
             LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= cg		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= Protein 	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
                    '''.format(nst)

        if inte == 'SD Algorithm':
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
            gmx make_ndx -f em.gro -o index2.ndx << EOF
            "Protein" | "Other"
            q
            EOF ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
            Protein_Other
            System
            EOF ''')
            os.system('''
            mv em_2.gro em.gro ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
            Protein_Other
            System
            EOF
            ''')
            os.system('''
                    mv em_2.pdb em.pdb ''')

        else:
            os.system(cmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system(cmd4)
            os.system('gmx grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
                        gmx make_ndx -f em.gro -o index2.ndx << EOF
                        "Protein" | "Other"
                        q
                        EOF ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                        Protein_Other
                        System
                        EOF ''')
            os.system('''
                        mv em_2.gro em.gro ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                        Protein_Other
                        System
                        EOF
                        ''')
            os.system('mv em_2.pdb em.pdb')


        stnvt = str(self.time_nvt.getvalue())
        temp = str(self.temperature.getvalue())
        t_st = str(self.time_st.getvalue())

        cmd5 = '''
        cat << EOF >| nvt.mdp
        title       = Protein-ligand complex NVT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = Protein LIG
        ; Bond parameters
        continuation    = no            ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        vdwtype             =  Cut-off
        ; Electrostatics
        coulombtype     = Cut-off       ; Particle Mesh Ewald for long-range electrostatics
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = no        ; no pressure coupling in NVT
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = yes       ; assign velocities from Maxwell distribution
        gen_temp    = {2}       ; temperature for Maxwell distribution
        gen_seed    = -1        ; generate a random seed
        EOF
        '''.format(stnvt, t_st, temp)

        os.system(cmd5)

        cmd8='pymol '+path+'/em.pdb'

        os.system(cmd8)

        if mbox.askyesno('View Complex', 'Is complex OK??'):
            pass
        else:
            mbox.showinfo('No', 'Process has been cancelled')
            sys.exit(0)

        os.system('echo 2|gmx genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
        os.system(r'''
        sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
        ''')
        os.system('mv trp2.top trp.top')

        os.system('cat << EOF > | queue.sh')

        os.system('gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr')
        os.system('gmx mdrun -v -deffnm nvt')


    def simulation_LIE_lig(self):
        try:
            os.chdir(path)
        except:
            pass
        os.system('babel Ligand.mol2 Ligand.sdf')
        os.system('babel Ligand.sdf Ligand.mol2')
        os.system("sed -i 's/ARG/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PHE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLN/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TYR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HIS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/CYS/LIG/g' Ligand.mol2")
        os.system("sed -i 's/MET/LIG/g' Ligand.mol2")
        os.system("sed -i 's/TRP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ASX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLX/LIG/g' Ligand.mol2")
        os.system("sed -i 's/GLU/LIG/g' Ligand.mol2")
        os.system("sed -i 's/PCA/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HYP/LIG/g' Ligand.mol2")
        os.system("sed -i 's/A1/LIG/g' Ligand.mol2")
        os.system("sed -i 's/UNK/LIG/g' Ligand.mol2")
        os.system("sed -i 's/ACE/LIG/g' Ligand.mol2")
        os.system("sed -i 's/FOR/LIG/g' Ligand.mol2")
        os.system("sed -i 's/HOH/LIG/g' Ligand.mol2")
        os.system("sed -i 's/DOD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/SO4/LIG/g' Ligand.mol2")
        os.system("sed -i 's/P04/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NAD/LIG/g' Ligand.mol2")
        os.system("sed -i 's/NDP/LIG/g' Ligand.mol2")
        usrm = self.b.getvar('var5')
        cmdd = 'acpype -i Ligand.mol2 {0}'.format(usrm)
        os.system(cmdd)
        os.system('cp Ligand.acpype/Ligand_GMX.gro Ligand.gro')

        ff = str(self.ff_menu.getcurselection())
        wt = str(self.wt_menu.getcurselection())

        if ff != 'oplsaa':
            os.system('cp Ligand.acpype/Ligand_GMX.itp Ligand.itp')
        else:
            os.system('cp Ligand.acpype/Ligand_GMX_OPLS.itp Ligand.itp')

        liecmd5 = '''
        cat << EOF >| topol.top
#include "{0}.ff/forcefield.itp"
#include "Ligand.itp"
#include "{0}.ff/{1}.itp"
#include "{0}.ff/ions.itp"
[ system ]
; Name
Ligand

[ molecules ]
; Compound        #mols
Ligand   1

        '''.format(ff,wt)

        os.system(liecmd5)

        bx = str(self.bx_menu.getcurselection())
        dst = str(self.dist.getvalue())

        liecmd1 = 'gmx editconf -bt {0} -f Ligand.gro -o trpb4solv.pdb -d {1}'.format(bx, dst)
        os.system(liecmd1)
        os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p topol.top')
        os.system('''
        cat << EOF >| em.mdp
; LINES STARTING WITH ';' ARE COMMENTS
                    title		= Minimization	; Title of run

                    ; Parameters describing what to do, when to stop and what to save
                    integrator	= steep		; Algorithm (steep = steepest descent minimization)
                    emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
                    emstep      = 0.01      ; Energy step size
                    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform
                    energygrps	= system	; Which energy group(s) to write to disk

                    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
                    cutoff-scheme   = Verlet
                    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
                    rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
                    coulombtype	    = PME		; Treatment of long range electrostatic interactions
                    rcoulomb	    = 1.0		; long range electrostatic cut-off
                    rvdw		    = 1.0		; long range Van der Waals cut-off
                    pbc             = xyz 		; Periodic Boundary Conditions
                    ''')
        os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p topol.top -o ion.tpr')

        if self.ion_menu.getcurselection() == 'Na':
            c_ion = '-np ' + self.ion_cc.getvalue()
        elif self.ion_menu.getcurselection() == 'Cl':
            c_ion = '-nn ' + self.ion_cc.getvalue()
        else:
            c_ion = '-conc ' + self.ion_cc.getvalue()
        liecmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p topol.top'.format(c_ion)
        os.system(liecmd2)

        lie_inte = self.min_menu.getcurselection()

        lie_nst = self.step.getvalue()

        liecmd3 = '''
            cat << EOF >| em_real.mdp
            ;LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= LIG 	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
            '''.format(lie_nst)
        liecmd4 = '''
            cat << EOF >| em_real.mdp
            ;LINES STARTING WITH ';' ARE COMMENTS
            title		= Minimization	; Title of run

            ; Parameters describing what to do, when to stop and what to save
            integrator	= cg		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
            emstep      = 0.01      ; Energy step size
            nsteps		= {0}	  	; Maximum number of (minimization) steps to perform
            energygrps	= LIG 	; Which energy group(s) to write to disk

            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; long range electrostatic cut-off
            rvdw		    = 1.0		; long range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions
            EOF
                    '''.format(lie_nst)

        if lie_inte == 'SD Algorithm':
            os.system(liecmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p topol.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro << EOF
            Other
            System
            EOF ''')
            os.system('''
            mv em_2.gro em.gro ''')
            os.system('''
            gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb << EOF
            Other
            System
            EOF
            ''')
            os.system('''
                    mv em_2.pdb em.pdb ''')

        else:
            os.system(liecmd3)
            os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p topol.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system(liecmd4)
            os.system('gmx grompp -f em_real.mdp -c em.gro -p topol.top -o em.tpr')
            os.system('gmx mdrun -v -deffnm em')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro << EOF
                        Other
                        System
                        EOF ''')
            os.system('''
                        mv em_2.gro em.gro ''')
            os.system('''
                        gmx trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb << EOF
                        Other
                        System
                        EOF
                        ''')
            os.system('mv em_2.pdb em.pdb')


        lie_stnvt = str(self.time_nvt.getvalue())
        lie_t_st = str(self.time_st.getvalue())
        lie_temp = str(self.temperature.getvalue())

        liecmd5 = '''
        cat << EOF >| nvt2.mdp
        title       = Ligand NVT equilibration
        define      = -DPOSRES  ; position restrain the protein and ligand
        ; Run parameters
        integrator  = md        ; leap-frog integrator
        nsteps      = {0}     ; 2 * 50000 = 100 ps
        dt          = {1}     ; 2 fs
        ; Output control
        nstxout     = 500       ; save coordinates every 1.0 ps
        nstvout     = 500       ; save velocities every 1.0 ps
        nstenergy   = 500       ; save energies every 1.0 ps
        nstlog      = 500       ; update log file every 1.0 ps
        energygrps  = LIG
        ; Bond parameters
        continuation    = no            ; first dynamics run
        constraint_algorithm = lincs    ; holonomic constraints
        constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
        lincs_iter      = 1             ; accuracy of LINCS
        lincs_order     = 4             ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme   = Verlet
        ns_type         = grid      ; search neighboring grid cells
        nstlist         = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb        = 1.4       ; short-range electrostatic cutoff (in nm)
        rvdw            = 1.4       ; short-range van der Waals cutoff (in nm)
        vdwtype             =  Cut-off
        ; Electrostatics
        coulombtype     = Cut-off       ; Particle Mesh Ewald for long-range electrostatics
        pme_order       = 4         ; cubic interpolation
        fourierspacing  = 0.16      ; grid spacing for FFT
        ; Temperature coupling
        tcoupl      = V-rescale                     ; modified Berendsen thermostat
        tc-grps     = LIG Water_Ion    ; two coupling groups - more accurate
        tau_t       = 0.1   0.1                     ; time constant, in ps
        ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
        ; Pressure coupling
        pcoupl      = no        ; no pressure coupling in NVT
        ; Periodic boundary conditions
        pbc         = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr    = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel     = yes       ; assign velocities from Maxwell distribution
        gen_temp    = {2}       ; temperature for Maxwell distribution
        gen_seed    = -1        ; generate a random seed
        EOF
        '''.format(lie_stnvt, lie_t_st, lie_temp)

        os.system(liecmd5)

        liecmd8='pymol '+path+'/em.pdb'

        os.system(liecmd8)

        if mbox.askyesno('View Complex', 'Is Ligand OK??'):
            pass
        else:
            mbox.showinfo('No', 'Process has been cancelled')
            sys.exit(0)
        os.system('''
                        gmx make_ndx -f em.gro -o index2.ndx << EOF
                        "Water" | "Ion"
                        q
                        EOF ''')
        os.system('gmx grompp -f nvt2.mdp -c em.pdb -p topol.top -o nvt2.tpr -n index2.ndx')
        os.system('gmx mdrun -v -deffnm nvt2')

    def lie_calculation(self):
        try:
                os.chdir(path)
        except:
            pass
        com1= '''echo 43 44 0 | gmx energy -f nvt.edr > saida1.txt'''
        os.system(com1)
        com2= '''echo 41 42 0 | gmx energy -f nvt2.edr > saida2.txt'''
        os.system(com2)

        df = pandas.read_table('saida1.txt', delim_whitespace=True, names=('A', 'B', 'C'))

        df.iloc[[6],[1]].to_csv(r'energy1.csv', header=None, index=None, sep=',', mode='a')
        df.iloc[[7],[1]].to_csv(r'energy1.csv', header=None, index=None, sep=',', mode='a')

        df1 = pandas.read_csv('energy1.csv', names=('A'))

        df1.astype(str)

        coul = df1.iloc[[0],[0]].values
        lj = df1.iloc[[1],[0]].values
        coull = re.sub('[\[\]]', '', numpy.array_str(coul))
        ljj = re.sub('[\[\]]', '', numpy.array_str(lj))


        df2 = pandas.read_table('saida2.txt', delim_whitespace=True, names=('A', 'B', 'C'))

        df2.iloc[[6],[1]].to_csv(r'energy2.csv', header=None, index=None, sep=',', mode='a')
        df2.iloc[[7],[1]].to_csv(r'energy2.csv', header=None, index=None, sep=',', mode='a')

        df2 = pandas.read_csv('energy2.csv', names=('A'))

        df2.astype(float)

        coul2 = df2.iloc[[0],[0]].values
        lj2 = df2.iloc[[1],[0]].values
        coull2 = re.sub('[\[\]]', '', numpy.array_str(coul2))
        ljj2 = re.sub('[\[\]]', '', numpy.array_str(lj2))

        os.system('gmx lie -f nvt.edr -o lie_comp.xvg -ligand LIG -Eqq {0} -Elj {1} > out1.txt'.format(coull,ljj))
        os.system('gmx lie -f nvt2.edr -o lie_lig.xvg -ligand LIG -Eqq {0} -Elj {1} >> out1.txt'.format(coull2,ljj2))

        df = pandas.read_table('out1.txt', delim_whitespace=True, names=('A', 'B','C','D','E'))

        df.iloc[[3],[2]].to_csv(r'energy3.csv', header=None, index=None, sep=',', mode='a')
        df.iloc[[7],[2]].to_csv(r'energy3.csv', header=None, index=None, sep=',', mode='a')

        df1 = pandas.read_csv('energy3.csv', names=('A'))

        df1.astype(float)

        a = df1.iloc[[0],[0]].values
        b = df1.iloc[[1],[0]].values
        self.c =float(a+b)
        self.lie_ener['text'] = str(self.c)



    def save_tprfile(self):
        os.system('chmod 777 queue.sh')
        os.system('./queue.sh')
        dr = str(self.save.getvalue())
        pj = str(self.project.getvalue())
        pj1 = pj
        md0 = '{0}.tpr'.format(pj1)
        shutil.copy('md.tpr', md0)
        shutil.copy2(md0, dr)
        mbox.showinfo('Finish', 'Job has finished')
        sys.exit(0)

    def run_simulation(self):
        os.system("echo 'gmx mdrun -v -deffnm md' >> queue.sh")
        os.system('chmod 777 queue.sh')
        os.system('./queue.sh')
        pj = str(self.project.getvalue())
        cmd8 = 'mkdir {0}'.format(pj)
        os.system(cmd8)
        md1 = '{0}.xtc'.format(pj)
        md2 = '{0}.edr'.format(pj)
        md3 = '{0}.gro'.format(pj)
        md4 = '{0}.tpr'.format(pj)
        try:
            os.chdir(path)
        except:
            pass
        find1=os.path.exists('Ligand.mol2')
        if find1 == False:

            os.system('''
                                             gmx trjconv -f md.gro -s md.tpr -pbc nojump -ur compact -center -o md_2.gro << EOF
                                             Protein
                                             System
                                             EOF ''')
            os.system('''
                                             mv md_2.gro md.gro ''')
            os.system('''
                                             gmx trjconv -f md.gro -s md.tpr -pbc mol -ur compact -center -o md_2.gro << EOF
                                             Protein
                                             System
                                             EOF
                                             ''')
            os.system('mv md_2.gro md.gro')
            os.system('''
                                                     gmx trjconv -f md.xtc -s md.tpr -pbc nojump -ur compact -center -o md_2.xtc << EOF
                                                     Protein
                                                     System
                                                     EOF ''')
            os.system('''
                                                     mv md_2.xtc md.xtc ''')
            os.system('''
                                                     gmx trjconv -f md.xtc -s md.tpr -pbc mol -ur compact -center -o md_2.xtc << EOF
                                                     Protein
                                                     System
                                                     EOF
                                                     ''')
            os.system('mv md_2.xtc md.xtc')
        else:
            os.system('''
                                             gmx make_ndx -f md.gro -o index3.ndx << EOF
                                             "Protein" | "Other"
                                             q
                                             EOF ''')
            os.system('''
                                             gmx trjconv -f md.gro -s md.tpr -pbc nojump -ur compact -center -o md_2.gro -n index3.ndx << EOF
                                             Protein_Other
                                             System
                                             EOF ''')
            os.system('''
                                             mv md_2.gro md.gro ''')
            os.system('''
                                             gmx trjconv -f md.gro -s md.tpr -pbc mol -ur compact -center -o md_2.gro -n index3.ndx << EOF
                                             Protein_Other
                                             System
                                             EOF
                                             ''')
            os.system('mv md_2.gro md.gro')
            os.system('''
                                                     gmx trjconv -f md.xtc -s md.tpr -pbc nojump -ur compact -center -o md_2.xtc -n index3.ndx << EOF
                                                     Protein_Other
                                                     System
                                                     EOF ''')
            os.system('''
                                                     mv md_2.xtc md.xtc ''')
            os.system('''
                                                     gmx trjconv -f md.xtc -s md.tpr -pbc mol -ur compact -center -o md_2.xtc -n index3.ndx << EOF
                                                     Protein_Other
                                                     System
                                                     EOF
                                                     ''')
            os.system('mv md_2.xtc md.xtc')
        shutil.copy('md.xtc', md1)
        shutil.copy('md.edr', md2)
        shutil.copy('md.gro', md3)
        shutil.copy('md.tpr', md4)
        dr = str(self.save.getvalue())
        shutil.copy2(md1, dr)
        shutil.copy2(md2, dr)
        shutil.copy2(md3, dr)
        shutil.copy2(md4, dr)
        mbox.showinfo('Finish', 'Job has finished')
        sys.exit(0)

    def opentprfile(self):
        try:
            self.tprfile = tkFileDialog.askopenfilename(initialdir=path2, title="Select TPR file",
                                                        filetypes=(("TPR files", "*.tpr"), ("all files", "*.*")))
            self.tpr1['text'] = self.tprfile
            shutil.copy(self.tprfile, path+'/md_an.tpr')
        except:
            self.tpr1['text'] = 'No select file'

    def openedrfile(self):
        try:
            self.edrfile = tkFileDialog.askopenfilename(initialdir=path2, title="Select EDR file",
                                                        filetypes=(("EDR files", "*.edr"), ("all files", "*.*")))
            self.edr1['text'] = self.edrfile
            shutil.copy(self.edrfile, path+'/md_an.edr')
        except:
            self.edr1['text'] = 'No select file'
            #

    def openxtcfile(self):
        try:
            self.xtcfile = tkFileDialog.askopenfilename(initialdir=path2, title="Select XTC file",
                                                         filetypes=(("XTC files", "*.xtc"), ("all files", "*.*")))
            self.xtc['text'] = self.xtcfile
            shutil.copy(self.xtcfile, path+'/md_an.xtc')
        except:
            self.xtc['text'] = 'No select file'

            #
    def RMSD(self):
        struct = str(self.structure_menu.getcurselection())
        analysis=str(self.analysis_menu.getcurselection())


        if analysis == 'RMSD':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                quit()
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                quit()
            else:
                pass
            command= ''' gmx rms -s md_an.tpr -f md_an.xtc -o rmsd.xvg -xvg none << EOF
                                    {0}
                                    {0}
                                    EOF '''.format(struct)
            os.system(command)
            pj = str(self.project.getvalue())
            an1 = 'rmsd_{0}.xvg'.format(pj)
            shutil.copy('rmsd.xvg', an1)
            directory = str(self.save.getvalue())
            shutil.copy2(an1, directory)
            t, rmsd = numpy.loadtxt("rmsd.xvg", unpack=True)
            fig = plt.figure(figsize=(10, 5))
            data = numpy.loadtxt('rmsd.xvg')
            min = ('%.5f' % (data[0:, 1].min()))
            max = ('%.5f' % (data[0:, 1].max()))
            mean = ('%.5f' % (data[0:, 1].mean()))
            ax = fig.add_subplot(111)
            plt.title(pj+' '+analysis+' '+struct, fontsize=20)
            fig.subplots_adjust(bottom=0.2)
            ax.set_xlabel("time $t$ (ps)")
            ax.set_ylabel(analysis+" (nm)")
            ax.fill_between(t, rmsd, color="blue", linestyle="-", alpha=0.1)
            ax.plot(t, rmsd, color="blue", linestyle="-")
            fig.patch.set_facecolor('white')
            fig.tight_layout()
            plt.show()
            data0 = ('min = ' + str(min) + '\nmax = ' + str(max) + '\nmean =' + str(mean))
            file = open('text.txt', 'w')
            text0 = """
LiGRO v 0.02 - Output of {0}
---------------------------------------------
            """.format(analysis)

            file.write(text0 + data0)
            file.close()
            command31 = 'gedit text.txt'
            os.system(command31)
            mbox.showinfo('Finish', 'Job has finished')
            sys.exit(0)


        elif analysis == 'RMSF':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                quit()
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                quit()
            else:
                pass
            command1 = ''' gmx rmsf -s md_an.tpr -f md_an.xtc -o rmsf.xvg -xvg none -res << EOF
                                                {0}
                                                EOF '''.format(struct)
            os.system(command1)
            pj = str(self.project.getvalue())
            an2 = 'rmsf_{0}.xvg'.format(pj)
            shutil.copy('rmsf.xvg', an2)
            directory = str(self.save.getvalue())
            shutil.copy2(an2, directory)
            t, rmsf = numpy.loadtxt("rmsf.xvg", unpack=True)
            fig = plt.figure(figsize=(10, 5))
            plt.title(pj+' '+analysis+' '+struct, fontsize=20)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(bottom=0.2)
            ax.set_xlabel("Residue")
            ax.set_ylabel(analysis + " (nm)")
            ax.fill_between(t, rmsf, color="blue", linestyle="-", alpha=0.1)
            ax.plot(t, rmsf, color="blue", linestyle="-")
            fig.patch.set_facecolor('white')
            fig.tight_layout()
            plt.show()
            mbox.showinfo('Finish', 'Job has finished')
            sys.exit(0)
        elif analysis == 'RG':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                quit()
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                quit()
            else:
                pass
            command2 = ''' gmx gyrate -s md_an.tpr -f md_an.xtc -o gyrate.xvg -xvg none << EOF
                                                {0}
                                                EOF '''.format(struct)
            os.system(command2)
            pj = str(self.project.getvalue())
            an3 = 'rg_{0}.xvg'.format(pj)
            shutil.copy('gyrate.xvg', an3)
            directory = str(self.save.getvalue())
            shutil.copy2(an3, directory)
            data = numpy.loadtxt('gyrate.xvg')
            fig = plt.figure(figsize=(10, 5))
            g1min = ('%.5f' % (data[0:, 1].min()))
            g2min = ('%.5f' % (data[0:, 2].min()))
            g3min = ('%.5f' % (data[0:, 3].min()))
            g4min = ('%.5f' % (data[0:, 4].min()))
            g1max = ('%.5f' % (data[0:, 1].max()))
            g2max = ('%.5f' % (data[0:, 2].max()))
            g3max = ('%.5f' % (data[0:, 3].max()))
            g4max = ('%.5f' % (data[0:, 3].max()))
            g1mean = ('%.5f' % (data[0:, 1].mean()))
            g2mean = ('%.5f' % (data[0:, 2].mean()))
            g3mean = ('%.5f' % (data[0:, 3].mean()))
            g4mean = ('%.5f' % (data[0:, 3].mean()))
            plt.title(pj+' '+analysis+' '+struct, fontsize=20)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(bottom=0.2)
            ax.set_xlabel("time $t$ (ps)")
            ax.set_ylabel(analysis + " (nm)")
            fig.patch.set_facecolor('white')
            fig.tight_layout()
            t = (data[0:, 0])
            g1 = (data[0:, 1])
            g2 = (data[0:, 2])
            g3 = (data[0:, 3])
            g4 = (data[0:, 4])
            plt.plot(t, g1)  # plotting t,a separately
            plt.plot(t, g2)  # plotting t,b separately
            plt.plot(t, g3)  # plotting t,c separately
            plt.plot(t, g4)  # plotting t,c separately
            plt.show()
            data1 = (
            'Rgmin = ' + str(g1min) + '\t\nRgmax = ' + str(g1max) + '\t\nRgmean =' + str(g1mean) + '\t\nRgxmin = ' + str(
                g2min)
            + '\t\nRgxmax = ' + str(g2max) + '\t\nRgxmean =' + str(g2mean) + '\t\nRgymin = ' + str(
                g3min) + '\t\nRgymax = ' + str(g3max) +
            '\t\nRgymean =' + str(g3mean) + '\t\nRgzmin = ' + str(g4min) + '\t\nRgzmax = ' + str(g4max) + '\t\nRgzmean =' + str(
                g4mean))
            file = open('text.txt', 'w')
            text = """
LiGRO v 0.02 - Output of {0}
---------------------------------------------
            """.format(analysis)

            file.write(text+data1)
            file.close()
            command31 = 'gedit text.txt'
            os.system(command31)
            mbox.showinfo('Finish', 'Job has finished')
            sys.exit(0)

        elif analysis == 'MSD':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                quit()
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                quit()
            else:
                pass
            command= ''' gmx msd -s md_an.tpr -f md_an.xtc -o msd.xvg -xvg none << EOF
                                    {0}
                                    EOF '''.format(struct)
            os.system(command)
            pj = str(self.project.getvalue())
            an4 = 'msd_{0}.xvg'.format(pj)
            shutil.copy('msd.xvg', an4)
            directory = str(self.save.getvalue())
            shutil.copy2(an4, directory)
            t, rmsd = numpy.loadtxt("msd.xvg", unpack=True)
            fig = plt.figure(figsize=(10, 5))
            data = numpy.loadtxt('msd.xvg')
            min = ('%.5f' % (data[0:, 1].min()))
            max = ('%.5f' % (data[0:, 1].max()))
            mean = ('%.5f' % (data[0:, 1].mean()))
            ax = fig.add_subplot(111)
            plt.title(pj+' '+analysis+' '+struct, fontsize=20)
            fig.subplots_adjust(bottom=0.2)
            ax.set_xlabel("time $t$ (ps)")
            ax.set_ylabel(analysis+" (nm^2)")
            ax.fill_between(t, rmsd, color="blue", linestyle="-", alpha=0.1)
            ax.plot(t, rmsd, color="blue", linestyle="-")
            fig.patch.set_facecolor('white')
            fig.tight_layout()
            plt.show()
            data0 = ('min = ' + str(min) + '\nmax = ' + str(max) + '\nmean =' + str(mean))
            file = open('text.txt', 'w')
            text0 = """
LiGRO v 0.02 - Output of {0}
---------------------------------------------
            """.format(analysis)

            file.write(text0 + data0)
            file.close()
            command31 = 'gedit text.txt'
            os.system(command31)
            mbox.showinfo('Finish', 'Job has finished')
            sys.exit(0)

        elif analysis == 'H_bond':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                quit()
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                quit()
            else:
                pass
            command= ''' gmx hbond -s md_an.tpr -f md_an.xtc -num hbond.xvg -xvg none << EOF
                                    Other
                                    Protein
                                    '''
            os.system(command)
            pj = str(self.project.getvalue())
            an6 = 'hbond_{0}.xvg'.format(pj)
            shutil.copy('hbond.xvg', an6)
            directory = str(self.save.getvalue())
            shutil.copy2(an6, directory)

            hbond = pylab.genfromtxt("hbond.xvg")
            x = hbond[:,0]
            y = hbond[:,1]
            fig = plt.figure(figsize=(10, 5))
            ax = fig.add_subplot(111)
            plt.title(pj+' '+analysis+' '+'H_bond', fontsize=20)
            fig.subplots_adjust(bottom=0.2)
            ax.set_xlabel("time $t$ (ps)")
            ax.set_ylabel(analysis+" (num)")
            fig = plt.figure()
            plt.bar(x, y, color = 'cadetblue', width = 100, edgecolor = "none")
            fig.patch.set_facecolor('white')
            fig.tight_layout()
            plt.show()

        elif analysis == 'LJSR-CoulSR IE':
            try:
                os.chdir(path)
            except:
                pass
            find=os.path.exists('md_an.edr')
            if find == False:
                mbox.showerror("Error", "EDR file not found.")
                quit()

            else:
                pass

            com1= '''echo 56 57 0 | gmx energy -f md_an.edr > saida.txt'''
            os.system(com1)


            df = pandas.read_table('saida.txt', delim_whitespace=True, names=('A', 'B', 'C'))

            df.iloc[[6],[1]].to_csv(r'energy.csv', header=None, index=None, sep=',', mode='a')
            df.iloc[[7],[1]].to_csv(r'energy.csv', header=None, index=None, sep=',', mode='a')

            df1 = pandas.read_csv('energy.csv', names=('A'))

            df1.astype(float)

            a = df1.iloc[[0],[0]].values
            b = df1.iloc[[1],[0]].values
            c =float(a+b)

            f = open('saida.txt', 'r')
            cont = f.readlines()
            cont.append('\n\nE int = <E LJ> + <E Coul> = ' + str(c) + ' kJ/mol')
            f.close()
            f = open('saida.txt', 'w')
            f.writelines(cont)
            f.close()


            with open('saida.txt') as infile, open('saida1.txt', 'w') as outfile:
                outfile.write('LiGRO v 0.02\n')
                outfile.write('Energy                      Average   Err.Est.       RMSD  Tot-Drift\n')
                copy = False
                for line in infile:
                    if line.strip() == "Energy                      Average   Err.Est.       RMSD  Tot-Drift":
                        copy = True
                    elif copy:
                        outfile.write(line)

            com2 = 'gedit saida1.txt'
            os.system(com2)

    def plip(self):
        try:
            os.chdir(path)
        except:
            pass
        find1=os.path.exists('md_an.tpr')
        find2=os.path.exists('md_an.xtc')
        if find1 == False:
            mbox.showerror("Error", "TPR or XTC file not found.")
            quit()
        elif find2 == False:
            mbox.showerror("Error", "TPR or XTC file not found.")
            quit()
        else:
            pass
        timeplip = str(self.ft.getvalue())
        command3 = '''gmx trjconv -s md_an.tpr -f md_an.xtc -o plip.pdb -b {0} -e {0} << EOF
                   non-Water
                   EOF'''.format(timeplip)
        command5 = 'gedit report.txt'
        os.system(command3)
        command4="plip -f plip.pdb -t -v"
        subprocess.call(['/bin/bash', '-i', '-c', command4])
        os.system(command5)
        mbox.showinfo('Finish', 'Job has finished')
        sys.exit(0)

        

######################################################################

# Create root window.

if __name__ == '__main__':
    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title(title)
    Tkinter.Label(root, text="LiGRO: Protein-Ligand Molecular Dynamics using ACPYPE, GROMACS and PLIP.").pack(fill='both',
                                                                                                        padx=1,
                                                                                                        pady=1)
    widget = GUI(root)
    root.geometry("1200x600+260+150")
    menubar = Tkinter.Menu(root)
    filemenu = Tkinter.Menu(menubar, tearoff=0)
    filemenu.add_separator()
    editmenu = Tkinter.Menu(menubar, tearoff=0)
    editmenu.add_separator()
    helpmenu = Tkinter.Menu(menubar, tearoff=0)


    def callback():
        if mbox.askyesno('Exit', 'Really quit?'):
            sys.exit(0)
        else:
            mbox.showinfo('No', 'Quit has been cancelled')


    exitButton = Tkinter.Button(root, text='Exit', command=callback)
    exitButton.pack()


    def help():
        webbrowser.open('https://www.ufrgs.br/lasomfarmacia/ligro/crbst_3.html', new=1, autoraise=True)


    def about():
        top = Tkinter.Toplevel()
        top.title('Ligro')
        window_txt = Tkinter.Text(top, width=100, height=10, font="Times 10")

        window_txt.pack()

        my_opening_message = """

        Autor: Luciano Porto Kagami
        Universidade Federal do Rio Grande do Sul
        Laboratorio de Sintese Organica Medicinal -LaSOM
        Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

                                                  """

        window_txt.insert(0.0, my_opening_message)
        button = Tkinter.Button(top, text="Dismiss", command=top.destroy)
        button.pack()
        top.geometry("+500+300")


    helpmenu.add_command(label="Help Index", command=help)
    helpmenu.add_command(label="About...", command=about)
    menubar.add_cascade(label="Help", menu=helpmenu)
    path4=path2+'/ligro'
    img = Tkinter.PhotoImage(file=path4+'/ic_launcher_.png')
    root.tk.call('wm', 'iconphoto', root._w, img)
    root.config(menu=menubar)
    root.mainloop()
