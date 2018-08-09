#!/usr/bin/python3
# -*- coding: utf-8 -*-
note ="""
---LiGRO - Version 0.3 ---

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

For PLIP please cite:

1. Salentin,S. et al. PLIP: fully automated protein-
ligand interaction profiler. Nucl. Acids Res. (1 July 2015) 43 (W1):
W443-W447. doi: 10.1093/nar/gkv315

For GROMACS please cite:

Abraham, M. J., Murtola, T., Schulz, R., Páll, S., Smith, J. C., Hess, B., & Lindahl, E. (2015). 
GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. 
SoftwareX, 1, 19-25.

For Pymol please cite:

Schrödinger, L. L. C. (2017). The PyMOL molecular graphics system, Version 1.8. 2015.

For this code plese cite:

Kagami, L. P., das Neves, G. M., da Silva, A. W. S., Caceres, R. A., Kawano, D. F., & Eifler-Lima, V. L. (2017). 
LiGRO: a graphical user interface for protein–ligand molecular dynamics. 
Journal of molecular modeling, 23(11), 304.

-------------------
"""
import os
import subprocess
import sys
import tkinter as Tkinter
from tkinter import constants as Tkconstants
from tkinter import filedialog as tkFileDialog
import Pmw
import shutil
import webbrowser
from tkinter import messagebox as mbox
import tempfile
import matplotlib.pyplot as plt
import numpy as np
from datetime import date
import re
import pylab
from biopandas.mol2 import PandasMol2
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
import time

sys.path[:0] = ['../../..']

path2 = os.environ.get('HOME')

try:
  path = tempfile.mkdtemp()
except:
    pass
def gromacs_flag(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

title = 'LiGRO: Version 0.3'

class GUI:
    def __init__(self, parent):
        #
        notebook = Pmw.NoteBook(parent)
        notebook.pack(fill='both', expand=2 ,padx=100, pady=10)

        #
        tab1 = notebook.add('Input Files')
        notebook.tab('Input Files').focus_set()

        #
        group_pdb = Pmw.Group(tab1, tag_text='Select PDB file:')
        group_pdb.pack(fill='both', expand=1, padx=100, pady=10)
        self.pdb1 = Tkinter.Label(group_pdb.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.pdb1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_pdb.interior(), text='').grid(row=0, column=2)
        pdb1_bt = Tkinter.Button(group_pdb.interior(), text='Browse...', command=self.openpdbfile)
        pdb1_bt.grid(row=0, column=3)

        group_cof = Pmw.Group(tab1, tag_text='Select Cofactor MOL2 file:')
        group_cof.pack(fill='both', expand=1, padx=100, pady=10)
        self.cof1 = Tkinter.Label(group_cof.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.cof1.grid(row=0, column=1, padx=10, pady=10)
        cof1_bt = Tkinter.Button(group_cof.interior(), text='Browse...', command=self.opencoffile)
        cof1_bt.grid(row=0, column=3)
        #
        group_mol2 = Pmw.Group(tab1, tag_text='Select MOL2 file')
        group_mol2.pack(fill='both', expand=1, padx=100, pady=10)
        self.mol21 = Tkinter.Label(group_mol2.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.mol21.grid(row=0, column=1, padx=10, pady=10)
        mol21_bt = Tkinter.Button(group_mol2.interior(), text='Browse...', command=self.openmol2file)
        mol21_bt.grid(row=0, column=3)
        #
        tab2 = notebook.add('System Preparation')
        group_pdb2gmx = Pmw.Group(tab2, tag_text='Protein Parameters')
        group_pdb2gmx.pack(fill='both', expand=1, padx=100, pady=10)
        self.ff_menu = Pmw.OptionMenu(group_pdb2gmx.interior(),
                                      labelpos='w',
                                      label_text='Select Force Field:',
                                      menubutton_textvariable=None,
                                      items=['amber03', 'amber94', 'amber99', 'amber99sb', 'amber99sb-ildn', 'amberGS',
                                             'oplsaa'],
                                      menubutton_width=15,
                                      )
        self.ff_menu.grid(row=1, column=0, padx=10, pady=10)

        self.metal = Pmw.OptionMenu(group_pdb2gmx.interior(),
                                    labelpos='w',
                                    label_text='Metal:',
                                    menubutton_textvariable=None,
                                    items=['None', 'CA','FE','MG', 'ZN'],
                                    menubutton_width=10,
                                    )
        self.metal.grid(row=1, column=1, padx=10, pady=10)
        self.hb = Tkinter.Checkbutton(group_pdb2gmx.interior(), text='Ignore hydrogen atoms', variable='var1',
                                      onvalue="-ignh", offvalue="")
        self.hb.select()
        self.hb.grid(row=3, column=0, padx=10, pady=10)
        self.hb = Tkinter.Checkbutton(group_pdb2gmx.interior(), text='Ignore hydrogen atoms', variable='var1',
                                      onvalue="-ignh", offvalue="")
        self.hb.select()
        self.hb.grid(row=3, column=0, padx=10, pady=10)
        #
        group_acpype = Pmw.Group(tab2, tag_text='Ligand Parameters (ACPYPE)')
        group_acpype.pack(fill='both', expand=1, padx=100, pady=10)
        self.cm_menu = Pmw.OptionMenu(group_acpype.interior(),
                                      labelpos='w',
                                      label_text='Charge method:',
                                      menubutton_textvariable=None,
                                      items=['bcc','gas', 'user'],
                                      menubutton_width=15,
                                      )
        self.cm_menu.grid(row=0, column=0, padx=1, pady=1)
        self.nc = Pmw.EntryField(group_acpype.interior(), labelpos='w', label_text='Net charge: ', value='0')
        self.nc.grid(row=0, column=1, padx=1, pady=1)
        self.mult = Pmw.EntryField(group_acpype.interior(), labelpos='w', label_text='Multiplicity: ', value='1',  validate = 'numeric')
        self.mult.grid(row=1, column=1, padx=1, pady=1)
        self.at_menu = Pmw.OptionMenu(group_acpype.interior(),
                                      labelpos='w',
                                      label_text='Atom type:        ',
                                      menubutton_textvariable=None,
                                      items=['gaff', 'amber','gaff2', 'amber2'],
                                      menubutton_width=15,
                                      )
        self.at_menu.grid(row=1, column=0, padx=1, pady=1)
        #
        group_solv = Pmw.Group(tab2, tag_text='Solvate')
        group_solv.pack(fill='both', expand=1, padx=100, pady=10)
        self.wt_menu = Pmw.OptionMenu(group_solv.interior(),
                                      labelpos='w',
                                      label_text='Select Water Model:',
                                      menubutton_textvariable=None,
                                      items=['spc', 'spce', 'tip3p'],
                                      menubutton_width=15,
                                      )
        self.wt_menu.grid(row=1, column=0, padx=10, pady=10)
        self.bx_menu = Pmw.OptionMenu(group_solv.interior(),
                                      labelpos='w',
                                      label_text='Select Box Type:     ',
                                      menubutton_textvariable=None,
                                      items=['triclinic', 'cubic', 'dodecahedron', 'octahedron'],
                                      menubutton_width=15,
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
                                       items=['Concentration (mol/liter)', 'Na (Number)', 'Cl (Number)'],
                                       menubutton_width=25,
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
                                       menubutton_width=15,
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
        group_options.pack(fill='both', expand=1, padx=10, pady=10)
        td = str(date.today())
        self.project = Pmw.EntryField(group_options.interior(), labelpos='w', value='MD_'+td, label_text='Project Name: ')
        self.project.grid(row=0, column=0, padx=20, pady=10)
        self.bkp = Tkinter.Checkbutton(group_options.interior(), text='Explain MD folder', variable='var2', onvalue=True, offvalue=False)
        self.bkp.grid(row=2, column=0, padx=10, pady=10)
        save_bt = Tkinter.Button(group_options.interior(), text='Save TPR file', command=self.choose_simulation_save_tpr_file)
        save_bt.grid(row=3, column=0)
        run_bt = Tkinter.Button(group_options.interior(), text='Run Dynamics', command=self.choose_run_simulation)
        run_bt.grid(row=3, column=1)
        self.save = Pmw.EntryField(group_options.interior(), labelpos='w', value=path2, label_text='Save directory:')
        self.save.grid(row=1, column=0, padx=20, pady=5)
        lie_group_min = Pmw.Group(tab4, tag_text='LIE Binding Energy Calculation')
        lie_group_min.pack(fill='both', expand=1, padx=10, pady=10)
        Tkinter.Label(lie_group_min.interior(), text='DG (kJ/mol)').grid(row=0, column=0)
        self.lie_ener= Tkinter.Label(lie_group_min.interior(), text='Empty', bg='white', padx=10, pady=5)
        self.lie_ener.grid(row=0, column=1, padx=40, pady=10)
        lie_run_bt = Tkinter.Button(lie_group_min.interior(), text='Run LIE calculation', command=self.mount_simulation_LIE_complex)
        lie_run_bt.grid(row=2, column=1)
        group_backup = Pmw.Group(tab4, tag_text='Backup (Only MD)')
        group_backup.pack(fill='both', expand=1, padx=10, pady=10)
        self.bkp_path = Tkinter.Label(group_backup.interior(), text='No select folder', bg='white', padx=100, pady=5)
        self.bkp_path.grid(row=1, column=0, padx=10, pady=10)
        bkp_bt = Tkinter.Button(group_backup.interior(), text='Browse...', command=self.open_backup_folder)
        bkp_bt.grid(row=1, column=2)
        run_bkp_bt = Tkinter.Button(group_backup.interior(), text='Run backup', command=self.run_bkp_file)
        run_bkp_bt.grid(row=2, column=0)
        #
        tab5 = notebook.add('Analysis')
        group_tpr = Pmw.Group(tab5, tag_text='Select TPR file:')
        group_tpr.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_tpr.interior(), text='                      ').grid(row=0, column=0)
        self.tpr1 = Tkinter.Label(group_tpr.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.tpr1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_tpr.interior(), text='                     ').grid(row=0, column=2)
        tpr1_bt = Tkinter.Button(group_tpr.interior(), text='Browse...', command=self.opentprfile)
        tpr1_bt.grid(row=0, column=3)
        #
        group_xtc = Pmw.Group(tab5, tag_text='Select XTC file')
        group_xtc.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_xtc.interior(), text='').grid(row=0, column=0)
        self.xtc = Tkinter.Label(group_xtc.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.xtc.grid(row=0, column=1, padx=100, pady=10)
        xtc1_bt = Tkinter.Button(group_xtc.interior(), text='Browse...', command=self.openxtcfile)
        xtc1_bt.grid(row=0, column=2)
        #
        group_edr = Pmw.Group(tab5, tag_text='Select EDR file:')
        group_edr.pack(fill='both', expand=1, padx=100, pady=10)
        Tkinter.Label(group_edr.interior(), text='                      ').grid(row=0, column=0)
        self.edr1 = Tkinter.Label(group_edr.interior(), text='No select file', bg='white', padx=100, pady=5)
        self.edr1.grid(row=0, column=1, padx=10, pady=10)
        Tkinter.Label(group_edr.interior(), text='                     ').grid(row=0, column=2)
        edr1_bt = Tkinter.Button(group_edr.interior(), text='Browse...', command=self.openedrfile)
        edr1_bt.grid(row=0, column=3)
        #
        group_analysis = Pmw.Group(tab5, tag_text='Analysis')
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
                                             menubutton_width=15,
                                             )
        self.analysis_menu.grid(row=0, column=2)
        run_analysis_bt = Tkinter.Button(group_analysis.interior(), text='Run', command=self.RMSD)
        run_analysis_bt.grid(row=0, column=4)
        Tkinter.Label(group_analysis.interior(),text='           ').grid(row=0, column=3)
        group_plip = Pmw.Group(tab5, tag_text='Protein-Ligand Interaction Profiler (PLIP) v1.4.2')
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

    def openmol2file(self):
        try:
            self.mol2file = tkFileDialog.askopenfilename(initialdir=path2, title="Select MOL2 file",
                                                         filetypes=(("MOL2 files", "*.mol2"), ("all files", "*.*")))
            self.mol21['text'] = self.mol2file
            shutil.copy(self.mol2file, path+'/Ligand.mol2')
        except:
            self.mol21['text'] = 'No select file'

    def opencoffile(self):
        try:
            self.coffile = tkFileDialog.askopenfilename(initialdir=path2, title="Select MOL2 file",
                                                         filetypes=(("MOL2 files", "*.mol2"), ("all files", "*.*")))
            self.cof1['text'] = self.coffile
            shutil.copy(self.coffile, path+'/Cofactor.mol2')
        except:
            self.cof1['text'] = 'No select file'

    def open_backup_folder(self):
        try:
            self.bkp_folder = tkFileDialog.askdirectory(initialdir=path2, title="Select Backup folder")
            self.bkp_path['text'] = self.bkp_folder
            self.path4 = self.bkp_folder
        
        except:
            self.bkp_path['text'] = 'No select folder' 

    def choose_simulation_save_tpr_file(self):

      dr = str(self.save.getvalue())
      pj = str(self.project.getvalue())
      pj1 = pj
      bk = self.bkp.getvar('var2')
      path3 = dr+'/'+pj
      if bk == True:
        if not os.path.exists(path3):
          os.makedirs(path3)
          try:
            shutil.copy(path+'/protein.pdb',path3)
          except:
            pass
          try:
            shutil.copy(path+'/Ligand.mol2',path3)
          except:
            pass
          try:
            shutil.copy(path+'/Cofactor.mol2',path3)
          except:
            pass
          try:
            os.chdir(path3)
          except:
            pass
        else:
          mbox.showinfo("INFO", "Please delete backup folder or change project name and try again")
          pass
      elif bk == False:
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
          pass
      elif find3 == False:
          mbox.showerror("Error", "PDB file not found.")
          pass

      elif find2 == True and find1 == False:
          mbox.showerror("Error", "It is not possible run Protein-Cofactor")
          pass

      else:
          mbox.showinfo("INFO", "Run Protein-Ligand simulation with Cofactor")
          self.mount_simulation_cof()
          self.save_tprfile()
          pass

    def choose_run_simulation(self):
      dr = str(self.save.getvalue())
      pj = str(self.project.getvalue())
      pj1 = pj
      bk = self.bkp.getvar('var2')
      path3 = dr+'/'+pj
      if bk == True:
        if not os.path.exists(path3):
          os.makedirs(path3)
          try:
            shutil.copy(path+'/protein.pdb',path3)
          except:
            pass
          try:
            shutil.copy(path+'/Ligand.mol2',path3)
          except:
            pass
          try:
            shutil.copy(path+'/Cofactor.mol2',path3)
          except:
            pass
          try:
            os.chdir(path3)
          except:
            pass
        else:
          mbox.showinfo("INFO", "Please delete backup folder or change project name and try again")
          quit()
      elif bk == False:
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
          pass
      elif find3 == False:
          mbox.showerror("Error", "PDB file not found.")
          pass

      elif find2 == True and find1 == False:
          mbox.showerror("Error", "It is not possible run Protein-Cofactor")
          pass

      else:
          mbox.showinfo("INFO", "Run Protein-Ligand simulation with Cofactor")
          self.mount_simulation_cof()
          self.run_simulation()
          pass

    def run_bkp_file(self):
        
        if self.path4 is not None:
          print(self.path4)
          find_tpr=os.path.exists(self.path4 + '/md.tpr')
          find_cpt=os.path.exists(self.path4 + '/md.cpt')
          if find_tpr == False:
            mbox.showerror("Error", "TPR file not found. Please try again.")
            pass
          elif find_cpt == False:
            mbox.showerror("Error", "CPT file not found. Please try again.")
            pass
          else:
            if gromacs_flag('mdrun'):
              cmd9 = 'mdrun -s {0}/md.tpr -cpi {0}/md.cpt -append no'.format(self.path4)
              os.system(cmd9)
            elif gromacs_flag('gmx'):
              cmd9 = 'gmx mdrun -s {0}/md.tpr -cpi {0}/md.cpt -append no'.format(self.path4)
              os.system(cmd9)

        else:
          mbox.showerror("Error", "Backup folder not found. Please try again.")
          pass


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
        ig = self.hb.getvar('var1')
        if gromacs_flag('mdrun'):
          cmd = 'pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
          out_pdb2gmx = os.system(cmd)
        elif gromacs_flag('gmx'):
          cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
          out_pdb2gmx = os.system(cmd)
          
        if out_pdb2gmx == 0:  
          try:
            pmol = PandasMol2().read_mol2('Ligand.mol2')
            subst = pmol.df.iloc[0]['subst_name']
            fh, abs_path = mkstemp()
            with fdopen(fh,'w') as new_file:
              with open('Ligand.mol2') as old_file:
                for line in old_file:
                  new_file.write(line.replace(subst, 'LIG'))
            #Remove original file
            remove('Ligand.mol2')
            #Move new file
            move(abs_path, 'Ligand.mol2')
          except:
            mbox.showerror("Error", "Mol2 file not recognized.. Please try again.")
            quit()

          cm = str(self.cm_menu.getcurselection())
          nc = str(self.nc.getvalue())
          mt = str(self.mult.getvalue())
          at = str(self.at_menu.getcurselection())
          cmdd = 'acpype -i Ligand.mol2 -c {0} -n {1} -m {2} -a {3}'.format(cm, nc, mt, at)
          acp = os.system(cmdd)
          if acp == 0:
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
            if gromacs_flag('mdrun'):
              cmd1 = 'editconf -bt {0} -f complex.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
              os.system(cmd1)
            elif gromacs_flag('gmx'):
              cmd1 = 'gmx editconf -bt {0} -f complex.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
              os.system(cmd1)
            if gromacs_flag('mdrun'):
              os.system('genbox -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
            elif gromacs_flag('gmx'):
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
            if gromacs_flag('mdrun'):
              os.system('grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')
            elif gromacs_flag('gmx'):
              os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')

            if self.ion_menu.getcurselection() == 'Na (Number)':
                c_ion = '-np ' + self.ion_cc.getvalue()
            elif self.ion_menu.getcurselection() == 'Cl (Number)':
                c_ion = '-nn ' + self.ion_cc.getvalue()
            else:
                c_ion = '-conc ' + self.ion_cc.getvalue()
            if gromacs_flag('mdrun'):
              cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              os.system(cmd2)
            elif gromacs_flag('gmx'):
              cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              
            os.system(cmd2)

            inte = self.min_menu.getcurselection()

            nst = self.step.getvalue()

            if gromacs_flag('mdrun'):
              cmd3 = '''
              cat << EOF >| em_real.mdp
              ; LINES STARTING WITH ';' ARE COMMENTS
              title   = Minimization  ; Title of run

              ; Parameters describing what to do, when to stop and what to save
              integrator  = steep   ; Algorithm (steep = steepest descent minimization)
              emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
              emstep      = 0.01      ; Energy step size
              nsteps    = {0}     ; Maximum number of (minimization) steps to perform
              energygrps  = Protein LIG ; Which energy group(s) to write to disk

              ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
              nstlist   = 1       ; Frequency to update the neighbor list and long range forces
              ns_type   = grid    ; Method to determine neighbor list (simple, grid)
              rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
              coulombtype = PME   ; Treatment of long range electrostatic interactions
              rcoulomb  = 1.0   ; long range electrostatic cut-off
              rvdw    = 1.0   ; long range Van der Waals cut-off
              pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
              '''.format(nst)
              cmd4 = '''
              cat << EOF >| em_real.mdp
              ; LINES STARTING WITH ';' ARE COMMENTS
              title   = Minimization  ; Title of run

              ; Parameters describing what to do, when to stop and what to save
              integrator  = cg    ; Algorithm (steep = steepest descent minimization)
              emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
              emstep      = 0.01      ; Energy step size
              nsteps    = {0}     ; Maximum number of (minimization) steps to perform
              energygrps  = Protein LIG ; Which energy group(s) to write to disk

              ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
              nstlist   = 1       ; Frequency to update the neighbor list and long range forces
              ns_type   = grid    ; Method to determine neighbor list (simple, grid)
              rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
              coulombtype = PME   ; Treatment of long range electrostatic interactions
              rcoulomb  = 1.0   ; long range electrostatic cut-off
              rvdw    = 1.0   ; long range Van der Waals cut-off
              pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
              '''.format(nst)

            elif gromacs_flag('gmx'):

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
            else:
              pass

            if inte == 'SD Algorithm':
                os.system(cmd3)
                if gromacs_flag('mdrun'):
                  os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                  os.system('mdrun -v -deffnm em')
                  os.system('''
                  make_ndx -f em.gro -o index2.ndx << EOF
                  "Protein" | "Other"
                  q
                  EOF ''')
                  os.system('''
                  trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                  Protein_Other
                  System
                  EOF ''')
                  os.system('''
                  mv em_2.gro em.gro ''')
                  os.system('''
                  trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                  Protein_Other
                  System
                  EOF
                  ''')
                  os.system('''
                          mv em_2.pdb em.pdb ''')
                elif gromacs_flag('gmx'):
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
              if gromacs_flag('mdrun'):
                os.system(cmd3)
                os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm em')
                os.system(cmd4)
                os.system('grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm em')
                os.system('''
                            make_ndx -f em.gro -o index2.ndx << EOF
                            "Protein" | "Other"
                            q
                            EOF ''')
                os.system('''
                            trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                            Protein_Other
                            System
                            EOF ''')
                os.system('''
                            mv em_2.gro em.gro ''')
                os.system('''
                            trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                            Protein_Other
                            System
                            EOF
                            ''')
                os.system('mv em_2.pdb em.pdb')

              elif gromacs_flag('gmx'):
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

              else:
                pass

            stnvt = str(self.time_nvt.getvalue())
            stnpt = str(self.time_npt.getvalue())
            stmd = str(self.time_md.getvalue())
            temp = str(self.temperature.getvalue())
            t_st = str(self.time_st.getvalue())

            if gromacs_flag('mdrun'):
              cmd5 = '''
              cat << EOF >| nvt.mdp
              title       = Protein-ligand complex NVT equilibration
              define      = -DPOSRES  ; position restrain the protein and ligand
              ; Run parameters
              integrator  = md        ; leap-frog integrator
              nsteps      = {0}     ; 2 * 50000 = 100 ps
              dt          = {1}     ; 2 fs
              ; Output control
              nstxout     = 100       ; save coordinates every 0.2 ps
              nstvout     = 100       ; save velocities every 0.2 ps
              nstenergy   = 100       ; save energies every 0.2 ps
              nstlog      = 100       ; update log file every 0.2 ps
              energygrps  = Protein LIG
              ; Bond parameters
              continuation    = no            ; first dynamics run
              constraint_algorithm = lincs    ; holonomic constraints
              constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
              lincs_iter      = 1             ; accuracy of LINCS
              lincs_order     = 4             ; also related to accuracy
              ; Neighborsearching
              ns_type     = grid      ; search neighboring grid cells
              nstlist     = 5         ; 10 fs
              rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
              rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
              rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
              ; Electrostatics
              coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
              pme_order       = 4         ; cubic interpolation
              fourierspacing  = 0.16      ; grid spacing for FFT
              ; Temperature coupling
              tcoupl      = V-rescale                     ; modified Berendsen thermostat
              tc-grps     = Protein non-Protein   ; two coupling groups - more accurate
              tau_t       = 0.1   0.1                    ; time constant, in ps
              ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
              ; Pressure coupling
              pcoupl      = no        ; no pressure coupling in NVT
              ; Periodic boundary conditions
              pbc         = xyz      ; 3-D PBC
              ; Dispersion correction
              DispCorr    = EnerPres  ; account for cut-off vdW scheme
              ; Velocity generation
              gen_vel     = yes       ; assign velocities from Maxwell distribution
              gen_temp    = {2}       ; temperature for Maxwell distribution
              gen_seed    = -1        ; generate a random seed'
              '''.format(stnvt, t_st, temp)

            elif gromacs_flag('gmx'):
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
            else:
              pass

            os.system(cmd5)

            if gromacs_flag('mdrun'):
              cmd6 = '''
              cat << EOF >| npt.mdp
              title       = Protein-ligand complex NPT equilibration
              define      = -DPOSRES  ; position restrain the protein and ligand
              ; Run parameters
              integrator  = md        ; leap-frog integrator
              nsteps      = {0}     ; 2 * 50000 = 100 ps
              dt          = {1}     ; 2 fs
              ; Output control
              nstxout     = 100       ; save coordinates every 0.2 ps
              nstvout     = 100       ; save velocities every 0.2 ps
              nstenergy   = 100       ; save energies every 0.2 ps
              nstlog      = 100       ; update log file every 0.2 ps
              energygrps  = Protein LIG
              ; Bond parameters
              continuation    = yes           ; first dynamics run
              constraint_algorithm = lincs    ; holonomic constraints
              constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
              lincs_iter      = 1             ; accuracy of LINCS
              lincs_order     = 4             ; also related to accuracy
              ; Neighborsearching
              ns_type     = grid      ; search neighboring grid cells
              nstlist     = 5         ; 10 fs
              rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
              rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
              rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
              ; Electrostatics
              coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
              pme_order       = 4         ; cubic interpolation
              fourierspacing  = 0.16      ; grid spacing for FFT
              ; Temperature coupling
              tcoupl      = V-rescale                     ; modified Berendsen thermostat
              tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
              tau_t       = 0.1   0.1                    ; time constant, in ps
              ref_t       = {2}  {2}                    ; reference temperature, one for each group, in K
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
              '''.format(stnpt, t_st, temp)


            elif gromacs_flag('gmx'):

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
            else:
              pass

            os.system(cmd6)


            if gromacs_flag('mdrun'):
              cmd7 = '''
              cat << EOF >| md.mdp
              title       = Protein-ligand complex MD simulation
              ; Run parameters
              integrator  = md        ; leap-frog integrator
              nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
              dt          = {1}     ; 2 fs
              ; Output control
              nstxout                  = 0                     ; [steps] freq to write coordinates to trajectory
              nstvout                  = 0                 ; [steps] freq to write velocities to trajectory
              nstfout                  = 0                     ; [steps] freq to write forces to trajectory
              nstlog                   = 100                   ; [steps] freq to write energies to log file
              nstenergy                = 500                   ; [steps] freq to write energies to energy file
              nstxtcout                = 500                   ; [steps] freq to write coordinates to xtc trajectory
              xtc_precision            = 1000                  ; [real] precision to write xtc trajectory
              xtc_grps                 = System                ; group(s) to write to xtc trajectory
              ; Bond parameters
              continuation    = yes           ; first dynamics run
              constraint_algorithm = lincs    ; holonomic constraints
              constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
              lincs_iter      = 1             ; accuracy of LINCS
              lincs_order     = 4             ; also related to accuracy
              ; Neighborsearching
              ns_type     = grid      ; search neighboring grid cells
              nstlist     = 5         ; 10 fs
              rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
              rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
              rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
              ; Electrostatics
              coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
              pme_order       = 4         ; cubic interpolation
              fourierspacing  = 0.16      ; grid spacing for FFT
              ; Temperature coupling
              tcoupl      = V-rescale                     ; modified Berendsen thermostat
              tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
              tau_t       = 0.1   0.1                    ; time constant, in ps
              ref_t       = {2}   {2}                   ; reference temperature, one for each group, in K
              ; Pressure coupling
              pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
              pcoupltype  = isotropic                     ; uniform scaling of box vectors
              tau_p       = 2.0                           ; time constant, in ps
              ref_p       = 1.0                           ; reference pressure, in bar
              compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
              ; Periodic boundary conditions
              pbc         = xyz      ; 3-D PBC
              ; Dispersion correction
              DispCorr    = EnerPres  ; account for cut-off vdW scheme
              ; Velocity generation
              gen_vel     = no        ; assign velocities from Maxwell distribution
              '''.format(stmd, t_st, temp)

            elif gromacs_flag('gmx'):
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
            else:
              pass

            os.system(cmd7)

            cmd8='pymol em.pdb'

            os.system(cmd8)

            if mbox.askyesno('View Complex', 'Is complex OK??'):
                pass
            else:
                mbox.showinfo('No', 'Process has been cancelled')
                quit()

            if gromacs_flag('mdrun'):
              os.system('echo 2|genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
              os.system(r'''
              sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
              ''')
              os.system('mv trp2.top trp.top')

              os.system('cat << EOF > | queue.sh')

              os.system("echo 'grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -maxwarn 1000' >> queue.sh")
              os.system("echo 'mdrun -v -deffnm nvt' >> queue.sh")
              os.system("echo 'grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -maxwarn 1000' >> queue.sh")
              os.system("echo 'mdrun -v -deffnm npt' >> queue.sh")
              os.system("echo 'grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -maxwarn 1000' >> queue.sh")

            elif gromacs_flag('gmx'):
              os.system('echo 2 | gmx genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
              os.system(r'''
              sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
              ''')
              os.system('mv trp2.top trp.top')

              os.system('cat << EOF > | queue.sh')

              os.system("echo 'gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -r em.gro -maxwarn 1000' >> queue.sh")
              os.system("echo 'gmx mdrun -v -deffnm nvt' >> queue.sh")
              os.system("echo 'gmx grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -r nvt.gro -maxwarn 1000' >> queue.sh")
              os.system("echo 'gmx mdrun -v -deffnm npt' >> queue.sh")
              os.system("echo 'gmx grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -r npt.gro -maxwarn 1000' >> queue.sh")
          else:
            mbox.showerror('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            quit()
        else:
          mbox.showerror('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

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
        ig = self.hb.getvar('var1')
        if gromacs_flag('mdrun'):
          cmd = 'pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        elif gromacs_flag('gmx'):
          cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        else:
          pass
        out_pdb2gmx = os.system(cmd)
        if out_pdb2gmx == 0:
          try:
            pmol = PandasMol2().read_mol2('Ligand.mol2')
            subst = pmol.df.iloc[0]['subst_name']
            fh, abs_path = mkstemp()
            with fdopen(fh,'w') as new_file:
              with open('Ligand.mol2') as old_file:
                for line in old_file:
                  new_file.write(line.replace(subst, 'LIG'))
            #Remove original file
            remove('Ligand.mol2')
            #Move new file
            move(abs_path, 'Ligand.mol2')
          except:
            mbox.showerror("Error", "Mol2 (ligand) file not recognized.. Please try again.")
            quit()

          
          try:
            pmol = PandasMol2().read_mol2('Cofactor.mol2')
            subst = pmol.df.iloc[0]['subst_name']
            fh, abs_path = mkstemp()
            with fdopen(fh,'w') as new_file:
              with open('Cofactor.mol2') as old_file:
                for line in old_file:
                  new_file.write(line.replace(subst, 'COF'))
            #Remove original file
            remove('Cofactor.mol2')
            #Move new file
            move(abs_path, 'Cofactor.mol2')
          except:
            mbox.showerror("Error", "Mol2 (cofactor) file not recognized.. Please try again.")
            quit()  
          
          cm = str(self.cm_menu.getcurselection())
          nc = str(self.nc.getvalue())
          mt = str(self.mult.getvalue())
          at = str(self.at_menu.getcurselection())
          cmdd0 = 'acpype -i Ligand.mol2 -c {0} -n {1} -m {2} -a {3}'.format(cm, nc, mt, at)
          cmdd = 'acpype -i Cofactor.mol2 -c {0} -n {1} -m {2} -a {3}'.format(cm, nc, mt, at)
          acp0 = os.system(cmdd0)
          acp1 = os.system(cmdd)
          
          if acp0 == 0:
            if acp1 == 0:
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

              if gromacs_flag('mdrun'):
                cmd1 = 'editconf -bt {0} -f complex2.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
              elif gromacs_flag('gmx'):
                cmd1 = 'gmx editconf -bt {0} -f complex2.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
              else:
                pass
              os.system(cmd1)
              if gromacs_flag('mdrun'):
                os.system('genbox -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
              elif gromacs_flag('gmx'):
                os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')  
              else:
                pass
              
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

              if self.ion_menu.getcurselection() == 'Na (Number)':
                  c_ion = '-np ' + self.ion_cc.getvalue()
              elif self.ion_menu.getcurselection() == 'Cl (Number)':
                  c_ion = '-nn ' + self.ion_cc.getvalue()
              else:
                  c_ion = '-conc ' + self.ion_cc.getvalue()
              cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              os.system(cmd2)

              inte = self.min_menu.getcurselection()

              nst = self.step.getvalue()

              if gromacs_flag('mdrun'):
                cmd3 = '''
                    cat << EOF >| em_real.mdp
                    ; LINES STARTING WITH ';' ARE COMMENTS
                    title   = Minimization  ; Title of run

                    ; Parameters describing what to do, when to stop and what to save
                    integrator  = steep   ; Algorithm (steep = steepest descent minimization)
                    emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
                    emstep      = 0.01      ; Energy step size
                    nsteps    = {0}     ; Maximum number of (minimization) steps to perform
                    energygrps  = Protein non-Protein ; Which energy group(s) to write to disk

                    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                    nstlist   = 1       ; Frequency to update the neighbor list and long range forces
                    ns_type   = grid    ; Method to determine neighbor list (simple, grid)
                    rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
                    coulombtype = PME   ; Treatment of long range electrostatic interactions
                    rcoulomb  = 1.0   ; long range electrostatic cut-off
                    rvdw    = 1.0   ; long range Van der Waals cut-off
                    pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
                    '''.format(nst)
                cmd4 = '''
                            cat << EOF >| em_real.mdp
                            ; LINES STARTING WITH ';' ARE COMMENTS
                            title   = Minimization  ; Title of run

                            ; Parameters describing what to do, when to stop and what to save
                            integrator  = cg    ; Algorithm (steep = steepest descent minimization)
                            emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
                            emstep      = 0.01      ; Energy step size
                            nsteps    = {0}     ; Maximum number of (minimization) steps to perform
                            energygrps  = Protein non-Protein ; Which energy group(s) to write to disk

                            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                            nstlist   = 1       ; Frequency to update the neighbor list and long range forces
                            ns_type   = grid    ; Method to determine neighbor list (simple, grid)
                            rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
                            coulombtype = PME   ; Treatment of long range electrostatic interactions
                            rcoulomb  = 1.0   ; long range electrostatic cut-off
                            rvdw    = 1.0   ; long range Van der Waals cut-off
                            pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
                            '''.format(nst)

              elif gromacs_flag('gmx'):
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
              else:
                pass

              if inte == 'SD Algorithm':
                  os.system(cmd3)
                  if gromacs_flag('mdrun'):
                    os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                    os.system('mdrun -v -deffnm em')
                    os.system('''
                    make_ndx -f em.gro -o index2.ndx << EOF
                    "Protein" | "Other"
                    q
                    EOF ''')
                    os.system('''
                    trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                    Protein_Other
                    System
                    EOF ''')
                    os.system('''
                    mv em_2.gro em.gro ''')
                    os.system('''
                    trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                    Protein_Other
                    System
                    EOF
                    ''')
                    os.system('''
                            mv em_2.pdb em.pdb ''')

                  elif gromacs_flag('gmx'):

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
                    pass

              else:
                  os.system(cmd3)
                  if gromacs_flag('mdrun'):
                    os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                    os.system('mdrun -v -deffnm em')
                    os.system(cmd4)
                    os.system('grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr -maxwarn 1000')
                    os.system('mdrun -v -deffnm em')
                    os.system('''
                                  make_ndx -f em.gro -o index2.ndx << EOF
                                  "Protein" | "Other"
                                  q
                                  EOF ''')
                    os.system('''
                                  trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                                  Protein_Other
                                  System
                                  EOF ''')
                    os.system('''
                                  mv em_2.gro em.gro ''')
                    os.system('''
                                  trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                                  Protein_Other
                                  System
                                  EOF
                                  ''')
                    os.system('mv em_2.pdb em.pdb')

                  elif gromacs_flag('gmx'):
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
                  else:
                      pass


              stnvt = str(self.time_nvt.getvalue())
              stnpt = str(self.time_npt.getvalue())
              stmd = str(self.time_md.getvalue())
              temp = str(self.temperature.getvalue())
              t_st = str(self.time_st.getvalue())

              if gromacs_flag('mdrun'):
                cmd5 = '''
                cat << EOF >| nvt.mdp
                title       = Protein-ligand complex NVT equilibration
                define      = -DPOSRES  ; position restrain the protein and ligand
                ; Run parameters
                integrator  = md        ; leap-frog integrator
                nsteps      = {0}     ; 2 * 50000 = 100 ps
                dt          = {1}     ; 2 fs
                ; Output control
                nstxout     = 100       ; save coordinates every 0.2 ps
                nstvout     = 100       ; save velocities every 0.2 ps
                nstenergy   = 100       ; save energies every 0.2 ps
                nstlog      = 100       ; update log file every 0.2 ps
                energygrps  = Protein LIG
                ; Bond parameters
                continuation    = no            ; first dynamics run
                constraint_algorithm = lincs    ; holonomic constraints
                constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
                lincs_iter      = 1             ; accuracy of LINCS
                lincs_order     = 4             ; also related to accuracy
                ; Neighborsearching
                ns_type     = grid      ; search neighboring grid cells
                nstlist     = 5         ; 10 fs
                rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
                rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
                rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
                ; Electrostatics
                coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
                pme_order       = 4         ; cubic interpolation
                fourierspacing  = 0.16      ; grid spacing for FFT
                ; Temperature coupling
                tcoupl      = V-rescale                     ; modified Berendsen thermostat
                tc-grps     = Protein non-Protein   ; two coupling groups - more accurate
                tau_t       = 0.1   0.1                    ; time constant, in ps
                ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
                ; Pressure coupling
                pcoupl      = no        ; no pressure coupling in NVT
                ; Periodic boundary conditions
                pbc         = xyz      ; 3-D PBC
                ; Dispersion correction
                DispCorr    = EnerPres  ; account for cut-off vdW scheme
                ; Velocity generation
                gen_vel     = yes       ; assign velocities from Maxwell distribution
                gen_temp    = {2}       ; temperature for Maxwell distribution
                gen_seed    = -1        ; generate a random seed'
                '''.format(stnvt, t_st, temp)


              elif gromacs_flag('gmx'):
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

              else:
                pass

              os.system(cmd5)

              if gromacs_flag('mdrun'):
                cmd6 = '''
                cat << EOF >| npt.mdp
                title       = Protein-ligand complex NPT equilibration
                define      = -DPOSRES  ; position restrain the protein and ligand
                ; Run parameters
                integrator  = md        ; leap-frog integrator
                nsteps      = {0}     ; 2 * 50000 = 100 ps
                dt          = {1}     ; 2 fs
                ; Output control
                nstxout     = 100       ; save coordinates every 0.2 ps
                nstvout     = 100       ; save velocities every 0.2 ps
                nstenergy   = 100       ; save energies every 0.2 ps
                nstlog      = 100       ; update log file every 0.2 ps
                energygrps  = Protein LIG
                ; Bond parameters
                continuation    = yes           ; first dynamics run
                constraint_algorithm = lincs    ; holonomic constraints
                constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
                lincs_iter      = 1             ; accuracy of LINCS
                lincs_order     = 4             ; also related to accuracy
                ; Neighborsearching
                ns_type     = grid      ; search neighboring grid cells
                nstlist     = 5         ; 10 fs
                rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
                rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
                rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
                ; Electrostatics
                coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
                pme_order       = 4         ; cubic interpolation
                fourierspacing  = 0.16      ; grid spacing for FFT
                ; Temperature coupling
                tcoupl      = V-rescale                     ; modified Berendsen thermostat
                tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
                tau_t       = 0.1   0.1                    ; time constant, in ps
                ref_t       = {2}  {2}                    ; reference temperature, one for each group, in K
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
                '''.format(stnpt, t_st, temp)

              elif gromacs_flag('gmx'):
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
              else:
                pass

              os.system(cmd6)

              if gromacs_flag('mdrun'):
                  cmd7 = '''
                  cat << EOF >| md.mdp
                  title       = Protein-ligand complex MD simulation
                  ; Run parameters
                  integrator  = md        ; leap-frog integrator
                  nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
                  dt          = {1}     ; 2 fs
                  ; Output control
                  nstxout                  = 0                     ; [steps] freq to write coordinates to trajectory
                  nstvout                  = 0                 ; [steps] freq to write velocities to trajectory
                  nstfout                  = 0                     ; [steps] freq to write forces to trajectory
                  nstlog                   = 100                   ; [steps] freq to write energies to log file
                  nstenergy                = 500                   ; [steps] freq to write energies to energy file
                  nstxtcout                = 500                   ; [steps] freq to write coordinates to xtc trajectory
                  xtc_precision            = 1000                  ; [real] precision to write xtc trajectory
                  xtc_grps                 = System                ; group(s) to write to xtc trajectory
                  ; Bond parameters
                  continuation    = yes           ; first dynamics run
                  constraint_algorithm = lincs    ; holonomic constraints
                  constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
                  lincs_iter      = 1             ; accuracy of LINCS
                  lincs_order     = 4             ; also related to accuracy
                  ; Neighborsearching
                  ns_type     = grid      ; search neighboring grid cells
                  nstlist     = 5         ; 10 fs
                  rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
                  rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
                  rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
                  ; Electrostatics
                  coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
                  pme_order       = 4         ; cubic interpolation
                  fourierspacing  = 0.16      ; grid spacing for FFT
                  ; Temperature coupling
                  tcoupl      = V-rescale                     ; modified Berendsen thermostat
                  tc-grps     = Protein non-Protein    ; two coupling groups - more accurate
                  tau_t       = 0.1   0.1                    ; time constant, in ps
                  ref_t       = {2}   {2}                   ; reference temperature, one for each group, in K
                  ; Pressure coupling
                  pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
                  pcoupltype  = isotropic                     ; uniform scaling of box vectors
                  tau_p       = 2.0                           ; time constant, in ps
                  ref_p       = 1.0                           ; reference pressure, in bar
                  compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
                  ; Periodic boundary conditions
                  pbc         = xyz      ; 3-D PBC
                  ; Dispersion correction
                  DispCorr    = EnerPres  ; account for cut-off vdW scheme
                  ; Velocity generation
                  gen_vel     = no        ; assign velocities from Maxwell distribution
                  '''.format(stmd, t_st, temp)
              elif gromacs_flag('gmx'):
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
              else:
                  pass

              os.system(cmd7)

              cmd8='pymol em.pdb'

              os.system(cmd8)

              if mbox.askyesno('View Complex', 'Is complex OK??'):
                  pass
              else:
                  mbox.showinfo('No', 'Process has been cancelled')
                  quit()
              if gromacs_flag('mdrun'):
                os.system('echo 2|genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
                os.system(r'''
                sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
                ''')
                os.system('mv trp2.top trp.top')

                os.system('cat << EOF > | queue.sh')

                os.system("echo 'grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -maxwarn 1000' >> queue.sh")
                os.system("echo  'mdrun -v -deffnm nvt' >> queue.sh")
                os.system("echo 'grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -maxwarn 1000' >> queue.sh")
                os.system("echo 'mdrun -v -deffnm npt' >> queue.sh")
                os.system("echo 'grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -maxwarn 1000' >> queue.sh")

              elif gromacs_flag('gmx'):
                os.system('echo 2 | gmx genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
                os.system(r'''
                sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
                ''')
                os.system('mv trp2.top trp.top')

                os.system('cat << EOF > | queue.sh')

                os.system("echo 'gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -maxwarn 1000 -r em.gro' >> queue.sh")
                os.system("echo 'gmx mdrun -v -deffnm nvt' >> queue.sh")
                os.system("echo 'gmx grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -maxwarn 1000 -r nvt.gro' >> queue.sh")
                os.system("echo 'gmx mdrun -v -deffnm npt' >> queue.sh")
                os.system("echo 'gmx grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -maxwarn 1000 -r npt.gro' >> queue.sh")
              else:
                pass

            else:
              mbox.showerror('Error', 'ACPYPE error in Cofactor. Process has been cancelled')
              quit()
          
          else:
            mbox.showerror('Error', 'ACPYPE error in Ligand. Process has been cancelled')
            quit()
        
        else:
          mbox.showerror('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

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
        ig = self.hb.getvar('var1')
        if gromacs_flag('mdrun'):
          cmd = 'pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        elif gromacs_flag('gmx'):
          cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)  
        else:
          pass
        out_pdb2gmx = os.system(cmd)
        if out_pdb2gmx == 0:
          bx = str(self.bx_menu.getcurselection())
          dst = str(self.dist.getvalue())
          if gromacs_flag('mdrun'):
            cmd1 = 'editconf -bt {0} -f trp.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
          elif gromacs_flag('gmx'):
            cmd1 = 'gmx editconf -bt {0} -f trp.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
          else:
            pass
          os.system(cmd1)
          if gromacs_flag('mdrun'):
            os.system('genbox -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
          elif gromacs_flag('gmx'):
            os.system('gmx solvate -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
          else:
            pass
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
          if gromacs_flag('mdrun'):
            os.system('grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')
          elif gromacs_flag('gmx'):
            os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000') 
          else:
            pass
          if self.ion_menu.getcurselection() == 'Na (Number)':
              c_ion = '-np ' + self.ion_cc.getvalue()
          elif self.ion_menu.getcurselection() == 'Cl (Number)':
              c_ion = '-nn ' + self.ion_cc.getvalue()
          else:
              c_ion = '-conc ' + self.ion_cc.getvalue()

          if gromacs_flag('mdrun'):
            cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
          elif gromacs_flag('gmx'):
            cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
          else:
            pass
          os.system(cmd2)

          inte = self.min_menu.getcurselection()

          nst = self.step.getvalue()

          if gromacs_flag('mdrun'):
            cmd3 = '''
              cat << EOF >| em_real.mdp
              ; LINES STARTING WITH ';' ARE COMMENTS
              title   = Minimization  ; Title of run

              ; Parameters describing what to do, when to stop and what to save
              integrator  = steep   ; Algorithm (steep = steepest descent minimization)
              emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
              emstep      = 0.01      ; Energy step size
              nsteps    = {0}     ; Maximum number of (minimization) steps to perform
              energygrps  = Protein   ; Which energy group(s) to write to disk

              ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
              nstlist   = 1       ; Frequency to update the neighbor list and long range forces
              ns_type   = grid    ; Method to determine neighbor list (simple, grid)
              rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
              coulombtype = PME   ; Treatment of long range electrostatic interactions
              rcoulomb  = 1.0   ; long range electrostatic cut-off
              rvdw    = 1.0   ; long range Van der Waals cut-off
              pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
              '''.format(nst)
            cmd4 = '''
                        cat << EOF >| em_real.mdp
                        ; LINES STARTING WITH ';' ARE COMMENTS
                        title   = Minimization  ; Title of run

                        ; Parameters describing what to do, when to stop and what to save
                        integrator  = cg    ; Algorithm (steep = steepest descent minimization)
                        emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
                        emstep      = 0.01      ; Energy step size
                        nsteps    = {0}     ; Maximum number of (minimization) steps to perform
                        energygrps  = Protein   ; Which energy group(s) to write to disk

                        ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                        nstlist   = 1       ; Frequency to update the neighbor list and long range forces
                        ns_type   = grid    ; Method to determine neighbor list (simple, grid)
                        rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
                        coulombtype = PME   ; Treatment of long range electrostatic interactions
                        rcoulomb  = 1.0   ; long range electrostatic cut-off
                        rvdw    = 1.0   ; long range Van der Waals cut-off
                        pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
                        '''.format(nst)

            if inte == 'SD Algorithm':
                os.system(cmd3)
                os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm em')
                os.system('''
                trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro << EOF
                Protein
                System
                EOF ''')
                os.system('''
                mv em_2.gro em.gro ''')
                os.system('''
                trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb << EOF
                Protein
                System
                EOF
                ''')
                os.system('''
                        mv em_2.pdb em.pdb ''')

            else:
                os.system(cmd3)
                os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm em')
                os.system(cmd4)
                os.system('grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm em')
                os.system('''
                            trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro << EOF
                            Protein
                            System
                            EOF ''')
                os.system('''
                            mv em_2.gro em.gro ''')
                os.system('''
                            trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb << EOF
                            Protein
                            System
                            EOF
                            ''')
                os.system('mv em_2.pdb em.pdb')

          elif gromacs_flag('gmx'):

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
                os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
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
                os.system('gmx grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                os.system('gmx mdrun -v -deffnm em')
                os.system(cmd4)
                os.system('gmx grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr -maxwarn 1000')
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
          else:
            pass

          stnvt = str(self.time_nvt.getvalue())
          stnpt = str(self.time_npt.getvalue())
          stmd = str(self.time_md.getvalue())
          temp = str(self.temperature.getvalue())
          t_st = str(self.time_st.getvalue())

          if gromacs_flag('mdrun'):
            cmd5 = '''
            cat << EOF >| nvt.mdp
            title       = Protein-ligand complex NVT equilibration
            define      = -DPOSRES  ; position restrain the protein and ligand
            ; Run parameters
            integrator  = md        ; leap-frog integrator
            nsteps      = {0}     ; 2 * 50000 = 100 ps
            dt          = {1}     ; 2 fs
            ; Output control
            nstxout     = 100       ; save coordinates every 0.2 ps
            nstvout     = 100       ; save velocities every 0.2 ps
            nstenergy   = 100       ; save energies every 0.2 ps
            nstlog      = 100       ; update log file every 0.2 ps
            energygrps  = System
            ; Bond parameters
            continuation    = no            ; first dynamics run
            constraint_algorithm = lincs    ; holonomic constraints
            constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
            lincs_iter      = 1             ; accuracy of LINCS
            lincs_order     = 4             ; also related to accuracy
            ; Neighborsearching
            ns_type     = grid      ; search neighboring grid cells
            nstlist     = 5         ; 10 fs
            rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
            rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
            rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
            ; Electrostatics
            coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
            pme_order       = 4         ; cubic interpolation
            fourierspacing  = 0.16      ; grid spacing for FFT
            ; Temperature coupling
            tcoupl      = V-rescale                     ; modified Berendsen thermostat
            tc-grps     = Protein Water_and_ions   ; two coupling groups - more accurate
            tau_t       = 0.1   0.1                    ; time constant, in ps
            ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
            ; Pressure coupling
            pcoupl      = no        ; no pressure coupling in NVT
            ; Periodic boundary conditions
            pbc         = xyz      ; 3-D PBC
            ; Dispersion correction
            DispCorr    = EnerPres  ; account for cut-off vdW scheme
            ; Velocity generation
            gen_vel     = yes       ; assign velocities from Maxwell distribution
            gen_temp    = {2}       ; temperature for Maxwell distribution
            gen_seed    = -1        ; generate a random seed'
            '''.format(stnvt, t_st, temp)

          elif gromacs_flag('gmx'):
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

          else:
            pass

          os.system(cmd5)

          

          if gromacs_flag('mdrun'):
            cmd6 = '''
            cat << EOF >| npt.mdp
            title       = Protein-ligand complex NPT equilibration
            define      = -DPOSRES  ; position restrain the protein and ligand
            ; Run parameters
            integrator  = md        ; leap-frog integrator
            nsteps      = {0}     ; 2 * 50000 = 100 ps
            dt          = {1}     ; 2 fs
            ; Output control
            nstxout     = 100       ; save coordinates every 0.2 ps
            nstvout     = 100       ; save velocities every 0.2 ps
            nstenergy   = 100       ; save energies every 0.2 ps
            nstlog      = 100       ; update log file every 0.2 ps
            energygrps  = System
            ; Bond parameters
            continuation    = yes           ; first dynamics run
            constraint_algorithm = lincs    ; holonomic constraints
            constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
            lincs_iter      = 1             ; accuracy of LINCS
            lincs_order     = 4             ; also related to accuracy
            ; Neighborsearching
            ns_type     = grid      ; search neighboring grid cells
            nstlist     = 5         ; 10 fs
            rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
            rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
            rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
            ; Electrostatics
            coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
            pme_order       = 4         ; cubic interpolation
            fourierspacing  = 0.16      ; grid spacing for FFT
            ; Temperature coupling
            tcoupl      = V-rescale                     ; modified Berendsen thermostat
            tc-grps     = Protein Water_and_ions    ; two coupling groups - more accurate
            tau_t       = 0.1   0.1                    ; time constant, in ps
            ref_t       = {2}  {2}                    ; reference temperature, one for each group, in K
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
            '''.format(stnpt, t_st, temp)

          elif gromacs_flag('gmx'):
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
          else:
            pass

          os.system(cmd6)

          if gromacs_flag('mdrun'):
            cmd7 = '''
            cat << EOF >| md.mdp
            title       = Protein-ligand complex MD simulation
            ; Run parameters
            integrator  = md        ; leap-frog integrator
            nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
            dt          = {1}     ; 2 fs
            ; Output control
            nstxout                  = 0                     ; [steps] freq to write coordinates to trajectory
            nstvout                  = 0                 ; [steps] freq to write velocities to trajectory
            nstfout                  = 0                     ; [steps] freq to write forces to trajectory
            nstlog                   = 100                   ; [steps] freq to write energies to log file
            nstenergy                = 500                   ; [steps] freq to write energies to energy file
            nstxtcout                = 500                   ; [steps] freq to write coordinates to xtc trajectory
            xtc_precision            = 1000                  ; [real] precision to write xtc trajectory
            xtc_grps                 = System                ; group(s) to write to xtc trajectory
            ; Bond parameters
            continuation    = yes           ; first dynamics run
            constraint_algorithm = lincs    ; holonomic constraints
            constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
            lincs_iter      = 1             ; accuracy of LINCS
            lincs_order     = 4             ; also related to accuracy
            ; Neighborsearching
            ns_type     = grid      ; search neighboring grid cells
            nstlist     = 5         ; 10 fs
            rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
            rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
            rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
            ; Electrostatics
            coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
            pme_order       = 4         ; cubic interpolation
            fourierspacing  = 0.16      ; grid spacing for FFT
            ; Temperature coupling
            tcoupl      = V-rescale                     ; modified Berendsen thermostat
            tc-grps     = Protein Water_and_ions    ; two coupling groups - more accurate
            tau_t       = 0.1   0.1                    ; time constant, in ps
            ref_t       = {2}   {2}                   ; reference temperature, one for each group, in K
            ; Pressure coupling
            pcoupl      = Parrinello-Rahman             ; pressure coupling is on for NPT
            pcoupltype  = isotropic                     ; uniform scaling of box vectors
            tau_p       = 2.0                           ; time constant, in ps
            ref_p       = 1.0                           ; reference pressure, in bar
            compressibility = 4.5e-5                    ; isothermal compressibility of water, bar^-1
            ; Periodic boundary conditions
            pbc         = xyz      ; 3-D PBC
            ; Dispersion correction
            DispCorr    = EnerPres  ; account for cut-off vdW scheme
            ; Velocity generation
            gen_vel     = no        ; assign velocities from Maxwell distribution
            '''.format(stmd, t_st, temp)

          elif gromacs_flag('gmx'):

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
          else:
            pass

          os.system(cmd7)

          cmd8='pymol em.pdb'

          os.system(cmd8)

          if mbox.askyesno('View Complex', 'Is protein OK??'):
              pass
          else:
              mbox.showinfo('No', 'Process has been cancelled')
              quit()

          os.system('cat << EOF > | queue.sh')

          if gromacs_flag('mdrun'):
            os.system("echo 'grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -maxwarn 1000' >> queue.sh")
            os.system("echo  'mdrun -v -deffnm nvt' >> queue.sh")
            os.system("echo 'grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -maxwarn 1000' >> queue.sh")
            os.system("echo 'mdrun -v -deffnm npt' >> queue.sh")
            os.system("echo 'grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -maxwarn 1000' >> queue.sh")

          elif gromacs_flag('gmx'):
            os.system("echo 'gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -r em.gro -maxwarn 1000' >> queue.sh")
            os.system("echo 'gmx mdrun -v -deffnm nvt' >> queue.sh")
            os.system("echo 'gmx grompp -f npt.mdp -c nvt.gro -p trp.top -o npt.tpr -r nvt.gro -maxwarn 1000' >> queue.sh")
            os.system("echo 'gmx mdrun -v -deffnm npt' >> queue.sh")
            os.system("echo 'gmx grompp -f md.mdp -c npt.gro -p trp.top -o md.tpr -r npt.gro -maxwarn 1000' >> queue.sh")
          else:
            pass
        else:
          mbox.showerror('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

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

        elif find3 == False:
            mbox.showerror("Error", "PDB file not found.")
            pass

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
        ig = self.hb.getvar('var1')
        if gromacs_flag('mdrun'):
          cmd = 'pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        elif gromacs_flag('gmx'):
          cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        else:
          pass
        out_pdb2gmx = os.system(cmd)
        if out_pdb2gmx == 0:
          try:
            pmol = PandasMol2().read_mol2('Ligand.mol2')
            subst = pmol.df.iloc[0]['subst_name']
            fh, abs_path = mkstemp()
            with fdopen(fh,'w') as new_file:
              with open('Ligand.mol2') as old_file:
                for line in old_file:
                  new_file.write(line.replace(subst, 'LIG'))
            #Remove original file
            remove('Ligand.mol2')
            #Move new file
            move(abs_path, 'Ligand.mol2')
          except:
            mbox.showerror("Error", "Mol2 file not recognized.. Please try again.")
            quit()
          
          cm = str(self.cm_menu.getcurselection())
          nc = str(self.nc.getvalue())
          mt = str(self.mult.getvalue())
          at = str(self.at_menu.getcurselection())
          cmdd = 'acpype -i Ligand.mol2 -c {0} -n {1} -m {2} -a {3}'.format(cm, nc, mt, at)
          acp2 = os.system(cmdd)
          if acp2 ==0:
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

            if gromacs_flag('mdrun'):
              cmd1 = 'gmx editconf -bt {0} -f complex.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
            elif gromacs_flag('gmx'):
              cmd1 = 'gmx editconf -bt {0} -f complex.pdb -o trpb4solv.pdb -d {1}'.format(bx, dst)
            else:
              pass
            os.system(cmd1)
            if gromacs_flag('mdrun'):
              os.system('genbox -cp trpb4solv.pdb -cs spc216.gro -o trpb4ion.pdb -p trp.top')
              os.system('''
                          cat << EOF >| em.mdp
                          ; LINES STARTING WITH ';' ARE COMMENTS
                          title   = Minimization  ; Title of run

                          ; Parameters describing what to do, when to stop and what to save
                          integrator  = steep   ; Algorithm (steep = steepest descent minimization)
                          emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
                          emstep      = 0.01      ; Energy step size
                          nsteps    = 50000     ; Maximum number of (minimization) steps to perform
                          energygrps  = system  ; Which energy group(s) to write to disk

                          ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
                          nstlist       = 1       ; Frequency to update the neighbor list and long range forces
                          cutoff-scheme   = Verlet
                          ns_type       = grid    ; Method to determine neighbor list (simple, grid)
                          rlist       = 1.0   ; Cut-off for making neighbor list (short range forces)
                          coulombtype     = PME   ; Treatment of long range electrostatic interactions
                          rcoulomb      = 1.0   ; long range electrostatic cut-off
                          rvdw        = 1.0   ; long range Van der Waals cut-off
                          pbc             = xyz     ; Periodic Boundary Conditions
                          ''')
              os.system('grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')

            elif gromacs_flag('gmx'):
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
            else:
              pass

            if self.ion_menu.getcurselection() == 'Na (Number)':
                c_ion = '-np ' + self.ion_cc.getvalue()
            elif self.ion_menu.getcurselection() == 'Cl (Number)':
                c_ion = '-nn ' + self.ion_cc.getvalue()
            else:
                c_ion = '-conc ' + self.ion_cc.getvalue()

            if gromacs_flag('mdrun'):
              cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
            elif gromacs_flag('gmx'):
              cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
            
            else:
              pass

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
            if gromacs_flag('mdrun'): 
              '''
               cat << EOF >| em_real.mdp
               ; LINES STARTING WITH ';' ARE COMMENTS
               title   = Minimization  ; Title of run

               ; Parameters describing what to do, when to stop and what to save
               integrator  = cg    ; Algorithm (steep = steepest descent minimization)
               emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
               emstep      = 0.01      ; Energy step size
               nsteps    = {0}     ; Maximum number of (minimization) steps to perform
               energygrps  = Protein LIG ; Which energy group(s) to write to disk
               ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
               nstlist   = 1       ; Frequency to update the neighbor list and long range forces
               ns_type   = grid    ; Method to determine neighbor list (simple, grid)
               rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
               coulombtype = PME   ; Treatment of long range electrostatic interactions
               rcoulomb  = 1.0   ; long range electrostatic cut-off
               rvdw    = 1.0   ; long range Van der Waals cut-off
               pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
               '''.format(nst)

              if inte == 'SD Algorithm':
                os.system(cmd3)
                os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm em')
                os.system('''
                make_ndx -f em.gro -o index2.ndx << EOF
                "Protein" | "Other"
                q
                EOF ''')
                os.system('''
                trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                Protein_Other
                System
                EOF ''')
                os.system('''
                mv em_2.gro em.gro ''')
                os.system('''
                trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                Protein_Other
                System
                EOF
                ''')
                
              else:
                  os.system(cmd3)
                  os.system('grompp -f em_real.mdp -c trpb4em.pdb -p trp.top -o em.tpr -maxwarn 1000')
                  os.system('mdrun -v -deffnm em')
                  os.system(cmd4)
                  os.system('grompp -f em_real.mdp -c em.gro -p trp.top -o em.tpr -maxwarn 1000')
                  os.system('mdrun -v -deffnm em')
                  os.system('''
                              make_ndx -f em.gro -o index2.ndx << EOF
                              "Protein" | "Other"
                              q
                              EOF ''')
                  os.system('''
                              trjconv -f em.gro -s em.tpr -pbc nojump -ur compact -center -o em_2.gro -n index2.ndx << EOF
                              Protein_Other
                              System
                              EOF ''')
                  os.system('''
                              mv em_2.gro em.gro ''')
                  os.system('''
                              trjconv -f em.gro -s em.tpr -pbc mol -ur compact -center -o em_2.pdb -n index2.ndx << EOF
                              Protein_Other
                              System
                              EOF
                              ''')
            elif gromacs_flag('gmx'):
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
            else:
              pass
                
            os.system('mv em_2.pdb em.pdb')


            stnvt = str(self.time_nvt.getvalue())
            temp = str(self.temperature.getvalue())
            t_st = str(self.time_st.getvalue())

            

            if gromacs_flag('mdrun'):
              cmd5 = '''
              cat << EOF >| nvt.mdp
              title       = Protein-ligand complex NVT equilibration
              define      = -DPOSRES  ; position restrain the protein and ligand
              ; Run parameters
              integrator  = md        ; leap-frog integrator
              nsteps      = {0}     ; 2 * 50000 = 100 ps
              dt          = {1}     ; 2 fs
              ; Output control
              nstxout     = 100       ; save coordinates every 0.2 ps
              nstvout     = 100       ; save velocities every 0.2 ps
              nstenergy   = 100       ; save energies every 0.2 ps
              nstlog      = 100       ; update log file every 0.2 ps
              energygrps  = Protein LIG
              ; Bond parameters
              continuation    = no            ; first dynamics run
              constraint_algorithm = lincs    ; holonomic constraints
              constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
              lincs_iter      = 1             ; accuracy of LINCS
              lincs_order     = 4             ; also related to accuracy
              ; Neighborsearching
              ns_type     = grid      ; search neighboring grid cells
              nstlist     = 5         ; 10 fs
              rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
              rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
              rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
              vdwtype             =  Cut-off
              ; Electrostatics
              coulombtype     = Cut-off       ; Particle Mesh Ewald for long-range electrostatics
              pme_order       = 4         ; cubic interpolation
              fourierspacing  = 0.16      ; grid spacing for FFT
              ; Temperature coupling
              tcoupl      = V-rescale                     ; modified Berendsen thermostat
              tc-grps     = Protein non-Protein   ; two coupling groups - more accurate
              tau_t       = 0.1   0.1                    ; time constant, in ps
              ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
              ; Pressure coupling
              pcoupl      = no        ; no pressure coupling in NVT
              ; Periodic boundary conditions
              pbc         = xyz      ; 3-D PBC
              ; Dispersion correction
              DispCorr    = EnerPres  ; account for cut-off vdW scheme
              ; Velocity generation
              gen_vel     = yes       ; assign velocities from Maxwell distribution
              gen_temp    = {2}       ; temperature for Maxwell distribution
              gen_seed    = -1        ; generate a random seed'
              '''.format(stnvt, t_st, temp)


            elif gromacs_flag('gmx'):
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

            else:
              pass

            os.system(cmd5)

            cmd8='pymol em.pdb'

            os.system(cmd8)

            if mbox.askyesno('View Complex', 'Is complex OK??'):
                pass
            else:
                mbox.showinfo('No', 'Process has been cancelled')
                quit()

            if gromacs_flag('mdrun'):
              os.system('echo 2|genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
              os.system(r'''
              sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
              ''')
              os.system('mv trp2.top trp.top')

              os.system('cat << EOF > | queue.sh')

              os.system('grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -maxwarn 1000')
              os.system('mdrun -v -deffnm nvt')


            elif gromacs_flag('gmx'):    
              os.system('echo 2|gmx genrestr -f Ligand.acpype/Ligand_GMX.gro -o posre_LIG.itp -fc 1000 1000 1000')
              os.system(r'''
              sed '/posre.itp/{p;s/.*/#endif \n\n; Ligand position restraints \n#ifdef POSRES \n#include "posre_LIG.itp"/;}' trp.top > trp2.top
              ''')
              os.system('mv trp2.top trp.top')

              os.system('cat << EOF > | queue.sh')

              os.system('gmx grompp -f nvt.mdp -c em.pdb -p trp.top -o nvt.tpr -r em.gro -maxwarn 1000')
              os.system('gmx mdrun -v -deffnm nvt')
            else:
              pass
          else:
            mbox.showerror('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            pass
        else:
          mbox.showerror('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

    def simulation_LIE_lig(self):
        try:
            os.chdir(path)
        except:
            pass
        try:
          pmol = PandasMol2().read_mol2('Ligand.mol2')
          subst = pmol.df.iloc[0]['subst_name']
          fh, abs_path = mkstemp()
          with fdopen(fh,'w') as new_file:
            with open('Ligand.mol2') as old_file:
              for line in old_file:
                new_file.write(line.replace(subst, 'LIG'))
          #Remove original file
          remove('Ligand.mol2')
          #Move new file
          move(abs_path, 'Ligand.mol2')
        except:
          mbox.showerror("Error", "Mol2 file not recognized.. Please try again.")
          quit()
        
        cm = str(self.cm_menu.getcurselection())
        nc = str(self.nc.getvalue())
        mt = str(self.mult.getvalue())
        at = str(self.at_menu.getcurselection())
        cmdd = 'acpype -i Ligand.mol2 -c {0} -n {1} -m {2} -a {3}'.format(cm, nc, mt, at)
        acp3 = os.system(cmdd)
        if acp3 == 0:
          os.system('cp Ligand.acpype/Ligand_GMX.gro Ligand.gro')

          ff = str(self.ff_menu.getcurselection())
          wt = str(self.wt_menu.getcurselection())

          if ff != 'oplsaa':
              os.system('cp Ligand.acpype/Ligand_GMX.itp Ligand.itp')
          else:
              os.system('cp Ligand.acpype/Ligand_GMX_OPLS.itp Ligand.itp')

          liecmd5 = '''
          cat << EOF >| topoll.top
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
          
          if gromacs_flag('mdrun'):
            liecmd1 = 'editconf -bt {0} -f Ligand.gro -o trpb4solvl.pdb -d {1}'.format(bx, dst)
            os.system(liecmd1)
            os.system('genbox -cp trpb4solvl.pdb -cs spc216.gro -o trpb4ionl.pdb -p topoll.top')

          elif gromacs_flag('gmx'):
            liecmd1 = 'gmx editconf -bt {0} -f Ligand.gro -o trpb4solvl.pdb -d {1}'.format(bx, dst)
            os.system(liecmd1)
            os.system('gmx solvate -cp trpb4solvl.pdb -cs spc216.gro -o trpb4ionl.pdb -p topoll.top')
          else:
            pass
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
          if gromacs_flag('mdrun'):
            os.system('grompp -f em.mdp -c trpb4ionl.pdb -p topoll.top -o ionl.tpr -maxwarn 1000')
          elif gromacs_flag('gmx'):
            os.system('gmx grompp -f em.mdp -c trpb4ionl.pdb -p topoll.top -o ionl.tpr -maxwarn 1000')
          else:
            pass  
          if self.ion_menu.getcurselection() == 'Na (Number)':
              c_ion = '-np ' + self.ion_cc.getvalue()
          elif self.ion_menu.getcurselection() == 'Cl (Number)':
              c_ion = '-nn ' + self.ion_cc.getvalue()
          else:
              c_ion = '-conc ' + self.ion_cc.getvalue()
          if gromacs_flag('mdrun'):
            liecmd2 = 'echo SOL|genion -s ionl.tpr -o trpb4eml.pdb -neutral {} -p topoll.top'.format(c_ion)
          elif gromacs_flag('gmx'):
            liecmd2 = 'echo SOL|gmx genion -s ionl.tpr -o trpb4eml.pdb -neutral {} -p topoll.top'.format(c_ion)
          else:
            pass  
          os.system(liecmd2)

          lie_inte = self.min_menu.getcurselection()

          lie_nst = self.step.getvalue()

          
          if gromacs_flag('mdrun'):
            liecmd3 = '''
            cat << EOF >| em_real.mdp
; LINES STARTING WITH ';' ARE COMMENTS
title   = Minimization  ; Title of run

; Parameters describing what to do, when to stop and what to save
integrator  = steep   ; Algorithm (steep = steepest descent minimization)
emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
emstep      = 0.01      ; Energy step size
nsteps    = {0}     ; Maximum number of (minimization) steps to perform
energygrps  = LIG ; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist   = 1       ; Frequency to update the neighbor list and long range forces
ns_type   = grid    ; Method to determine neighbor list (simple, grid)
rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
coulombtype = PME   ; Treatment of long range electrostatic interactions
rcoulomb  = 1.0   ; long range electrostatic cut-off
rvdw    = 1.0   ; long range Van der Waals cut-off
pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
            '''.format(lie_nst)
            liecmd4 = '''
cat << EOF >| em_real.mdp
; LINES STARTING WITH ';' ARE COMMENTS
title   = Minimization  ; Title of run

; Parameters describing what to do, when to stop and what to save
integrator  = cg    ; Algorithm (steep = steepest descent minimization)
emtol   = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
emstep      = 0.01      ; Energy step size
nsteps    = {0}     ; Maximum number of (minimization) steps to perform
energygrps  = LIG ; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist   = 1       ; Frequency to update the neighbor list and long range forces
ns_type   = grid    ; Method to determine neighbor list (simple, grid)
rlist   = 1.0   ; Cut-off for making neighbor list (short range forces)
coulombtype = PME   ; Treatment of long range electrostatic interactions
rcoulomb  = 1.0   ; long range electrostatic cut-off
rvdw    = 1.0   ; long range Van der Waals cut-off
pbc       = xyz     ; Periodic Boundary Conditions (yes/no)
                    '''.format(lie_nst)

            if lie_inte == 'SD Algorithm':
                os.system(liecmd3)
                os.system('grompp -f em_real.mdp -c trpb4eml.pdb -p topoll.top -o eml.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm eml')
                os.system('''
                trjconv -f eml.gro -s eml.tpr -pbc nojump -ur compact -center -o em_2l.gro << EOF
                Other
                System
                EOF ''')
                os.system('''
                mv em_2l.gro eml.gro ''')
                os.system('''
                trjconv -f eml.gro -s eml.tpr -pbc mol -ur compact -center -o em_2l.pdb << EOF
                Other
                System
                EOF
                ''')
                
            else:
                os.system(liecmd3)
                os.system('grompp -f em_real.mdp -c trpb4eml.pdb -p topoll.top -o eml.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm eml')
                os.system(liecmd4)
                os.system('grompp -f em_real.mdp -c eml.gro -p topoll.top -o eml.tpr -maxwarn 1000')
                os.system('mdrun -v -deffnm eml')
                os.system('''
                            trjconv -f eml.gro -s eml.tpr -pbc nojump -ur compact -center -o em_2l.gro << EOF
                            Other
                            System
                            EOF ''')
                os.system('''
                            mv em_2l.gro eml.gro ''')
                os.system('''
                            trjconv -f eml.gro -s eml.tpr -pbc mol -ur compact -center -o em_2l.pdb << EOF
                            Other
                            System
                            EOF
                            ''')

          elif gromacs_flag('gmx'):      

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
                  os.system('gmx grompp -f em_real.mdp -c trpb4eml.pdb -p topoll.top -o eml.tpr -maxwarn 1000')
                  os.system('gmx mdrun -v -deffnm eml')
                  os.system('''
                  gmx trjconv -f eml.gro -s eml.tpr -pbc nojump -ur compact -center -o em_2l.gro << EOF
                  Other
                  System
                  EOF ''')
                  os.system('''
                  mv em_2l.gro eml.gro ''')
                  os.system('''
                  gmx trjconv -f eml.gro -s eml.tpr -pbc mol -ur compact -center -o em_2l.pdb << EOF
                  Other
                  System
                  EOF
                  ''')
                  
              else:
                  os.system(liecmd3)
                  os.system('gmx grompp -f em_real.mdp -c trpb4eml.pdb -p topoll.top -o eml.tpr -maxwarn 1000')
                  os.system('gmx mdrun -v -deffnm eml')
                  os.system(liecmd4)
                  os.system('gmx grompp -f em_real.mdp -c eml.gro -p topoll.top -o eml.tpr -maxwarn 1000')
                  os.system('gmx mdrun -v -deffnm eml')
                  os.system('''
                              gmx trjconv -f eml.gro -s eml.tpr -pbc nojump -ur compact -center -o em_2l.gro << EOF
                              Other
                              System
                              EOF ''')
                  os.system('''
                              mv em_2l.gro eml.gro ''')
                  os.system('''
                              gmx trjconv -f eml.gro -s eml.tpr -pbc mol -ur compact -center -o em_2l.pdb << EOF
                              Other
                              System
                              EOF
                              ''')
          
          else:
            pass
          
          os.system('mv em_2l.pdb eml.pdb')


          lie_stnvt = str(self.time_nvt.getvalue())
          lie_t_st = str(self.time_st.getvalue())
          lie_temp = str(self.temperature.getvalue())

          
          if gromacs_flag('mdrun'):
            liecmd5 = '''
            cat << EOF >| nvt2.mdp
            title       = Ligand NVT equilibration
            define      = -DPOSRES  ; position restrain the protein and ligand
            ; Run parameters
            integrator  = md        ; leap-frog integrator
            nsteps      = {0}     ; 2 * 50000 = 100 ps
            dt          = {1}     ; 2 fs
            ; Output control
            nstxout     = 100       ; save coordinates every 0.2 ps
            nstvout     = 100       ; save velocities every 0.2 ps
            nstenergy   = 100       ; save energies every 0.2 ps
            nstlog      = 100       ; update log file every 0.2 ps
            energygrps  = LIG
            ; Bond parameters
            continuation    = no            ; first dynamics run
            constraint_algorithm = lincs    ; holonomic constraints
            constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
            lincs_iter      = 1             ; accuracy of LINCS
            lincs_order     = 4             ; also related to accuracy
            ; Neighborsearching
            ns_type     = grid      ; search neighboring grid cells
            nstlist     = 5         ; 10 fs
            rlist       = 0.9       ; short-range neighborlist cutoff (in nm)
            rcoulomb    = 0.9       ; short-range electrostatic cutoff (in nm)
            rvdw        = 1.4       ; short-range van der Waals cutoff (in nm)
            vdwtype             =  Cut-off
            ; Electrostatics
            coulombtype     = Cut-off   ;
            pme_order       = 4         ; cubic interpolation
            fourierspacing  = 0.16      ; grid spacing for FFT
            ; Temperature coupling
            tcoupl      = V-rescale                     ; modified Berendsen thermostat
            tc-grps     = LIG   Water_Ion   ; two coupling groups - more accurate
            tau_t       = 0.1   0.1                    ; time constant, in ps
            ref_t       = {2}   {2}                    ; reference temperature, one for each group, in K
            ; Pressure coupling
            pcoupl      = no        ; no pressure coupling in NVT
            ; Periodic boundary conditions
            pbc         = xyz      ; 3-D PBC
            ; Dispersion correction
            DispCorr    = EnerPres  ; account for cut-off vdW scheme
            ; Velocity generation
            gen_vel     = yes       ; assign velocities from Maxwell distribution
            gen_temp    = {2}       ; temperature for Maxwell distribution
            gen_seed    = -1        ; generate a random seed'
            '''.format(lie_stnvt, lie_t_st, lie_temp) 
          elif gromacs_flag('gmx'):  
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
          else:
            pass  
          os.system(liecmd5)

          liecmd8='pymol eml.pdb'

          os.system(liecmd8)

          if mbox.askyesno('View Complex', 'Is Ligand OK??'):
              pass
          else:
              mbox.showinfo('No', 'Process has been cancelled')
              pass

          if gromacs_flag('mdrun'):
            os.system('''
                        make_ndx -f eml.gro -o index2l.ndx << EOF
                        "Water" | "Ion"
                        q
                        EOF ''')
            os.system('grompp -f nvt2.mdp -c eml.pdb -p topoll.top -o nvt2l.tpr -n index2l.ndx -maxwarn 1000')
            os.system('mdrun -v -deffnm nvt2l')    
          
          elif gromacs_flag('gmx'):

            os.system('''
                            gmx make_ndx -f eml.gro -o index2l.ndx << EOF
                            "Water" | "Ion"
                            q
                            EOF ''')
            os.system('gmx grompp -f nvt2.mdp -c eml.pdb -p topoll.top -o nvt2l.tpr -n index2l.ndx -r eml.gro -maxwarn 1000')
            os.system('gmx mdrun -v -deffnm nvt2l')
          else:
            pass  
        else:
          mbox.showerror('Error', 'ACPYPE error (Ligand). Process has been cancelled')
          quit()

    def lie_calculation(self):
        try:
            os.chdir(path)
        except:
            pass
        if gromacs_flag('mdrun'):
          com1= '''echo 43 44 0 | g_energy -f nvt.edr > out1.txt'''
        elif gromacs_flag('gmx'):
          com1= '''echo 43 44 0 | gmx energy -f nvt.edr > out1.txt'''
        else:
          pass
        os.system(com1)
        if gromacs_flag('mdrun'):
          com2= '''echo 41 42 0 | g_energy -f nvt2l.edr > out2.txt'''
        elif gromacs_flag('gmx'):
          com2= '''echo 41 42 0 | gmx energy -f nvt2l.edr > out2.txt'''
        else:
          pass
        os.system(com2)
        try:
          with open('out1.txt', 'r') as f:
            lines = f.readlines()
            for line in lines:
              if re.search(r'Coul-14:Protein-LIG', line):
                match0 = re.split(r'\s+', line)

              elif re.search(r'LJ-14:Protein-LIG', line):
                match1 = re.split(r'\s+', line)
                pass
        except UnboundLocalError:
          mbox.showerror('Error', 'Please try to use other simulation type box (Dodecahedron).')
          quit()


        coul = match0[1]
        lj = match1[1]
        
        with open('out1.txt', 'r') as f:
          lines = f.readlines()
          for line in lines:
            if re.search(r'Coul-14:Protein-LIG', line):
              match02 = re.split(r'\s+', line)

            elif re.search(r'LJ-14:Protein-LIG', line):
              match12 = re.split(r'\s+', line)
              pass

        coul2 = match02[1]
        lj2 = match12[1]

        if gromacs_flag('mdrun'):
          os.system('g_lie -f nvt.edr -o lie_comp.xvg -ligand LIG -Eqq {0} -Elj {1} > out_1.txt'.format(coul,lj))
          os.system('g_lie -f nvt2l.edr -o lie_lig.xvg -ligand LIG -Eqq {0} -Elj {1} >> out_2.txt'.format(coul2,lj2))
        elif gromacs_flag('gmx'):
          os.system('gmx lie -f nvt.edr -o lie_comp.xvg -ligand LIG -Eqq {0} -Elj {1} > out_1.txt'.format(coul,lj))
          os.system('gmx lie -f nvt2l.edr -o lie_lig.xvg -ligand LIG -Eqq {0} -Elj {1} >> out_2.txt'.format(coul2,lj2))
        else:
          pass
        with open('out_1.txt', 'r') as f:
          lines = f.readlines()
          for line in lines:
            if re.search(r'DGbind =', line):
              match001 = re.split(r'\s+', line)

        with open('out_2.txt', 'r') as f:
          lines = f.readlines()
          for line in lines:
            if re.search(r'DGbind =', line):
              match002 = re.split(r'\s+', line)
        a = float(match001[2])
        b = float(match002[2])
        c =float(a+b)
        y = format(c, '.2f')
        self.lie_ener['text'] = str(y)


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
        return None

    def run_simulation(self):

      if gromacs_flag('mdrun'):
        os.system("echo 'mdrun -v -deffnm md -cpt 1' >> queue.sh")
      
      elif gromacs_flag('gmx'):
        os.system("echo 'gmx mdrun -v -deffnm md -cpt 1' >> queue.sh")
      else:
        pass
      os.system('chmod 777 queue.sh')
      os.system('./queue.sh')
      pj = str(self.project.getvalue())
      
      md1 = '{0}.xtc'.format(pj)
      md2 = '{0}.edr'.format(pj)
      md3 = '{0}.gro'.format(pj)
      md4 = '{0}.tpr'.format(pj)
      try:
          os.chdir(path)
      except:
          pass
      find1=os.path.exists('Ligand.mol2')
     
      if gromacs_flag('mdrun'):
        if find1 == False:
          os.system('''
                                             trjconv -f md.gro -s md.tpr -pbc nojump -ur compact -center -o md_2.gro << EOF
                                             Protein
                                             System
                                             EOF ''')
          os.system('''
                                             mv md_2.gro md.gro ''')
          os.system('''
                                             trjconv -f md.gro -s md.tpr -pbc mol -ur compact -center -o md_2.gro << EOF
                                             Protein
                                             System
                                             EOF
                                             ''')
          os.system('mv md_2.gro md.gro')
          os.system('''
                                                     trjconv -f md.xtc -s md.tpr -pbc nojump -ur compact -center -o md_2.xtc << EOF
                                                     Protein
                                                     System
                                                     EOF ''')
          os.system('''
                                                     mv md_2.xtc md.xtc ''')
          os.system('''
                                                     trjconv -f md.xtc -s md.tpr -pbc mol -ur compact -center -o md_2.xtc << EOF
                                                     Protein
                                                     System
                                                     EOF
                                                     ''')
          os.system('mv md_2.xtc md.xtc')
        else:
          os.system('''
                                             make_ndx -f md.gro -o index3.ndx << EOF
                                             "Protein" | "Other"
                                             q
                                             EOF ''')
          os.system('''
                                             trjconv -f md.gro -s md.tpr -pbc nojump -ur compact -center -o md_2.gro -n index3.ndx << EOF
                                             Protein_Other
                                             System
                                             EOF ''')
          os.system('''
                                             mv md_2.gro md.gro ''')
          os.system('''
                                             trjconv -f md.gro -s md.tpr -pbc mol -ur compact -center -o md_2.gro -n index3.ndx << EOF
                                             Protein_Other
                                             System
                                             EOF
                                             ''')
          os.system('mv md_2.gro md.gro')
          os.system('''
                                                     trjconv -f md.xtc -s md.tpr -pbc nojump -ur compact -center -o md_2.xtc -n index3.ndx << EOF
                                                     Protein_Other
                                                     System
                                                     EOF ''')
          os.system('''
                                                     mv md_2.xtc md.xtc ''')
          os.system('''
                                                     trjconv -f md.xtc -s md.tpr -pbc mol -ur compact -center -o md_2.xtc -n index3.ndx << EOF
                                                     Protein_Other
                                                     System
                                                     EOF
                                                     ''')
      elif gromacs_flag('gmx'):
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
      else:
        pass
      os.system('mv md_2.xtc md.xtc')
      shutil.copy('md.xtc', md1)
      shutil.copy('md.edr', md2)
      shutil.copy('md.gro', md3)
      shutil.copy('md.tpr', md4)
      dr = str(self.save.getvalue())
      results_path = dr + '/'+ pj + '_results'
      if not os.path.exists(results_path):
        os.makedirs(results_path)
      else:
        pass
      try:
        shutil.copy2(md1, results_path)
        shutil.copy2(md2, results_path)
        shutil.copy2(md3, results_path)
        shutil.copy2(md4, results_path)
        mbox.showinfo('Finish', 'Job has finished')
        pass
      except:
        mbox.showerror("Error", " It is not possible copy MD files.")
        quit()

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
          print('Please wait load file process is finished...')
          
          self.xtc['text'] = self.xtcfile
          
          shutil.copy(self.xtcfile, path+'/md_an.xtc')
          print('OK')
                    
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

            if gromacs_flag('mdrun'):
              command= ''' g_rms -s md_an.tpr -f md_an.xtc -o rmsd.xvg -xvg 
               << EOF
                                    {0}
                                    {0}
                                    EOF '''.format(struct)
            elif gromacs_flag('gmx'):    
              command= ''' gmx rms -s md_an.tpr -f md_an.xtc -o rmsd.xvg -xvg none << EOF
                                    {0}
                                    {0}
                                    EOF '''.format(struct)
            else:
              pass                        
            os.system(command)
            pj = str(self.project.getvalue())
            an1 = 'rmsd_{0}.xvg'.format(pj)
            shutil.copy('rmsd.xvg', an1)
            directory = str(self.save.getvalue())
            shutil.copy2(an1, directory)
            t, rmsd = np.loadtxt("rmsd.xvg", unpack=True)
            fig = plt.figure(figsize=(10, 5))
            data = np.loadtxt('rmsd.xvg')
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
            text0 = """
LiGRO v 0.03 - Output of {0}
---------------------------------------------
            """.format(analysis)
            try:
              with open('rmsd_out.txt', 'w') as rmsd_infile:
                rmsd_out = tkFileDialog.asksaveasfile(mode='w', initialdir=path2, filetypes =(("Text File", "*.txt"),("All Files","*.*")),
    title = "Save RMSD statistic file.", initialfile='ligro_RMSD')
                rmsd_out.write(text0+data0)
                rmsd_out.close()
                pass
            except:
              pass

        elif analysis == 'RMSF':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                pass
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                pass
            else:
                pass
            if gromacs_flag('mdrun'):    
              command1 = ''' g_rmsf -s md_an.tpr -f md_an.xtc -o rmsf.xvg -xvg none -res << EOF
                                                {0}
                                                EOF '''.format(struct)
            elif gromacs_flag('gmx'):    
              command1 = ''' gmx rmsf -s md_an.tpr -f md_an.xtc -o rmsf.xvg -xvg none -res << EOF
                                                {0}
                                                EOF '''.format(struct)
            else:
              pass
            os.system(command1)
            pj = str(self.project.getvalue())
            an2 = 'rmsf_{0}.xvg'.format(pj)
            shutil.copy('rmsf.xvg', an2)
            directory = str(self.save.getvalue())
            shutil.copy2(an2, directory)
            t, rmsf = np.loadtxt("rmsf.xvg", unpack=True)
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
            pass
        elif analysis == 'RG':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                pass
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                pass
            else:
                pass
            if gromacs_flag('mdrun'):    
              command2 = ''' g_gyrate -s md_an.tpr -f md_an.xtc -o gyrate.xvg -xvg none << EOF
                                                {0}
                                                EOF '''.format(struct)    
            elif gromacs_flag('gmx'):    
              command2 = ''' gmx gyrate -s md_an.tpr -f md_an.xtc -o gyrate.xvg -xvg none << EOF
                                                {0}
                                                EOF '''.format(struct)
            else:
              pass
            os.system(command2)
            pj = str(self.project.getvalue())
            an3 = 'rg_{0}.xvg'.format(pj)
            shutil.copy('gyrate.xvg', an3)
            directory = str(self.save.getvalue())
            shutil.copy2(an3, directory)
            data = np.loadtxt('gyrate.xvg')
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
            # plotting separately
            plt.plot([t, g1], label='Rg')
            plt.plot([t, g2], label='RgX')
            plt.plot([t, g3], label='RgY')
            plt.plot([t, g4], label='RgZ')
            plt.show()
            data1 = (
            'Rgmin = ' + str(g1min) + '\t\nRgmax = ' + str(g1max) + '\t\nRgmean =' + str(g1mean) + '\t\nRgxmin = ' + str(
                g2min)
            + '\t\nRgxmax = ' + str(g2max) + '\t\nRgxmean =' + str(g2mean) + '\t\nRgymin = ' + str(
                g3min) + '\t\nRgymax = ' + str(g3max) +
            '\t\nRgymean =' + str(g3mean) + '\t\nRgzmin = ' + str(g4min) + '\t\nRgzmax = ' + str(g4max) + '\t\nRgzmean =' + str(
                g4mean))
            text = """
LiGRO v 0.03 - Output of {0}
---------------------------------------------
            """.format(analysis)
            try:
              with open('rg_out.txt', 'w') as rg_infile:
                rg_out = tkFileDialog.asksaveasfile(mode='w', initialdir=path2, filetypes =(("Text File", "*.txt"),("All Files","*.*")),
    title = "Save RG statistic file.", initialfile='ligro_RG')
                rg_out.write(text+data1)
                rg_out.close()
                pass
            except:
              pass
            
        elif analysis == 'MSD':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                pass
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                pass
            else:
                pass
            if gromacs_flag('mdrun'):    
              command= ''' g_msd -s md_an.tpr -f md_an.xtc -o msd.xvg -xvg none << EOF
                                    {0}
                                    EOF '''.format(struct)    
            elif gromacs_flag('gmx'):    
              command= ''' gmx msd -s md_an.tpr -f md_an.xtc -o msd.xvg -xvg none << EOF
                                    {0}
                                    EOF '''.format(struct)
            else:
              pass
            os.system(command)
            pj = str(self.project.getvalue())
            an4 = 'msd_{0}.xvg'.format(pj)
            shutil.copy('msd.xvg', an4)
            directory = str(self.save.getvalue())
            shutil.copy2(an4, directory)
            t, rmsd = np.loadtxt("msd.xvg", unpack=True)
            fig = plt.figure(figsize=(10, 5))
            data = np.loadtxt('msd.xvg')
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
            text0 = """
LiGRO v 0.03 - Output of {0}
---------------------------------------------
            """.format(analysis)
            try:
              with open('msd_out.txt', 'w') as msd_infile:
                msd_out = tkFileDialog.asksaveasfile(mode='w', initialdir=path2, filetypes =(("Text File", "*.txt"),("All Files","*.*")),
    title = "Save MSD statistic file.", initialfile='ligro_MSD')
                msd_out.write(text0 + data0)
                msd_out.close()
                pass
            except:
              pass

        elif analysis == 'H_bond':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                mbox.showerror("Error", "TPR file not found.")
                pass
            elif find2 == False:
                mbox.showerror("Error", "XTC file not found.")
                pass
            else:
                pass
            if gromacs_flag('mdrun'):    
              command= ''' g_hbond -s md_an.tpr -f md_an.xtc -num hbond.xvg -xvg none << EOF
                                    Other
                                    Protein
                                    '''    
            elif gromacs_flag('gmx'):    
              command= ''' gmx hbond -s md_an.tpr -f md_an.xtc -num hbond.xvg -xvg none << EOF
                                    Other
                                    Protein
                                    '''
            else:
              pass
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
                pass

            else:
                pass

            if gromacs_flag('mdrun'):
              com1= '''echo 56 57 0 | g_energy -f md_an.edr > out.txt'''
            elif gromacs_flag('gmx'):
              com1= '''echo 56 57 0 | gmx energy -f md_an.edr > out.txt'''
            else:
              pass
            os.system(com1)


            with open('out.txt', 'r') as f:
              lines = f.readlines()
              for line in lines:
                if re.search(r'Coul-SR:', line):
                  match0 = re.split(r'\s+', line)

                elif re.search(r'LJ-SR:', line):
                  match1 = re.split(r'\s+', line)
                  pass

            coul = match0[1]
            lj = match1[1]
            a = float(coul)
            b = float(lj)
            c =float(a+b)
            y = format(c, '.2f')
            f = open('out.txt', 'r')
            cont = f.readlines()
            cont.append('\n\nE int = <E LJ> + <E Coul> = ' + str(y) + ' kJ/mol')
            f.close()
            f = open('out.txt', 'w')
            f.writelines(cont)
            f.close()


            try:
              with open('out.txt') as infile, tkFileDialog.asksaveasfile(mode='w', initialdir=path2, filetypes =(("Text File", "*.txt"),("All Files","*.*")),
  title = "Save LJSR-CoulSR IE txt file.", initialfile='LJSR-CoulSR_IE') as outfile:
                  outfile.write('LiGRO v 0.03\n')
                  outfile.write('Energy                      Average   Err.Est.       RMSD  Tot-Drift\n')
                  copy = False
                  for line in infile:
                      if line.strip() == "Energy                      Average   Err.Est.       RMSD  Tot-Drift":
                          copy = True
                      elif copy:
                          outfile.write(line)
            except:
              pass

    def plip(self):
        try:
            os.chdir(path)
        except:
            pass
        find1=os.path.exists('md_an.tpr')
        find2=os.path.exists('md_an.xtc')
        if find1 == False:
            mbox.showerror("Error", "TPR or XTC file not found.")
            pass
        elif find2 == False:
            mbox.showerror("Error", "TPR or XTC file not found.")
            pass
        else:
            pass
        timeplip = str(self.ft.getvalue())
        
        if gromacs_flag('mdrun'):
          command3 = '''trjconv -s md_an.tpr -f md_an.xtc -o plip.pdb -b {0} -e {0} << EOF
                   non-Water
                   EOF'''.format(timeplip)

        elif gromacs_flag('gmx'):
          command3 = '''gmx trjconv -s md_an.tpr -f md_an.xtc -o plip.pdb -b {0} -e {0} << EOF
                   non-Water
                   EOF'''.format(timeplip)
        else:
          pass
        os.system(command3)
        command4="plip -f plip.pdb -v -t --name plip_out"
        subprocess.call(['/bin/bash', '-i', '-c', command4])
        with open('plip_out.txt') as plip_infile:
          plip_lines = plip_infile.read()
          plip_out = tkFileDialog.asksaveasfile(mode='w', initialdir=path2, filetypes =(("Text File", "*.txt"),("All Files","*.*")),
  title = "Save PLIP report file.", initialfile='report')
          plip_out.write(plip_lines)
          plip_out.close()
          pass

        

######################################################################

# Create root window.

def main_init():
  if gromacs_flag('mdrun'):
    pass
  elif gromacs_flag('gmx'):
    pass
  else:
    print('Please wait to install GROMACS package')
    os.system('conda install -c bioconda gromacs')
    pass

  print(note)
  root = Tkinter.Tk()
  Pmw.initialise(root)
  root.option_add("*Font", "helvetica 10 bold")
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
          sys.exit()
      else:
          mbox.showinfo('No', 'Quit has been cancelled')


  exitButton = Tkinter.Button(root, text='Exit', command=callback)
  exitButton.pack()


  def help():
      webbrowser.open('https://www.ufrgs.br/lasomfarmacia/node/8', new=1, autoraise=True)


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
  img = Tkinter.PhotoImage(file = r'ic_launcher_.png')
  root.tk.call('wm', 'iconphoto', root._w, img)
  root.config(menu=menubar)
  root.mainloop()

if __name__ == '__main__':
  main_init()
