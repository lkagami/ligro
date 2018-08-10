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
from PIL import ImageTk
import zlib, base64

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
  ICO = zlib.decompress(base64.b64decode('''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'''))
  img = Tkinter.PhotoImage(data=ICO)
  root.tk.call('wm', 'iconphoto', root._w, img)
  root.config(menu=menubar)
  root.mainloop()

if __name__ == '__main__':
  main_init()