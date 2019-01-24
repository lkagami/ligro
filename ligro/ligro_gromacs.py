#!/usr/bin/python3
# -*- coding: utf-8 -*-
note ="""
---LiGRO - Version 0.7 ---

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
from plip.modules import report
plipver = report.__version__

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

title = 'LiGRO: Version 0.7'

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
        group_plip = Pmw.Group(tab5, tag_text='Protein-Ligand Interaction Profiler (PLIP) v{0}'.format(plipver))
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
                        EOF
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
              EOF
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
              EOF
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
              EOF
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
              EOF
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
              EOF
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
                          EOF
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
                    EOF
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
                            EOF
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
                EOF
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
                EOF
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
                  EOF
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
                      EOF
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
              EOF
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
                        EOF
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
            EOF
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
            EOF
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
            EOF
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
                        EOF
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
              EOF
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
              EOF
              '''.format(nst)

            elif gromacs_flag('gmx'):

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
                nstlist       = 1       ; Frequency to update the neighbor list and long range forces
                cutoff-scheme   = Verlet
                ns_type       = grid    ; Method to determine neighbor list (simple, grid)
                rlist       = 1.0   ; Cut-off for making neighbor list (short range forces)
                coulombtype     = PME   ; Treatment of long range electrostatic interactions
                rcoulomb      = 1.0   ; long range electrostatic cut-off
                rvdw        = 1.0   ; long range Van der Waals cut-off
                pbc           = xyz     ; Periodic Boundary Conditions
                EOF
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
                nstlist       = 1       ; Frequency to update the neighbor list and long range forces
                cutoff-scheme   = Verlet
                ns_type       = grid    ; Method to determine neighbor list (simple, grid)
                rlist       = 1.0   ; Cut-off for making neighbor list (short range forces)
                coulombtype     = PME   ; Treatment of long range electrostatic interactions
                rcoulomb      = 1.0   ; long range electrostatic cut-off
                rvdw        = 1.0   ; long range Van der Waals cut-off
                pbc           = xyz     ; Periodic Boundary Conditions
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
              EOF
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
              EOF
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
              coulombtype     = Cut-off       ;
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
              EOF
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
              nstenergy           = 500      ; save energies every 1.0 ps
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
              coulombtype     = Cut-off       ;
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
              os.system("echo 'mdrun -v -deffnm md' >> queue.sh")
              os.system('chmod 777 queue.sh')
              os.system('./queue.sh')

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
              os.system("echo 'gmx mdrun -v -deffnm md' >> queue.sh")
              os.system('chmod 777 queue.sh')
              os.system('./queue.sh')
          else:
            mbox.showerror('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            quit()
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
                      EOF
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
EOF
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
EOF
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
          os.system('''
                            make_ndx -f eml.gro -o index_lie.ndx << EOF
                            "Water" | "Ion"
                            q
                            EOF ''')


          stnvt = str(self.time_nvt.getvalue())
          stnpt = str(self.time_npt.getvalue())
          stmd = str(self.time_md.getvalue())
          temp = str(self.temperature.getvalue())
          t_st = str(self.time_st.getvalue())

          
          if gromacs_flag('mdrun'):
              os.system('''
                            make_ndx -f eml.gro -o index_lie.ndx << EOF
                            "Water" | "Ion"
                            q
                            EOF ''')
              cmd5 = '''
              cat << EOF >| nvt2.mdp
              title       = Ligand NVT
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
              ; Electrostatics
              coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
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
              EOF
              '''.format(stnvt, t_st, temp)

          elif gromacs_flag('gmx'):
             os.system('''
                            gmx make_ndx -f eml.gro -o index_lie.ndx << EOF
                            "Water" | "Ion"
                            q
                            EOF ''')
             cmd5 = '''cat << EOF >| nvt2.mdp
             title       = Ligand NVT
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
             ; Electrostatics
             coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
             pme_order       = 4         ; cubic interpolation
             fourierspacing  = 0.16      ; grid spacing for FFT
             ; Temperature coupling
             tcoupl      = V-rescale                     ; modified Berendsen thermostat
             tc-grps     = LIG   Water_Ion; two coupling groups - more accurate
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
             cat << EOF >| npt2.mdp
             title       = Ligand NPT
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
             tc-grps     = LIG   Water_Ion    ; two coupling groups - more accurate
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
             EOF
             '''.format(stnpt, t_st, temp)

          elif gromacs_flag('gmx'):
             cmd6 = '''
             cat << EOF >| npt2.mdp
             title       = Ligand NPT
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
             tc-grps     = LIG   Water_Ion   ; two coupling groups - more accurate
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
             cat << EOF >| md2.mdp
             title       = Ligand MD
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
             coulombtype     = Cut-off       ;
             ; Temperature coupling
             tcoupl      = V-rescale                     ; modified Berendsen thermostat
             tc-grps     = LIG Water_Ion    ; two coupling groups - more accurate
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
             EOF
             '''.format(stmd, t_st, temp)

          elif gromacs_flag('gmx'):
             cmd7 = '''
             cat << EOF >| md2.mdp
             title       = Ligand MD
             ; Run parameters
             integrator  = md        ; leap-frog integrator
             nsteps      = {0}    ; 2 * 500000 = 1000 ps (1 ns)
             dt          = {1}     ; 2 fs
             ; Output control
             nstxout             = 0         ; suppress .trr output
             nstvout             = 0         ; suppress .trr output
             nstenergy           = 500      ; save energies every 1.0 ps
             nstlog              = 5000      ; update log file every 10.0 ps
             nstxout-compressed  = 5000      ; write .xtc trajectory every 10.0 ps
             compressed-x-grps   = System
             energygrps          = LIG
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
             coulombtype     = Cut-off       ;
             ; Temperature coupling
             tcoupl      = V-rescale                     ; modified Berendsen thermostat
             tc-grps     = LIG Water_Ion    ; two coupling groups - more accurate
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

          cmd8='pymol eml.pdb'

          os.system(cmd8)

          if mbox.askyesno('View Ligand', 'Is solvated ligand OK??'):
            pass
          else:
               mbox.showinfo('No', 'Process has been cancelled')
               quit()

          if gromacs_flag('mdrun'):
             os.system('cat << EOF > | queue2.sh')
             os.system("echo 'grompp -f nvt2.mdp -c eml.pdb -p topoll.top -o nvt2.tpr -n index_lie.ndx -maxwarn 1000' >> queue2.sh")
             os.system("echo 'mdrun -v -deffnm nvt2' >> queue2.sh")
             os.system("echo 'grompp -f npt2.mdp -c nvt2.gro -p topoll.top -o npt2.tpr -n index_lie.ndx -maxwarn 1000' >> queue2.sh")
             os.system("echo 'mdrun -v -deffnm npt2' >> queue2.sh")
             os.system("echo 'grompp -f md2.mdp -c npt2.gro -p topoll.top -o md2.tpr -n index_lie.ndx -maxwarn 1000' >> queue2.sh")
             os.system("echo 'mdrun -v -deffnm md2' >> queue2.sh")
             os.system('chmod 777 queue2.sh')
             os.system('./queue2.sh')

          elif gromacs_flag('gmx'):
             os.system('cat << EOF > | queue2.sh')
             os.system("echo 'gmx grompp -f nvt2.mdp -c eml.pdb -p topoll.top -o nvt2.tpr -n index_lie.ndx -maxwarn 1000' >> queue2.sh")
             os.system("echo 'gmx mdrun -v -deffnm nvt2' >> queue2.sh")
             os.system("echo 'gmx grompp -f npt2.mdp -c nvt2.gro -p topoll.top -o npt2.tpr -n index_lie.ndx -maxwarn 1000' >> queue2.sh")
             os.system("echo 'gmx mdrun -v -deffnm npt2' >> queue2.sh")
             os.system("echo 'gmx grompp -f md2.mdp -c npt2.gro -p topoll.top -o md2.tpr -n index_lie.ndx -maxwarn 1000' >> queue2.sh")
             os.system("echo 'gmx mdrun -v -deffnm md2' >> queue2.sh")
             os.system('chmod 777 queue2.sh')
             os.system('./queue2.sh')
          else:
            mbox.showerror('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            quit()
        else:
          mbox.showerror('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

    def lie_calculation(self):
        try:
            os.chdir(path)
        except:
            pass
        if gromacs_flag('mdrun'):
          com1= '''echo 49 50 0 | g_energy -f md2.edr > out1.txt'''
        elif gromacs_flag('gmx'):
          com1= '''echo 50 51 0 | gmx energy -f md2.edr > out1.txt'''
        else:
          pass
        os.system(com1)
        if gromacs_flag('mdrun'):
          try:
            with open('out1.txt', 'r') as f:
              lines = f.readlines()
              for line in lines:
                if re.search(r'Coulomb \(SR\)', line):
                  print('Coulomb')
                  match0 = re.split(r'\s+', line)
                  pass
                elif re.search(r'LJ \(SR\)', line):
                  print('LJ')
                  match1 = re.split(r'\s+', line)
                  pass
                
          except UnboundLocalError:
            mbox.showerror('Error', 'Please try to use other simulation type box.')
            quit()
        elif gromacs_flag('gmx'):
          try:
            with open('out1.txt', 'r') as f:
              lines = f.readlines()
              for line in lines:
                if re.search(r'Coul-SR:', line):
                  match0 = re.split(r'\s+', line)

                elif re.search(r'LJ-SR:', line):
                  match1 = re.split(r'\s+', line)
                  pass
          except UnboundLocalError:
            mbox.showerror('Error', 'Please try to use other simulation type box (Dodecahedron).')
            quit()
        else:
          pass
        if gromacs_flag('mdrun'):  
          coul = match0[5]
          lj = match1[5]
        elif gromacs_flag('gmx'):
          coul = match0[4]
          lj = match1[4]
        else:
          pass

        if coul is not None:
          pass
        else:
          mbox.showerror('Error', 'Please try to use other simulation type box (Dodecahedron).')
          quit()

        if gromacs_flag('mdrun'):
          os.system('g_lie -f md.edr -o lie_lig.xvg -ligand LIG -Eqq {0} -Elj {1} >> out.txt'.format(coul,lj))
        elif gromacs_flag('gmx'):
          os.system('gmx lie -f md.edr -o lie_lig.xvg -ligand LIG -Eqq {0} -Elj {1} >> out.txt'.format(coul,lj))
        else:
          pass
        with open('out.txt', 'r') as f:
          lines = f.readlines()
          for line in lines:
            if re.search(r'DGbind =', line):
              match001 = re.split(r'\s+', line)

        dr = str(self.save.getvalue())
        pj = str(self.project.getvalue())
        pj1 = pj+'_LIE'
        lie0 = '{0}.txt'.format(pj1)
        shutil.copy('out.txt', lie0)
        shutil.copy2(lie0, dr)

        a = float(match001[2])
        y = format(a, '.2f')
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
LiGRO v 0.7 - Output of {0}
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
LiGRO v 0.7 - Output of {0}
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
LiGRO v 0.7 - Output of {0}
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
              com1= '''echo 9 7 0 | g_energy -f md_an.edr > out.txt'''
            elif gromacs_flag('gmx'):
              com1= '''echo 9 7 0 | gmx energy -f md_an.edr > out.txt'''
            else:
              pass
            os.system(com1)


            with open('out.txt', 'r') as f:
              lines = f.readlines()
              for line in lines:
                if re.search(r'Coulomb', line):
                  match0 = re.split(r'\s+\s+', line)

                elif re.search(r'LJ', line):
                  match1 = re.split(r'\s+\s+', line)
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
                  outfile.write('LiGRO v 0.7\n')
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
  ICO = base64.decodestring(b'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\n')
  img = Tkinter.PhotoImage(data=ICO)
  root.tk.call('wm', 'iconphoto', root._w, img)
  root.config(menu=menubar)
  root.mainloop()

if __name__ == '__main__':
  main_init()