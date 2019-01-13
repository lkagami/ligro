# -*- coding: utf-8 -*-

# need install libgl1-mesa-glx
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
import os
import subprocess
import sys
import shutil
import webbrowser
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
#from PIL import ImageTk
import zlib, base64
#from plip.modules import report
#plipver = report.__version__

sys.path[:0] = ['../../..']

path2 = os.environ.get('HOME')

td = str(date.today())

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

title = 'LiGRO: Version 1.0'

class Ui_MainWindo(object):
    def setupUi(self, MainWindo):
        MainWindo.setObjectName("MainWindo")
        MainWindo.resize(800, 545)
        self.centralwidget = QtWidgets.QWidget(MainWindo)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 60, 781, 441))
        self.tabWidget.setAccessibleName("")
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.groupBox_3 = QtWidgets.QGroupBox(self.tab)
        self.groupBox_3.setGeometry(QtCore.QRect(20, 220, 711, 80))
        self.groupBox_3.setObjectName("groupBox_3")
        self.browse_lig_mol2 = QtWidgets.QPushButton(self.groupBox_3)
        self.browse_lig_mol2.setGeometry(QtCore.QRect(550, 40, 75, 23))
        self.browse_lig_mol2.setObjectName("browse_lig_mol2")
        self.browse_lig_mol2.clicked.connect(self.openconfmol2file)
        self.sel_lig_mol2 = QtWidgets.QLineEdit(self.groupBox_3)
        self.sel_lig_mol2.setGeometry(QtCore.QRect(30, 40, 481, 21))
        self.sel_lig_mol2.setObjectName("sel_lig_mol2")
        self.sel_lig_mol2.setText('No select file')
        self.groupBox_2 = QtWidgets.QGroupBox(self.tab)
        self.groupBox_2.setGeometry(QtCore.QRect(20, 120, 711, 81))
        self.groupBox_2.setObjectName("groupBox_2")
        self.browse_cof_mol2 = QtWidgets.QPushButton(self.groupBox_2)
        self.browse_cof_mol2.setGeometry(QtCore.QRect(550, 40, 75, 23))
        self.browse_cof_mol2.setObjectName("browse_cof_mol2")
        self.browse_cof_mol2.clicked.connect(self.openconfmol2file)
        self.sel_cof_mol2 = QtWidgets.QLineEdit(self.groupBox_2)
        self.sel_cof_mol2.setGeometry(QtCore.QRect(30, 40, 481, 21))
        self.sel_cof_mol2.setObjectName("sel_cof_mol2")
        self.sel_cof_mol2.setText('No select file')
        self.groupBox = QtWidgets.QGroupBox(self.tab)
        self.groupBox.setGeometry(QtCore.QRect(20, 30, 711, 80))
        self.groupBox.setObjectName("groupBox")
        self.sel_pdb = QtWidgets.QLineEdit(self.groupBox)
        self.sel_pdb.setGeometry(QtCore.QRect(30, 40, 481, 21))
        self.sel_pdb.setObjectName("sel_pdb")
        self.sel_pdb.setText('No select file')
        self.browse_pdb = QtWidgets.QPushButton(self.groupBox)
        self.browse_pdb.setGeometry(QtCore.QRect(550, 40, 75, 23))
        self.browse_pdb.setObjectName("browse_pdb")
        self.browse_pdb.clicked.connect(self.openpdbfile)
        self.tabWidget.addTab(self.tab, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.groupBox_4 = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox_4.setGeometry(QtCore.QRect(10, 100, 731, 80))
        self.groupBox_4.setObjectName("groupBox_4")
        self.label_4 = QtWidgets.QLabel(self.groupBox_4)
        self.label_4.setGeometry(QtCore.QRect(10, 30, 81, 16))
        self.label_4.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.label_4.setObjectName("label_4")
        self.csel_ch_met = QtWidgets.QComboBox(self.groupBox_4)
        self.csel_ch_met.setGeometry(QtCore.QRect(110, 30, 69, 22))
        self.csel_ch_met.setObjectName("csel_ch_met")
        self.csel_ch_met.addItem("")
        self.csel_ch_met.addItem("")
        self.csel_ch_met.addItem("")
        self.label_5 = QtWidgets.QLabel(self.groupBox_4)
        self.label_5.setGeometry(QtCore.QRect(190, 30, 61, 16))
        self.label_5.setObjectName("label_5")
        self.sel_net_ch = QtWidgets.QSpinBox(self.groupBox_4)
        self.sel_net_ch.setGeometry(QtCore.QRect(260, 30, 42, 22))
        self.sel_net_ch.setObjectName("sel_net_ch")
        self.label_6 = QtWidgets.QLabel(self.groupBox_4)
        self.label_6.setGeometry(QtCore.QRect(320, 30, 61, 16))
        self.label_6.setObjectName("label_6")
        self.sel_at_type = QtWidgets.QComboBox(self.groupBox_4)
        self.sel_at_type.setGeometry(QtCore.QRect(380, 30, 69, 22))
        self.sel_at_type.setObjectName("sel_at_type")
        self.sel_at_type.addItem("")
        self.sel_at_type.addItem("")
        self.sel_at_type.addItem("")
        self.sel_at_type.addItem("")
        self.label_7 = QtWidgets.QLabel(self.groupBox_4)
        self.label_7.setGeometry(QtCore.QRect(470, 30, 61, 16))
        self.label_7.setObjectName("label_7")
        self.sel_mult = QtWidgets.QSpinBox(self.groupBox_4)
        self.sel_mult.setGeometry(QtCore.QRect(530, 30, 42, 22))
        self.sel_mult.setProperty("value", 1)
        self.sel_mult.setObjectName("sel_mult")
        self.groupBox_7 = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox_7.setGeometry(QtCore.QRect(10, 10, 731, 80))
        self.groupBox_7.setObjectName("groupBox_7")
        self.label_2 = QtWidgets.QLabel(self.groupBox_7)
        self.label_2.setGeometry(QtCore.QRect(10, 30, 101, 16))
        self.label_2.setObjectName("label_2")
        self.sel_ff = QtWidgets.QComboBox(self.groupBox_7)
        self.sel_ff.setGeometry(QtCore.QRect(110, 30, 111, 22))
        self.sel_ff.setAccessibleName("")
        self.sel_ff.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.sel_ff.setObjectName("sel_ff")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.sel_ff.addItem("")
        self.label_3 = QtWidgets.QLabel(self.groupBox_7)
        self.label_3.setGeometry(QtCore.QRect(320, 30, 47, 13))
        self.label_3.setObjectName("label_3")
        self.sel_metal = QtWidgets.QComboBox(self.groupBox_7)
        self.sel_metal.setGeometry(QtCore.QRect(370, 30, 51, 22))
        self.sel_metal.setObjectName("sel_metal")
        self.sel_metal.addItem("")
        self.sel_metal.addItem("")
        self.sel_metal.addItem("")
        self.sel_metal.addItem("")
        self.sel_metal.addItem("")
        self.ig_hyd = QtWidgets.QCheckBox(self.groupBox_7)
        self.ig_hyd.setGeometry(QtCore.QRect(520, 30, 151, 20))
        self.ig_hyd.setChecked(True)
        self.ig_hyd.setObjectName("ig_hyd")
        self.groupBox_5 = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox_5.setGeometry(QtCore.QRect(10, 210, 731, 80))
        self.groupBox_5.setObjectName("groupBox_5")
        self.label_8 = QtWidgets.QLabel(self.groupBox_5)
        self.label_8.setGeometry(QtCore.QRect(10, 30, 101, 16))
        self.label_8.setObjectName("label_8")
        self.comboBox_5 = QtWidgets.QComboBox(self.groupBox_5)
        self.comboBox_5.setGeometry(QtCore.QRect(360, 30, 91, 22))
        self.comboBox_5.setObjectName("comboBox_5")
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.sel_box = QtWidgets.QLabel(self.groupBox_5)
        self.sel_box.setGeometry(QtCore.QRect(270, 30, 91, 16))
        self.sel_box.setObjectName("sel_box")
        self.sel_water = QtWidgets.QComboBox(self.groupBox_5)
        self.sel_water.setGeometry(QtCore.QRect(110, 30, 69, 22))
        self.sel_water.setObjectName("sel_water")
        self.sel_water.addItem("")
        self.sel_water.addItem("")
        self.sel_water.addItem("")
        self.label_10 = QtWidgets.QLabel(self.groupBox_5)
        self.label_10.setGeometry(QtCore.QRect(530, 30, 47, 13))
        self.label_10.setObjectName("label_10")
        self.sel_dist = QtWidgets.QDoubleSpinBox(self.groupBox_5)
        self.sel_dist.setGeometry(QtCore.QRect(590, 30, 62, 22))
        self.sel_dist.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.sel_dist.setSingleStep(0.01)
        self.sel_dist.setProperty("value", 1.0)
        self.sel_dist.setObjectName("sel_dist")
        self.groupBox_6 = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox_6.setGeometry(QtCore.QRect(10, 310, 731, 80))
        self.groupBox_6.setObjectName("groupBox_6")
        self.sel_ioniz = QtWidgets.QComboBox(self.groupBox_6)
        self.sel_ioniz.setGeometry(QtCore.QRect(10, 30, 141, 22))
        self.sel_ioniz.setObjectName("sel_ioniz")
        self.sel_ioniz.addItem("")
        self.sel_ioniz.addItem("")
        self.sel_ioniz.addItem("")
        self.sel_conc = QtWidgets.QDoubleSpinBox(self.groupBox_6)
        self.sel_conc.setGeometry(QtCore.QRect(170, 30, 62, 22))
        self.sel_conc.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.sel_conc.setProperty("value", 0.15)
        self.sel_conc.setObjectName("sel_conc")
        self.tabWidget.addTab(self.tab_3, "")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.groupBox_8 = QtWidgets.QGroupBox(self.tab_4)
        self.groupBox_8.setGeometry(QtCore.QRect(10, 20, 241, 80))
        self.groupBox_8.setObjectName("groupBox_8")
        self.label_11 = QtWidgets.QLabel(self.groupBox_8)
        self.label_11.setGeometry(QtCore.QRect(10, 20, 151, 16))
        self.label_11.setObjectName("label_11")
        self.label_12 = QtWidgets.QLabel(self.groupBox_8)
        self.label_12.setGeometry(QtCore.QRect(10, 50, 91, 16))
        self.label_12.setObjectName("label_12")
        self.sel_int_step = QtWidgets.QDoubleSpinBox(self.groupBox_8)
        self.sel_int_step.setGeometry(QtCore.QRect(160, 20, 62, 22))
        self.sel_int_step.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.sel_int_step.setDecimals(3)
        self.sel_int_step.setSingleStep(0.001)
        self.sel_int_step.setProperty("value", 0.002)
        self.sel_int_step.setObjectName("sel_int_step")
        self.sel_temp = QtWidgets.QDoubleSpinBox(self.groupBox_8)
        self.sel_temp.setGeometry(QtCore.QRect(160, 50, 62, 22))
        self.sel_temp.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.sel_temp.setDecimals(2)
        self.sel_temp.setMaximum(400.0)
        self.sel_temp.setProperty("value", 310.15)
        self.sel_temp.setObjectName("sel_temp")
        self.groupBox_11 = QtWidgets.QGroupBox(self.tab_4)
        self.groupBox_11.setGeometry(QtCore.QRect(410, 20, 131, 80))
        self.groupBox_11.setObjectName("groupBox_11")
        self.sel_nvt_step = QtWidgets.QSpinBox(self.groupBox_11)
        self.sel_nvt_step.setGeometry(QtCore.QRect(50, 30, 71, 22))
        self.sel_nvt_step.setMinimum(100)
        self.sel_nvt_step.setMaximum(999999999)
        self.sel_nvt_step.setSingleStep(1000)
        self.sel_nvt_step.setProperty("value", 5000000)
        self.sel_nvt_step.setObjectName("sel_nvt_step")
        self.label_14 = QtWidgets.QLabel(self.groupBox_11)
        self.label_14.setGeometry(QtCore.QRect(10, 30, 47, 13))
        self.label_14.setObjectName("label_14")
        self.groupBox_12 = QtWidgets.QGroupBox(self.tab_4)
        self.groupBox_12.setGeometry(QtCore.QRect(569, 20, 141, 80))
        self.groupBox_12.setObjectName("groupBox_12")
        self.label_15 = QtWidgets.QLabel(self.groupBox_12)
        self.label_15.setGeometry(QtCore.QRect(10, 30, 47, 13))
        self.label_15.setObjectName("label_15")
        self.sel_npt_step = QtWidgets.QSpinBox(self.groupBox_12)
        self.sel_npt_step.setGeometry(QtCore.QRect(50, 30, 71, 22))
        self.sel_npt_step.setMinimum(100)
        self.sel_npt_step.setMaximum(999999999)
        self.sel_npt_step.setSingleStep(1000)
        self.sel_npt_step.setProperty("value", 5000000)
        self.sel_npt_step.setObjectName("sel_npt_step")
        self.groupBox_13 = QtWidgets.QGroupBox(self.tab_4)
        self.groupBox_13.setGeometry(QtCore.QRect(10, 120, 141, 80))
        self.groupBox_13.setObjectName("groupBox_13")
        self.label_16 = QtWidgets.QLabel(self.groupBox_13)
        self.label_16.setGeometry(QtCore.QRect(20, 30, 47, 13))
        self.label_16.setObjectName("label_16")
        self.sel_md_step = QtWidgets.QSpinBox(self.groupBox_13)
        self.sel_md_step.setGeometry(QtCore.QRect(60, 30, 71, 22))
        self.sel_md_step.setMinimum(100)
        self.sel_md_step.setMaximum(999999999)
        self.sel_md_step.setSingleStep(1000)
        self.sel_md_step.setProperty("value", 50000000)
        self.sel_md_step.setObjectName("sel_md_step")
        self.groupBox_10 = QtWidgets.QGroupBox(self.tab_4)
        self.groupBox_10.setGeometry(QtCore.QRect(270, 20, 131, 80))
        self.groupBox_10.setObjectName("groupBox_10")
        self.label_13 = QtWidgets.QLabel(self.groupBox_10)
        self.label_13.setGeometry(QtCore.QRect(10, 20, 47, 13))
        self.label_13.setObjectName("label_13")
        self.sel_min_step = QtWidgets.QSpinBox(self.groupBox_10)
        self.sel_min_step.setGeometry(QtCore.QRect(50, 20, 61, 22))
        self.sel_min_step.setMinimum(100)
        self.sel_min_step.setMaximum(999999999)
        self.sel_min_step.setSingleStep(100)
        self.sel_min_step.setProperty("value", 1000)
        self.sel_min_step.setObjectName("sel_min_step")
        self.sel_min_alg = QtWidgets.QComboBox(self.groupBox_10)
        self.sel_min_alg.setGeometry(QtCore.QRect(8, 50, 101, 22))
        self.sel_min_alg.setObjectName("sel_min_alg")
        self.sel_min_alg.addItem("")
        self.sel_min_alg.addItem("")
        self.tabWidget.addTab(self.tab_4, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.groupBox_9 = QtWidgets.QGroupBox(self.tab_5)
        self.groupBox_9.setGeometry(QtCore.QRect(10, 20, 761, 80))
        self.groupBox_9.setObjectName("groupBox_9")
        self.label_17 = QtWidgets.QLabel(self.groupBox_9)
        self.label_17.setGeometry(QtCore.QRect(10, 20, 71, 16))
        self.label_17.setObjectName("label_17")
        self.sel_proj_name = QtWidgets.QLineEdit(self.groupBox_9)
        self.sel_proj_name.setGeometry(QtCore.QRect(90, 20, 181, 20))
        self.sel_proj_name.setObjectName("sel_proj_name")
        self.sel_proj_name.setText('MD_'+td)
        self.sel_sav_direc = QtWidgets.QLineEdit(self.groupBox_9)
        self.sel_sav_direc.setGeometry(QtCore.QRect(370, 20, 361, 20))
        self.sel_sav_direc.setObjectName("sel_sav_direc")
        self.sel_sav_direc.setText(path2)
        self.label_18 = QtWidgets.QLabel(self.groupBox_9)
        self.label_18.setGeometry(QtCore.QRect(290, 20, 81, 16))
        self.label_18.setObjectName("label_18")
        self.sel_exp_fol = QtWidgets.QCheckBox(self.groupBox_9)
        self.sel_exp_fol.setGeometry(QtCore.QRect(10, 50, 111, 17))
        self.sel_exp_fol.setChecked(False)
        self.sel_exp_fol.setObjectName("sel_exp_fol")
        self.save_tpr = QtWidgets.QPushButton(self.groupBox_9)
        self.save_tpr.setGeometry(QtCore.QRect(290, 50, 75, 23))
        self.save_tpr.setObjectName("save_tpr")
        self.run_md = QtWidgets.QPushButton(self.groupBox_9)
        self.run_md.setGeometry(QtCore.QRect(660, 50, 75, 23))
        self.run_md.setObjectName("run_md")
        self.groupBox_14 = QtWidgets.QGroupBox(self.tab_5)
        self.groupBox_14.setGeometry(QtCore.QRect(10, 110, 281, 80))
        self.groupBox_14.setObjectName("groupBox_14")
        self.label_19 = QtWidgets.QLabel(self.groupBox_14)
        self.label_19.setGeometry(QtCore.QRect(10, 30, 61, 16))
        self.label_19.setObjectName("label_19")
        self.lie_val = QtWidgets.QLineEdit(self.groupBox_14)
        self.lie_val.setGeometry(QtCore.QRect(90, 30, 61, 20))
        self.lie_val.setObjectName("lie_val")
        self.run_lie = QtWidgets.QPushButton(self.groupBox_14)
        self.run_lie.setGeometry(QtCore.QRect(160, 30, 111, 23))
        self.run_lie.setObjectName("run_lie")
        self.groupBox_15 = QtWidgets.QGroupBox(self.tab_5)
        self.groupBox_15.setGeometry(QtCore.QRect(300, 110, 471, 80))
        self.groupBox_15.setObjectName("groupBox_15")
        self.sel_bkp_fol = QtWidgets.QPushButton(self.groupBox_15)
        self.sel_bkp_fol.setGeometry(QtCore.QRect(370, 30, 75, 23))
        self.sel_bkp_fol.setObjectName("sel_bkp_fol")
        self.bkp_fold = QtWidgets.QLineEdit(self.groupBox_15)
        self.bkp_fold.setGeometry(QtCore.QRect(10, 30, 341, 21))
        self.bkp_fold.setObjectName("bkp_fold")
        self.groupBox_16 = QtWidgets.QGroupBox(self.tab_5)
        self.groupBox_16.setGeometry(QtCore.QRect(10, 200, 761, 201))
        self.groupBox_16.setObjectName("groupBox_16")
        self.md_out = QtWidgets.QPlainTextEdit(self.groupBox_16)
        self.md_out.setGeometry(QtCore.QRect(0, 20, 761, 171))
        self.md_out.setObjectName("md_out")
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.groupBox_17 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_17.setGeometry(QtCore.QRect(10, 10, 751, 111))
        self.groupBox_17.setObjectName("groupBox_17")
        self.sel_tpr = QtWidgets.QLineEdit(self.groupBox_17)
        self.sel_tpr.setGeometry(QtCore.QRect(130, 20, 481, 21))
        self.sel_tpr.setObjectName("sel_tpr")
        self.browse_tpr = QtWidgets.QPushButton(self.groupBox_17)
        self.browse_tpr.setGeometry(QtCore.QRect(650, 20, 75, 23))
        self.browse_tpr.setObjectName("browse_tpr")
        self.browse_xtc = QtWidgets.QPushButton(self.groupBox_17)
        self.browse_xtc.setGeometry(QtCore.QRect(650, 50, 75, 23))
        self.browse_xtc.setObjectName("browse_xtc")
        self.sel_xtc = QtWidgets.QLineEdit(self.groupBox_17)
        self.sel_xtc.setGeometry(QtCore.QRect(130, 50, 481, 21))
        self.sel_xtc.setObjectName("sel_xtc")
        self.sel_edr = QtWidgets.QLineEdit(self.groupBox_17)
        self.sel_edr.setGeometry(QtCore.QRect(130, 80, 481, 21))
        self.sel_edr.setObjectName("sel_edr")
        self.browse_edr = QtWidgets.QPushButton(self.groupBox_17)
        self.browse_edr.setGeometry(QtCore.QRect(650, 80, 75, 23))
        self.browse_edr.setObjectName("browse_edr")
        self.label_20 = QtWidgets.QLabel(self.groupBox_17)
        self.label_20.setGeometry(QtCore.QRect(16, 20, 81, 20))
        self.label_20.setObjectName("label_20")
        self.label_21 = QtWidgets.QLabel(self.groupBox_17)
        self.label_21.setGeometry(QtCore.QRect(16, 50, 81, 20))
        self.label_21.setObjectName("label_21")
        self.label_22 = QtWidgets.QLabel(self.groupBox_17)
        self.label_22.setGeometry(QtCore.QRect(16, 80, 81, 20))
        self.label_22.setObjectName("label_22")
        self.groupBox_18 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_18.setGeometry(QtCore.QRect(10, 130, 751, 80))
        self.groupBox_18.setObjectName("groupBox_18")
        self.label_23 = QtWidgets.QLabel(self.groupBox_18)
        self.label_23.setGeometry(QtCore.QRect(10, 30, 81, 16))
        self.label_23.setObjectName("label_23")
        self.sel_struc = QtWidgets.QComboBox(self.groupBox_18)
        self.sel_struc.setGeometry(QtCore.QRect(100, 30, 101, 22))
        self.sel_struc.setObjectName("sel_struc")
        self.sel_struc.addItem("")
        self.sel_struc.addItem("")
        self.sel_struc.addItem("")
        self.sel_struc.addItem("")
        self.sel_struc.addItem("")
        self.label_24 = QtWidgets.QLabel(self.groupBox_18)
        self.label_24.setGeometry(QtCore.QRect(220, 30, 81, 16))
        self.label_24.setObjectName("label_24")
        self.sel_analy = QtWidgets.QComboBox(self.groupBox_18)
        self.sel_analy.setGeometry(QtCore.QRect(310, 30, 101, 22))
        self.sel_analy.setObjectName("sel_analy")
        self.sel_analy.addItem("")
        self.sel_analy.addItem("")
        self.sel_analy.addItem("")
        self.sel_analy.addItem("")
        self.sel_analy.addItem("")
        self.sel_analy.addItem("")
        self.run_analy = QtWidgets.QPushButton(self.groupBox_18)
        self.run_analy.setGeometry(QtCore.QRect(460, 30, 75, 23))
        self.run_analy.setObjectName("run_analy")
        self.groupBox_19 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_19.setGeometry(QtCore.QRect(10, 220, 751, 51))
        self.groupBox_19.setObjectName("groupBox_19")
        self.label_25 = QtWidgets.QLabel(self.groupBox_19)
        self.label_25.setGeometry(QtCore.QRect(20, 20, 81, 16))
        self.label_25.setObjectName("label_25")
        self.sel_plip_time = QtWidgets.QLineEdit(self.groupBox_19)
        self.sel_plip_time.setGeometry(QtCore.QRect(100, 20, 113, 20))
        self.sel_plip_time.setObjectName("sel_plip_time")
        self.run_plip = QtWidgets.QPushButton(self.groupBox_19)
        self.run_plip.setGeometry(QtCore.QRect(230, 20, 75, 23))
        self.run_plip.setObjectName("run_plip")
        self.groupBox_20 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_20.setGeometry(QtCore.QRect(10, 280, 751, 121))
        self.groupBox_20.setObjectName("groupBox_20")
        self.analy_out = QtWidgets.QPlainTextEdit(self.groupBox_20)
        self.analy_out.setGeometry(QtCore.QRect(0, 20, 751, 101))
        self.analy_out.setObjectName("analy_out")
        self.tabWidget.addTab(self.tab_2, "")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(90, 10, 631, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setGeometry(QtCore.QRect(10, 40, 781, 16))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        MainWindo.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindo)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        self.menu_Help = QtWidgets.QMenu(self.menubar)
        self.menu_Help.setObjectName("menu_Help")
        self.menuHelp = QtWidgets.QMenu(self.menu_Help)
        self.menuHelp.setObjectName("menuHelp")
        MainWindo.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindo)
        self.statusbar.setObjectName("statusbar")
        MainWindo.setStatusBar(self.statusbar)
        self.actionExit = QtWidgets.QAction(MainWindo)
        self.actionExit.setObjectName("actionExit")
        self.actionAbout = QtWidgets.QAction(MainWindo)
        self.actionAbout.setObjectName("actionAbout")
        self.actionHelp_2 = QtWidgets.QAction(MainWindo)
        self.actionHelp_2.setObjectName("actionHelp_2")
        self.menuHelp.addAction(self.actionAbout)
        self.actionAbout.triggered.connect(self.about)
        self.menuHelp.addAction(self.actionHelp_2)
        self.actionHelp_2.triggered.connect(self.help)
        self.menu_Help.addAction(self.actionExit)
        self.actionExit.triggered.connect(self.exit)
        self.menu_Help.addAction(self.menuHelp.menuAction())
        self.menubar.addAction(self.menu_Help.menuAction())

        self.retranslateUi(MainWindo)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindo)

    def retranslateUi(self, MainWindo):
        _translate = QtCore.QCoreApplication.translate
        MainWindo.setWindowTitle(_translate("MainWindo", title))
        self.groupBox_3.setTitle(_translate("MainWindo", "Select Ligand MOL2 File"))
        self.browse_lig_mol2.setText(_translate("MainWindo", "Browse"))
        self.groupBox_2.setTitle(_translate("MainWindo", "Select Cofactor MOL2 File"))
        self.browse_cof_mol2.setText(_translate("MainWindo", "Browse"))
        self.groupBox.setTitle(_translate("MainWindo", "Select PDB Files"))
        self.browse_pdb.setText(_translate("MainWindo", "Browse"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindo", "Input Files"))
        self.groupBox_4.setTitle(_translate("MainWindo", "Ligand / Cofactor Parameters (ACPYPE)"))
        self.label_4.setText(_translate("MainWindo", "Charge Method:"))
        self.csel_ch_met.setItemText(0, _translate("MainWindo", "bcc"))
        self.csel_ch_met.setItemText(1, _translate("MainWindo", "gas"))
        self.csel_ch_met.setItemText(2, _translate("MainWindo", "user"))
        self.label_5.setText(_translate("MainWindo", "Net Charge:"))
        self.label_6.setText(_translate("MainWindo", "Atom type:"))
        self.sel_at_type.setItemText(0, _translate("MainWindo", "gaff"))
        self.sel_at_type.setItemText(1, _translate("MainWindo", "amber"))
        self.sel_at_type.setItemText(2, _translate("MainWindo", "gaff2"))
        self.sel_at_type.setItemText(3, _translate("MainWindo", "amber2"))
        self.label_7.setText(_translate("MainWindo", "Multiplicity:"))
        self.groupBox_7.setTitle(_translate("MainWindo", "Protein Parameters"))
        self.label_2.setText(_translate("MainWindo", "Select Force Field:"))
        self.sel_ff.setCurrentText(_translate("MainWindo", "AMBER03"))
        self.sel_ff.setItemText(0, _translate("MainWindo", "AMBER03"))
        self.sel_ff.setItemText(1, _translate("MainWindo", "AMBER94"))
        self.sel_ff.setItemText(2, _translate("MainWindo", "AMBER96"))
        self.sel_ff.setItemText(3, _translate("MainWindo", "AMBER99"))
        self.sel_ff.setItemText(4, _translate("MainWindo", "AMBER99SB"))
        self.sel_ff.setItemText(5, _translate("MainWindo", "AMBER99SB-ILDN"))
        self.sel_ff.setItemText(6, _translate("MainWindo", "AMBERGS"))
        self.sel_ff.setItemText(7, _translate("MainWindo", "OPLS-AA/L"))
        self.label_3.setText(_translate("MainWindo", "Metal:"))
        self.sel_metal.setItemText(0, _translate("MainWindo", "None"))
        self.sel_metal.setItemText(1, _translate("MainWindo", "CA"))
        self.sel_metal.setItemText(2, _translate("MainWindo", "FE"))
        self.sel_metal.setItemText(3, _translate("MainWindo", "MG"))
        self.sel_metal.setItemText(4, _translate("MainWindo", "ZN"))
        self.ig_hyd.setText(_translate("MainWindo", "Ignore Hydrogen Atoms"))
        self.groupBox_5.setTitle(_translate("MainWindo", "Solvate"))
        self.label_8.setText(_translate("MainWindo", "Select Water Model:"))
        self.comboBox_5.setItemText(0, _translate("MainWindo", "triclinic"))
        self.comboBox_5.setItemText(1, _translate("MainWindo", "cubic"))
        self.comboBox_5.setItemText(2, _translate("MainWindo", "dodecahedron"))
        self.comboBox_5.setItemText(3, _translate("MainWindo", "octahedron"))
        self.sel_box.setText(_translate("MainWindo", "Select Box Type:"))
        self.sel_water.setItemText(0, _translate("MainWindo", "spc"))
        self.sel_water.setItemText(1, _translate("MainWindo", "spce"))
        self.sel_water.setItemText(2, _translate("MainWindo", "tip3p"))
        self.label_10.setText(_translate("MainWindo", "Distance:"))
        self.groupBox_6.setTitle(_translate("MainWindo", "Neutralization / Ionize"))
        self.sel_ioniz.setItemText(0, _translate("MainWindo", "Concentration (mol/liter)"))
        self.sel_ioniz.setItemText(1, _translate("MainWindo", "Na (Number)"))
        self.sel_ioniz.setItemText(2, _translate("MainWindo", "Cl (Number)"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindo", "System Parametrization"))
        self.groupBox_8.setTitle(_translate("MainWindo", "General Options"))
        self.label_11.setText(_translate("MainWindo", "Time step for integration (ps):"))
        self.label_12.setText(_translate("MainWindo", "Temperature (K):"))
        self.groupBox_11.setTitle(_translate("MainWindo", "NVT"))
        self.label_14.setText(_translate("MainWindo", "Steps:"))
        self.groupBox_12.setTitle(_translate("MainWindo", "NPT"))
        self.label_15.setText(_translate("MainWindo", "Steps:"))
        self.groupBox_13.setTitle(_translate("MainWindo", "Molecular Dynamics"))
        self.label_16.setText(_translate("MainWindo", "Steps:"))
        self.groupBox_10.setTitle(_translate("MainWindo", "Minimization"))
        self.label_13.setText(_translate("MainWindo", "Steps:"))
        self.sel_min_alg.setItemText(0, _translate("MainWindo", "SD Algorithm"))
        self.sel_min_alg.setItemText(1, _translate("MainWindo", "SD + CG"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("MainWindo", "MDP Parameters"))
        self.groupBox_9.setTitle(_translate("MainWindo", "Options"))
        self.label_17.setText(_translate("MainWindo", "Project Name:"))
        self.label_18.setText(_translate("MainWindo", "Save Directory:"))
        self.sel_exp_fol.setText(_translate("MainWindo", "Explain MD folder"))
        self.save_tpr.setText(_translate("MainWindo", "Save TPR File"))
        self.run_md.setText(_translate("MainWindo", "Run Dynamics"))
        self.groupBox_14.setTitle(_translate("MainWindo", "LIE Energy Calculation"))
        self.label_19.setText(_translate("MainWindo", "DG (kJ/mol)"))
        self.lie_val.setText(_translate("MainWindo", "Empty"))
        self.run_lie.setText(_translate("MainWindo", "Run LIE Calculation"))
        self.groupBox_15.setTitle(_translate("MainWindo", "Backup (Only MD)"))
        self.sel_bkp_fol.setText(_translate("MainWindo", "Browse"))
        self.bkp_fold.setText(_translate("MainWindo", "No select folder"))
        self.groupBox_16.setTitle(_translate("MainWindo", "Output"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), _translate("MainWindo", "Run MD"))
        self.groupBox_17.setTitle(_translate("MainWindo", "Select Input Files"))
        self.browse_tpr.setText(_translate("MainWindo", "Browse"))
        self.browse_xtc.setText(_translate("MainWindo", "Browse"))
        self.browse_edr.setText(_translate("MainWindo", "Browse"))
        self.label_20.setText(_translate("MainWindo", "Select TPR file:"))
        self.label_21.setText(_translate("MainWindo", "Select XTC file:"))
        self.label_22.setText(_translate("MainWindo", "Select EDR file:"))
        self.groupBox_18.setTitle(_translate("MainWindo", "Analysis"))
        self.label_23.setText(_translate("MainWindo", "Select structure:"))
        self.sel_struc.setItemText(0, _translate("MainWindo", "Protein"))
        self.sel_struc.setItemText(1, _translate("MainWindo", "Protein-H"))
        self.sel_struc.setItemText(2, _translate("MainWindo", "C-alpha"))
        self.sel_struc.setItemText(3, _translate("MainWindo", "Backbone"))
        self.sel_struc.setItemText(4, _translate("MainWindo", "LIG"))
        self.label_24.setText(_translate("MainWindo", "Select Analysis:"))
        self.sel_analy.setItemText(0, _translate("MainWindo", "RMSD"))
        self.sel_analy.setItemText(1, _translate("MainWindo", "RMSF"))
        self.sel_analy.setItemText(2, _translate("MainWindo", "RG"))
        self.sel_analy.setItemText(3, _translate("MainWindo", "MSD"))
        self.sel_analy.setItemText(4, _translate("MainWindo", "H_bond"))
        self.sel_analy.setItemText(5, _translate("MainWindo", "LJSR-CoulSR IE"))
        self.run_analy.setText(_translate("MainWindo", "Run"))
        self.groupBox_19.setTitle(_translate("MainWindo", "Protein-Ligand Interaction Profiler (PLIP) v"))
        self.label_25.setText(_translate("MainWindo", "Frame time (ps):"))
        self.sel_plip_time.setText(_translate("MainWindo", "100"))
        self.run_plip.setText(_translate("MainWindo", "Run"))
        self.groupBox_20.setTitle(_translate("MainWindo", "Output"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindo", "Analysis"))
        self.label.setText(_translate("MainWindo", "LiGRO: a graphical user interface for protein–ligand molecular dynamics."))
        self.menu_Help.setTitle(_translate("MainWindo", "&File"))
        self.menuHelp.setTitle(_translate("MainWindo", "Help"))
        self.actionExit.setText(_translate("MainWindo", "Exit"))
        self.actionAbout.setText(_translate("MainWindo", "About"))
        self.actionHelp_2.setText(_translate("MainWindo", "Help"))

    
    def openpdbfile(self):
    	try:
    		pdb_dialog = QtWidgets.QFileDialog()
    		pdb_file = pdb_dialog.getOpenFileName(None, "Select PDB file", path2, "PDB files (*.pdb)")
    		shutil.copy(self.pdbfile, path+'/protein.pdb')
    		self.sel_pdb.setText(pdb_file[0])
    	except:
    		self.sel_pdb.setText('No select file')

    def openligmol2file(self):
    	try:
    		ligmol2_dialog = QtWidgets.QFileDialog()
    		ligmol2_file = ligmol2_dialog.getOpenFileName(None, "Select MOL2 file", path2, "MOL2 files (*.mol2)")
    		shutil.copy(self.mol2file, path+'/Ligand.mol2')
    		self.sel_lig_mol2.setText(ligmol2_file[0])
    	except:
    		self.sel_lig_mol2.setText('No select file')

    def openconfmol2file(self):
    	try:
    		confmol2_dialog = QtWidgets.QFileDialog()
    		confmol2_file = confmol2_dialog.getOpenFileName(None, "Select MOL2 file", path2, "MOL2 files (*.mol2)")
    		shutil.copy(self.coffile, path+'/Cofactor.mol2')
    		self.sel_cof_mol2.setText(confmol2_file[0])
    	except:
    		self.sel_cof_mol2.setText('No select file')



    def open_backup_folder(self):
        try:
            self.bkp_folder = QtWidgets.QFileDialog.getExistingDirectory(None, 
                                                         "Select Backup folder", 
                                                         path2, 
                                                         QtWidgets.QFileDialog.ShowDirsOnly)
            self.sel_bkp_fol.setText(self.bkp_folder)
            self.path4 = self.bkp_folder
        
        except:
        	self.bkp_path['text'] = 'No select folder' 
    
    def showdialog(msgtitle,msgtxt):
    	mb = QMessageBox()
    	mb.setIcon(QMessageBox.Information)
    	mb.setWindowTitle(msgtitle)
    	mb.setText(msgtxt)
    	mb.setStandardButtons(QMessageBox.Ok)
    	mb.exec_()

    def choose_simulation_save_tpr_file(self):

      dr = str(self.sel_sav_direc.text())
      pj = str(self.sel_proj_name.text())
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
          showdialog("INFO", "Please delete backup folder or change project name and try again")
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
          showdialog("INFO", "Run Protein simulation")
          self.mount_simulation_prot()
          self.save_tprfile()
      elif find2 == False:
          showdialog("INFO","Run Protein-Ligand simulation")
          self.mount_simulation_lig()
          self.save_tprfile()
          pass
      elif find3 == False:
          showdialog("Error", "PDB file not found.")
          pass

      elif find2 == True and find1 == False:
          showdialog("Error", "It is not possible run Protein-Cofactor")
          pass

      else:
          showdialog("INFO", "Run Protein-Ligand simulation with Cofactor")
          self.mount_simulation_cof()
          self.save_tprfile()
          pass

    def choose_run_simulation(self):
      dr = str(self.sel_sav_direc.text())
      pj = str(self.sel_sav_direc.text())
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
          showdialog("INFO", "Please delete backup folder or change project name and try again")
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
          showdialog("INFO", "Run Protein simulation")
          self.mount_simulation_prot()
          self.run_simulation()
      elif find2 == False:
          showdialog("INFO","Run Protein-Ligand simulation")
          self.mount_simulation_lig()
          self.run_simulation()
          pass
      elif find3 == False:
          showdialog("Error", "PDB file not found.")
          pass

      elif find2 == True and find1 == False:
          showdialog("Error", "It is not possible run Protein-Cofactor")
          pass

      else:
          showdialog("INFO", "Run Protein-Ligand simulation with Cofactor")
          self.mount_simulation_cof()
          self.run_simulation()
          pass

    def run_bkp_file(self):
        
        if self.path4 is not None:
          print(self.path4)
          find_tpr=os.path.exists(self.path4 + '/md.tpr')
          find_cpt=os.path.exists(self.path4 + '/md.cpt')
          if find_tpr == False:
            showdialog("Error", "TPR file not found. Please try again.")
            pass
          elif find_cpt == False:
            showdialog("Error", "CPT file not found. Please try again.")
            pass
          else:
            if gromacs_flag('mdrun'):
              cmd9 = 'mdrun -s {0}/md.tpr -cpi {0}/md.cpt -append no'.format(self.path4)
              os.system(cmd9)
            elif gromacs_flag('gmx'):
              cmd9 = 'gmx mdrun -s {0}/md.tpr -cpi {0}/md.cpt -append no'.format(self.path4)
              os.system(cmd9)

        else:
          showdialog("Error", "Backup folder not found. Please try again.")
          pass


    def mount_simulation_lig(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.sel_metal.currentText())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.sel_ff.currentText())
        wt = str(self.sel_water.currentText())

        if self.ig_hyd.isChecked():
            ig = '-ignh'
        else:
            ig = None

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
            showdialog("Error", "Mol2 file not recognized.. Please try again.")
            quit()

          cm = str(self.csel_ch_met.currentText())
          nc = str(self.sel_net_ch.text())
          mt = str(self.sel_mult.text())
          at = str(self.sel_at_type.currentText())
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
            bx = str(self.sel_box.currentText())
            dst = str(self.sel_dist.text())
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

            if self.sel_ioniz.currentText() == 'Na (Number)':
                c_ion = '-np ' + self.sel_conc.text()
            elif self.sel_ioniz.currentText() == 'Cl (Number)':
                c_ion = '-nn ' + self.sel_conc.text()
            else:
                c_ion = '-conc ' + self.sel_conc.text()
            if gromacs_flag('mdrun'):
              cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)              
            elif gromacs_flag('gmx'):
              cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              
            os.system(cmd2)

            inte = self.sel_min_alg.currentText()

            nst = self.sel_min_step.text()

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

            stnvt = str(self.sel_nvt_step.text())
            stnpt = str(self.sel_npt_step.text())
            stmd = str(self.sel_md_step.text())
            temp = str(self.sel_temp.text())
            t_st = str(self.sel_int_step.text())

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

            reply = QMessageBox.question(None, 'View Complex','Is complex OK?', QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.Yes:
                pass
            else:
                showdialog('No', 'Process has been cancelled')
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
            showdialog('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            quit()
        else:
          showdialog('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

    def mount_simulation_cof(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.sel_metal.currentText())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.sel_ff.currentText())
        wt = str(self.sel_water.currentText())

        if self.ig_hyd.isChecked():
            ig = '-ignh'
        else:
            ig = None

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
            showdialog("Error", "Mol2 (ligand) file not recognized.. Please try again.")
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
            showdialog("Error", "Mol2 (cofactor) file not recognized.. Please try again.")
            quit()  
          
          cm = str(self.csel_ch_met.currentText())
          nc = str(self.sel_net_ch.text())
          mt = str(self.sel_mult.text())
          at = str(self.sel_at_type.currentText())
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

              bx = str(self.sel_box.currentText())
              dst = str(self.sel_dist.text())

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

              if self.sel_ioniz.currentText() == 'Na (Number)':
                c_ion = '-np ' + self.sel_conc.text()
              elif self.sel_ioniz.currentText() == 'Cl (Number)':
                c_ion = '-nn ' + self.sel_conc.text()
              else:
                c_ion = '-conc ' + self.sel_conc.text()

              if gromacs_flag('mdrun'):
                cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              elif gromacs_flag('gmx'):
                cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)

              os.system(cmd2)

              inte = self.sel_min_alg.currentText()
              nst = self.sel_min_step.text()

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


              stnvt = str(self.sel_nvt_step.text())
              stnpt = str(self.sel_npt_step.text())
              stmd = str(self.sel_md_step.text())
              temp = str(self.sel_temp.text())
              t_st = str(self.sel_int_step.text())

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

              reply = QMessageBox.question(None, 'View Complex','Is complex OK?', QMessageBox.Yes, QMessageBox.No)
              if reply == QMessageBox.Yes:
                pass
              else:
                showdialog('No', 'Process has been cancelled')
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
              showdialog('Error', 'ACPYPE error in Cofactor. Process has been cancelled')
              quit()
          
          else:
            showdialog('Error', 'ACPYPE error in Ligand. Process has been cancelled')
            quit()
        
        else:
          showdialog('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

    def mount_simulation_prot(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.sel_metal.currentText())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.sel_ff.currentText())
        wt = str(self.sel_water.currentText())

        if self.ig_hyd.isChecked():
            ig = '-ignh'
        else:
            ig = None

        if gromacs_flag('mdrun'):
          cmd = 'pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)
        elif gromacs_flag('gmx'):
          cmd = 'gmx pdb2gmx -ff {0} -f protein_clean.pdb -o trp.pdb -p trp.top -water {1} {2}'.format(ff, wt, ig)  
        else:
          pass

        out_pdb2gmx = os.system(cmd)

        if out_pdb2gmx == 0:
          bx = str(self.bx_menu.currentText())
          dst = str(self.dist.text())
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
          if self.sel_ioniz.currentText() == 'Na (Number)':
                c_ion = '-np ' + self.sel_conc.text()
          elif self.sel_ioniz.currentText() == 'Cl (Number)':
                c_ion = '-nn ' + self.sel_conc.text()
          else:
                c_ion = '-conc ' + self.sel_conc.text()

          if gromacs_flag('mdrun'):
            cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
          elif gromacs_flag('gmx'):
            cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
          else:
            pass

          os.system(cmd2)

          inte = self.sel_min_alg.currentText()
          nst = self.sel_min_step.text()

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

          stnvt = str(self.sel_nvt_step.text())
          stnpt = str(self.sel_npt_step.text())
          stmd = str(self.sel_md_step.text())
          temp = str(self.sel_temp.text())
          t_st = str(self.sel_int_step.text())

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

          reply = QMessageBox.question(None, 'View Complex','Is complex OK?', QMessageBox.Yes, QMessageBox.No)
          if reply == QMessageBox.Yes:
            pass
          else:
            showdialog('No', 'Process has been cancelled')
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
          showdialog('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
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
            showdialog("INFO","Ligand not found")
            quit()
        elif find2 == False:
            showdialog("INFO","Run LIE calculation")
            self.simulation_LIE_complex()
            self.simulation_LIE_lig()
            self.lie_calculation()
            showdialog("INFO","LIE calculation is finished")

        elif find3 == False:
            showdialog("Error", "PDB file not found.")
            pass

        else:
            showdialog("Error", "LIE calculation can not be performed with Cofactor.")
            quit()

    def simulation_LIE_complex(self):
        os.system("grep 'ATOM ' protein.pdb > protein_clean.pdb")
        mt = str(self.sel_metal.currentText())
        cmd0 = "grep {0} protein.pdb >> protein_clean.pdb".format(mt)

        if mt == 'None':
            pass
        else:
            os.system(cmd0)

        ff = str(self.sel_ff.currentText())
        wt = str(self.sel_water.currentText())

        if self.ig_hyd.isChecked():
            ig = '-ignh'
        else:
            ig = None

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
            showdialog("Error", "Mol2 file not recognized.. Please try again.")
            quit()

          cm = str(self.csel_ch_met.currentText())
          nc = str(self.sel_net_ch.text())
          mt = str(self.sel_mult.text())
          at = str(self.sel_at_type.currentText())
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
            bx = str(self.sel_box.currentText())
            dst = str(self.sel_dist.text())
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
                        ''')
            if gromacs_flag('mdrun'):
              os.system('grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')
            elif gromacs_flag('gmx'):
              os.system('gmx grompp -f em.mdp -c trpb4ion.pdb -p trp.top -o ion.tpr -maxwarn 1000')

            if self.sel_ioniz.currentText() == 'Na (Number)':
                c_ion = '-np ' + self.sel_conc.text()
            elif self.sel_ioniz.currentText() == 'Cl (Number)':
                c_ion = '-nn ' + self.sel_conc.text()
            else:
                c_ion = '-conc ' + self.sel_conc.text()
            if gromacs_flag('mdrun'):
              cmd2 = 'echo SOL|genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              
            elif gromacs_flag('gmx'):
              cmd2 = 'echo SOL|gmx genion -s ion.tpr -o trpb4em.pdb -neutral {} -p trp.top'.format(c_ion)
              
            os.system(cmd2)

            inte = self.sel_min_alg.currentText()

            nst = self.sel_min_step.text()

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

            stnvt = str(self.sel_nvt_step.text())
            stnpt = str(self.sel_npt_step.text())
            stmd = str(self.sel_md_step.text())
            temp = str(self.sel_temp.text())
            t_st = str(self.sel_int_step.text())

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

            reply = QMessageBox.question(None, 'View Complex','Is complex OK?', QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.Yes:
                pass
            else:
                showdialog('No', 'Process has been cancelled')
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
            showdialog('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            quit()
        else:
          showdialog('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
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
          showdialog("Error", "Mol2 file not recognized.. Please try again.")
          quit()
        
        cm = str(self.csel_ch_met.currentText())
        nc = str(self.sel_net_ch.text())
        mt = str(self.sel_mult.text())
        at = str(self.sel_at_type.currentText())
        cmdd = 'acpype -i Ligand.mol2 -c {0} -n {1} -m {2} -a {3}'.format(cm, nc, mt, at)
        acp3 = os.system(cmdd)
        if acp3 == 0:
          os.system('cp Ligand.acpype/Ligand_GMX.gro Ligand.gro')

          ff = str(self.sel_ff.currentText())
          wt = str(self.sel_water.currentText())

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

          bx = str(self.sel_box.currentText())
          dst = str(self.sel_dist.text())
          
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
          if self.sel_ioniz.currentText()  == 'Na (Number)':
              c_ion = '-np ' + self.sel_conc.text()
          elif self.sel_ioniz.currentText()  == 'Cl (Number)':
              c_ion = '-nn ' + self.sel_conc.text()
          else:
              c_ion = '-conc ' + self.sel_conc.text()

          if gromacs_flag('mdrun'):
            liecmd2 = 'echo SOL|genion -s ionl.tpr -o trpb4eml.pdb -neutral {} -p topoll.top'.format(c_ion)
          elif gromacs_flag('gmx'):
            liecmd2 = 'echo SOL|gmx genion -s ionl.tpr -o trpb4eml.pdb -neutral {} -p topoll.top'.format(c_ion)
          else:
            pass  
          os.system(liecmd2)

          lie_inte = self.sel_min_alg.currentText()

          lie_nst = self.sel_min_step.text()

          
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
          os.system('''
                            make_ndx -f eml.gro -o index_lie.ndx << EOF
                            "Water" | "Ion"
                            q
                            EOF ''')


          stnvt = str(self.sel_nvt_step.text())
          stnpt = str(self.sel_npt_step.text())
          stmd = str(self.sel_md_step.text())
          temp = str(self.sel_temp.text())
          t_st = str(self.sel_int_step.text())

          
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

          reply = QMessageBox.question(None, 'View Complex','Is complex OK?', QMessageBox.Yes, QMessageBox.No)
          if reply == QMessageBox.Yes:
            pass
          else:
            showdialog('No', 'Process has been cancelled')
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
            showdialog('Error', 'ACPYPE error (Ligand). Process has been cancelled')
            quit()
        else:
          showdialog('Error', 'PDB2GMX error (Receptor). Process has been cancelled')
          quit()

    def lie_calculation(self):
        try:
            os.chdir(path)
        except:
            pass
        if gromacs_flag('mdrun'):
          com1= '''echo 9 6 0 | g_energy -f md2.edr > out1.txt'''
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
            showdialog('Error', 'Please try to use other simulation type box.')
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
            showdialog('Error', 'Please try to use other simulation type box (Dodecahedron).')
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
          showdialog('Error', 'Please try to use other simulation type box (Dodecahedron).')
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

        dr = str(self.sel_sav_direc.text())
        pj = str(self.sel_proj_name.text())
        pj1 = pj+'_LIE'
        lie0 = '{0}.txt'.format(pj1)
        shutil.copy('out.txt', lie0)
        shutil.copy2(lie0, dr)

        a = float(match001[2])
        y = format(a, '.2f')
        self.lie_val.setText(str(y))


    def save_tprfile(self):
        os.system('chmod 777 queue.sh')
        os.system('./queue.sh')
        dr = str(self.sel_sav_direc.text())
        pj = str(self.sel_proj_name.text())
        pj1 = pj
        md0 = '{0}.tpr'.format(pj1)
        shutil.copy('md.tpr', md0)
        shutil.copy2(md0, dr)
        showdialog('Finish', 'Job has finished')
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
      pj = str(self.sel_proj_name.text())
      
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
      dr = str(self.save.text())
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
        showdialog('Finish', 'Job has finished')
        pass
      except:
        showdialog("Error", " It is not possible copy MD files.")
        quit()

    def opentprfile(self):
        try:
            tpr_dialog = QtWidgets.QFileDialog()
            tpr_file = tpr_dialog.getOpenFileName(None, "Select TPR file", path2, "TPR files (*.tpr)")
            shutil.copy(tpr_file[0], path+'/md_an.tpr')
            self.sel_tpr.setText(tpr_file[0])
        except:
            self.sel_tpr.setText('No select file')

    def openedrfile(self):
        try:
            edr_dialog = QtWidgets.QFileDialog()
            edr_file = edr_dialog.getOpenFileName(None, "Select EDR file", path2, "EDR files (*.pdb)")
            shutil.copy(edr_file[0], path+'/md_an.edr')
            self.sel_edr.setText(edr_file[0])
        except:
            self.sel_edr.setText('No select file')

    def openxtcfile(self):
        try:
            xtc_dialog = QtWidgets.QFileDialog()
            xtc_file = xtc_dialog.getOpenFileName(None, "Select PDB file", path2, "PDB files (*.pdb)")
            shutil.copy(xtc_file[0], path+'/md_an.xtc')
            self.sel_xtc.setText(xtc_file[0])
        except:
            self.sel_xtc.setText('No select file')

            #
    def RMSD(self):
        
        struct = str(self.sel_struc.currentText())
        analysis=str(self.sel_analy.currentText())


        if analysis == 'RMSD':
            try:
                os.chdir(path)
            except:
                pass
            find1=os.path.exists('md_an.tpr')
            find2=os.path.exists('md_an.xtc')
            if find1 == False:
                showdialog("Error", "TPR file not found.")
                quit()
            elif find2 == False:
                showdialog("Error", "XTC file not found.")
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
            pj = str(self.sel_proj_name.text())
            an1 = 'rmsd_{0}.xvg'.format(pj)
            shutil.copy('rmsd.xvg', an1)
            directory = str(self.sel_sav_direc.text())
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
LiGRO v 1.0 - Output of {0}
---------------------------------------------
            """.format(analysis)
            try:
              with open('rmsd_out.txt', 'w') as rmsd_infile:
                rmsd_out, _ = QFileDialog.getSaveFileName(self,"Save RMSD statistic file", os.path.join(path2+'ligro_rmsd_out.txt'), 'Text file(*.txt)')
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
                showdialog("Error", "TPR file not found.")
                pass
            elif find2 == False:
                showdialog("Error", "XTC file not found.")
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
            pj = str(self.project.text())
            an2 = 'rmsf_{0}.xvg'.format(pj)
            shutil.copy('rmsf.xvg', an2)
            directory = str(self.save.text())
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
                showdialog("Error", "TPR file not found.")
                pass
            elif find2 == False:
                showdialog("Error", "XTC file not found.")
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
            pj = str(self.project.text())
            an3 = 'rg_{0}.xvg'.format(pj)
            shutil.copy('gyrate.xvg', an3)
            directory = str(self.save.text())
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
LiGRO v 1.0 - Output of {0}
---------------------------------------------
            """.format(analysis)
            try:
              with open('rg_out.txt', 'w') as rg_infile:
                rg_out, _ = QFileDialog.getSaveFileName(self,"Save RG statistic file", os.path.join(path2+'ligro_rg_out.txt'), 'Text file(*.txt)')
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
                showdialog("Error", "TPR file not found.")
                pass
            elif find2 == False:
                showdialog("Error", "XTC file not found.")
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
            pj = str(self.project.text())
            an4 = 'msd_{0}.xvg'.format(pj)
            shutil.copy('msd.xvg', an4)
            directory = str(self.save.text())
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
LiGRO v 1.0 - Output of {0}
---------------------------------------------
            """.format(analysis)
            try:
              with open('msd_out.txt', 'w') as msd_infile:
                msd_out, _ = QFileDialog.getSaveFileName(self,"Save MSD statistic file", os.path.join(path2+'ligro_msd_out.txt'), 'Text file(*.txt)')
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
                showdialog("Error", "TPR file not found.")
                pass
            elif find2 == False:
                showdialog("Error", "XTC file not found.")
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
            pj = str(self.sel_proj_name.text())
            an6 = 'hbond_{0}.xvg'.format(pj)
            shutil.copy('hbond.xvg', an6)
            directory = str(self.save.text())
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
                showdialog("Error", "EDR file not found.")
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
              with open('out.txt') as infile, QFileDialog.getSaveFileName(self,"Save LJSR-CoulSR IE txt", os.path.join(path2+'ligro_LJSR-CoulSR_IE_out.txt'), 'Text file(*.txt)') as outfile:
                outfile.write('LiGRO v 1.0\n')
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
            showdialog("Error", "TPR or XTC file not found.")
            pass
        elif find2 == False:
            showdialog("Error", "TPR or XTC file not found.")
            pass
        else:
            pass
        timeplip = str(self.sel_plip_time.text())
        
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
          plip_out, _ = QFileDialog.getSaveFileName(self,"Save PLIP reort file", os.path.join(path2+'ligro_plip_report_out.txt'), 'Text file(*.txt)')
          plip_out.write(plip_lines)
          plip_out.close()
          pass

    def about(self):
        msgBox = QMessageBox()
        msgBox.setWindowTitle('About')
        msgBox.setText('''
LiGRO - Version 1.0
This software is available to you under the terms of the GPL-3. See ~/ligro/LICENCE for more informations. Software is created and maintained by Laboratorio de Sintese Organica Medicinal-LaSOM at Universidade Federal do Rio Grande do Sul.

Contributors:
MSc. Luciano Porto Kagami

luciano_dot_kagami_at_ufrgs_dot_br

Universidade Federal do Rio Grande do Sul Laboratorio de Sintese Organica Medicinal -LaSOM

Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

''')
        msgBox.exec_()

    def help(self):
        webbrowser.open('https://www.ufrgs.br/lasomfarmacia/node/8', new=1, autoraise=True)


    def exit(self):
        reply = QMessageBox.question(None, 'Exit','Do you really want to quit?', QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.Yes:
            quit()
        else:
            None
          
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindo = QtWidgets.QMainWindow()
    ui = Ui_MainWindo()
    ui.setupUi(MainWindo)
    MainWindo.show()
    sys.exit(app.exec_())