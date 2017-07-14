#!/usr/bin/python
"""
---LiGRO - version 0.2---

This software is available to you under the terms of the GPL-3. See ~/ligro/LICENCE for more informations.
               Software is created and maintained by Laboratorio de Sintese Organica Medicinal-LaSOM at
               Universidade Federal do Rio Grande do Sul.
               Contributors:
               Luciano Porto Kagami
               luciano_dot_kagami_at_ufrgs_dot_br
               Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

               For Antechamber, please cite:

               1.  Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and
               bond type perception in molecular mechanical calculations". Journal of
               Molecular Graphics and Modelling , 25, 2006, 247260.
               2.  Wang, J., Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A.
               "Development and testing of a general AMBER force field". Journal of
               Computational Chemistry, 25, 2004, 1157-1174.

               For ACPYPE, please cite:

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

               Salentin,S. et al. PLIP: fully automated protein-
               ligand interaction profiler. Nucl. Acids Res. (1 July 2015) 43 (W1):
               W443-W447. doi: 10.1093/nar/gkv315

               For GROMACS please cite:

               For this code plese cite:
               Manuscript to be submitted.

-------------------
"""
title = 'LiGRO: Version 0.2 - Install'

print('''
---LiGRO - version 0.2---

    LiGRO INSTALLER

-------------------------

''')

a = input("Select 1 for Ubuntu distro or select 2 for OpenSuse distro:       " )

if a == 1:
    dt = 'apt-get'
elif a == 2:
    dt = 'zypper'
else:
    print("Plese, choose a valid number")
    quit()

import sys
import os
command000="sudo {0} install python-Pmw python-tk python-pandas ".format(dt)
os.system(command000)

sys.path[:0] = ['../../..']
import Tkinter, Tkconstants, tkFileDialog
from Tkinter import *
import Pmw
import shutil
import webbrowser
import tkMessageBox as mbox
import subprocess

path2 = os.environ.get('HOME')
os.system('mkdir '+path2+'/ligro_temp')

try:
    path=path2 + '/ligro_temp'

except:
    pass


# Create root window.

if __name__ == '__main__':

    def openamberfile():
        try:
            amberfile = tkFileDialog.askopenfilename(initialdir=path2, title="Select AMBER TAR.BZ2 file",
                                                        filetypes=(("TAR.BZ2 files", "*.tar.bz2"), ("all files", "*.*")))
            amb1['text'] = amberfile
            shutil.copy(amberfile, path)
        except:
            amb1['text'] = 'No select file'


    def callback():
        if mbox.askyesno('Cancel', 'Really cancel?'):
            sys.exit(0)
        else:
            None


    def help():
        webbrowser.open('https://www.ufrgs.br/lasomfarmacia/', new=1, autoraise=True)


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

    def is_tool(name):
        try:
            devnull = open(os.devnull)
            subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                return False
        return True


    def install_ambertools():
            find1 = os.path.exists('AmberTools16.tar.bz2')
            if find1 == False:
                mbox.showerror("Error", "AmberTools16.tar.bz2 file not found.")
                quit()
            else:
                os.system(
                    'sudo {0} install make csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev patch python-tk python-matplotlib'.format(
                        dt))
                os.system('tar xvfj AmberTools16.tar.bz2')
            os.system('mv {0}/amber16 {1}'.format(path, path2))
            os.putenv('AMBERHOME', path2 + '/amber16')
            path3 = path2 + '/amber16'
            os.chdir(path3)
            os.system('./configure gnu')
            os.system('source amber.sh')
            os.system('make install')
            os.chdir(path2)
            with open('.bashrc', 'a') as file:
                command40 = "source {0}/amber16/amber.sh\n".format(path2)
                file.write(command40)
            file.close()

    def install_acpype():
            command6 = 'git clone https://github.com/t-/acpype.git'
            os.system(command6)
            path4 = path2 + '/acpype'
            os.chdir(path4)
            command7 = 'sudo ln -s $PWD/acpype.py /usr/local/bin/acpype'
            os.system(command7)

    def install_plip():
            command8='sudo {0} install python-openbabel'.format(dt)
            os.system(command8)
            command9='git clone https://github.com/ssalentin/plip.git ~/pliptool'
            os.system(command9)
            os.chdir(path2)
            with open('.bashrc', 'a') as file:
                command40 = "alias plip='{0}/pliptool/plip/plipcmd\n'".format(path2)
                file.write(command40)
            file.close()

    def install_gromacs():
            command10 = 'sudo {0} install gromacs'.format(dt)
            os.system(command10)

    def install_ligro4x():
        command10='sudo {0} install python-matplotlib python-numpy pymol gedit'.format(dt)
        os.system(command10)
        os.chdir(path2)
        with open('.bashrc', 'a') as file:
            command40 = "alias ligro='python {0}/ligro/ligro_gromacs_v4x.py'\n".format(path2)
            file.write(command40)
        file.close()


    def install_ligro5x():
        command10='sudo {0} install python-matplotlib python-numpy pymol gedit'.format(dt)
        os.system(command10)
        os.chdir(path2)
        with open('.bashrc', 'a') as file:
            command40 = "alias ligro='python {0}/ligro/ligro_gromacs_v5x.py'\n".format(path2)
            file.write(command40)
        file.close()

    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title(title)
    Tkinter.Label(root, text="LiGRO: Protein-Ligand Molecular Dynamics using ACPYPE, GROMACS and PLIP.").grid(row=0, column=1)
    widget = root
    root.geometry("890x400+260+150")
    menubar = Tkinter.Menu(root)
    filemenu = Tkinter.Menu(menubar, tearoff=0)
    filemenu.add_separator()
    editmenu = Tkinter.Menu(menubar, tearoff=0)
    editmenu.add_separator()
    helpmenu = Tkinter.Menu(menubar, tearoff=0)

    w = Pmw.Group(root, tag_text='Install Packeges')
    w.grid(row=2, column=1)
    amb1 = Tkinter.Label(w.interior(), text='No select file', bg='white', padx=100, pady=5)
    amb1.grid(row=3, column=1, padx=10, pady=10)

    var1 = BooleanVar()
    var2 = BooleanVar()
    var3 = BooleanVar()
    var4 = BooleanVar()

    b1=Tkinter.Checkbutton(w.interior(), text='   AmberTools', variable=var1).grid(row=3, column=0, padx=10, pady=10)
    b2=Tkinter.Checkbutton(w.interior(), text='         ACPYPE', variable=var2).grid(row=4, column=0, padx=10, pady=10)
    b3=Tkinter.Checkbutton(w.interior(), text='             PLIP', variable=var3).grid(row=5, column=0, padx=10, pady=10)
    b3=Tkinter.Checkbutton(w.interior(), text='     GROMACS', variable=var4).grid(row=6, column=0, padx=10, pady=10)
    brw1_bt = Tkinter.Button(w.interior(), text='Browse...', command=openamberfile)
    brw1_bt.grid(row=3, column=2)
    scrollbar = Scrollbar(root)
    scrollbar.grid(row=1, column=3, sticky=W)
    window_txt = Text(root, width=106, height=8, font="Courier 10", yscrollcommand=scrollbar.set)
    window_txt.grid(row=1, column=0, columnspan=3, sticky=W)
    scrollbar.config(command=window_txt.yview)

    def install_all():
        command5 = 'sudo {0} install git'.format(dt)
        os.system(command5)
        command9 = 'git clone https://github.com/lkagami/ligro.git ~/ligro'
        os.system(command9)
        if var1.get() == True:
           install_ambertools()
        elif var1.get() == 0:
            pass

        if var2.get() == True:
           install_acpype()
        elif var2.get() == 0:
            pass

        if var3.get() == True:
           install_plip()
        elif var3.get() == 0:
            pass

        if var4.get() == True:
           install_gromacs()
        elif var4.get() == 0:
            pass
        if is_tool('gmx') == False:
            pass
        else:
            install_ligro5x()
        if is_tool('mdrun') == False:
            pass
        else:
            install_ligro4x()
        os.system('sudo rm -r {0}'.format(path))
        os.chdir(path2)
        os.system('. ~/.bashrc')
        mbox.showinfo('Finish', 'Installation has finished, type ligro for to execute')
        sys.exit(0)
    my_opening_message0 = """
               This software is available to you under the terms of the GPL-3. See ~/ligro/LICENCE for more informations.
               Software is created and maintained by Laboratorio de Sintese Organica Medicinal-LaSOM at
               Universidade Federal do Rio Grande do Sul.
               Contributors:
               Luciano Porto Kagami
               luciano_dot_kagami_at_ufrgs_dot_br
               Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

               For Antechamber, please cite:

               1.  Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and
               bond type perception in molecular mechanical calculations". Journal of
               Molecular Graphics and Modelling , 25, 2006, 247260.
               2.  Wang, J., Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A.
               "Development and testing of a general AMBER force field". Journal of
               Computational Chemistry, 25, 2004, 1157-1174.

               For ACPYPE, please cite:

               SOUSA DA SILVA, A. W.; VRANKEN, W. F.; LAUE, E. D. ACPYPE - AnteChamber
               PYthon Parser Interface. Manuscript to be submitted.

               Alan Wilter S. da Silva, D.Sc. - CCPN Research Associate
               Department of Biochemistry, University of Cambridge.
               80 Tennis Court Road, Cambridge CB2 1GA, UK.
               >>http://www.bio.cam.ac.uk/~awd28<<

               alanwilter _at_ gmail _dot_ com

               For PLIP please cite:

               Salentin,S. et al. PLIP: fully automated protein-
               ligand interaction profiler. Nucl. Acids Res. (1 July 2015) 43 (W1):
               W443-W447. doi: 10.1093/nar/gkv315

               For GROMACS please cite:

               For this code plese cite:
               """

    window_txt.insert(0.0, my_opening_message0)

    helpmenu.add_command(label="Help Index", command=help)
    helpmenu.add_command(label="About...", command=about)
    menubar.add_cascade(label="Help", menu=helpmenu)
    CamcelButton = Tkinter.Button(root, text='Cancel', command=callback).grid(row=6, column=0)
    nextButton = Tkinter.Button(root, text='Next', command=install_all).grid(row=6, column=2)
    root.config(menu=menubar)
    root.mainloop()
