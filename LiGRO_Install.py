#!/usr/bin/python
"""
---LiGRO - Trial version---

Luciano Porto Kagami
Universidade Federal do Rio Grande do Sul
Laboratorio de Sintese Organica Medicinal -LaSOM
Av. Ipiranga, 2752 - Azenha, Porto Alegre - RS, 90610-000 - Brazil

-------------------
"""
title = 'LiGRO: Version 0.1 - Install'
import sys

sys.path[:0] = ['../../..']
import Tkinter, Tkconstants, tkFileDialog
from Tkinter import *
import Pmw
import os
import shutil
import webbrowser
import tkMessageBox as mbox
import tempfile
import subprocess


try:
    path = tempfile.mkdtemp()

except:
    pass
path2 = os.environ.get('HOME')

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
        try:
            os.chdir(path)
        except:
            pass
        # verificar se AMBERTools esta instalado...
        if is_tool('antechamber') == False:
            find1 = os.path.exists('AmberTools16.tar.bz2')
            if find1 == False:
                mbox.showerror("Error", "AmberTools16.tar.bz2 file not found.")
                quit()
            else:
                dist = str(distro_menu.getcurselection())
                if dist == 'Ubuntu':
                    dt = 'apt-get'
                elif dist == 'OpenSuse':
                    dt = 'zypper'
                command0='sudo {0} install make csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev patch python-tk python-matplotlib'.format(dt)
                os.system(command0)
                command1='tar xvfj AmberTools16.tar.bz2'
                os.system(command1)
                shutil.copy2(path+'/amber16', path2)
                os.putenv('AMBERHOME',path2+'/amber16')
                path3=path2+'/amber16'
                os.chdir(path3)
                command3='source amber.sh'
                os.system(command3)
                command4='make install'
                os.system(command4)
                os.chdir(path)
                bashrc = open('.basrc', 'a')
                bashrc.write('source amber.sh')
                bashrc.close()
        else: pass


    def install_acpype():
        dist = str(distro_menu.getcurselection())
        if dist == 'Ubuntu':
            dt = 'apt-get'
        elif dist == 'OpenSuse':
            dt = 'zypper'
        if is_tool('acpype') == False:
            command5='sudo {0} install git'.format(dt)
            os.system(command5)
            command6='git clone https://github.com/t-/acpype.git'
            os.system(command6)
            path4=path2+'/acpype-master'
            os.chdir(path4)
            command7='ln -s $PWD/acpype.py /usr/local/bin/acpype'
            os.system(command7)
        else: pass

    def install_plip():
        dist = str(distro_menu.getcurselection())
        if dist == 'Ubuntu':
            dt = 'apt-get'
        elif dist == 'OpenSuse':
            dt = 'zypper'
        if is_tool('plip') == False:
            command8='sudo {0} install python-openbabel git'.format(dt)
            os.system(command8)
            command9='git clone https://github.com/ssalentin/plip.git ~/pliptool'
            os.system(command9)
            os.chdir(path)
            command13 = '''echo "alias plip='{0}/pliptool/plip/plipcmd'" >> .basrc'''.format(path2)
            os.system(command13)
        else: pass


    def install_ligro():
        dist = str(distro_menu.getcurselection())
        if dist == 'Ubuntu':
            dt = 'apt-get'
        elif dist == 'OpenSuse':
            dt = 'zypper'
        command10='sudo {0} install python-webbroser python-tempfile python-matplotlib python-numpy pymol gedit'.format(dt)
        os.system(command10)
        os.chdir(path)
        if is_tool('pdb2gmx') == True:
            command14 = '''echo "alias ligro='{0}/ligro/ligro_gromacs_v4x.py'" >> .basrc'''.format(path2)
            os.system(command14)
        elif is_tool('gmx'):
            command15 = '''echo "alias ligro='{0}/ligro/ligro_gromacs_v5x.py'" >> .basrc'''.format(path2)
            os.system(command15)

    def install_gromacs():
        dist = str(distro_menu.getcurselection())
        if dist == 'Ubuntu':
            dt = 'apt-get'
        elif dist == 'OpenSuse':
            dt = 'zypper'
        if is_tool('gmx') or is_tool('pdb2gmx') == False:
            command12='sudo {0} install gromacs'.format(dt)
            os.system(command12)
        else: pass

    def install_all():
        install_ambertools()
        install_acpype()
        install_gromacs()
        #install_ligro()
        mbox.showinfo('Finish', 'Installation has finished, type ligro for to execute')
        sys.exit(0)


    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title(title)
    Tkinter.Label(root, text="LiGRO: Protein-Ligand Molecular Dynamics using ACPYPE, GROMACS and PLIP.").grid(row=0, column=1)
    widget = root
    root.geometry("890x300+260+150")
    menubar = Tkinter.Menu(root)
    filemenu = Tkinter.Menu(menubar, tearoff=0)
    filemenu.add_separator()
    editmenu = Tkinter.Menu(menubar, tearoff=0)
    editmenu.add_separator()
    helpmenu = Tkinter.Menu(menubar, tearoff=0)

    Tkinter.Label(root, text='AmberTools16.tar.bz2 file:').grid(row=3, column=0)
    amb1 = Tkinter.Label(root, text='No select file', bg='white', padx=100, pady=5)
    amb1.grid(row=3, column=1, padx=10, pady=10)
    brw1_bt = Tkinter.Button(root, text='Browse...', command=openamberfile)
    brw1_bt.grid(row=3, column=2)
    distro_menu = Pmw.OptionMenu(root,
                                    labelpos='w',
                                    label_text='Select your distro:',
                                    menubutton_textvariable=None,
                                    items=['Ubuntu', 'OpenSuse'],
                                    menubutton_width=10,
                                    )
    distro_menu.grid(row=4, column=1)
    scrollbar = Scrollbar(root)
    scrollbar.grid(row=1, column=3, sticky=W)
    window_txt = Text(root, width=106, height=8, font="Courier 10", yscrollcommand=scrollbar.set)
    window_txt.grid(row=1, column=0, columnspan=3, sticky=W)
    scrollbar.config(command=window_txt.yview)

    my_opening_message0 = """
               This software is available to you under the terms of the GPL-3.
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
    CamcelButton = Tkinter.Button(root, text='Cancel', command=callback).grid(row=4, column=0)
    nextButton = Tkinter.Button(root, text='Next', command=install_all).grid(row=4, column=2)
    img = PhotoImage(file='ic_launcher_.png')
    root.tk.call('wm', 'iconphoto', root._w, img)
    root.config(menu=menubar)
    root.mainloop()