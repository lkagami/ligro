#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import numpy as np
from Cython.Build import cythonize
from setuptools import setup, Extension, find_packages

extensions = [
    Extension("vismol.utils.matrix_operations",
              sources=["src/vismol/utils/matrix_operations.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.PDBFiles",
              sources=["src/vismol/utils/PDBFiles.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.GROFiles",
              sources=["src/vismol/utils/GROFiles.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.c_distances",
              sources=["src/vismol/utils/c_distances.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.selectors",
              sources=["src/vismol/utils/selectors.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.mesh_decimation",
              sources=["src/vismol/utils/mesh_decimation.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.mesh_smoothing",
              sources=["src/vismol/utils/mesh_smoothing.pyx"],
              include_dirs=[np.get_include()]),
    #Extension("vismol.utils.cartoon",
    #          sources=["src/vismol/utils/cartoon.pyx"],
    #          include_dirs=[np.get_include()]),
]

setup(
    python_requires=">=3.7",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={
        "vismol.gui.gtk_widgets": ["*.glade"],
        "vismol.libgl.fonts": ["*.ttf"],
        # "gui": ["gtk_widgets/*.glade"],
        # "libgl": ["fonts/*.ttf"],
    },
    scripts=["scripts/vismol_simple.py"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: Creative Commons Attribution-NonCommercial 4.0",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research"
        ],
    ext_modules=cythonize(extensions, compiler_directives={"language_level" : "3"}),
    )
