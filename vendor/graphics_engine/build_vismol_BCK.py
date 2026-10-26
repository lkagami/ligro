#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("vismol.utils.matrix_operations",
              sources=["vismol/utils/matrix_operations.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.PDBFiles",
              sources=["vismol/utils/PDBFiles.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.GROFiles",
              sources=["vismol/utils/GROFiles.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.c_distances",
              sources=["vismol/utils/c_distances.pyx"],
              include_dirs=[np.get_include()]),
    Extension("vismol.utils.cartoon",
              sources=["vismol/utils/cartoon.pyx"],
              include_dirs=[np.get_include()]),
]

setup(
    name="vismol",
    ext_modules=cythonize(extensions, compiler_directives={"language_level" : "3"}),
)
