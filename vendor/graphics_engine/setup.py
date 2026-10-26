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

with open("LICENSE", "r") as f:
    license = f.read()

setup(
    name="vismol",
    version=__import__("src").__version__,
    description="Python molecular viewer library",
    url="https://github.com/casebor/graphics_engines",
    project_urls={
        "Source": "https://github.com/casebor/graphics_engines"
        },
    python_requires=">=3.7",
    author=__import__("src").__author__,
    author_email=__import__("src").__mail__,
    license=license,
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={
        "vismol.gui.gtk_widgets": ["*.glade"],
        "vismol.libgl.fonts": ["*.ttf"],
        # "gui": ["gtk_widgets/*.glade"],
        # "libgl": ["fonts/*.ttf"],
    },
    scripts=["scripts/vismol_simple.py"],
    install_requires=[
        "numpy>=1.20.0",
        "Cython>=0.29.12",
        "PyOpenGL>=3.1.0",
        "PyQt5>=5.15"
        ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: Creative Commons Attribution-NonCommercial 4.0",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research"
        ],
    ext_modules=cythonize(extensions, compiler_directives={"language_level" : "3"}),
    )
