#!/usr/bin/env python
#
# File: setup.py
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zongxiao He"
__date__ = "11/5/15"


from distutils.core import setup

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

import sys, os

try:
    import argparse, scipy
except ImportError:
    sys.exit('The program requires Python 2.7 + SciPy 0.16 or Anaconda 2.3.')
#
setup(name = 'RVGDT',
    version='1.0',
    description='Rare Variant Extension of the Generalized Disequilibrium Test',
    author='Zongxiao He',
    author_email='zongxiah@bcm.edu',
    url='XXX',
    package_dir={'RVGDT':'source'},
    packages=['RVGDT'],
    scripts = ['rvgdt'],
    cmdclass = {'build_py': build_py}
    )