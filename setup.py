from platform import version
from setuptools import setup, Extension

added_sources = []
added_sources += ['res/wam.c','res/ddg.c','res/lnorm.c','res/jacobi.c']
added_sources += ['res/spk.c','res/fit.c']
added_sources += ['res/eigenpap.c','res/tools.c']
toload = ['spkmeansmodule.c'] + added_sources
setup(name = 'spkmeans',
version='1.0',
description='c-api for spkmeans final project',
ext_modules=[Extension('spkmeans', sources = toload)])
