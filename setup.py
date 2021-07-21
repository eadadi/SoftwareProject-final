from platform import version
from setuptools import setup, Extension

added_sources = ['res\ddg.c','res\jacobi.c','res\lnorm.c','res\spk.c','res\wam.c']
toload = ['spkmeansmodule.c'] + added_sources
setup(name = 'spkmeansm',
version='1.0',
description='c-api for spkmeans final project',
ext_modules=[Extension('spkmeans', sources = toload)])
