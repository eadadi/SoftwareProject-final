from platform import version
from setuptools import setup, Extension

setup(name = 'spkmeansmodule',
version='1.0',
description='c-api for spkmeans final project',
ext_modules=[Extension('spkmeansmodule', sources=['spkmeans.c'])])
