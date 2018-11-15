from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = cythonize("csg/**/*.pyx")
for e in extensions:
  e.include_dirs.append(numpy.get_include())

with open("requirements.txt","rt") as f:
  requirements = f.read().split()

setup(
  name = "csg",
  version = "0.2.1",
  author_email = "",
  author = "Center for Statistical Genetics",
  packages = [
    "csg",
    "csg/genetics",
    "csg/genome",
    "csg/notify",
    "csg/plotting",
    "csg/statistics",
    "csg/intervaltree",
    "csg/pedigree",
    "csg/gwas"
  ],
  classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.5",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
  ],
  install_requires = requirements,
  ext_modules = extensions,
  zip_safe = False,
  scripts = [
    "csg/gwas/Harmonize/Harmonize",
    "csg/gwas/TransformEffect",
    "csg/gwas/Filter",
    "csg/gwas/AnnotateRs",
    "csg/gwas/SignificantLoci",
    "csg/gwas/GCcorrect",
    "csg/gwas/IndependentHits",
    "csg/gwas/InflationFromFile",
    "csg/pedigree/trios/TrioPuller"
  ],
)
