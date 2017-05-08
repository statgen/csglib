
from distutils.core import setup
from distutils.extension import Extension

def build_ext(*args, **kwargs):
  from Cython.Distutils import build_ext
  return build_ext(*args, **kwargs)

setup(
  name = "csg",
  version = "0.2",
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
    "csg/pedigree"
  ],
  cmdclass = {"build_ext" : build_ext},
  ext_modules = [
    Extension(
      "csg.statistics.invnorm",
      sources = ["csg/statistics/invnorm.pyx"],
    ),
    Extension(
      "csg.genetics.cyalleles",
      sources = ["csg/genetics/cyalleles.pyx"],
    )
  ],
  classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.5',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ],
  setup_requires = [
    'setuptools>=18.0',
    'cython',
    'numpy',
  ],
)


# -*- python-indent-offset: 2 -*-
