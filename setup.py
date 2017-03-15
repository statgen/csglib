from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

setup(
  name = "csg",
  version = "0.2",
  author_email = "welchr@umich.edu",
  author = "Ryan Welch",
  packages = [
    "csg",
    "csg/genetics",
    "csg/genome",
    "csg/notify",
    "csg/plotting",
    "csg/statistics"
  ],
  cmdclass = {"build_ext" : build_ext},
  ext_modules = cythonize([
    Extension(
      "csg.statistics.invnorm",
      sources = ["csg/statistics/invnorm.pyx"],
    ),
    Extension(
      "csg.genetics.cyalleles",
      sources = ["csg/genetics/cyalleles.pyx"],
    )
  ]),
  classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.5',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ]
)

