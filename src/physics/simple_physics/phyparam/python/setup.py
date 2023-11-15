from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup( name = "Wrapper to simple physics package",
       ext_modules = cythonize([ Extension("phyparam",["phyparam.pyx"], 
                                           libraries=['phyparam'],
                                           library_dirs=['./lib'])] )
       )
