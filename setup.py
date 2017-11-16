from distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.core import setup as npsetup

wdt = Extension(name='_wdt', sources=['fortran/mheap.f90', 'fortran/wdt_module.f90'])

if __name__ == "__main__":
    npsetup(name='FORTRAN modules',
            description="FORTRAN module for the weighted distance transform",
            author="Omar Richardson",
            ext_modules=[wdt])
