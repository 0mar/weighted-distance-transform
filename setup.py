from distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.core import setup as npsetup

prop_dist = Extension(name='wdt_fortran', sources=['fortran/propagate_distance.f90'])

if __name__ == "__main__":
    npsetup(name='FORTRAN modules',
            description="FORTRAN module for the weighted distance transform script",
            author="Omar Richardson",
            ext_modules=[prop_dist])
