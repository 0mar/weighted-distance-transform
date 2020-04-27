from distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.core import setup as npsetup

wdt = Extension(name='_wdt', sources=['fortran/mheap.f90', 'fortran/wdt_module.f90'])

with open("README.md","r") as f:
    long_description = f.read()

setuptools.setup(
        name="weighted distance transform",
        version="1.0.0",
        author="Omar Richardson",
        author_email="omsrichardson@gmail.com",
        description="Python library for computing weighted distance transforms",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/0mar/weighted-distance-transform",
        packages=['wdt'],
        classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
if __name__ == "__main__":
    npsetup(name='FORTRAN modules',
            description="FORTRAN module for the weighted distance transform",
            author="Omar Richardson",
            ext_modules=[wdt])
