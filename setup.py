import sys
import os
from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("preprocess", ["./src/preprocess.py"]),Extension("mapreads",["./src/mapreads.py"]),Extension("findorfs",["./src/findorfs.py"]),Extension("findscforfs",["./src/findscforfs.py"]),Extension("scaffold",["./src/scaffold.py"])]

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "metAMOS",
    version = "0.35",
    author = "Todd J. Treangen",
    author_email = "treangen+metAMOS@gmail.com",
    description = ("A metagenomic assembly all-in-one pipeline."),
    license = "GPLv3+",
    keywords = "metagenomic assembly bioinformatics",
    url = "http://github.com/treangen/metAMOS/wiki",
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
    scripts=['src/initPipeline.py', 'src/runPipeline.py'],
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext}
)
