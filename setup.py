import sys
import os
from setuptools import setup

if (sys.version_info[0] < 2) or (sys.version_info[0] == 2 and sys.version_info[1] < 6): 
  print "Python version is %s. metAMOS requires at least 2.6"%(sys.version)
  sys.exit(1)

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "metAMOS",
    version = "0.3",
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
    scripts=['src/createProject.py', 'src/runPipeline.py'],
)
