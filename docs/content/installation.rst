############
Installation
############

Before your start
=================
The most common cause of a failed run is a missing package or dependency. We provide two paths to simplify the task of
downloading and manual install of all required dependencies: an install script (INSTALL.py) and a frozen binary.
See below sections for further details.

Prerequisites
==============
MetAMOS has several dependencies/prerequisites, the large majority of which are automatically downloaded
and installed when running INSTALL.py (see next section). In addition, several dependences/prerequisites 
are not installed by INSTALL.py and must be available on your system:

    * perl (5.8.8+)
    * python (2.7.3+)
    * R (2.11.1+ with PNG support)
    * gcc (4.7+ for full functionality)
    * curl 
    * wget

Here is a list of currently supported Operating Systems:

1. Mac OSX (10.7 or newer)
2. Linux 64-bit (tested on CentOS, Fedora, RedHat, OpenSUSE and Ubuntu)

And here is a current list of required packages/libs installed by metAMOS:

1. Perl 
    * File::Copy::Link
    * Statistics::Descriptive 
    * Time::Piece
    * XML::Parser
    * XML::Simple
    
2. Python
    * Cython
    * matplolib
    * NumPy
    * psutil
    * PySam
    * setuptools

3. Other
    * Boost
    * cmake
    * Jellyfish
    * SparseHash

Automated installation
======================

MetAMOS contains an automated installation script which installs
MetAMOS along with required Python dependencies, third party software
and necessary data files. If you encounter issues during installation, you can 
try installing the required dependencies manually and re-running INSTALL.py. If you 
continue to encounter issues, plase provide the output from INSTALL.py as a new `issue <https://github.com/marbl/metAMOS/issues?state=open>`. 

To download the software release package:

.. code-block:: bash

    $ wget https://github.com/marbl/metAMOS/archive/v1.5rc2.zip

If you see a certificate not trusted error, you can add the following option to wget:

.. code-block:: bash

    $ --no-check-certificate

And if ``wget`` not available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl -L https://github.com/marbl/metAMOS/archive/v1.5rc2.zip > v1.5rc2.zip

You can also browse the https://github.com/marbl/MetAMOS/tree/v1.5rc2
and click on Downloads. 

Once downloaded, extract to unpack:

.. code-block:: bash

    $ unzip v1.5rc2.zip

Change to MetAMOS directory:

.. code-block:: bash

    $ cd metAMOS-v1.5rc2

Once inside the MetAMOS directory, run:

.. code-block:: bash

    $ python INSTALL.py core

This will download and install the external dependencies which may
take minutes or hours to download depending on your connection speed.
metAMOS supports workflows to install subsets of tools for faster installation.
By default only the core dependencies are installed. 

To install iMetAMOS run:

.. code-block:: bash

    $ python INSTALL.py iMetAMOS


Also, you can run:

.. code-block:: bash

    $ python INSTALL.py -h

to get a listing of available workflows and programs. You can specify either
workflows or programs as arguments to INSTALL.py. For example, to install the
core workflow plus PhyloSift, run:

.. code-block:: bash

    $ python INSTALL.py core phylosift


To install the programs which are part of the optional workflow run:

.. code-block:: bash

    $ python INSTALL.py optional


If all dependencies are downloaded (including optional/deprecated ones), this will take
quite awhile to complete (plan on a few hours to 2 days).

Running the test suite
===========================
MetAMOS comes with a comprehensive `test suite <testsuite.html>`_ to make sure that installation has succeeded
on your system. To run a quick test and very installation succeeded run: 

.. code-block:: bash

    $ cd ./Test
    $ ./run_pipeline_test.sh

