############
Installation
############

Automated installation
======================

MetAMOS contains an automated installation script which installs
MetAMOS along with required Python dependencies, third party software
and necessary data files. 

To download the software release package:

.. code-block:: bash

    $ wget https://github.com/marbl/metAMOS/archive/Release1.5.zip

If you see a certificate not trusted error, you can add the following option to wget:

.. code-block:: bash

    $ --no-check-certificate

And if ``wget`` not available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl https://github.com/marbl/metAMOS/archive/Release1.5.zip > Release1.5.zip

You can also browse the https://github.com/marbl/MetAMOS/tree/Release1.5
and click on Downloads. 

Once downloaded, extract to unpack:

.. code-block:: bash

    $ unzip Release1.5.zip

Change to MetAMOS directory:

.. code-block:: bash

    $ cd metAMOS-Release1.5

Alternatively, you can download the INSTALL.py directly from:

.. code-block:: bash

    $ wget https://raw2.github.com/marbl/metAMOS/master/INSTALL.py

Again, if ``wget`` isn't available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl https://raw2.github.com/marbl/metAMOS/master/INSTALL.py > INSTALL.py

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
MetAMOS comes with a comprehensive test suite to make sure that installation has succeeded
on your system. 

.. code-block:: bash

    $ bash run_master_test.sh
