############
Single binary
############

MetAMOS PyInstaller single file binary 
======================

In attempt to further simplify the MetAMOS installation process, we are happy to announce the availability of a 'frozen' MetAMOS binary for Linux-x68_64 platforms. Along with this binary comes a significantly reduced list of prerequisites:

* Java 1.6 (or newer)
* Perl 5.8.8 (or newer)
* 64-bit *nix OS or Mac OSX 10.7+ (you may need to install MacPorts for full functionality)

**Disclaimer: This is provided as-is, has limited support and reduced functionality/features compared to installing from source.**

First, select your flavor (DBs below are required but provided separately):

Linux 64-bit: 

.. code-block:: bash
    
    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/MA_fb_v1.5rc1_linux.tar.gz
    $ tar -xf MA_fb_v1.5rc1_linux.tar.gz

OSX 64-bit: 

.. code-block:: bash

    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/MA_fb_v1.5rc1_OSX.tar.gz
    $ tar -xf MA_fb_v1.5rc1_OSX.tar.gz

Then add the toppings.

ALL DBS: 

.. code-block:: bash

    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/allDBs.tar.gz
    $ cd [$METAMOS_BIN_INSTALL_DIR]
    $ tar -xf ftp://ftp.cbcb.umd.edu/pub/data/treangen/allDBs.tar.gz -C [$METAMOS_BIN_INSTALL_DIR]/DB

LIGHT DBS: 

.. code-block:: bash

    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/minDBs.tar.gz
    $ cd [$METAMOS_BIN_INSTALL_DIR]
    $ tar -xf ftp://ftp.cbcb.umd.edu/pub/data/treangen/minDBs.tar.gz -C	[$METAMOS_BIN_INSTALL_DIR]/DB

Finally, run a quick test:

.. code-block:: bash

    $ cd ./Test
    $ ./run_test.sh


This will take a moment to extract and run. You will get a "No DBs found ERROR!" if you do not download any DBs. The DB dir needs to be placed inside of the frozen binary install dir. We recommend you try the Light version first, and if you need the extended DBs (for BLAST, FCP, QUAST, etc) grab the "ALL DBS" tarball. This should be a download once and only once operation. In addition, your existing DBs could work (assuming they are same format, etc). Further details on the expected DBs on the readthedocs page.

**Note: please use caution! this binaries eat up disk space quickly. Please ensure you have ample free space (100GB+) before download & use. 

