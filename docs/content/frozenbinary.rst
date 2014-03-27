############
Single binary
############

MetAMOS PyInstaller single file binary 
======================

In attempt to further simplify the MetAMOS installation process, we are happy to announce the availability of a 'frozen' MetAMOS binary for Linux-x68_64 platforms. Along with this binary comes a significantly reduced list of prerequisites:

* Java 1.6 (or newer)
* Perl 5.8.8
	* Newer versions of Perl are not backwards compatible so if you have Perl 5.10 (or newer) you may need to install the following packages for the frozen binary to work:
		* Statistics::Descriptive
		* Bio::Seq
		* Time::Piece
		* XML::Simple
		* Storable
		* XML::Parser
		* File::Spec
* R 2.11.1 (or newer)
* 64-bit \*nix OS or Mac OSX 10.7+ (you may need to install MacPorts for full functionality)

**Disclaimer: The frozen binary is provided as-is and has limited support and reduced functionality/features compared to installing from source. If you encouter issues with a frozen binary, please try installing the latest release.**

First, select your flavor (DBs below are required but provided separately):

Linux 64-bit: 

.. code-block:: bash
    
    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/MA_fb_v1.5rc1_linux.tar.gz
    $ tar -xf MA_fb_v1.5rc1_linux.tar.gz
    $ chmod u+rwx metAMOS_v1.5rc1_binary

OSX 64-bit: 

.. code-block:: bash

    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/MA_fb_v1.5rc1_OSX.tar.gz
    $ tar -xf MA_fb_v1.5rc1_OSX.tar.gz
    $ chmod u+rwx metAMOS/*

Then add the toppings. The light DB is recommended if you are testing/getting started with metAMOS installation. If you are planning to do analysis, the full DB download is recommended. The full DB adds support for:
	* FCP classifier
	* BLAST databases required for BLAST-based classification
	* RefSeq database for required to recruit references for validation

When using the miniature DB, some features will be automatically disabled. Only Kraken can be used as the classifier with its miniature database and QUAST cannot be used as no reference will be available for recruitment. This should be a download once and only once operation. Updated frozen binaries will be backwards compatible with a previous DB download. Further details on the expected DBs on the readthedocs page.

ALL DBS: 

.. code-block:: bash

    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/allDBs.tar.gz
    $ tar -xf ftp://ftp.cbcb.umd.edu/pub/data/treangen/allDBs.tar.gz -C [$METAMOS_BIN_INSTALL_DIR]/

LIGHT DBS: 

.. code-block:: bash

    $ wget ftp://ftp.cbcb.umd.edu/pub/data/treangen/minDBs.tar.gz
    $ tar -xf ftp://ftp.cbcb.umd.edu/pub/data/treangen/minDBs.tar.gz -C	[$METAMOS_BIN_INSTALL_DIR]/

Finally, run a quick test:

.. code-block:: bash

    $ cd ./Test
    $ ./run_pipeline_test.sh

The frozen binary is actually a collection of programs that extracts/runs/cleans up automatically using `PyInstaller <http://www.pyinstaller.org/>`_. By default, PyInstaller will use the following directories to extract into:

 * The directory named by the TMPDIR environment variable.
 * The directory named by the TEMP environment variable.
 * The directory named by the TMP environment variable.

If your system is missing all of the above, does not have sufficient space, or is missing write-premissions, runPipeline will not be able to extract itself and will report: INTERNAL ERROR: cannot create temporary directory!. The extracted runPipeline requires at least 4GB of free temporary disk space. You will get a "No DBs found ERROR!" if you do not download any DBs. The DB dir needs to be placed inside of the frozen binary install dir. 

**Note: please use caution! this binaries eat up disk space quickly. Please ensure you have ample free space (100GB+) before download & use.** 

