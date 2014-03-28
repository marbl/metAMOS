iMetAMOS 
########

What is iMetAMOS
------------------
iMetAMOS is an extension of metAMOS to isolate genome assembly. It is a `workflow <workflows.html>`_ which, by default, uses multiple assemblers and validation tools to select the best assembly for a given sample. Effectively, this is equivalent to GAGE-in-a-box or ensemble assembly. iMetAMOS is included in the `frozen binary <frozenbinary.html>`_.

If you have used iMetAMOS for analyzing your, please cite (in addition to the individual software component citations listed in main output):
   Koren S, Treangen TJ, Hill CM, Pop M, Phillippy AM
   Automated ensemble assembly and validation of microbial genomes
   bioRxivdoi: 10.1101/002469

Please also consider citing the original metAMOS publication:
    Treangen TJ\*, Koren S\*, Sommer DD, Liu B, Astrovskaya I, Ondov B, Darling AE, Phillippy AM, Pop M.
    MetAMOS: a modular and open source metagenomic assembly and analysis pipeline.
    Genome Biol. 2013 Jan 15;14(1):R2. PMID: 23320958.

*Indicates both authors contributed equally to this work


To install iMetAMOS without using a frozen binary, run:

.. code-block:: bash

    $ curl -L https://github.com/marbl/metAMOS/archive/v1.5rc1 > v1.5rc1.zip

    $ unzip v1.5rc1.zip

    $ cd metAMOS-v1.5rc1

    $ python INSTALL.py iMetAMOS

To enable iMetAMOS, specify it as an option to initPipelien using the -W flag. Below is a simple example of running of iMetAMOS to assemble an SRA dataset::

    initPipeline -q -1 SRR987657 -d projectDir -W iMetAMOS
    runPipeline -d projectDir -p 16

