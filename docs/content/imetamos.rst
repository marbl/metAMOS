iMetAMOS 
########

What is iMetAMOS
------------------
iMetAMOS is an extension of metAMOS to isolate genome assembly. It is a `workflow <workflows.html>`_ which, by default, uses multiple assemblers and validation tools to select the best assembly for a given sample. Effectively, this is equivalent to GAGE-in-a-box or ensemble assembly. iMetAMOS is included in the `frozen binary <frozenbinary.html>`_.

To install iMetAMOS without using a frozen binary, run:

.. code-block:: bash

    $ curl -L https://github.com/marbl/metAMOS/archive/v1.5rc1 > v1.5rc1.zip

    $ unzip v1.5rc1.zip

    $ cd metAMOS-v1.5rc1

    $ python INSTALL.py iMetAMOS

To enable iMetAMOS, specify it as an option to initPipelien using the -W flag. Below is a simple example of running of iMetAMOS to assemble an SRA dataset::

    initPipeline -q -1 SRR987657 -d projectDir -W iMetAMOS
    runPipeline -d projectDir -p 16

