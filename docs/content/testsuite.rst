############
Test suite
############

Test scripts and sanity checks
===============================

We have developed a set of scripts for testing the various features of
MetAMOS. All of these regression test scripts are available inside the
/Test directory and include all necessary datasets to run them. Here
is a brief listing of the test scripts we currently include:

Test initPipeline
-----------------

./Test/test_create.sh


Test runPipeline
----------------

./Test/run_test.sh

Test Preprocess filtration of non-interleaved fastq files
---------------------------------------------------------

./Test/test_filter_noninterleaved_fastq.sh

Test iMetAMOS
-------------

./Test/test_ima.sh

Test SRA download
-----------------

./Test/test_sra.sh

Test Newbler (if available)
---------------------------

./Test/test_newbler.sh

Test CA (fasta)
---------------

./Test/test_ca_fasta.sh

Test CA (fastq)
---------------

./Test/test_ca.sh

Test SOAPdenovo
---------------

./Test/test_soap.sh

Test MetaVelvet
---------------

./Test/test_metavelvet.sh

Test SparseAssembler
--------------------

./Test/test_sparse.sh

Test Velvet
-----------

./Test/test_velvet.sh

Test FCP
--------

./Test/test_fcp.sh

Test Spades
-----------

./Test/test_spades.sh

Test BLAST
----------

./Test/test_blast.sh

