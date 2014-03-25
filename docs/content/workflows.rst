############
Workflows
############

Workflows (and common use cases) 
===============================

What is a workflow?
------------------
Good question! A workflow is a text-file that specified command-line options and input sequences
required to run metAMOS. A workflow may optionally inherit options/data from other
workflows. A workflow may also be immutable if the parameters should not be modifiable
by a user. 

Example workflow
----------------

An example workflow::

    inherit:            isolate
    modify:             True
    command:            -q -u -r -v -I -c kraken -p 16 -a spades,velvet-sc,abyss,ray,edena,sga,masurca,soap,soap2,velvet -t metamos -n FunctionalAnnotation -f Postprocess -z phylum 
    asmcontigs:         /Users/skoren/Personal/Research/metAMOS/Test/test.asm,ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Candidatus_Carsonella_ruddii_uid58773/NC_008512.fna 
    lib1format:         fasta
    lib1mated:          True
    lib1innie:          True
    lib1interleaved:	True
    lib1f1:             /Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe_filt.fna.gz,2000,5000,3500,500

Available options
-----------------

The available options are:

*   inherit   - any other workflows to inherit from. In this case, the workflow inherits options from the isolate workflow
*   modify    - whether users are allowed to specify command-line parameters at runtime. If false, command-line options are ignored
*   command   - command-line options to specify for runPipeline
*   asmcontigs - optional, pre-assembled contigs to include in analysis. Can be remote file. Multiple files can be separated using commas.
*   lib#format - input type for lib #. Can be fasta/fastq/sff
*   lib#mated - whether the library is mated or not
*   lib#innie - whether the mates are in the innie (Illumina paired-end) format or not (Illumina mate-pair)
*   lib#interleaved - whether the input sequences are in a single file or in two separate files
*   lib#f1 	    - the name of the input file, along with library min, max, mean, stdev

An arbitrary number of libraries may be specified in the above format. The below example shows an unmated library::

    lib1format:	        fasta
    lib1mated:		False
    lib1innie:		False
    lib1interleaved:	False
    lib1frg:		/Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe_filt.fna.gz


as well as a non-interleaved library::


    lib1format:		fasta
    lib1mated:      	True
    lib1innie:      	True
    lib1interleaved:    False
    lib1f1: 		/Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe_1.fna.gz,2000,5000,3500,500
    lib1f2: 		/Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe.2.fna.gz,2000,5000,3500,500

Sharing your favorite MA workflows with others
----------------------------------------------
Workflows may be shared between users, as long as the input files are accessible (i.e. they are on a remote server or the systems share a file system). Workflow files should be placed in the metAMOS/workflows directory or the working directory where MetAMOS is launched.

