############
Output
############

.. raw:: html

    <style> .red {background-color:red} </style>
    <style> .yellow {background-color:yellow} </style>
    <style> .green {background-color:lightgreen} </style>

Full listing of expected output files
===============================

MetAMOS generates an interactive web page once a run successfully completes::

     http://www.cbcb.umd.edu/~sergek/imetamos/gageb/Postprocess/out/html/summary.html

This includes summary statistics and taxonomic information based on Krona [1].
The easiest way to interact with the results is through the web interface. The web interface has
been tested in several browsers. The currently known issues are:

.. role:: yellow
.. role:: red
.. role:: green
==================  ==============  ====================================
Browser             Version         Issues
==================  ==============  ====================================
:yellow:`Chrome`     33.0.1750.152   MapReads and Assemble steps do not show
:green:`Safari`      6.1.2           None
:yellow:`Firefox`    28              QUAST reports do not show for Validate
:red:`IE`            9               Not Tested
==================  ==============  ====================================

The Postprocess/out directory contains the results of the analysis. By default, 
metAMOS uses the prefix "proba" (Galician for test). Thus, files will have the name "proba".*.

abundance.krona.html
--------------------
    
Krona [1] plot of abundances using the tool selected for abundance (MetaPhyler [2] by default)

annotate.krona.html
-------------------

Krona [1] plot of abundances using the tool selected for classification (Kraken [3] by default)


asm.scores
----------

    Validation scores for each assembly/kmer combination run. Header contains information on scores generated

best.asm
--------

    The name of the assembly/kmer combination that was selected as the best

<taxonomy>.classified
---------------------

Subdirectory containing each level of the selected taxonomy (class by default) and the contigs/reads/orfs belonging to each

<taxonomy>.original.annots
--------------------------

    Tab-delimited taxonomic level assignments for each contig/unassembled read. Class IDs correspond to NCBI taxonomy IDs.

<taxonomy>.original.reads.annots
--------------------------------

Tab-delimited taxonomic level assignments as above, where contigs are replaced with their constituent sequences.

<taxonomy>.propagated.annots
----------------------------

Tab-delimited file as above after assembly graph-based propagation of assignments to contigs.

<taxonomy>.propagated.reads.annots
----------------------------------

Tab-delimited file as above after propagation and having contigs replaced with their constituent reads.

html (directory)
----------------

    HTML output from the pipeline. summary.html contains an interactive results view.

proba.bnk
---------

    AMOS bank format of the assembly that can be visualized using Hawkeye.

proba.classify.txt 	
------------------

    The raw output of the abundances using the tool selected for abundance estimations (MetaPhyler [2] by default)

proba.ctg.cnt	  
---------------    	  

    The number of sequences mapped to each assembly contig	

proba.ctg.cvg	  	    	   
-------------

    The coverage of each assembly contig

proba.ctg.fa	    	 
------------

    The assembled contigs

proba.hits			
----------

    The raw output of the contig/unassembled reads classifications using the selected tool (Kraken [3]) by default.

proba.lib1.contig.reads 
-----------------------

    The per-library assignment of sequences to contigs

proba.lib1.unaligned.fasta   
--------------------------

    The per-library unassembled sequences

proba.scf.fa				
------------

    The assembled scaffolds

proba.motifs.fa		
---------------

    The motifs within scaffolds identified by Bambus 2

proba.orf.faa
-------------

    The protein sequences of identified open reading frames (ORFs) in the assembly and unassembled reads

proba.orf.fna
-------------

    The fasta sequences of identified open reading frames (ORFs) in the assembly and unassembled reads

proba.scf.orf.faa
-----------------

    The protein sequences of identified open reading frames (ORFs) in the scaffolds

proba.scf.orf.fna
-----------------

    The protein sequences of identified open reading frames (ORFs) in the scaffolds
    
ref.fasta			
---------

    The recruited reference genome used for validation (iMetAMOS only)

ref.name	  	    	   
--------

    The name of the recruited reference genome (iMetAMOS only)

Additional details for each step are available under <STEP NAME>/out. This includes the raw
output (as well as any intermediate files) of any tools run during that step. For example, 
Annotate/out/proba.prokka includes the full Prokka annotation output. 
Assemble/out/abyss*/ contains the intermediate files output by ABySS. Additionally, 
since MetAMOS stores all of its results in an AMOS bank, the assemblies 
can be visualized with Hawkeye.

[1] Ondov BD, Bergman NH, Phillippy AM.. Interactive
metagenomic visualization in a Web browser. BMC Bioinformatics. 2011
Sep 30;12:385.  PMID: 21961884

[2] Liu B, Gibbons T, Ghodsi M, Treangen T, Pop M. Accurate and fast estimation of taxonomic profiles from metagenomic shotgun sequences. BMC Genomics. 2011;12 Suppl 2:S4. Epub 2011 Jul 27.

[3] Wood DE, Salzberg SL. Rapid phylogenetic sequence classification through repeated exact alignment. In preparation.
