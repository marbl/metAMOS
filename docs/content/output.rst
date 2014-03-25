############
Output
############

Directory layout/description
============================

All of the step mentioned below have the following directory structure:: 

    [STEP]/in  -> required input
    [STEP]/out -> generated output

We will now describe in detail the functionality of each step, along with the expected input & output.

Preprocess
-------------

Required step?
^^^^^^^

* Yes

Software currently supported
^^^^^^^

* ea-utils (code.google.com/p/ea-utils)
* FastQC (bioinformatics.babraham.ac.uk)
* KmerGenie (Chikhi et al 2014)

What it does
^^^^^^^

* Quality control
* Read filtering 
* Read trimming
* Sanity checks on fasta/q files
* Conversion to required formats

Expected input
^^^^^^
* Raw reads

Expected output
^^^^^^
* Cleaned reads
* Quality report
* Converted files

Assemble
-------------

Required step?
^^^^^^^

* No

Software currently supported
^^^^^^^

* ABySS (Simpson et al 2009)
* CABOG (Miller et al 2008)
* IDBA-UD (Peng et al 2012)
* MaSuRCA (Zimin et al 2013)
* MetaVelvet (Namiki et al 2011)
* Mira (Chevreux et al 1999)
* RayMeta (Boisvert et al 2012)
* SGA (Simpson et al 2012)
* SOAPdenovo2 (Luo et al 2012)
* SPAdes (Bankevich et al 2012)
* SparseAssembler (Ye et al 2012)
* Velvet (Zerbino et al 2008)
* Velvet-SC (Chitsaz et al 2011)

What it does
^^^^^^^

* Construct assembly (no scaffolds)

Expected input
^^^^^^

* Cleaned reads

Expected output
^^^^^^

* Unitigs
* Contigs
* Singletons
* Degenerates/Surrogates

FindORFs
-------------

Required step?
^^^^^^^

* No

Software currently supported
^^^^^^^

* FragGeneScan
* MetaGeneMark

What it does
^^^^^^^

* Finds/predicts ORFs in contigs

Expected input
^^^^^^

* Assembled contigs in fasta format (>300bp) 

Expected output
^^^^^^

* ORFs in multi-fasta format (FAA,FNA)


Validate
-------------

Required step?
^^^^^^^

* No

Software currently supported
^^^^^^^

* ALE (Clark et al 2013)
* CGAL (Rahman et al 2013)
* FRCbam (Vezzi et al 2013)
* FreeBayes (Garrison et al 2012)
* LAP (Ghodsi et al 2013)
* QUAST (Gurevich et al 2013)
* REAPR (Hunt et al 2013)


What it does
^^^^^^^

* Checks assembly correctness using intrinsic quality metrics

Expected input
^^^^^^

* Assembled contigs in fasta format

Expected output
^^^^^^

* List of errors
* Poorly assembled regions 
* Assembly quality metrics

FindRepeats (deprecated)
-------------

This step was initially added to help speed up Bambus 2 repeat identification step; optimizations to Bambus 2 have made this speed-up unnecessary. Step is turned off by default.

Required step?
^^^^^^^

* No


Software currently supported
^^^^^^^

* Repeatoire

What it does
^^^^^^^

* Find contigs (or parts of contigs) that appear to be repetitive and flag for further steps.

Expected input
^^^^^^

* Assembled contigs in fasta format

Expected output
^^^^^^

* List of contigs likely to be repeats


Abundance (deprecated)
-------------

This step was created to estimate taxonomic abundance of a give metagenomic sample 

Required step?
^^^^^^^

* No


Software currently supported
^^^^^^^

* MetaPhyler (Liu et al 2011)

What it does
^^^^^^^

* Find contigs (or parts of contigs) that appear to be repetitive and flag for further steps.

Expected input
^^^^^^

* Assembled contigs in fasta format

Expected output
^^^^^^

* List of contigs likely to be repeats

FunctionalAnnotation
-------------

Required step?
^^^^^^^

* No


Software currently supported
^^^^^^^

* Prokka (Seemann, 2013)
* BLAST

What it does
^^^^^^^

* Assigns functional annotation to ORFs

Expected input
^^^^^^

* ORFs in multi-fasta format (FAA,FNA)

Expected output
^^^^^^

* Text file containing functional labels for ORFs

Scaffold
-------------

Required step?
^^^^^^^

* Yes

Software currently supported
^^^^^^^

* Bambus2 (Koren, 2011)

What it does
^^^^^^^

* Link together contigs using mate-pairs.  Also identify variant patterns.

Expected input
^^^^^^

* Assembled contigs in fasta format

Expected output
^^^^^^

* scaffolds in agp format
* scaffolds in fasta format
* motifs/variants 
* longer contigs in fasta format

Postprocess
-------------

Required step?
^^^^^^^

* Yes


Software currently supported
^^^^^^^

* Krona (Ondov, 2010)

What it does
^^^^^^^

* Generates summary reports
* Collates output
* Generates combined HTML page

Expected input
^^^^^^

* Majority of the aforementioned outputs

Expected output
^^^^^^

* HTML summary file
* Output directory tree


Full listing of expected output files
===============================

MetAMOS generates an interactive web page once a run successfully completes::

     http://www.cbcb.umd.edu/~sergek/imetamos/gageb/Postprocess/out/html/summary.html

This includes summary statistics and taxonomic information based on Krona [1].
The easiest way to interact with the results is through the web interface.
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
