###########################
MetAMOS directory structure
############################

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

Annotate (a.k.a Taxonomic classifcation)
-------------

Required step?
^^^^^^^

* Yes

Software currently supported
^^^^^^^

* FCP
* Kraken
* Phylosift

What it does
^^^^^^^

* Labels contigs with taxonomic id

Expected input
^^^^^^

* Multi-fasta file of contigs

Expected output
^^^^^^

* Text file containing contig id to taxonomic id 1-to-1 mapping

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

Propagate
-------------

Required step?
^^^^^^^

* Yes

Software currently supported
^^^^^^^

* NA

What it does
^^^^^^^

* Propagate taxonomic labels along scaffolds

Expected input
^^^^^^

* Scaffolds in agp format
* Contig taxonomic labels

Expected output
^^^^^^

* contig taxonomic labels

FindScaffoldORFs
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

* Find ORFs in scaffolds, mainly serves as an extra validation step after Scaffold.

Expected input
^^^^^^

* Scaffolds in agp format

Expected output
^^^^^^

* Multi-fasta file of ORFs as fna,faa

Classify (a.k.a Bin)
-------------

Required step?
^^^^^^^

* Yes

Software currently supported
^^^^^^^

* NA


What it does
^^^^^^^

* Bins contigs/scaffold by taxonomic label

Expected input
^^^^^^

* Multi-fasta file of contigs
* Multi-fasta file of scaffolds

Expected output
^^^^^^

* Binned out contigs/scaffolds by directory


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
