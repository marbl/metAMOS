############
Experimental: TweetAssembler v0.1b
############

Introduction
===============

TweetAssembler is a twitter-based interface to an isolate genome assembly server powered by iMetAMOS::

	Automated ensemble assembly and validation of microbial genomes.
	Sergey Koren, Todd J Treangen, Christopher M Hill, Mihai Pop, Adam M Phillippy
	BMC Bioinformatics 15:126, 2014.
	http://www.biomedcentral.com/1471-2105/15/126/abstract


Why Twitter?
==============

- Good question! The main Raison d'être of TweetAssembler is to highlight the utility of iMetAMOS; just point it to your reads (no other params required!) and it will preprocess, tune, assemble, validate and create an HTML report of the results. This enables the submit command to be readily constructed in fewer than 140 chars. 

Limitations
===============

Before proceeding, its important to higlight a few important points:

- The server behind TweetAssembler is only able to assemble a couple of requests (at best) per day. Specs are: 32GB RAM & 32GhZ of compute... be gentle!
- There exists limitations on the size of the input data. i.e. MiSeq ok, HiSeq not ok. 
- TweetAssembler is nothing more than a tweet-based interface to an iMetAMOS webserver.
- Given the limited resources, job queue management is disabled. You will only be able to run a job if no other jobs are active; your only indication that your job was accepted is the confirmation tweet (see below). 
- Twitter has a maximum # of tweets allowed per day (1000), as well hourly limits. If TweetAssembler goes over any of these limits it will be deactivated for approx. 1 hour, potentially longer.
- No guarantees on preservation of output! Assemblies & associated output can & will be deleted regularly.

Quick Start
===============

1) First, issue a request to follow @imetamos:

.. image:: f0.png

2) Next, contact the developers to get your twitter account added to the `allowed accounts` list:

- Todd J Treangen (treangen@gmail.com)
- Sergey Koren (sergekoren@gmail.com)

3) Once approved, compose a tweet to @imetamos using the following syntax:

.. code-block:: bash

    @imetamos [fastq_pair_1] [fastq_pair_2] [#ASSEMBLE] [id]

- Think of an automated #icanhazpdf but for genome assemblies (#icanhazasm). 
- Currently reads need to be in non-interleaved fastq format.
- id simply needs to be a job-unique integer to avoid duplicate tweets in case you have to submit your job multiple times before it runs. 
- You should notice that no parameters are required (except for the input data). In practice this works thanks to several software packages, e.g. kmergenie (http://kmergenie.bx.psu.edu/), and an ensemble assembly approach (powered by several assembly and assembly validation tools). 

.. image:: f1.png

4) You should immediately receive a response tweet similar to:

.. image:: f2.png

5) Then simply wait for the confirmation tweet that the job was successful. 

.. image:: f3.png

6) Upon completion, you will be able to view the HTML report :

.. image:: f4.png

7) and download your assembly:

.. image:: f5.png

8) Suggestions & comments welcome! 


Viewing Output
===================

Your output will be located at http://www.traingene.com/tweetasm/P_[TIMESTAMP]/out/html/summary.html.

- Example output: http://www.traingene.com/tweetasm/P_2014_02_11_142926937305/out/html/summary.html

To save assembly:

.. code-block:: bash

     wget http://www.traingene.com/tweetasm/P_[TIMESTAMP]/out/proba.ctg.fa 

Supported Software
====================

Last but not least, we would like to acknowledge all of the wonderful software that provides the firepower behind TweetAssembler and iMetAMOS:

[Preprocess]

- ea-utils (code.google.com/p/ea-utils)
- FastQC (bioinformatics.babraham.ac.uk)
- KmerGenie (Chikhi et al 2014)

[Assemble]

- ABySS (Simpson et al 2009)
- CABOG (Miller et al 2008)
- IDBA-UD (Peng et al 2012)
- MaSuRCA (Zimin et al 2013) 
- MetaVelvet (Namiki et al 2011)
- Mira (Chevreux et al 1999)
- RayMeta (Boisvert et al 2012) 
- SGA (Simpson et al 2012)
- SOAPdenovo2 (Luo et al 2012)
- SPAdes (Bankevich et al 2012)
- SparseAssembler (Ye et al 2012)
- Velvet (Zerbino et al 2008)
- Velvet-SC (Chitsaz et al 2011)

[MapReads]

- Bowtie (Langmead  et al 2009) 
- Bowtie2 (Langmead  et al 2012) 

[Validate]

- ALE (Clark et al 2013)
- CGAL (Rahman et al 2013)
- FRCbam (Vezzi et al 2013)
- FreeBayes (Garrison et al 2012)
- LAP (Ghodsi et al 2013)
- QUAST (Gurevich et al 2013)
- REAPR (Hunt et al 2013)

[FindORFS/Annotate]

- Prokka (Seemann, 2013)

thanks!
