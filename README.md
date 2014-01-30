# MetAMOS v1.5rc1 "Praline Brownie" README
Last updated: January 24th 2014
***********************************************************************************
We are happy to announce version 1.5rc1, a.k.a. Praline Brownie 

## Table of Contents ##

- [NEWS](#news)
- [MetAMOS single file binary](#metamos-single-file-binary)
    - [Download](http://www.cbcb.umd.edu/confcour/temp/metAMOS_binary_noblastdbs.tar.gz)
- [SUMMARY](#summary)
    - [HARDWARE REQUIREMENTS](#a-hardware-requirements)
         - 32GB of RAM
    - [SOFTWARE REQUIREMENTS](#b-software-requirements)
        - git
        - gcc
        - automake
        - python-tools
        - python-devel
        - zlib-devel
        - numpy
        - freetype, freetype-devel
        - libpng-devel
        - matplotlib
        - curl
    - [INSTALLING METAMOS](#c-installing-metamos)
        - python INSTALL.py core
    - [QUICK START](#d-quick-start)
        - initPipeline -q -1 fastq.1 -2 fastq.2 -d projectDir -W core
        - runPipeline -d projectDir -p 16
    - [WORKFLOWS](#e-workflows)
        - core
        - imetamos
    - [GENERIC TOOL](#f-generic-tools)
    - [TEST SUITE](#g-test-suite)
    - [EXAMPLE OUTPUT](#h-example-output)
    - [CONTACT](#i-contact)
    - [CITE](#j-cite)
        - PMID: 23320958
    - [ISSUES](#k-issues)
- [FIN](#end)

## NEWS
       1. iMetAMOS now available.
       2. Updated 64-bit frozen Linux binary now available (fixed FastQC issue)
       3. Generic framework for adding assemblers/classifiers
       4. SRA run identifiers supported.
       5. Remote input files, as well as compressed input files are supported.
       6. Numerous bug fixes

*on deck: kmher, and Viritas

-----------------------------------------------------------------------------------
[TOC](#table-of-contents)
## MetAMOS single file binary
In attempt to further simplify the MetAMOS installation process, we are happy to announce the availability of a 'frozen' MetAMOS binary for Linux-x68_64 platforms. Along with this binary comes a significantly reduced list of prerequisites:

* Java 1.6 (or newer)
* Perl 5.8.8 (or newer)
* 64-bit *nix OS

You can download a frozen binary with all BLAST DBs [here](http://www.cbcb.umd.edu/confcour/temp/metAMOS_binary.tar.gz) 

*Note: this is a big file to download and will take awhile to finish (~13GB compressed, 27GB unzipped)

Alternatively, you can download a frozen binary with NO BLAST DBs [here](http://www.cbcb.umd.edu/confcour/temp/metAMOS_binary_noblastdbs.tar.gz) 

*Note: ~2GB compressed, 4GB unzipped

The revised installation procedure is simply:

1. Download tarball
2. Extract
3. Run a test script or two

Mac OSX equivalent will be soon appearing [here](http://www.cbcb.umd.edu/confcour/temp/metAMOS_binary_noblastdbs_osx.tar.gz) 

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
## SUMMARY 
        A) HARDWARE REQUIREMENTS
        B) SOFTWARE REQUIREMENTS
        C) INSTALLING MetAMOS
        D) QUICK START
        E) WORKFLOWS
        F) GENERIC TOOL
        G) TEST SUITE
        H) EXAMPLE OUTPUT
        I) CONTACT
        J) CITE
        K) ISSUES

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### A) HARDWARE REQUIREMENTS 

MetAMOS was designed to work on any standard 64bit Linux
environment. To use MetAMOS for tutorial/teaching purposes, a minimum
of 8 GB RAM is required. To get started on real data sets a minimum of
32 GB of RAM is recommended, and anywhere from 64-1000 GB may be
necessary for larger datasets. In our experience, for most 50-100
million read datasets, 64 GB is a good place to start (128 GB of memory now
available on High Memory Instance at Amazon Elastic Compute Cloud ). 

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### B) SOFTWARE REQUIREMENTS 

If you are using the frozen binary, you can skip this section. The MetAMOS frozen
binary includes dependencies and only requires perl 5.8+ and java 1.6+.

The main prerequisite software for installing/running MetAMOS is python 2.6+, perl 5.8+,
and java 1.6+. Depending on your platform/Linux distribution, you might also need to download and 
install the following BEFORE running INSTALL.py:

1. git
2. gcc
3. automake
4. python-tools
5. python-devel
6. zlib-devel
7. numpy
8. freetype, freetype-devel
9. libpng-devel
10. matplotlib
11. curl

Additional software will be downloaded by MetAMOS as needed.
Additionally, there is some software that MetAMOS can
incorporate into its pipeline that we are not allowed to distribute,
such as MetaGeneMark and Newbler. To get a license to use MetaGeneMark, please
visit: http://exon.gatech.edu/license_download.cgi. Once the tool is installed,
add it to your PATH variable and MetAMOS will then enable its use in the pipeline.

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### C) INSTALLING MetAMOS 

To download the software release package, go [here](https://github.com/treangen/metAMOS/archive/Release1.5.zip).
You can also browse the [repository](https://github.com/treangen/MetAMOS/tree/Release1.5)
and click on Downloads. Once downloaded, simply unpack the files and
open the MetAMOS directory. Once inside the MetAMOS directory, run:
```
python INSTALL.py
```
This will download and install the external dependencies which may 
take minutes or hours to download depending on your connection speed. 
metAMOS supports workflows to install subsets of tools for faster installation.
By default only the core dependencies are installed. To install iMetAMOS run
```
python INSTALL.py iMetAMOS
```
You can run:
```
python INSTALL.py -h
```
to get a listing of available workflows and programs. You can specify either
workflows or programs as arguments to INSTALL.py. For example, to install the 
core workflow plus PhyloSift, run
```
python INSTALL.py core phylosift
```
To install the programs which are part of the optional workflow run
```
python INSTALL.py optional
```
If all dependencies are downloaded (including optional/deprecated ones), this will take 
quite awhile to complete (plan on a few hours to 2 days).

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### D) QUICK START 

Before you get started using MetAMOS/iMetAMOS a brief review of its design will
help clarify its intended use. MetAMOS gas two main components:

1. initPipeline
2. runPipeline

Below is a simple example of running of iMetAMOS to assemble an SRA dataset:
```
initPipeline -q -1 SRR987657 -d projectDir -W iMetAMOS
runPipeline -d projectDir -p 16
```

The first component, initPipeline, is for creating new projects and
also initializing sequence libraries. Currently interleaved &
non-interleaved fasta, fastq, and SFF files are supported. Input
files can be compressed (bzip2, gzip) and can reside on remote
servers (in this case the full URL must be specified). SRA run identifiers
are also supported. 

The file-type flags (-f, -q, and -s) must be specified before the file.
Once specified, they remain in effect until a different file type is specified.

usage: initPipeline -f/-q -1 file.fastq.1 -2 file.fastq.2 -d projectDir -i 300:500 

options: -s -c -q, -f, -1, -2, -d, -m, -i

* -1: either non-paired file of reads or first file in pair, can be list of multiple separated by a comma
* -2: second paired read file, can be list of multiple separated by a comma
* -c:  fasta file containing contigs
* -d: output project directory (required)
* -f: boolean, reads are in fasta format (default is fastq)
* -h: display help message
* -i: insert size of library, can be list separated by commas for multiple libraries
* -l: SFF linker type
* -m: interleaved file of paired reads
* -o: reads are in outtie orientation (default innie)
* -q: boolean, reads are in fastq format (default is fastq)
* -s/--sff: boolean, reads are in SFF format (default is fastq)
* -W: string: workflow name (-W iMetAMOS will run iMetAMOS). A workflow can specify parameters as well as data. A workflow can be immutable in which case any command-line parameters will not be used. Otherwise, command-line parameters take priority over workflow defaults.

For example, to input a:

(non-interleaved fastq, single library)
```
initPipeline -q -1 file.fastq.1 -2 file.fastq.2 -d projectDir -i 300:500
```
(non-interleaved fasta, single library)
```
initPipeline -f -1 file.fastq.1 -2 file.fastq.2 -d projectDir -i 300:500
```
(interleaved fastq, single library)
```
initPipeline -q -m file.fastq.12  -d projectDir -i 300:500
```
(interleaved fastq, multiple libraries)
```
initPipeline -q -m file.fastq.12,file2.fastq.12  -d projectDir -i 300:500,1000:2000
```
(interleaved fastq, multiple libraries, existing assembly)
```
initPipeline -q -m file.fastq.12,file2.fastq.12 -c file.contig.fa -d projectDir -i 300:500,1000:2000
```
(non-interleaved remote fastq, single library)
```
initPipeline -q -1 ftp://ftp.cbcb.umd.edu/pub/data/metamos/gage-b-rb.miseq.1.fastq.gz -2 ftp://ftp.cbcb.umd.edu/pub/data/metamos/gage-b-rb.miseq.2.fastq.gz -d projectDir -i 300:500
```
(unpaired SRA run using iMetAMOS)
```
initPipeline 1 <SRA RUN ID> -d projectDir -W iMetAMOS
```
(paired-end SRA run using iMetAMOS)
```
initPipeline -m <SRA RUN ID> -d projectDir -i 300:500 -W imetAMOS
```
The second component, runPipeline, takes a project directory as
input and runs the following steps by default:

1. Preprocess
2. Assemble
3. FindORFs
4. Validate
5. FindRepeats
6. Abundance
7. Annotate
8. FunctionalAnnotation
9. Scaffold
10. Propagate 
11. FindScaffoldORFs
12. Classify 
13. Postprocess

usage info:

usage: runPipeline [options] -d projectdir

*   -h = <bool>:   print help [this message]
*   -j = <bool>:   just output all of the programs and citations then exit (default = NO)
*   -v = <bool>:   verbose output? (default = NO)
*   -d = <string>: directory created by initPipeline (REQUIRED)

[options]: [pipeline_opts] [misc_opts]

[pipeline_opts]: options that affect the pipeline execution

Pipeline consists of the following steps:

Preprocess, Assemble, FindORFS, MapReads, Validate, Abundance, Annotate,
  Scaffold, Propagate, Classify, Postprocess

Each of these steps can be referred to by the following options:

*   -f = <string>: force this step to be run (default = NONE)
*   -s = <string>: start at this step in the pipeline (default = Preprocess)
*   -e = <string>: end at this step in the pipeline (default = Postprocess)
*   -n = <string>: step to skip in pipeline (default=NONE)

For each step you can fine-tune the execution as follows

[Preprocess]

*   -t = <string>:   enable filter of input reads (default = metAMOS, options = metAMOS, EA-UTILS, PBcR for PacBio sequences)
*   -q = <bool>:   produce FastQC quality report for reads with quality information (fastq or sff)? (default = NO)

[Assemble]

*   -a = <string>: genome assembler to use (default = SOAPdenovo).
    	  This can also be a comma-separated list of assembler (for example: soap,velvet)
		  in this case, all selected assemblers will be run and the best selected for subsequent analysis
*   -k = <kmer size>: k-mer size to be used for assembly (default = auto-selected).
*         This can also be a comma-separated list of kmers to use
*   -o = <int>:    min overlap length

[MapReads]

*   -m = <string>: read mapper to use? (default = bowtie)
*   -i = <bool>:   save bowtie (i)ndex (default = NO)
*   -b = <bool>:   create library specific per bp coverage of assembled contigs (default = NO)

[FindORFS]

*   -g = <string>: gene caller to use (default=FragGeneScan)
*   -l = <int>:    min contig length to use for ORF call (default = 300)
*   -x = <int>:    min contig coverage to use for ORF call (default = 3X)

[Validate]

*   -X = <string>: comma-separated list of validators to run on the assembly. (default = lap, supported = reapr,orf,lap,ale,quast,frcbam,freebayes,cgal,n50)
*   -S = <string>: comma-separated list of scores to use to select the winning assembly. By default, all validation tools specified by -X will be run. For each score, an optional weight can be specified as SCORE:WEIGHT. For example, LAP:1,CGAL:2 (supported = all,lap,ale,cgal,snp,frcbam,orf,reapr,n50)

[Annotate]

*   -c = <string>: classifier to use for annotation (default = FCP)
*   -u = <bool>:   annotate unassembled reads (default = NO)

[Classify]

*   -z = <string>: taxonomic level to categorize at (default = class)

[misc_opts]: Miscellaneous options

*   -r = <bool>:   retain the AMOS bank  (default = NO)
*   -p = <int>:    number of threads to use (be greedy!) (default=1)
*   -4 = <bool>:   454 data (default = NO)

For example, to enable read filtering:
```
-t
```
and to enable meta-IDBA as the assembler:
```
-a metaidba
```
And to use PhyloSift to annotate:
```
-c phylosift
```
Any single step in the pipeline can be skipped by passing the
following parameter to runPipeline:
```
-n,--skipsteps=Step1,..
```
MetAMOS reruns steps based on timestamp information, so if the input
files for a step in the pipeline hasn't changed since the last run, it
will be skipped automatically. However, you can forcefully run any step
in the pipeline by passing the following parameter to runPipeline:
```
-f,--force=Step1,..
```
MetAMOS stores a summary of the input libraries in pipeline.ini 
in the working directory. The pipeline.conf file stores the list 
of programs available to MetAMOS. Finally, pipeline.run stores the 
selected parameters and programs for the current run. MetAMOS also stores 
detailed logs of all commands executed by the pipeline in Logs/COMMANDS.log 
and a log for each step of the pipeline in Logs/<STEP NAME>.log

Upon completion, all of the final results will be stored in the
Postprocess/out directory. A component, create_summary.py, takes
this directory as input and as output, generates an HTML page with
with summary statistics and a few plots. An optional component, create_plots.py,
takes one or multiple Postprocess/out directories as input and generates
comparative plots.


----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### E) WORKFLOWS

A workflow is a text-file that specified command-line options and input sequences
required to run metAMOS. A workflow may optionally inherit options/data from other
workflows. A workflow may also be immutable if the parameters should not be modifiable
by a user. An example workflow is below:
```
inherit:isolate
modify:True
command: -q -u -r -v -I -c kraken -p 16 -a spades,velvet-sc,abyss,ray,edena,sga,masurca,soap,soap2,velvet -t metamos -n FunctionalAnnotation -f Postprocess -z phylum
asmcontigs:	/Users/skoren/Personal/Research/metAMOS/Test/test.asm,ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Candidatus_Carsonella_ruddii_uid58773/NC_008512.fna
lib1format:	fasta
lib1mated:	True
lib1innie:	True
lib1interleaved:	True
lib1f1:	/Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe_filt.fna.gz,2000,5000,3500,500
```
The available options are:

*   inherit	- any other workflows to inherit from. In this case, the workflow inherits options from the isolate workflow
*	modify	- whether users are allowed to specify command-line parameters at runtime. If false, command-line options are ignored
*	command	- command-line options to specify for runPipeline
*	asmcontigs - optional, pre-assembled contigs to include in analysis. Can be remote file. Multiple files can be separated using commas.
*	lib#format - input type for lib #. Can be fasta/fastq/sff
*	lib#mated - whether the library is mated or not
*	lib#innie - whether the mates are in the innie (Illumina paired-end) format or not (Illumina mate-pair)
*	lib#interleaved - whether the input sequences are in a single file or in two separate files
*	lib#f1 	- the name of the input file, along with library min, max, mean, stdev

An arbitrary number of libraries may be specified in the above format. The below example shows an unmated library:
```
lib1format:	fasta
lib1mated:	False
lib1innie:	False
lib1interleaved:	False
lib1frg:	/Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe_filt.fna.gz
```
as well as a non-interleaved library:
```
lib1format:     fasta
lib1mated:      True
lib1innie:      True
lib1interleaved:        False
lib1f1: /Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe_1.fna.gz,2000,5000,3500,500
lib1f2: /Users/skoren/Personal/Research/metAMOS/Test/carsonella_pe.2.fna.gz,2000,5000,3500,500
```

Workflows may be shared between users, as long as the input files are accessible (i.e. they are on a remote server or the systems share a file system). Workflow files should be placed in the metAMOS/workflows directory or the working directory where MetAMOS is launched.

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### F) GENERIC TOOLS

MetAMOS allows new tools to be added to the ASSEMBLE and ANNOTATE steps without requiring code changes. The addition of a tool is a three-step process. 

###### 1) Add the tool name under metAMOS/Utilities/\<STEPNAME\>.generic. For example. if you want to add a new assembler, you would modify ASSEMBLE.generic. This file contains one tool name per line. The tool name is arbitrary text and will be used by MetAMOS to look up detailed configuration. The current ASSEMBLE.generic looks like:
```
% cat Utilities/config/ASSEMBLE.generic
abyss
sga
spades
ray
masurca
mira
edena
idba-ud
```
You can add multiple versions of an assembler. In this documentation, we will add SOAPdenovo v1.05 in addition to the above tools. First, we will add soap_v105 to the end of ASSEMBLE.generic:
```
% cat Utilities/config/ASSEMBLE.generic
abyss
sga
spades
ray
masurca
mira
edena
idba-ud
soap_v105
```
###### 2) Write a configuration file for the tool. The configuration file specifies input requirements for the program as well as a name, output, and executable location. Within configuration files, several keywords may be specified that are updated at runtime. The list of currently supported keywords can be found at the end of this section. In the above example, MetAMOS would expect a file named soap_v105.spec.

Below is an example configuration file used for Ray:
```
% cat Utilities/config/ray.spec 
[CONFIG]
maxlibs 1
input FASTQ
name Ray
output [PREFIX]_ray/Contigs.fasta
scaffoldOutput [PREFIX]_ray/Scaffolds.fasta
location cpp/[MACHINE]/Ray/bin
threads -n
paired_interleaved -i [FIRST]
paired -p [FIRST] [SECOND]
commands rm -rf [RUNDIR]/ray && \
		[MPI] [THREADS] Ray \
		-o [RUNDIR]/[PREFIX]_ray [INPUT]
unpaired -s [FIRST]
[Ray]
k	[KMER]
```
The [CONFIG] section is the generic configuration section, you can specify step-specific configuration later on. Here, most properties of where the tool is located, what its output is, and what input it requires is specified:

 * input - the type of input (FASTQ in this case)
 * name - the full name of the tool you want to report later on. This can be arbitrary text.
 * output - where the output contigs from the tool are. For assemblers, this is contigs. [PREFIX] is a keyword for the MetAMOS prefix for the assembly when it is run. This is assumed to be relative to the MetAMOS run directory.
 * scaffoldOutput - where the output scaffolds from the tool are, if available.
 * backupOutput - some assemblers fail to generate their final output on some datasets. In this case, this can specify preliminary contig output which will only be used if the main output is not available.
 * location - path to the executable. This is relative to metAMOS/Utilities. You can specify [MACHINE] to substitute your machine type into the executable path (i.e. Linux-x86_64). The user path will be searched if the tool is not found in the specified location
 * threads - the parameter to pass number of threads to use for the program, if available
 * paired - how to pass paired-end (assumed innie) interleaved data (FIRST refers to left mates, SECOND to right)
 * paired_interleaved - how to pass paired-end (assumed innie) non-interleaved. FIRST refers to the interleaved file.
 * mated - how to pass mate-pair data (assumed outtie) non-interleaved data (FIRST refers to left mates, SECOND to right)
 * mated_interleaved - how to pass mate-pair data (assumed outtie) interleaved mates
 * unpaired - how to pass fragment data to the program. FIRST refers to the unmated file.
 * commands - an arbitrary list of commands to run to execute the tool. Multiple lines are supported with the \ character. Multiple commands can be specified using &&. In the above example, rm -rf will run first followed by Ray. Common useful keywords are:
    * [PREFIX] - the prefix to use for output 
    * [RUNDIR] where the program is running
    * [KMER] - the selected k-mer to use for assembly
    * [MEM] - available memory
    * [THREADS] - the threads parameter and number of threads requested by the user
    * [INPUT] - the formatted input based on the libraries provided to metAMOS

The [Ray] section is a step-specific configuration. This is based on the executable names used in commands above. By default the parameters will be passed with prefixed - so here Ray will be run with -k [KMER]

Some assemblers (SOAPdenovo, MaSuRCA, etc) require an input configuration file rather than taking parameters on the command line. In this case, we need both a spec and template file (soap_v105.spec and soap_v105.template) which will get updated at runtime and passed to the assembler. The [CONFIG] section then includes a config option which specifies the template and the keyword [INPUT] will pass the configuration file rather than library information. 

Below is an example spec file for SOAPdenovo that requires a template and spec file:
```
% cat Utilities/config/soap_v105.spec
[CONFIG]
input FASTQ
name soap_v105
threads -p
output [PREFIX]/[PREFIX].asm.contig
location cpp/[MACHINE]/SOAPdenovo_1.05/
scaffoldOutput [PREFIX]/[PREFIX].asm.scafSeq
config config/soap_v105.template
mated rank=[LIB]\navg_ins=[MEAN]\nreverse_seq=1\nasm_flags=2\nq1=[FIRST]\nq2=[SECOND]
paired rank=[LIB]\navg_ins=[MEAN]\nreverse_seq=0\nasm_flags=3\nq1=[FIRST]\nq2=[SECOND]
unpaired rank=[LIB]\navg_ins=0\nq=[FIRST]
commands rm -rf [PREFIX] && \
         mkdir [PREFIX] && \
	 SOAPdenovo all -s [INPUT] -o [PREFIX]/[PREFIX].asm -K [KMER] [THREADS]
```
```
% cat Utilities/config/soap_v105.template
#maximal read length
max_rd_len=150
[LIB]
[INPUT]
```
Here, the config template is specified (again relative to metAMOS/Utilities) and the [INPUT] keyword will be replaced by the library information at run time.

###### 3) Add a citation to the tool under metAMOS/doc/citations.rst
Citations are tab-delimited and specify the lower-case tool alias, full tool-name, and citation information. For example:
```
soap_v105	SOAPdenovo v1.05	Li Y, Hu Y, Bolund L, Wang J: State of the art de novo assembly of human genomes from massively parallel sequencing data.Human genomics 2010, 4:271-277.
```
The citation will be automatically printed by MetAMOS whenever a run uses the specified tool. 

###### 4) For ANNOTATE tools, we also need a way to convert the output to Krona. By default, MetAMOS will look for an Import\<toolName\>.pl script. If one is not found, it will rely on a generic import which will assumed a tab-delimited format:
```
contig/readID	NCBI Taxonomy ID
```
The currently supported list of keywords:

*   MEM - max memory limit
*   LIB - library identifier (i.e. 1, 2, 3, etc)
*   INPUT - replace with input to the program (a collection of input files or libraries depending on the step or a configuration file)
*   MACHINE - replaced with Linux-x86_64, Darwin-x86_64, etc
*   FIRST - replaced with left mates in mated read or interleaved or unpaired reads otherwise
*   SECOND - replaced with right mates, in paired non-interleaved libs
*   ORIENTATION - replaced with the word innie or outtie
*   ORIENTATION_FIGURE - replaced with ---> <--- or <--- ---> for pe and mp, respectively
*   MEAN - replaced with library mean
*   SD - replaced with library standard dev
*   MIN - replaced with library min
*   MAX - replaced with library max
*   THREADS - replaced with thread parameter specified and requested number of threads
*   KMER - the kmer requested
*   OFFSET - the phred offset (33/64) of the input files
*   PREFIX - the desired prefix for the program output
*   DB - the location of the MetAMOS DBs (i.e. Utilities/DB)
*   RUNDIR - the location where the program is running (i.e. MetAMOS run directory)
*   LOCATION - the location where the program executable lives
*   TECHNOLOGY - the type of sequencing data (454, Illumina, etc) 

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### G) Test suite

We have developed a set of scripts for testing the various features of
MetAMOS. All of these regression test scripts are available inside the
/Test directory and include all necessary datasets to run them. Here
is a brief listing of the test scripts we currently include:

*Test initPipeline
./Test/test_create.sh

*Vanilla test
./Test/run_test.sh

*Test PhlyoSift
./Test/test_amphora.sh

*Test Minimus
./Test/test_minimus.sh

*Test Preprocess filtration of non-interleaved fastq files
./Test/test_filter_noninterleaved_fastq.sh

*Test iMetAMOS
./Test/test_ima.sh

*Test Newbler (if available)
./Test/test_newbler.sh

*Test CA (fasta)
./Test/test_ca_fasta.sh

*Test CA (fastq)
./Test/test_ca.sh

*Test SOAPdenovo
./Test/test_soap.sh

*Test MetaVelvet
./Test/test_metavelvet.sh

*Test SparseAssembler
./Test/test_sparse.sh

*Test Velvet
./Test/test_velvet.sh

*Test FCP
./Test/test_fcp.sh

*Test Spades
./Test/test_spades.sh

*Test BLAST
./Test/test_blast.sh

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### H) Example output

MetAMOS generates an interactive web page once a run successfully completes:
http://treangen.github.io/metAMOS/example/html/summary.html

This includes summary statistics and taxonomic information based on Krona [1].
The easiest way to interact with the results is through the web interface.
The Postprocess/out directory contains the results of the analysis. By default, 
metAMOS uses the prefix "proba" (Galician for test). Thus, files will have the name "proba".*.

* abundance.krona.html
    
    Krona [1] plot of abundances using the tool selected for abundance (MetaPhyler [2] by default)

* annotate.krona.html	

    Krona [1] plot of abundances using the tool selected for classification (Kraken [3] by default)

* asm.scores		

    Validation scores for each assembly/kmer combination run. Header contains information on scores generated

* best.asm		

    The name of the assembly/kmer combination that was selected as the best

* <taxonomy>.classified	

    Subdirectory containing each level of the selected taxonomy (class by default) and the contigs/reads/orfs belonging to each

* <taxonomy>.original.annots	

    Tab-delimited taxonomic level assignments for each contig/unassembled read. Class IDs correspond to NCBI taxonomy IDs.

* <taxonomy>.original.reads.annots	

    Tab-delimited taxonomic level assignments as above, where contigs are replaced with their constituent sequences.

* <taxonomy>.propagated.annots		

    Tab-delimited file as above after assembly graph-based propagation of assignments to contigs.

* <taxonomy>.propagated.reads.annots	

    Tab-delimited file as above after propagation and having contigs replaced with their constituent reads.

* html				

    HTML output from the pipeline. summary.html contains an interactive results view.

* proba.bnk			

    AMOS bank format of the assembly that can be visualized using Hawkeye.

* proba.classify.txt		

    The raw output of the abundances using the tool selected for abundance estimations (MetaPhyler [2] by default)

* proba.ctg.cnt			

    The number of sequences mapped to each assembly contig	

* proba.ctg.cvg			

    The coverage of each assembly contig

* proba.ctg.fa			

    The assembled contigs

* proba.hits			

    The raw output of the contig/unassembled reads classifications using the selected tool (Kraken [3]) by default.

* proba.lib1.contig.reads	

    The per-library assignment of sequences to contigs

* proba.lib1.unaligned.fasta	

    The per-library unassembled sequences

* proba.scf.fa			

    The assembled scaffolds

* proba.motifs.fa    		

    The motifs within scaffolds identified by Bambus 2

* proba.orf.faa

    The protein sequences of identified open reading frames (ORFs) in the assembly and unassembled reads

* proba.orf.fna

    The fasta sequences of identified open reading frames (ORFs) in the assembly and unassembled reads

* proba.scf.orf.faa

    The protein sequences of identified open reading frames (ORFs) in the scaffolds

* proba.scf.orf.fna

    The protein sequences of identified open reading frames (ORFs) in the scaffolds
    
* ref.fasta			

    The recruited reference genome used for validation (iMetAMOS only)

* ref.name			

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

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### I) CONTACT

If you encounter any problems/bugs, please check the known issues pages:
https://github.com/treangen/MetAMOS/issues?direction=desc&sort=created&state=open
to see if it has already been documented.

If not, please report the issue either using the contact information below or 
by submitting a new issue online. Please include information on your run,
any output produced by runPipeline, as well as the pipeline.* files and the 
Log/<LAST_STEP> file (if not too large).

Who to contact to report bugs, forward complaints, feature requests:

Todd Treangen: treangen@gmail.com
Sergey Koren: sergek@umd.edu

----------------------------------------------------------------------------------
[TOC](#table-of-contents)
### J) CITE

Treangen TJ\*, Koren S\*, Sommer DD, Liu B, Astrovskaya I, Ondov B,
Darling AE, Phillippy AM, Pop M.  MetAMOS: a modular and open source
metagenomic assembly and analysis pipeline. Genome Biol. 2013 Jan
15;14(1):R2. PMID: 23320958.

url: http://genomebiology.com/content/pdf/gb-2013-14-1-r2.pdf

*Indicates both authors contributed equally to this work

[TOC](#table-of-contents)
### K) ISSUES

Here is a [link](https://github.com/treangen/metAMOS/issues?state=open) to known, open issues.

## FIN
