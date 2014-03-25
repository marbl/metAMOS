############
Quick Start
############

Getting started
===============

Before you get started using MetAMOS/iMetAMOS a brief review of its design will
help clarify its intended use. MetAMOS gas two main components:

1. initPipeline
2. runPipeline

Below is a simple example of running of iMetAMOS to assemble an SRA dataset:
```
initPipeline -q -1 SRR987657 -d projectDir -W iMetAMOS
runPipeline -d projectDir -p 16
```

initPipeline
===============

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

runPipeline
===============

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



Preprocess
----------

*   -t = <string>:   enable filter of input reads (default = metAMOS, options = metAMOS, EA-UTILS, PBcR for PacBio sequences)
*   -q = <bool>:   produce FastQC quality report for reads with quality information (fastq or sff)? (default = NO)


Assemble
--------
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
