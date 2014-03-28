
Generic tools (or plug-in framework)
===================================

Description
------------
MetAMOS allows new tools to be added to the ASSEMBLE and ANNOTATE steps without requiring code changes. 

Contributing to metAMOS
~~~~~~~~~~~~~~~~~~~~~~
If you add an assembler or classifier that you believe will benefit the community, please post the required spec file and citation either as a new `issue <https://github.com/marbl/metAMOS/issues?state=open>`_ through a `pull <https://github.com/marbl/metAMOS/pulls>`_ request. 

How-to-use
----------

The addition of a tool is a three (or four) step process; we will now review the required four steps.

Step 1: Add the tool name to spec file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add the tool name under metAMOS/Utilities/\<STEPNAME\>.generic. For example. if you want to add a new assembler, you would modify ASSEMBLE.generic. This file contains one tool name per line. The tool name is arbitrary text and will be used by MetAMOS to look up detailed configuration. The current ASSEMBLE.generic looks like::

    >cat Utilities/config/ASSEMBLE.generic

    abyss
    sga
    spades
    ray
    masurca
    mira
    edena
    idba-ud


Note: You can add multiple versions of an assembler. In this documentation, we will add SOAPdenovo v1.05 in addition to the above tools. First, we will add soap_v105 to the end of ASSEMBLE.generic::


    > cat Utilities/config/ASSEMBLE.generic
    abyss
    sga
    spades
    ray
    masurca
    mira
    edena
    idba-ud
    soap_v105

Step 2: Write a configuration file for the tool. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The configuration file specifies input requirements for the program as well as a name, output, and executable location. Within configuration files, several keywords may be specified that are updated at runtime. The list of currently supported keywords can be found at the end of this section. In the above example, MetAMOS would expect a file named soap_v105.spec.

Below is an example configuration file used for Ray::

    > cat Utilities/config/ray.spec 

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
    commands rm -rf [RUNDIR]/ray && [MPI] [THREADS] Ray -o [RUNDIR]/[PREFIX]_ray [INPUT]
    unpaired -s [FIRST]
   [Ray]
   k	
   [KMER]


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
 * commands - a list of commands to run to execute the tool. Multiple lines are supported with the \ character. Multiple commands can be specified using &&. A limited set of system commands is supported (mkdir, mv, bash, ln, rm, cp, ls, echo). Other system commands are currently not supported. In the above example, rm -rf will run first followed by Ray. Common useful keywords are:
    * [PREFIX] - the prefix to use for output 
    * [RUNDIR] where the program is running
    * [KMER] - the selected k-mer to use for assembly
    * [MEM] - available memory
    * [THREADS] - the threads parameter and number of threads requested by the user
    * [INPUT] - the formatted input based on the libraries provided to metAMOS

The [Ray] section is a step-specific configuration. This is based on the executable names used in commands above. By default the parameters will be passed with prefixed - so here Ray will be run with -k [KMER]

Some assemblers (SOAPdenovo, MaSuRCA, etc) require an input configuration file rather than taking parameters on the command line. In this case, we need both a spec and template file (soap_v105.spec and soap_v105.template) which will get updated at runtime and passed to the assembler. The [CONFIG] section then includes a config option which specifies the template and the keyword [INPUT] will pass the configuration file rather than library information. 

Below is an example spec file for SOAPdenovo that requires a template and spec file::

    >cat Utilities/config/soap_v105.spec
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
    commands rm -rf [PREFIX] && mkdir [PREFIX] && SOAPdenovo all -s [INPUT] -o [PREFIX]/[PREFIX].asm -K [KMER] [THREADS]

    >cat Utilities/config/soap_v105.template
    #maximal read length
    max_rd_len=150
    [LIB]
    [INPUT]

Here, the config template is specified (again relative to metAMOS/Utilities) and the [INPUT] keyword will be replaced by the library information at run time.

Step 3: Add a citation to the tool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Citations are tab-delimited and specify the lower-case tool alias, full tool-name, and citation information. For example::

    soap_v105	SOAPdenovo v1.05	Li Y, Hu Y, Bolund L, Wang J: State of the art de novo assembly of human genomes from massively parallel sequencing data.Human genomics 2010, 4:271-277.

The citation will be automatically printed by MetAMOS whenever a run uses the specified tool. 

Step 4: (ANNOTATE) convert to Krona input format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For ANNOTATE tools, we also need a way to convert the output to Krona. By default, MetAMOS will look for an Import\<toolName\>.pl script. If one is not found, it will rely on a generic import which will assumed a tab-delimited format::

    contig/readID	NCBI Taxonomy ID

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
