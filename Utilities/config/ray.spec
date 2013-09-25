[CONFIG]
maxlibs 1
input FASTQ
name Ray
output ray/Contigs.fasta
scaffoldOutput ray/Scaffolds.fasta
location cpp/[MACHINE]/Ray/bin
threads -n
#disable threading, not yet tested
threadsSupport NONE
paired_interleaved -i [FIRST]
paired -p [FIRST] [SECOND]
commands rm -rf [RUNDIR]/ray && \
		[MPI] [THREADS] Ray \
		-o [RUNDIR]/ray [INPUT]
unpaired -s [FIRST]
[Ray]
k	[KMER]
