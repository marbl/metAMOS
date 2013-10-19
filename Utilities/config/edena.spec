# how to specify multiple commands
# lots of special keywords being used now, get away from this
[CONFIG]
input FASTQ
name Edena
output [PREFIX]_contigs.fasta
location cpp/[MACHINE]/edena/bin
threads -nThreads
paired -paired [FIRST] [SECOND]
mated -matePairs [FIRST]  [SECOND]
commands \
	 edena -p [PREFIX] [THREADS] [INPUT] -minOverlap [KMER] && \
         edena -p [PREFIX] -e [PREFIX].ovl
unpaired singleEnd [FIRST]
