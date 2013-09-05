# how to specify multiple commands
# lots of special keywords being used now, get away from this
[CONFIG]
input FASTQ
name SPAdes
output contigs.fasta
location cpp/[MACHINE]/spades/bin
threads -t
paired_interleaved --pe[LIB]-12 [FIRST]
paired --pe[LIB]-1 [FIRST] --pe[LIB]-2 [SECOND]
commands spades.py \
		-o ./ -m [MEM] [THREADS] [INPUT]
unpaired --pe[LIB]-s [FIRST]
[spades.py]
k	21,33,[KMER]
-careful
-phred-offset [OFFSET]
