# how to specify multiple commands
# lots of special keywords being used now, get away from this
[CONFIG]
input FASTQ
name SPAdes
output contigs.fasta
location cpp/MACHINE/spades/bin
threads -t
paired_interleaved --peLIB-12 FIRST
paired --peLIB-1 FIRST --peLIB-2 SECOND
commands spades.py \
		-o ./ -m MEM THREADS INPUT
unpaired --peLIB-s FIRST
[spades.py]
k	21,33,KMER
-careful
-phred-offset OFFSET
