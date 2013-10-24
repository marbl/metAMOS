[CONFIG]
input FASTA
required PAIRED
name IDBA-UD
output [PREFIX].asm/contig.fa
location cpp/[MACHINE]/idba/bin
threads --num_threads
paired_interleaved --read [FIRST]
commands idba_ud \
		--out [PREFIX].asm [THREADS] [INPUT]
[idba_ud]
-mink 21
-maxk [KMER]
-pre_correction
