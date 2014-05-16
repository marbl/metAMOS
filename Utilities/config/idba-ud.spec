[CONFIG]
input FASTA
required PAIRED
name IDBA-UD
backupOutput [PREFIX].asm/contig.fa
scaffoldOutput [PREFIX].asm/scaffold.fa
output [PREFIX].asm/scaffold.fa
location cpp/[MACHINE]/idba/bin
threads --num_threads
paired_interleaved --read [FIRST]
commands idba_ud \
		--out [PREFIX].asm [THREADS] [INPUT]
[idba_ud]
-mink 21
-maxk [READLEN]
-pre_correction
