[CONFIG]
input FASTQ
name ABySS
output [PREFIX]-contigs.fa
scaffoldOutput [PREFIX]-scaffolds.fa
location cpp/[MACHINE]/abyss/bin
threads np=
paired lib[LIB]='[FIRST] [SECOND]'
mated lib[LIB]='[FIRST] [SECOND]'

commands abyss-pe k=[KMER] name=[PREFIX] lib='[PELIST]' mp='[MPLIST]' \
		[THREADS] [INPUT]
required PAIRED
