[CONFIG]
input FASTQ
name MaSuRCA
output CA/10-gapclose/genome.ctg.fasta
scaffoldOutput cA/10-gapclose/genome.scf.fasta
location cpp/[MACHINE]/MaSuRCA/bin
config config/masurca.template
paired PE=p[LIB] [MEAN] [SD] [FIRST] [SECOND]
mated JUMP=s[LIB] [MEAN] [SD] [FIRST] [SECOND]
required PAIRED
commands rm -rf CA && \
 	 rm -rf work1 && \
 	 rm -rf super* && \
 	 rm -rf guillaumeKUnitigsAtLeast32bases_all.* && \
	 rm -rf k_u_* && \
	 runSRCA.pl [INPUT] && \
	 bash assemble.sh
