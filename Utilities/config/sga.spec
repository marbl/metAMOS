[CONFIG]
maxlibs 1
input FASTQ
name SGA
output [PREFIX]-contigs.fa
scaffoldOutput [PREFIX]-scaffolds.fa
location cpp/[MACHINE]/sga/bin
threads -t
threadsSupport	LINUX
paired --pe-mode 1 [FIRST] [SECOND]
paired_interleaved --pe-mode 2 [FIRST]
unpaired --pe-mode 0 [FIRST]
commands sga preprocess [INPUT] -o [PREFIX].pp.fastq && \
         sga index -a ropebwt [THREADS] --no-reverse [PREFIX].pp.fastq && \
         sga correct -k 21 --learn [THREADS] -o [PREFIX].ec.fastq [PREFIX].pp.fastq && \
         sga index -a ropebwt [THREADS] [PREFIX].ec.fastq && \
         sga filter -x 2 [THREADS] [PREFIX].ec.fastq && \
         sga overlap -m [KMER] [THREADS] [PREFIX].ec.filter.pass.fa && \
         sga assemble -m [KMER] --min-branch-length 500 -o [PREFIX] [PREFIX].ec.filter.pass.asqg.gz
# 	 skip scaffolds in sgs, it requires parts of abyss to be in system path as well as pysam
#         bwa index [PREFIX]-contigs.fa && \
#         bwa aln [THREADS] [PREFIX]-contigs.fa [FIRST] > [PREFIX].1.sai && \
#         bwa aln [THREADS] [PREFIX]-contigs.fa [SECOND] > [PREFIX].2.sai && \
#         bwa sampe [PREFIX]-contigs.fa [PREFIX].1.sai [PREFIX].2.sai [FIRST] [SECOND] > [PREFIX].sam && \
#         samtools view -Sb [PREFIX].sam > [PREFIX].bam && \
#         sga-bam2de.pl -n 10 -m 200 --prefix [PREFIX] [PREFIX].sam && \
#         sga-astat.py -m 200 [PREFIX].bam > [PREFIX].astat && \
#         sga scaffold -m 200 -a [PREFIX].asta -o [PREFIX].scaf --pe [PREFIX].de [PREFIX]-contigs.fa && \
#         sga scaffold2fasta --use-overlap --write-unplaced -m 200 [PREFIX]-graph.asqg.gz -o [PREFIX]-scaffolds.fa [PREFIX].scaf
