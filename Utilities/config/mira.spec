[CONFIG]
input FASTQ
name MIRA
output mira_assembly/mira_d_results/mira_out.unpadded.fasta 
scaffoldOutput mira_assembly/mira_d_results/mira_out.padded.fasta
location cpp/[MACHINE]/mira/bin
config config/mira.template
paired readgroup = pe[LIB]\n data = [FIRST] [SECOND]\n technology = [TECHNOLOGY]\ntemplatesize = [MIN] [MAX]\nsegmentplacement = [ORIENTATION_FIGURE] \n segmentnaming = solexa
mated readgroup = mp[LIB]\n data = [FIRST] [SECOND]\n technology = [TECHNOLOGY]\ntemplatesize = [MIN] [MAX]\nsegmentplacement = [ORIENTATION_FIGURE] \n segmentnaming = solexa
unpaired readgroup = up[LIB]\n data = [FIRST]\n technology = [TECHNOLOGY]
solexa -AL:mo=[KMER] -CO:fnicpst=yes
sanger -AL:mo=[KMER] -CO:fnicpst=yes
454 -AL:mo=[KMER] -CO:fnicpst=yes
commands rm -rf mira_assembly && \
         mira [INPUT]
