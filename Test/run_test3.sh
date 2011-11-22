#/bin/sh
#../initPipeline -f -m MA_test2.filt2.fna -d test2 -i 1300:1700
#../runPipeline -c metaphyler -p 8 -d test2 -k 31 -f Scaffold,Propagate,Classify -v
#now test contigs
#mv ./test2/Assemble/out/proba.asm.contig test.asm
../initPipeline -f -m MA_test2.filt2.fna -d test3 -i 1300:1700 -c test.asm
../runPipeline -c metaphyler -p 8 -d test3 -k 31 -f Scaffold,Propagate,Classify -v
