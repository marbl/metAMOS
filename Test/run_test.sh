#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test1  -i 500:3500
../runPipeline -a metaidba -c phylosift -g fraggenescan -p 15 -d test1 -k 45 -f Assemble,MapReads,FindORFS,FindScaffoldORFS,Abundance,Classify,Annotate
