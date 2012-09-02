#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test1  -i 500:3500
../runPipeline -a soap -c fcp -g fraggenescan -p 15 -d test1 -k 45 -t -f Assemble,MapReads,FindORFS,FindScaffoldORFS,Abundance,Classify,Annotate
