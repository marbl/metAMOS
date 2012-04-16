#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test_fcp  -i 500:3500
../runPipeline -a soap -c fcp -g fraggenescan -p 15 -d test_fcp -k 45 -f Assemble,MapReads,FindORFS,FindScaffoldORFS,Abundance,Classify,Annotate
