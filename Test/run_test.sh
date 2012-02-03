#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test1  -i 150:450
../runPipeline -a metaidba -c amphora2 -g fraggenescan -p 15 -d test1 -k 45 -f FindORFS,Abundance,Classify,Annotate
