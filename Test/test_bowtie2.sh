#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test_bowtie2  -i 500:3500
../runPipeline -a soap -m bowtie2 -c fcp -g fraggenescan -p 15 -d test_bowtie2 -k 25 -t -f Assemble,MapReads,FindORFS,FindScaffoldORFS,Abundance,Classify,Annotate
