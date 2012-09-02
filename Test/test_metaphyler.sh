#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test_metaphyler -i 500:3500
../runPipeline -a soap -m bowtie -c metaphyler -g fraggenescan -p 15 -d test_metaphyler -k 25 -t -f Assemble,MapReads,FindORFS,FindScaffoldORFS,Abundance,Classify,Annotate -n FindRepeats
