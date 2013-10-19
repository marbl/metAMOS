#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test_fast  -i 500:3500
../dist/pdfposter/runPipeline  -a soap -c fcp -g fraggenescan -y -p 15 -d test_fast -k 55 -t -f Assemble,MapReads,FindORFS,FindScaffoldORFS,Abundance,Classify,Annotate
