../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_blast -i 3000:4000 
../runPipeline -r -v -c blast -p 8 -a soap -d test_blast -k 55 -f Preprocess,Assemble,Annotate,Postprocess -z phylum
