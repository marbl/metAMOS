#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_phymm -i 3000:4000 
../runPipeline -r -v -c phymm -p 8 -a soap -d test_phymm -k 55 -t -f Annotate,Postprocess -z phylum
