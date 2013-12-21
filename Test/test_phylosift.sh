#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_phylosift -i 3000:4000 
../runPipeline -r -v -c phylosift -p 8 -a soap -d test_phylosift -k 55 -f Preprocess -z phylum
