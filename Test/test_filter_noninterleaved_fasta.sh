#bzip2 -d test4.sff.bz2
../initPipeline -f -m carsonella_pe_filt.fna -d test_filter -i 3000:4000 
../runPipeline -r -v -c kraken -p 8 -a soap -d test_filter -k 55 -t metamos -f FunctionalAnnotation,Postprocess -n Classify,Propagate -z phylum
