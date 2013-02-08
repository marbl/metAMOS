#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_ns.fna -d test_filter -i 3000:4000 
../runPipeline -r -v -c fcp -p 8 -a soap -d test_filter -k 55 -t -f FunctionalAnnotation,Postprocess -n Classify,Propagate -z phylum
