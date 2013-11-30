#bzip2 -d test4.sff.bz2
../initPipeline -f -m combined.tiny.fa -d test_soap_mg_tiny -i 200:500 
../runPipeline -c kraken -p 16 -d test_soap_mg_tiny -k 25 -t metamos -y -f Preprocess,Abundance
