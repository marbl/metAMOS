#bzip2 -d test4.sff.bz2
../initPipeline -f -m combined.fa -d test_soap_mg -i 200:500 
../runPipeline -c fcp -p 16 -d test_soap_mg -k 25 -t metamos -f Preprocess,Abundance
