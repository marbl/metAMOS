#bzip2 -d test4.sff.bz2
../createProject -i 3000:4000 -sm test4.sff -d test_newbler 
../runPipeline -c amphora -a newbler -p 4 -d test_newbler -k 55 -t -f Preprocess
