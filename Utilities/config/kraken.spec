[CONFIG]
input CONTIGS
name Kraken
output [PREFIX].hits
location cpp/[MACHINE]/kraken/bin
threads --threads
commands kraken [THREADS] -db [DB]/kraken --output [OUTPUT] [INPUT] 
[kraken]
-preload
