
system("./SOAPdenovo all -M 3 -s config.txt -o test");

# step by step
system("./SOAPdenovo pregraph -s config.txt -K 23 -o idmates");

system("./SOAPdenovo contig -M 3  -g idmates");
#system("./SOAPdenovo map -s config.txt -g test");
