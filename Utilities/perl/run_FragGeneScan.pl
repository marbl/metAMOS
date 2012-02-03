#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $genome_file = "";
my $FGS_result = "";
my $FGS_whole = -1;
my $FGS_train_file = "";
my $command;
my $program = $0;
my $dir = substr($0, 0, length($0)-19);
my $train_file;


GetOptions(
           'genome=s' => \$genome_file,
           'out=s' => \$FGS_result,
           'complete=s' => \$FGS_whole,
           'train=s' => \$FGS_train_file,
           );

if (length($genome_file)==0){
    print "ERROR: An input genome file was not specified.\n";
    print_usage();
    exit;
}elsif (! -e $genome_file){
    print "ERROR: The input genome file [$genome_file] does not exist.\n";
    print_usage();
    exit;
}

if (length($FGS_result) == 0 ){
    print "ERROR: An output file name was not specified.\n";
    print_usage();
    exit;
}

unless ($FGS_whole eq "1" || $FGS_whole eq "0"){
    print "ERROR: An incorrect value for the option -complete was entered.\n";
    print_usage();
    exit;
}  

$train_file = $dir."train/".$FGS_train_file;
if (length($FGS_train_file)==0){
    print  "ERROR: A file for model parameters was not specified.\n";
    print_usage();
    exit;
}elsif (! -e $train_file){
    print  "ERROR: The file for model parameter [$train_file] does not exist.\n";
    print_usage();
    exit;
}


$command = $dir."FragGeneScan";
$command .= " -s ".$genome_file;
$command .= " -o ".$FGS_result;
$command .= " -w ".$FGS_whole ;
$command .= " -t ".$FGS_train_file;
system($command); 

if ($FGS_whole eq "1"){
    system($dir."post_process.pl -genome=".$genome_file." -pre=".$FGS_result." -post=".$FGS_result.".out");
    system("rm ".$FGS_result);

}else{
    system("mv ".$FGS_result." ".$FGS_result.".out");
}


sub print_usage{

    print "USAGE: ./run_FragGeneScan.pl -genome=[seq_file_name] -out=[output_file_name] -complete=[1 or 0] -train=[train_file_name]\n";
    print "       [seq_file_name]:    sequence file name including the full path\n";
    print "       [output_file_name]: output file name including the full path\n";
    print "       [1 or 0]:           1 if the sequence file has complete genomic sequences\n";
    print "                           0 if the sequence file has short sequence reads\n";
    print "       [train_file_name]:  file name that contains model parameters; this file should be in the \"train\" directory\n";
    print "                           Note that four files containing model parameters already exist in the \"train\" directory\n"; 
    print "                           [complete] for complete genomic sequences or short sequence reads without sequencing error\n";
    print "                           [sanger_5] for Sanger sequencing reads with about 0.5% error rate\n";
    print "                           [sanger_10] for Sanger sequencing reads with about 1% error rate\n";
    print "                           [454_10] for 454 pyrosequencing reads with about 1% error rate\n";
    print "                           [454_30] for 454 pyrosequencing reads with about 3% error rate\n";
    print "                           [illumina_5] for Illumina sequencing reads with about 0.5% error rate\n";
    print "                           [illumina_10] for Illumina sequencing reads with about 1% error rate\n\n";
}
