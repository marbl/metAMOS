#!/usr/bin/perl
use strict;
use warnings;

my $dir = "./out/";
my $input = $ARGV[0];#"./in/s12.fa";
my $output = $ARGV[1];#"./out/s12.blastn";
my $db = $ARGV[2];#'./in/markers.fna';
my $command = "blastall -p tblastn -a 20 -F F -m 8 -i $input -d $db -e 0.0001 -b 20 -v 20 -o $output";
system("$command");
