#!/usr/bin/perl

$filenameIn = $ARGV[0];
$filename1 = $ARGV[1];
$filename2 = $ARGV[2];

open FILEA, "< $filenameIn";

open OUTFILE1, "> $filename1";
open OUTFILE2, "> $filename2";

my $lineA = "";
while(defined $lineA) {
        if ($lineA eq "") { 
           $lineA = <FILEA>;
        }
        print OUTFILE1 $lineA;
        $lineA = <FILEA>;
        while (defined $lineA && $lineA !~ m/>/) {
                print OUTFILE1 $lineA;
                $lineA = <FILEA>;
        }

        print OUTFILE2 $lineA;
        $lineA = <FILEA>;
        while (defined $lineA && $lineA !~ m/>/) {
                print OUTFILE2 $lineA;
                $lineA = <FILEA>;
        }
}
