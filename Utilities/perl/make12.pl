#!/usr/bin/perl

my $master = $ARGV[0];
my $file1 = $ARGV[1];
my $file2 = $ARGV[2];
my $prefix = $ARGV[3];

print STDERR "make12 master in1 in2 outprefix\n";

open(M, $master) || die ("Cannot open $master: $!\n");
open(IN1, $file1) || die ("Cannot open $file1: $!\n");
open(IN2, $file2) || die ("Cannot open $file2: $!\n");
open(BOTH1, ">$prefix.1.fastq") || die ("Cannot open $prefix.1.fastq: $!\n");
open(BOTH2, ">$prefix.2.fastq") || die ("Cannot open $prefix.2.fastq: $!\n");
open(UN, ">$prefix.un.fastq") || die ("Cannot open $prefix.un.fastq: $!\n");


$line1 = <IN1>;
$line2 = <IN2>;
while (<M>){
    my $head = $_;
    chomp $head;
#    if ($head !~ /^@/){
#	die ("Wrong line: $head");
#    }
#    $head =~ s/@//;
#    $head =~ s/(.*)\/\d$/\1/; # get rid of sequence number
#    <M>;<M>;<M>; # skip rest of header
    $head1 = $line1;
    chomp $head1;
    $head1 =~ s/^@//;
    $head1 =~ s/(.*)\/\d$/\1/;
    $head2 = $line2;
    chomp $head2;
    $head2 =~ s/^@//;
    $head2 =~ s/(.*)\/\d$/\1/;

    if ($head1 eq $head && $head2 eq $head) {# going into both
	print BOTH1 $line1;
        print BOTH2 $line2;
        $line1=<IN1>; print BOTH1 $line1;
        $line1=<IN1>; print BOTH1 $line1;
        $line1=<IN1>; print BOTH1 $line1;
        $line1=<IN1>;
        $line2=<IN2>; print BOTH2 $line2; 
        $line2=<IN2>; print BOTH2 $line2; 
        $line2=<IN2>; print BOTH2 $line2; 
        $line2=<IN2>;
    } elsif ($head1 eq $head) { # only going into 1
        print UN $line1;
        $line1=<IN1>; print UN $line1;
        $line1=<IN1>; print UN $line1;
        $line1=<IN1>; print UN $line1;
        $line1=<IN1>;
    } elsif ($head2 eq $head) { # only going into 2
        print UN $line2;
        $line2 = <IN2>; print UN $line2;
        $line2 = <IN2>; print UN $line2;
        $line2 = <IN2>; print UN $line2;
        $line2=<IN2>;
    } else {
        print "$head $head1 $head2 do not match\n";
    }
}
close(IN1); close(IN2); close(M); 
close(BOTH1); close(BOTH2); close(UN);

