#!/usr/bin/perl

$filenameIn = $ARGV[0];
$filename1 = $ARGV[1];
$filename2 = $ARGV[2];

open $FILEA, "< $filenameIn";

open $OUTFILE1, "> $filename1";
open $OUTFILE2, "> $filename2";

while(<$FILEA>) {
	print $OUTFILE1 $_;
	$_ = <$FILEA>;
	print $OUTFILE1 $_; 
	$_ = <$FILEA>;
	print $OUTFILE1 $_; 
	$_ = <$FILEA>;
	print $OUTFILE1 $_; 

	$_ = <$FILEA>;
	print $OUTFILE2 $_; 
	$_ = <$FILEA>;
	print $OUTFILE2 $_;
	$_ = <$FILEA>;
	print $OUTFILE2 $_;
	$_ = <$FILEA>;
	print $OUTFILE2 $_;
}
