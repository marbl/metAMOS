#!/usr/bin/perl


#####################################################
#
# Program: MetaPyler 
#          Taxonomic profiling for metagenomic sequences
#          
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Last Update: Fri Jul 16 15:59:20 EDT 2010
#
#####################################################


use strict;
use warnings;




#---------------------------------------------------#
# command line parameters
#---------------------------------------------------#
my $blastfile = "";       # input blastx file
my $prefix = "";          # output files prefix
my $contigCovfile = "";
my $outdir = ".";
my $rundir = ".";
if (scalar @ARGV == 5) {
    $blastfile = $ARGV[0];
    $prefix = $ARGV[1];
    $contigCovfile = $ARGV[2];
    $outdir = $ARGV[3];
    $rundir = $ARGV[4];
} else {
    Usage();
}

my $dir = "$rundir/config/";
#---------------------------------------------------#


#---------------------------------------------------#
# depth of coverage for contigs
#---------------------------------------------------#
my %contigcov = ();
open(FH, "$contigCovfile") or die("Could not open file $contigCovfile\n");
foreach my $line (<FH>) {
    chomp $line;
    
    my ($id, $cov) = split("\t", $line);
    $contigcov{$id} = $cov;
}
close FH;

#---------------------------------------------------#


#---------------------------------------------------#
# taxonomy information for reference marker genes 
#---------------------------------------------------#
# $gtax{gene id}{tax level} = tax name
my %gtax = ();   

# taxonomic info. for reference marker genes
my $gtaxfile = "$dir/markerstax.tab";
open(FH, $gtaxfile) or die("Could not open file $gtaxfile\n");
foreach my $line (<FH>) {
    $line =~ /^(\S+)\t(\S+)\t(.+?)\t(\w+)$/;
    my ($acc, $tid, $tname, $lev) = ($1, $2, $3, $4);
    $gtax{$acc}{$lev} = $tname;
}
#---------------------------------------------------#


#---------------------------------------------------#
# load classifiers
#---------------------------------------------------#

# $g2m{$gene} = marker name
my %g2m = ();

# $abund{tax level}{tax name} = number of reads mapped
my %abund = ();

# total number of reads mapped
my $numMapped = 0;

# classifier for each gene
my %gcut = ();

# classifier file
open(FH, "$dir/lr.tab") or die("Could not open file\n");
foreach my $line (<FH>) {
    chomp $line;
    
    # @values: $a, $b, $minlen, $level, $tid, $tname, $m
    my ($acc, @values) = split("\t", $line);
    my $a = $values[0];
    my $b = $values[1];
    my $lev = $values[3];
    $gcut{$acc}{$lev} = [$a, $b];
    $g2m{$acc} = $values[6];
}
close FH;
#----------------------------------------------------#


#---------------------------------------------------#
# taxonomic classification for query reads
#---------------------------------------------------#

# classify at these tax levels
my @tlevs = ("genus", "family", "order", "class", "phylum");

# $hits{query id} = [hit ids]
my %hits = ();

# $hitsbit{query id} = [hit bit scores]
my %hitsbit = ();

# $hitspct{query id} = [hit percentages]
my %hitspct = ();

# $hitshspl{query id} = [hsp lengths]
my %hitshspl = ();

# open input blastx file in -m 8 format
open(FH, $blastfile);
foreach my $line (<FH>) {
    chomp $line;
    # all lines should match
    unless ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+).*\s(\S+)$/) {
	print "This line does not match:\n";
	print $line;
	exit;
    }
    my ($qacc, $hacc, $pct, $hspl, $bits) = ($1, $2, $3, $4, $5);

    if (exists $hits{$qacc}) {

	my $prebits = $hitsbit{$qacc} -> [0];
	my $prepct = $hitspct{$qacc} -> [0];

	# equally good
	if ($bits == $prebits) {
	}
	else {
	    next;
	}
    }

    # ignore if HSP length < 60bp
    if ($hspl < 20) { next;}

    # ignore if bit score <= 35
    if ($bits <= 35) { next;}

    if ($pct <= 80) { next;}

    # ignore if previous hits have been used
    # and % identity of current hit < 98%
    
    # store blast information
    push(@{$hits{$qacc}}, $hacc);
    push(@{$hitspct{$qacc}}, $pct);
    push(@{$hitsbit{$qacc}}, $bits);
    push(@{$hitshspl{$qacc}}, $hspl);
}
close FH;
#---------------------------------------------------#


#---------------------------------------------------#
# classify query reads
#---------------------------------------------------#

# output classification file
open(OUTPUT, ">$outdir/$prefix.classify.txt");
select OUTPUT;

foreach my $rid (keys %hits) {

    
    # accession number of the best hit
    my $hacc = $hits{$rid} -> [0];

    # this contig does not have coverage information;
    # then ignore it.
    unless (exists $contigcov{$rid}) {
	#print STDOUT "$rid\n";
	next;
    }
    my $cov = $contigcov{$rid};
    #print STDOUT "$cov\n";

    # bit score of the best hit
    my $bestbit = $hitsbit{$rid} -> [0];

    # HSP length of the best hit
    my $hspl = $hitshspl{$rid} -> [0];

    # classify at each tax level
    foreach my $lev (@tlevs) {

	my $assign = 0;
	
	##################################
	# if top hits are equally good
	# then see if their taxonomic labels
	# agree with each other
	my $num = scalar @{$hits{$rid}};
	my @tids = ();
	for (my $i = 0; $i < $num; $i++) {
	    
	    my $refid = $hits{$rid} -> [$i];
	    my $pct = $hitspct{$rid} -> [$i];
	    if ($pct < 97) { next;}
	    if (exists $gtax{$refid}{"$lev"}) {
		push(@tids, $gtax{$refid}{$lev});
	    }
	}

	my $num2 = scalar @tids;
	for (my $i = 1; $i < $num2; $i++) {
	    my $preid = $tids[$i - 1];
	    my $id = $tids[$i];

	    # agree
	    if ($preid eq $id) {
		$assign = 1;
	    }

	    # do not agree
	    else {
		$assign = 2;
		last;
	    }
	}
	##################################
	
	# only one hit >= 98%
	# then classify to it
	if ($num2 == 1) {
	    $assign = 1;
	}

	if ($assign == 1 && $num ==1 && $hacc =~ /^tm7/) {
		$numMapped++;#=($numMapped*$cov);
		
		print "$rid\t$hacc\tNA\t";
		
		print "TM7[genus]\tTM7[family]\tTM7[order]\tTM7[class]\tTM7[phylum]\n";
#		$abund{"genus"}{"TM7[genus]"}++;
#		$abund{"family"}{"TM7[family]"}++;
#		$abund{"order"}{"TM7[order]"}++;
#		$abund{"class"}{"TM7[class]"}++;
#		$abund{"phylum"}{"TM7[phylum]"}++;
		$abund{"genus"}{"TM7[genus]"} += $cov*$hspl;
		$abund{"family"}{"TM7[family]"} += $cov*$hspl;
		$abund{"order"}{"TM7[order]"} += $cov*$hspl;
		$abund{"class"}{"TM7[class]"} += $cov*$hspl;
		$abund{"phylum"}{"TM7[phylum]"} += $cov*$hspl;
		next;
	}
	
	if ($num == 2) {
	    my $refid1 = $hits{$rid} -> [0];
	    my $refid2 = $hits{$rid} -> [1];
	    if (exists $gtax{$refid1}{"$lev"} && exists $gtax{$refid2}{$lev}) {
		if ($gtax{$refid1}{$lev} ne $gtax{$refid2}{$lev}) {
		    $assign = 2;
		}
	    }
	}
	
	# top hits do not agree
	# do not classify
	if ($assign == 2) { next;}
	

	# use bit score to classify
	if ($assign != 1) {
	    if (exists $gcut{$hacc}{"$lev"}) {
		my $a = $gcut{$hacc}{"$lev"} -> [0];
		my $b = $gcut{$hacc}{"$lev"} -> [1];
		my $cutvalue = $a + $b * $hspl;
		if ($bestbit >= $cutvalue) {
		    $assign = 1;
		}
	    }
	}

	# read is classified at $lev level
	if ($assign == 1) {

	    unless (exists $g2m{$hacc}) {
		next;
	    }

	    if ($g2m{$hacc} eq "rpoB") { next;}
	    if ($g2m{$hacc} eq "tsf") { next;}

	    # total number of reads mapped
	    #$numMapped++;
	    $numMapped+=($cov*$hspl);	    
	    print "$rid\t$hacc\t$g2m{$hacc}\t";
	    my $levtag = 0;


	    if ($lev eq "genus") {
		if (exists $gtax{$hacc}{$lev} && $gtax{$hacc}{$lev} eq "cyanobacteria") {
		    print "synechococcus\tNA\tchroococcales\tNA\tcyanobacteria\n";
		    $abund{"genus"}{"synechococcus"} += $cov*$hspl;
		    $abund{"family"}{"NA"} += $cov*$hspl;
		    $abund{"order"}{"chroococcales"} += $cov*$hspl;
		    $abund{"class"}{"NA"} += $cov*$hspl;
		    $abund{"phylum"}{"cyanobacteria"} += $cov*$hspl;
		    next;
		}
	    }
	    
	    
	    foreach my $tlev (@tlevs) {

		# classified at this level
		if ($tlev eq $lev) {
		    $levtag = 1;
		}

		# print out taxonomic information
		if ($levtag == 1) {

		    # has name at this tax level
		    if (exists $gtax{$hacc}{$tlev}) {
			my $tname = $gtax{$hacc}{$tlev};
			print "$tname\t";
			$abund{$tlev}{$tname} += $cov*$hspl;
		    }

		    # does not have name at this tax level
		    else {
			print "NA\t";
			$abund{$tlev}{"NA"} += $cov*$hspl;
		    }
		}

		# not classified at this level
		else {
		    print "NA\t";
		    $abund{$tlev}{"NA"} += $cov*$hspl;
		}
	    }
	    
	    print "\n";
	    last;
	}
    }
}
close OUTPUT;
select STDOUT;


#---------------------------------------------------#
# taxonomic profile at each level
#---------------------------------------------------#
# output taxonomic profile file
open(PCT, ">$outdir/$prefix.taxprof.pct.txt");
open(COUNT, ">$outdir/$prefix.taxprof.count.txt");
foreach my $lev (@tlevs) {
    
    print PCT ">$lev\n";
    print COUNT ">$lev\n";

    # sort by abundances from high to low
    foreach my $tname (sort {$abund{$lev}{$b} <=> $abund{$lev}{$a}} keys %{$abund{$lev}}) {
	my $abund = $abund{$lev}{$tname};

	# % abundance
	my $pct = $abund / $numMapped;

	# if % abundance <= 0.0005, ignore
	# if <= 2 reads classified at this tax level, ignore
	if ($tname eq "NA") {
#	    if ($lev eq "phylum") { next;}
	}
	else {
#	    if ($pct <= 0.001) { next;}
#	    if ($abund <= 1) { last;}
	}
	
	print PCT "$tname\t$pct\n";
	print COUNT "$tname\t$abund\n";
    }
    print PCT "\n";
    print COUNT "\n";
}
close PCT;
close COUNT;
#---------------------------------------------------#


exit;


sub Usage {
    die("
Usage:
       perl metaphyler.pl <BLAST file> <output prefix>

Output:
       \"prefix\".classify
           column 1: query read id
           column 2: reference gene id
           column 3: reference gene name
           column 4: taxonomic label at genus level
           column 5: taxonomic label at family level
           column 6: taxonomic label at order level
           column 7: taxonomic label at class level
           column 8: taxonomic label at phylum level

       \"prefix\".taxprof
           >taxonomic_level
           taxonomic_name percent_abundance

Notes:
       1, The input file is the BLASTP or BLASTX file (-m 8 format)
          of the metagenomic sequences against the reference marker genes.

Example:
       perl metaphyler.pl ./test/test.blastx test

");
}
