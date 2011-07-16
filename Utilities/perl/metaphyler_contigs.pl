#!/usr/bin/perl

#####################################################
#
# Program: MetaPyler Version 2.12
#
#          Estimate bacterial composition from
#          metagenomic sequences. 
#          
#
# Contact: Bo Liu, boliu@umiacs.umd.edu
#
# Last Update: Mon Feb 14 10:58:43 EST 2011
#
#####################################################

use strict;
use warnings;


#---------------------------------------------------#
# read command line parameters
#---------------------------------------------------#
my $blastfile = "";       # input blastx file
my $prefix = "";          # output files prefix
my $contigCov = "";
if (scalar @ARGV == 3) {
    $blastfile = $ARGV[0];
    $prefix = $ARGV[1];
    $contigCov = $ARGV[2];
} else {
    Usage();
}
#---------------------------------------------------#


#---------------------------------------------------#
# data files
#---------------------------------------------------#
my $dir = "/fs/szasmg3/boliu/blPrograms/pipeline/metaphyler/data/";
my $gtaxfile = "$dir/id2tax.tab";
my $id2mfile = "$dir/id2m.tab";
my $cutofffile = "$dir/cutoff.tab";
my $lenfile = "$dir/markerslen.tab";
#---------------------------------------------------#


#---------------------------------------------------#
# taxonomic information for reference marker genes 
#---------------------------------------------------#
# $gtax{reference gene id}{tax level} = tax name
my %gtax = ();   
open(FH, $gtaxfile) or die("Could not open file $gtaxfile\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($acc, $s, $g, $f, $o, $c, $p) = split("\t", $line);
    $gtax{$acc}{"genus"} = $g;
    $gtax{$acc}{"family"} = $f;
    $gtax{$acc}{"order"} = $o;
    $gtax{$acc}{"class"} = $c;
    $gtax{$acc}{"phylum"} = $p;
}
#---------------------------------------------------#



#---------------------------------------------------#
# marker gene names for reference genes
#---------------------------------------------------#

# $g2m{reference gene id} = marker gene name
my %g2m = ();
open(FH, "$id2mfile") or die("Could not open file $id2mfile\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($acc, $m) = split("\t", $line);
    $g2m{$acc} = $m;
}
close FH;
#---------------------------------------------------#


#---------------------------------------------------#
# read in classification cutoffs
#---------------------------------------------------#
# $cuts{reference gene id}{tax lev}{nucleotide position} = cutoff
my %cuts = ();
open(FH, "$cutofffile") or die("Could not open file $cutofffile\n");
my $pos = 0;
my $id = "";
foreach my $line (<FH>) {
    if ($line =~ /^\>\>(\d+)/) {
	$pos = $1;
    }
    elsif ($line =~ /^\>(\S+)/) {
	$id = $1;
    }
    elsif ($line =~ /^(\S+)\t(\S+)/) {
	$cuts{$id}{$1}{$pos} = $2;
    }
}
close FH;
#----------------------------------------------------#


#---------------------------------------------------#
# read in BLASTN mapping information
#---------------------------------------------------#

# taxonomic levels
my @tlevs = ("genus", "family", "order", "class", "phylum");

# $hitsloc{$query seq id} = center HSP locus in reference of the best hit
my %hitsloc = ();

# $hits{query seq id} = [best hit ids]
my %hits = ();


# $hitspct{query seq id} = %similarity of the best hit
my %hitspct = ();

#  $bestbits{query seq id} = bit score of the best hit
my %bestbits = ();

# $readhspl{query seq id} = HSP length of the best hit
my %readhspl = ();

open(FH, $blastfile) or die("Could not open file $blastfile\n");
foreach my $line (<FH>) {

    chomp $line;
    my ($qacc, $hacc, $pct, $hspl, $mis, $insert, $qs, $qe, $s, $e, $evalue, $bits) = split("\t", $line);

    # this is not the top 1 hit
    if (exists $bestbits{$qacc}) {

	# if %similarity < 98, then only use the best hit
	# even if this hit is as good as the best one
	if ($pct < 98) { next;}

	# test if this hit is as good as the best hit
	if ($bits != $bestbits{$qacc}) { next;}
    }

    # this is the top 1 hit
    else {
	$bestbits{$qacc} = $bits;
    }

    # ignore if HSP length < 60bp
    if ($hspl < 60) { next;}

    # ignore if %similarity < 80
    if ($pct <= 90) { next;}

    # this blast hit passes all previous criteria
    # store its info.
    push(@{$hits{$qacc}}, $hacc);
    $hitspct{$qacc} = $pct;
    $readhspl{$qacc} = $hspl;
    unless (exists $hitsloc{$qacc}) {
	$hitsloc{$qacc} = int (($e + $s) / 2);
    }
}
close FH;
#---------------------------------------------------#


#---------------------------------------------------#
# classify query sequences
#---------------------------------------------------#
# $read2tax{query seq id}{tax lev} = tax name
my %read2tax = ();

# $read2m{query seq id} = marker gene name
my %read2m = ();

# $read2tag{query seq id} = classification approach
my %read2tag = ();

foreach my $rid (keys %hits) {

    # %similarity of the best hit
    my $bestpct = $hitspct{$rid};

    # only one best hit
    if (scalar @{$hits{$rid}} == 1) {

	# best hit ref id
	my $refid = $hits{$rid} -> [0];

	# if %similarity of the single best hit > 98%
	# classify it directly
	if ($bestpct >= 98) {
	    foreach my $tlev (@tlevs) {
		$read2tax{$rid}{$tlev} = $gtax{$refid}{$tlev};
		$read2m{$rid} = $g2m{$refid};
		$read2tag{$rid} = "best";
	    }
	}

	# if %similarity of the single best hit < 98%
	# classify it using learned cutoff values
	else {

	    # tracking if this query seq is classified
	    my $tag = 2;

	    # classified at this tax lev
	    my $levnum = 0;

	    # locus of the HSP
	    my $loc = $hitsloc{$rid};

	    for ($levnum = 0; $levnum < 5; $levnum++) {

		my $tlev = $tlevs[$levnum];

		# exists classification cutoffs
		if (exists $cuts{$refid}{$tlev}) {

		    # using cutoffs from nearest neighborhood within 30bp
		    for (my $i = 0; $i <= 30; $i++) {
			my $rightloc = $i + $loc;
			my $leftloc = $loc - $i;
			if (exists $cuts{$refid}{$tlev}{$rightloc}) {
			    if ($bestpct >= $cuts{$refid}{$tlev}{$rightloc}) {
				$tag = 1;
			    }
			    else {
				$tag = 0;
			    }
			    last;
			}

			if (exists $cuts{$refid}{$tlev}{$leftloc}) {
			    if ($bestpct >= $cuts{$refid}{$tlev}{$leftloc}) {
				$tag = 1;
			    }
			    else {
				$tag = 0;
			    }
			    last;
			}
		    }
		}

		# if classified then go out
		if ($tag == 1) { last;}
	    }

	    # classified
	    if ($tag == 1) {
		$read2tag{$rid} = "cutoff";    
		$read2m{$rid} = $g2m{$refid};

		for (my $i = 0; $i < 5; $i++) {
		    my $tlev = $tlevs[$i];

		    # classified at tax level higher than i
		    if ($i < $levnum) {
			$read2tax{$rid}{$tlev} = "NA";
		    }
		    else {
			$read2tax{$rid}{$tlev} = $gtax{$refid}{$tlev};
		    }
		}
	    }
#	    elsif ($bestpct >= 90) {
	    # no cutoff available
	    elsif ($tag == 2) {
		$read2tag{$rid} = "bestphylum";    
		$read2m{$rid} = $g2m{$refid};
		for (my $i = 0; $i < 4; $i++) {
		    my $tlev = $tlevs[$i];
		    $read2tax{$rid}{$tlev} = "NA";
		}
		$read2tax{$rid}{"phylum"} = $gtax{$refid}{"phylum"};
	    }
	}
    } # if (scalar @{$hits{$rid}} == 1)

    # more than 1 best hits
    # then classify using consensus approach
    else {

	my $firstrefid = $hits{$rid} -> [0];
	$read2tag{$rid} = "consensus";    
	$read2m{$rid} = $g2m{$firstrefid};

	# match tax names at each level
	foreach my $tlev (@tlevs) {
	    my $pretname = $gtax{$firstrefid}{$tlev};
	    my $tag = 0;
	    foreach my $refid (@{$hits{$rid}}) {
		my $tname = $gtax{$refid}{$tlev};

		# if one is NA then, say they agree
		if ($tname eq "NA") { next;}
		if ($pretname eq "NA") {
		    $pretname = $tname;
		    next;
		}

		# do not agree
		if ($pretname ne $tname) {
		    $tag = 1;
		}
	    }
	    
	    if ($tag == 1) {
		$read2tax{$rid}{$tlev} = "NA";
	    }
	    else {
		$read2tax{$rid}{$tlev} = $pretname;
	    }
	}
    } # end of more than 1 best hits
}
#---------------------------------------------------#


#---------------------------------------------------#
# Classification refinement 1:
# Using consensus of the top best hits,
# query can not be classified at the genus level.
# Then, we classify the query to one of the best hits
# randomly using weights.
#---------------------------------------------------#
# $refid2count{ref id} = number of queries mapped
my %refid2count = ();
foreach my $rid (keys %read2tax) {

    my $genus = $read2tax{$rid}{"genus"};

    # this is the case to be improved, so ignore
    if ($genus eq "NA" && $read2tag{$rid} eq "consensus") { next;}
    my $refid = $hits{$rid} -> [0];
    $refid2count{$refid}++;
}

foreach my $rid (keys %read2tax) {
    
    my $genus = $read2tax{$rid}{"genus"};

    # improve this classification
    if ($genus eq "NA" && $read2tag{$rid} eq "consensus") {

	# store weighted tax names
	my %abund = ();
	foreach my $refid (@{$hits{$rid}}) {
	    if (exists $refid2count{$refid}) {
		$abund{$refid} = $refid2count{$refid};
	    }
	}

	# if all features do not have weights,
	# then randomly assign using equal weights.
	if (scalar keys %abund == 0) {
	    foreach my $refid (@{$hits{$rid}}) {
		$abund{$refid} = 1;
	    }
	}
	
	my $randomid = randomSelect(\%abund);
	foreach my $tlev (@tlevs) {
	    $read2tax{$rid}{$tlev} = $gtax{$randomid}{$tlev};
	}
	$read2tag{$rid} = "random";
    }
}
#---------------------------------------------------#


#---------------------------------------------------#
# Classification refinement 2:
#---------------------------------------------------#
my %refid2tax = ();
foreach my $rid (keys %read2tax) {
    foreach my $tlev (@tlevs) {
	if ($read2tax{$rid}{$tlev} ne "NA") {
	    foreach my $refid (@{$hits{$rid}}) {
		$refid2tax{$refid}{$tlev} = $read2tax{$rid}{$tlev};
	    }
	}
    }
}

foreach my $rid (keys %read2tax) {
    my $genus = $read2tax{$rid}{"genus"};
    my $refid = $hits{$rid} -> [0];
    
    if ($genus eq "NA" && ($read2tag{$rid} eq "cutoff" || $read2tag{$rid} eq "bestphylum")) {
	foreach my $tlev (@tlevs) {
	    if (exists $refid2tax{$refid}{$tlev}) {
		$read2tax{$rid}{$tlev} = $refid2tax{$refid}{$tlev};
	    }
	}
    }
}
#---------------------------------------------------#


#---------------------------------------------------#
# length of ref genes
#---------------------------------------------------#
my %len = ();
open(FH, "$lenfile") or die("Could not open file $lenfile\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($len, $id) = split("\t", $line);
    $len{$id} = $len;
}
close FH;
#---------------------------------------------------#


#---------------------------------------------------#
# output classification results
#---------------------------------------------------#
# $taxcounts{tax lev} {tax name} = number of reads
my %taxcounts = ();

# $dcov{tax lev} {tax name}{marker gene name} = depth of coverage
my %dcov = ();

my %genus2pct = ();

# $ctgcov{contig id} = [read coverage/depth]
my %ctgcov = ();

foreach my $rid (keys %read2tax) {
    $ctgcov{$rid} = 1;
}
open(CC, $contigCov) or die("Could not open file $contigCov\n");
foreach my $line (<CC>) {

    chomp $line;
    my ($cid, $ccov) = split("\t", $line);
    $ctgcov{$cid} = $ccov;


open(OUTPUT, ">$prefix.classify.tab");
print OUTPUT "\@contig ID\t\@gene name\t\@best BLAST hit\t\@percent similarity\t\@classification\t\@genus\t\@family\t\@order\t\@class\t\@phylum\n";
foreach my $rid (keys %read2tax) {
    my $refid = $hits{$rid} -> [0];
    my $m = $read2m{$rid};
    my $cov = $readhspl{$rid} / $len{$refid};

    $cov = $cov * $ctgcov{$rid};

    my $pct = $hitspct{$rid};
    my $tag = $read2tag{$rid};
    print OUTPUT "$rid\t$m\t$refid\t$pct\t$tag\t";
    
    my $prename = "NA";
    my $prelev = "phylum";
    my @tnames = ();
    for (my $i = 4; $i >= 0; $i--) {
	my $tlev = $tlevs[$i];
	my $tname = $read2tax{$rid}{$tlev};

	# replace name NA with higher tax names
	if ($tname eq "NA") {
	    $tname = $prename . "{$prelev}";
	}
	else {
	    $prelev = $tlevs[$i];
	    $prename = $tname;
	}
	$taxcounts{$tlev}{$tname}++;
	$dcov{$tlev}{$tname}{$m} += $cov;
	push(@tnames, $tname);
    }
    push(@{$genus2pct{$tnames[4]}}, $pct);
    print OUTPUT join("\t", reverse(@tnames)), "\n";
}
close OUTPUT;

#---------------------------------------------------#


#---------------------------------------------------#
# output taxonomic profile at each level 
#---------------------------------------------------#
# average similarity at the genus level
my %genus2avepct = ();
foreach my $genus (keys %genus2pct) {
    my $sum = 0;
    foreach my $pct (@{$genus2pct{$genus}}) {
	$sum += $pct;
    }
    $genus2avepct{$genus} = $sum / scalar @{$genus2pct{$genus}};
}

foreach my $tlev (@tlevs) {
    open(OUTPUT, ">$prefix.$tlev.tab");
    if ($tlev eq "genus") {
	print OUTPUT "\@$tlev\t\@% abundance\t\@depth of coverage\t\@number of reads\t\@similarity with reference\n";
    }
    else {
	print OUTPUT "\@$tlev\t\@% abundance\t\@depth of coverage\t\@number of reads\n";
    }
    my %tcov = ();
    my $covsum = 0;
    foreach my $tax (keys %{$dcov{$tlev}}) {

	my @depths = values %{$dcov{$tlev}{$tax}};
	my $num = scalar @depths;
	
	my $cov = median(@depths) * $num / 29;

	$tcov{$tax} = $cov;
	$covsum += $cov;
    }
    foreach my $tax (sort {$tcov{$b} <=> $tcov{$a}} keys %tcov) {
	my $cov = $tcov{$tax};
	my $pct = int ($cov * 10000/ $covsum) / 100;

	# if % abundance < 0.01%, ignore
#	if ($pct == 0) { last;}
	print OUTPUT "$tax\t$pct\t";
	print OUTPUT int ($cov * 100) / 100, "\t";
	print OUTPUT $taxcounts{$tlev}{$tax}, "\t";
	if ($tlev eq "genus") {
	    print OUTPUT int ($genus2avepct{$tax} * 100) / 100, "\n";
	}
	else {
	    print OUTPUT "\n";
	}
    }
    close OUTPUT;
}
#---------------------------------------------------#


exit;


sub Usage {
    die("
USAGE:
       perl metaphyler_contigs.pl <BLAST file> <Output prefix> <CONTIG COVERAGE file>

INPUT:
       BLAST file:
           BLASTN mappings of query sequences against
           MetaPhyler reference marker genes.

       Output prefix:
           6 output files will be generated:\"prefix\".classify.tab,
           \"prefix\".genus.tab, \"prefix\".family.tab,
           \"prefix\".order.tab, \"prefix\".class.tab, \"prefix\".phylum.tab
OUTPUT:
       \"prefix\".classify.tab
           column 1: query sequence id
           column 2: phylogenetic marker gene name
           column 3: best reference gene hit
           column 4: % similarity with best hit
           column 5: classification rule
           column 6: taxonomic label at genus level
           column 7: taxonomic label at family level
           column 8: taxonomic label at order level
           column 9: taxonomic label at class level
           column 10: taxonomic label at phylum level

       \"prefix\".genus|family|order|class|phylum.tab
           column 1: taxonomic clade name
           column 2: % relative abundances
           column 3: depth of coverage of genomes
           column 4: number of sequences binned to this clade
           column 5: similarity with reference genes (only available at the genus level)

NOTES:
       1, Here we assume the query sequences have been BLASTNed against the
          reference marker genes, which is used as the input of this program.
          To run MetaPhyler directly for a set of query sequences, use program
          metaphylerall.pl

EXAMPLE:
       perl metaphyler.pl ./test/test.blastn test

");
}


# given a set of features and associated weights
# randomly return a features
sub randomSelect {
    my ($hash) = @_;
    my @feas = sort keys %{$hash};
    my $sum = 0;
    foreach my $value (values %{$hash}) {
	$sum += $value;
    }
    my $rnum = int(rand($sum));
    my $j = 0;
    foreach my $fea (@feas) {
	$j += $hash -> {$fea};
	if ($rnum < $j) {
	    return $fea;
	}
    }
}


# given an array of numbers
# return the median
sub median {
    my (@values) = @_;
    @values = sort {$a <=> $b} @values;
    my $num = scalar @values;
    if ($num == 1) {
	return $values[0];
    }
    if ($num % 2 == 0) {
	return ($values[$num / 2 - 1] + $values[$num / 2]) / 2;
    }
    else {
	return $values[($num - 1) / 2];
    }
}
}
