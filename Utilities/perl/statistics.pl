#!/usr/bin/perl

# use lib "...";  # where AMOS::ParseFasta lives if not in standard location

use AMOS::ParseFasta;
use Getopt::Long;
use Statistics::Descriptive;
use strict;


my $prefix = $ARGV[0];

my $nohead = 0;
my $justhead = undef;
my $limit = undef;
my $N50Base = undef;
my $hmp = undef;
my $help = undef;
GetOptions(
    'h|help' => \$help, # print help
    'nohead' => \$nohead, # don't print a header
    'justhead' => \$justhead, # just print the header
    'limit=i' => \$limit, #ignore contigs less than
    'hmp'=>\$hmp, # if provided assumes HMP file setup
    'n50base=i' => \$N50Base); # use this value when computing N50

if (defined $help) {
print STDERR "Usage: statistics.pl [options] PREFIX > statistics.csv\n";
print STDERR q~
Available options:
-h|--help    -- print this message
--hmp        -- assume input is in HMP format (PREFIX.contigs.fa, PREFIX.scaffolds.fa)
--nohead     -- do not print header line
--justhead   -- just print header line
--n50base NN -- use NN as the base genome size for computing the N50 value
--limit NN   -- only report results for contigs >= NN bp
PREFIX       -- if --hmp selected results are reported on the files 
                PREFIX.contigs.fa and PREFIX.scaffolds.fa
                otherwise PREFIX is the name of a fasta file containing the
                contigs or scaffolds for which you want statistics computed

NOTE: files can be compressed with gzip or bzip2 and the program automatically
decompresses them.
~;
exit(0);
}

if ($nohead == 0){
    print "File\tNumber\tTotal Size\tMin Size\tMax Size\tAverage Size\tMedian Size\tN50\tSize @ 1Mbp\tNumber @ 1Mbp\tSize @ 2Mbp\tNumber @ 2Mbp\tSize @ 4Mbp\tNumber @ 4Mbp\tSize @ 10Mbp\tNumber @ 10Mbp\n";
  if (defined $justhead){exit(0);}
}

my @files = ();
if (defined $hmp){
 foreach my $suffix('contigs.fa', 'scaffolds.fa') {
   push @files, "$ARGV[0].$suffix";
 }
} else {
  push @files, $ARGV[0];
}
foreach my $file (@files)
{
    # print filename
    print "$file\t";

    if ($file =~ /\.gz$/){ # gzipped file
       open(IN, "gzip -dc $file |") || die ("Cannot open $file\n");
    } elsif ($file =~ /\.bz2$/) { #bzipped file
       open(IN, "bzip2 -dc $file |") || die ("Cannot open $file\n");
    } else {
       open(IN, "$file") || die ("Cannot open $file\n");
    }
    my $fr = new AMOS::ParseFasta(\*IN);
    if (! defined $fr){ print STDERR ("Cannot parse file\n"); exit 0;}

    my @sizes = ();
    my $stats = new Statistics::Descriptive::Full;

    while (my ($head, $body) = $fr->getRecord()){
        if (defined $limit && length($body) < $limit) { next;}
	push @sizes, length($body);
	$stats->add_data(length($body));
    }
    close(IN);

    print $stats->count(), "\t"; # count 
    print $stats->sum(), "\t"; # total size
    print $stats->min(), "\t"; # min size
    print $stats->max(), "\t"; # max size
    print sprintf("%.2f\t",$stats->mean()); # mean size
    print $stats->median(), "\t"; # median size
    
    # now for the N* statistics
    if (! defined $N50Base){
	$N50Base = $stats->sum() / 2;
    }
    my $total = 0;
    my @sizes = sort {$b <=> $a} @sizes;

    my $n50 = undef;
    my $n1m = undef;
    my $n2m = undef;
    my $n4m = undef;
    my $n10m = undef;
    my $s1m = undef;
    my $s2m = undef;
    my $s4m = undef;
    my $s10m = undef;

    for (my $i = 0; $i <= $#sizes; $i++){
	$total += $sizes[$i];
	if ($total >= $N50Base && ! defined $n50){
	    $n50 = $sizes[$i];
	}
	if ($total >= 1000000 && ! defined $n1m){
	    $n1m = $i + 1;
	    $s1m = $sizes[$i];
	}
	if ($total >= 2000000 && ! defined $n2m){
	    $n2m = $i + 1;
	    $s2m = $sizes[$i];
	}
	if ($total >= 4000000 && ! defined $n4m){
	    $n4m = $i + 1;
	    $s4m = $sizes[$i];
	}
	if ($total >= 10000000 && ! defined $n10m){
	    $n10m = $i + 1;
	    $s10m = $sizes[$i];
	}
    }

    print "$n50\t$s1m\t$n1m\t$s2m\t$n2m\t$s4m\t$n4m\t$s10m\t$n10m\n";
}
