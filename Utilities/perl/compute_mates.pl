#!/usr/bin/perl

# Computes mate-pair statistics from Bowtie output

use Getopt::Long;
use Statistics::Descriptive;
use strict;

# Location of bowtie
my $bowtie = `which bowtie`;
chomp $bowtie;

print STDERR "Usage: compute_mates.pl -d[atabase] BOWTIEDB [-trim TT] [-threads NN] [-limit NN] mate1.fastq[.gz] mate2.fastq[.gz]\n";
print STDERR "BOWTIEDB - location of bowtie database\n";
print STDERR "--trim TT - exclude TT fraction of outliers\n";
print STDERR "--threads NN - use NN threads when running bowtie\n";
print STDERR "--limit NN - only align first NN reads in each file\n";
print STDERR "\nNOTE: if the output files called *.bout exist bowtie is not run\n";

my $trim = undef;
my $Database = undef;
my $threads = undef;
my $limit = undef;

GetOptions(
    "d|database=s" => \$Database, # database
    "trim=f" => \$trim,
    "threads=i" => \$threads,
    "limit=i" => \$limit
    );

if (! defined $Database){
    die ("Please provide a bowtie database with -d option\n");
}

my $mate1 = $ARGV[0];
my $mate2 = $ARGV[1];

my $bowtie_param = "-v 2 --suppress 6,7,8";
if (defined $threads) {$bowtie_param .= " -p $threads";}
if (defined $limit) {$bowtie_param .= " -u $limit";}
$bowtie_param .= " $Database";

my $files = "";
# run bowtie first
foreach my $file ($mate1, $mate2){
    my $todo = $file;
    $todo =~ /(.*)\.\S+/;
    my $prefix = $1;
    if ($file =~ /(.*)\.gz/){
	$todo = $1;
	$todo =~ /(.*)\.\S+/;
	$prefix = $1;
	system("gzip -dc $file > $todo") unless (-f "$prefix.bout");
    }
    if ($file =~ /(.*)\.bz2/){
	$todo = $1;
	$todo =~ /(.*)\.\S+/;
	$prefix = $1;
	system("bzip2 -dc $file > $todo") unless (-f "$prefix.bout");
    }

    $files .= "$prefix.bout ";

    if (-f "$prefix.bout"){
	print STDERR "WARNING $prefix.bout found! Bowtie will not run.  Remove $prefix.bout if you want bowtie to execute!\n";
    } else {
	print STDERR "Running $bowtie $bowtie_param $todo | sed \'s/\\/[12]//\' > $prefix.bout\n";
	system("$bowtie $bowtie_param $todo | sed \'s/\\/[12]//\' > $prefix.bout");
    }
    if ($todo ne $file){
	unlink $todo;
    }
}

# read the paired alignments
print STDERR "Running join -1 1 -2 1 $files\n";
open(ALIGNS, "join -1 1 -2 1 $files |") || die ("Cannot run join -1 1 -2 1 $files: $!\n");

my $stats = new Statistics::Descriptive::Full;
while (<ALIGNS>){
    my @fields = split(/\s+/, $_);
    if ($fields[1] eq $fields[5]) {next;} # both reads in same orientation
    if ($fields[2] ne $fields[6]) {next;} # match different reference
    my $min1 = $fields[3];
    my $max1 = $fields[3] + length($fields[4]);
    my $min2 = $fields[7];
    my $max2 = $fields[7] + length($fields[8]);
    my $min = ($min1 < $min2) ? $min1 : $min2;
    my $max = ($max1 > $max2) ? $max1 : $max2;
print "Min: $min, Max: $max\n";
    my $size = $max - $min;
    $stats->add_data($size);
}
close(ALIGNS);

# unlink $files;

my $mean;
my $stdev;
if (defined $trim){
    $mean = $stats->trimmed_mean($trim);
} else {
    $mean = $stats->mean();
}
$stdev = $stats->standard_deviation();

print "Number\tMean\t10%Mean\tStandard deviation\n";
print $stats->count(), "$mean\t", $stats->trimmed_mean(0.1), "\t$stdev\n";
