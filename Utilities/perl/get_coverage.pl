#!/usr/bin/perl

# Estimates singletons from assemblies 

use Getopt::Long;
use Statistics::Descriptive;
use strict;
use AMOS::ParseFasta;

# Location of bowtie
my $bowtie = `which bowtie`;
my $bowtiebuild = `which bowtie-build`;
chomp $bowtie;
chomp $bowtiebuild;

print STDERR "Usage: get_coverage.pl -reads READDIR -assembly ASMDIR [-nobowtie] [--threads NN]\n";
print STDERR "READDIR - location of raw reads\n";
print STDERR "ASMDIR - location of assembly output (in PGA format)\n";
print STDERR "--threads NN - use NN threads when running bowtie\n";
print STDERR "-nobowtie  - do not run bowtie (assumes  the output files called *.bout exist)\n";
print STDERR "\nNOTE: assumes data are in the PGA format as stored at the HMP DACC, i.e. the assembly directory contains a file named PREFIX.contigs.fa\n";

my $readdir = undef;
my $asmdir = undef;
my $threads = undef;
my $trim = 25;
my $nobowtie = undef;

GetOptions(
    "reads=s" => \$readdir, # database
    "assembly=s" => \$asmdir, 
    "trim=f" => \$trim,
    "threads=i" => \$threads,
    "nobowtie" => \$nobowtie,
    );

if (! defined $readdir){
    die ("Please provide a directory where the reads are located\n");
}

if (! defined $asmdir) {
    die ("Please provide a directory where the assembly is located\n");
}


my $bowtie_param = "-v 1 -M 2 --suppress 1,2,5,6,7,8";
if (defined $threads) {$bowtie_param .= " -p $threads";}

opendir(ASM, $asmdir) || die ("Cannot open $asmdir: $!\n");
my @files = grep {/^.*\.contigs.fa$/} readdir(ASM);
closedir(ASM);

# build bowtie index
if (! -e "$asmdir/IDX.1.ebwt" && ! defined $nobowtie) {system "$bowtie-build $asmdir/$files[0] $asmdir/IDX"};

opendir(RDS, $readdir) || die ("Cannot open $readdir: $!\n");
my @rds = grep {/^.*\.fastq.*/} readdir(RDS);
closedir(RDS);
if (defined $nobowtie){goto NOBO;}
foreach my $file (@rds){
    # first trim to 25bp
    $file =~ /(.*\.denovo_duplicates_marked.trimmed\.\w+).*/;
    my $prefix = $1;
    
    open(TRIM, ">$readdir/$prefix.trim.fastq") || die ("Cannot open $readdir/$prefix.trim.fastq: $!\n");

    if ($file =~ /.*\.gz/){
	open(IN, "gzcat $readdir/$file |") || die ("Cannot open $readdir/$file: $!\n"); 
    } elsif ($file =~ /.*\.bz2/){
	open(IN, "bzip2 -dc $readdir/$file |") || die ("Cannot open $readdir/$file: $!\n");
    } else {
	open(IN, "$readdir/$file") || die ("Cannot open $readdir/$file: $!\n");
    }

    while(<IN>){
	if ($. % 2 == 0){
	    print TRIM substr($_, 0, 25), "\n";
	} else {
	    print TRIM;
	}
    }
    close(IN);
    close(TRIM);
    
    # run bowtie
    if (! -f "$readdir/$prefix.bout") {system("$bowtie $bowtie_param $asmdir/IDX $readdir/$prefix.trim.fastq > $readdir/$prefix.bout")};
#    unlink("$readdir/$prefix.trim.fastq");
}
NOBO:
# now the bowtie results are there, simply report for each contig the # of
# reads that map to it


    
my %contigsize; 
my %contigreads; 
open(CON, "$asmdir/$files[0]") || die ("cannot open contigs $files[0]: $!\n");
my $pf = new AMOS::ParseFasta(\*CON);
print STDERR "Reading contigs from $asmdir/$files[0]\n";
my $ncontig = 0;
while (my ($head, $data) = $pf->getRecord()){
	$contigsize{$head} = length($data);
	$ncontig++;
}
close(CON);
print STDERR "Got $ncontig contigs\n";

my %done;
foreach my $file (@rds){
    # first trim to 25bp
    $file =~ /(.*\.denovo_duplicates_marked.trimmed\.\w+).*/;
    my $prefix = $1;
    if (! exists $done{"$readdir/$prefix.bout"}){
	open(IN, "$readdir/$prefix.bout") || die ("Cannot open bowtie out: $prefix.bout $!\n");
	print STDERR "Getting alignment information from $readdir/$prefix.bout\n";
	$ncontig = 0;
	while (<IN>){
	    my @fields = split('\t', $_);
#	print STDERR "Got $fields[0]\n";
	    $contigreads{$fields[0]}++;
	    $ncontig++;
	}
	close(IN);
	$done{"$readdir/$prefix.bout"} = 1;
	print STDERR "Got $ncontig alignments\n";
    }
}

print STDERR "Writing results to $asmdir/$files[0].cov\n";
open(OUT, ">$asmdir/$files[0].cov") || die ("Cannot open $files[0].cov: $!\n");
$ncontig = 0;
while (my ($name, $cnt) = each %contigreads){
  print OUT "$name\t$contigsize{$name}\t$cnt\t", sprintf("%.2f\n",$cnt * 100/$contigsize{$name});
  $ncontig++;
}
close(OUT);
print "Wrote $ncontig contigs\n";
exit(0);
