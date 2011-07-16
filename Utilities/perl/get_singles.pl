#!/usr/bin/perl

# Estimates singletons from assemblies 

use Getopt::Long;
use strict;

# Location of bowtie
my $resultDir = "./test/bowtie";
my $bowtie = `which bowtie`;
my $bowtiebuild = `which bowtie-build`;
chomp $bowtie;
chomp $bowtiebuild;
chomp $resultDir;

$resultDir .= "/";

print STDERR "Usage: get_singles.pl -reads READDIR -assembly ASMDIR [--threads NN]\n";
print STDERR "READDIR - location of raw reads\n";
print STDERR "ASMDIR - location of assembly output (in PGA format)\n";
print STDERR "--threads NN - use NN threads when running bowtie\n";
print STDERR "\nNOTE: if the output files called *.bout exist bowtie is not run\n";
print STDERR "\nNOTE: assumes data are in the PGA format as stored at the HMP DACC, i.e. the assembly directory contains a file named PREFIX.contigs.fa\n";

my $readdir = undef;
my $asmdir = undef;
my $threads = undef;
my $trim = 25;

GetOptions(
    "reads=s" => \$readdir, # database
    "assembly=s" => \$asmdir, 
    "trim=f" => \$trim,
    "threads=i" => \$threads
    );

if (! defined $readdir){
    die ("Please provide a directory where the reads are located\n");
}

if (! defined $asmdir) {
    die ("Please provide a directory where the assembly is located\n");
}


my $bowtie_param = "-v 1 -M 2";
if (defined $threads) {$bowtie_param .= " -p $threads";}

opendir(ASM, $asmdir) || die ("Cannot open $asmdir: $!\n");
#my @files = grep {/^.*\.contigs.fa/} readdir(ASM);
#closedir(ASM);
#print "$files[0]";
# build bowtie index
system "$bowtie-build $asmdir/contig.fa $resultDir/IDX";

opendir(RDS, $readdir) || die ("Cannot open $readdir: $!\n");
my @rds = grep {/^.*\.fastq.*/} readdir(RDS);
closedir(RDS);

foreach my $file (@rds){
    print "$file";
    # first trim to 25bp
    #$file =~ /(.*\.\.\w+).*/
    #$file =~ s/\.[^.]+$//;
    #pass prefix as param?
    my $prefix = "proba";
    print "$prefix";
    open(TRIM, ">$resultDir/$prefix.trim.fastq") || die ("Cannot open $resultDir/$prefix.trim.fastq: $!\n");

    if ($file =~ /.*\.gz/){
	open(IN, "gzcat $readdir/$file |") || die ("Cannot open $readdir/$file: $!\n"); 
    } elsif ($file =~ /.*\.bz2/){
	open(IN, "bzip2 -dc $readdir/$file |") || die ("Cannot open $readdir/$file: $!\n");
    } elsif ($file =~ /.*\.1/ || $file =~ /.*\.2/) {
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
    system("$bowtie $bowtie_param $resultDir/IDX $resultDir/$prefix.trim.fastq >& $resultDir/$prefix.bout");
    #unlink("$resultDir/$prefix.trim.fastq");
}

print stderr "DONE RUNNING BOWTIE, BUILDING UNMAPPED READS\n";
die;
# now the bowtie results are there, simply select from the original files those
# which don't match 
my %files;
my $prefix;
foreach my $file (@rds){
    # first trim to 25bp
    $file =~ /(.*\.)\.(\w+).*/;
    $prefix = $1;
    my $num = $2;
    my $fh;

    if ($file =~ /.*\.gz/){
	open($fh, "gzcat $readdir/$file |") || die ("Cannot open $readdir/$file: $!\n"); 
    } elsif ($file =~ /.*\.bz2/){
	open($fh, "bzip2 -dc $readdir/$file |") || die ("Cannot open $readdir/$file: $!\n");
    } else {
	open($fh, "$readdir/$file") || die ("Cannot open $readdir/$file: $!\n");
    }

    print "Got file: $file, prefix: $prefix, num: $num\n";

    $files{$num} = $fh;
}

my %one;
my %two;
my %singles;
open(ONE, "$resultDir/$prefix.1.bout") || die ("Cannot open $resultDir/$prefix.1.bout: $!\n");
while (<ONE>){
    chomp;
    $_ =~ /^(\S+)\/\d/;
    $one{$1} = 1;
}
print "Got ", scalar keys %one, " elements in one\n";
my @blah = keys %one;
print "First: $blah[0]\n";
close(ONE);

open(TWO, "$resultDir/$prefix.2.bout") || die ("Cannot open $resultDir/$prefix.2.bout: $!\n");
while (<TWO>){
    chomp;
    $_ =~ /^(\S+)\/\d/;
    $two{$1} = 1;
}
print "Got ", scalar keys %two, " elements in two\n";
my @blah = keys %two;
print "First: $blah[0]\n";
close(TWO);
open(SINGLES, "$resultDir/$prefix.singleton.bout") || die ("Cannot open $resultDir/$prefix.singleton.bout: $!\n");
while (<SINGLES>){
    chomp;
    $_ =~ /^(\S+)\/\d/;
    $singles{$1} = 1;
}
print "Got ", scalar keys %singles, " elements in singles\n";
my @blah = keys %singles;
print "First: $blah[0]\n";
close(SINGLES);

open(OUTONE, ">$resultDir/$prefix.unplaced.1.fastq") || die ("cannot open $asmdir/$prefix.unplaced.1.fastq:$!\n");
open(OUTTWO, ">$resultDir/$prefix.unplaced.2.fastq") || die ("cannot open $asmdir/$prefix.unplaced.2.fastq:$!\n");
open(OUTS, ">$resultDir/$prefix.unplaced.singleton.fastq") || die ("cannot open $asmdir/$prefix.unplaced.singleton.fastq:$!\n");

my $f1 = $files{1};
my $f2 = $files{2};
my $fs = $files{singleton};
while (<$f1>){
    if ($_ =~ /^@(\S+)\/\d/){
	my $head1 = $1;
	if (exists $one{$head1} && exists $two{$head1}){
	    # skip both sequences
	    print "got both: $head1\n";
	    <$f1>;<$f1>;<$f1>; # skip 4 lines
	    <$f2>;<$f2>;<$f2>;<$f2>; # skip 4 lines
	} elsif (exists $one{$head1}){
	    print "got one: $head1\n";
	    # skip first sequence
	    <$f1>;<$f1>;<$f1>; # skip 4 lines
	    # print second sequence in singletons
	    $_ = <$f2>;print OUTS;
	    $_ = <$f2>;print OUTS;
	    $_ = <$f2>;print OUTS;
	    $_ = <$f2>;print OUTS;
	} elsif (exists $two{$head1}){
	    print "got two: $head1\n";
	    # print first sequence in singletons
	    print OUTS;
	    $_ = <$f1>;print OUTS;
	    $_ = <$f1>;print OUTS;
	    $_ = <$f1>;print OUTS;
            # skip second
	    <$f2>;<$f2>;<$f2>;<$f2>; # skip 4 lines
	} else {
	    print "got none: $head1\n";
	    # print both in corresponding files
	    print OUTONE;
	    $_ = <$f1>;print OUTONE;
	    $_ = <$f1>;print OUTONE;
	    $_ = <$f1>;print OUTONE;
	    $_ = <$f2>;print OUTTWO;
	    $_ = <$f2>;print OUTTWO;
	    $_ = <$f2>;print OUTTWO;
	    $_ = <$f2>;print OUTTWO;
	}
    } else {
	die ("Weird line: $_");
    }
}
while (<$fs>){
    if ($_ =~ /^@(\S+)\/\d/){
	my $head1 = $1;
	if (exists $singles{$head1}){
	    print "got single: $head1\n";
	    # skip sequence
	    <$fs>;
	    <$fs>;
	    <$fs>;
	} else {
	    print "did not get single: $head1\n";
	    print OUTS;
	    $_ = <$fs>;print OUTS;
	    $_ = <$fs>;print OUTS;
	    $_ = <$fs>;print OUTS;	    
	}
    } else {
	die ("Weird line: $_");
    }
}

close(OUTONE); close(OUTTWO); close(OUTS);
close($files{1}); close($files{2}); close($files{singleton});
exit(0);
