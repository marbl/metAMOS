#! /usr/bin/perl

# Copyright © 2011, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.


use strict;

# get the path of this script; dependencies are relative
#
my $scriptPath;
BEGIN
{
	use Cwd 'abs_path';
	abs_path($0) =~ /(.*)\//;
	$scriptPath = $1;
}
#use lib "$scriptPath/../lib";
use lib "$scriptPath/";

use Getopt::Long;
use Krona;

my $totalMag;
my $outFile = 'report.krona.html';
my $include;
my $random;
my $colorIdentity;
my $colorBitScore;
my $combine;
my $local;
my $verbose;

GetOptions(
	'o=s' => \$outFile,
	'i'   => \$include,
	'r'   => \$random,
	'p'   => \$colorIdentity,
	'b'   => \$colorBitScore,
	'c'   => \$combine,
	'l'   => \$local,
	'v'   => \$verbose
	);

if
(
	@ARGV < 1
)
{
 	print '

Description:
   Infers taxonomic abundance from Amphora results.  By default, each
   Amphora result file will create a separate dataset named after the file
   (see -c).

Usage:

ktImportBLAST [options] \
   phmmer_output_1[:magnitude_file_1] \
   phmmer_output_2[:magnitude_file_2] \
   ...

Input:

   blast_output     File containing BLAST results in tabular format ("Hit table
	                (text)" when downloading from NCBI).  If running BLAST
                    locally, subject IDs in the local database must contain GI
                    numbers in "gi|12345" format.
   
   [magnitue_file]  Optional file listing query IDs with magnitudes, separated
                    by tabs.  The can be used to account for read length or
                    contig depth to obtain a more accurate representation of
                    abundance.  By default, query sequences without specified
                    magnitudes will be assigned a magnitude of 1.  Magnitude
                    files for Newbler or Celera Assembler assemblies can be
                    created with getContigMagnitudesNewbler.pl or
                    getContigMagnitudesCA.pl
Options:

   [-o <string>]  Output file name.  Default is blast.krona.html.

   [-i]           Include queries with no hits in the total magnitude of the
                  chart.  If running BLAST locally, the output must have comment
                  lines ("-m 9") to use this option.

   [-r]           Break ties for the top hit randomly.  Default is to use the
                  lowest common ancestor of all ties for the top hit.

   [-p]           Color taxa by average percent identity instead of average log
                  e-value.

   [-b]           Color taxa by average bit score instead of average log
                  e-value.

   [-c]           Combine input files into single dataset.

   [-l]           Create a local chart, which does not require an internet
                  connection to view (but will only work on this computer).

   [-v]           Verbose.

';
	exit;
}

my %tree;

# load taxonomy

print "Loading taxonomy...\n";
loadTaxonomy();

# parse BLAST results

my $set = 0;
my @datasetNames;
my $zeroEVal;

foreach my $input (@ARGV)
{
	my $totalMagnitude;
	
	my ($fileName, $magFile) = split /:/, $input;
	
	$fileName =~ /([^\/]+)\./;
	
	if ( ! $combine )
	{
		push @datasetNames, $1;
	}
	
	my %magnitudes;
	
	# load magnitudes
	
	if ( defined $magFile )
	{
		print "Loading magnitudes for $fileName...\n";
		
		open MAG, "<$magFile" or die $!;
		
		while ( my $line = <MAG> )
		{
			chomp $line;
			my ( $id, $mag ) = split /\t/, $line;
			$magnitudes{$id} = $mag;
		}
		
		close MAG;
	}
	
	print "Importing $fileName...\n";
	
	open BLAST, "<$fileName";
	
	my $topScore;
	my $ties;
	my $taxID;
	my $score = 0;
	
	while ( 1 )
	{
		my $line = <BLAST>;
		
		chomp $line;
                #print "$line";
		if ( $line =~ /^#/ )
		{
		my
		(
			$taxID,
                        $blank1,
			$taxaName,
			$n1,
			$n2,
			$n3,
			$n4,
			$n5,
			$n6,
			$n7

		) = split /\s+/, $line; #split /\t/, $line;
                if ( defined $taxID )
		{
			# add the chosen hit
			my $magnitude = 1;
			add($set, \%tree, $taxID, $magnitude, $score);
			$ties = 1;
		}
		
		if ( ! defined $taxID )
		{
			last; # EOF
		}
		
		if ( $colorIdentity )
		{
			$score = 0;
		}
		elsif ( $colorBitScore )
		{
			$score = 0;
		}
		else
		{
			if ( $eVal > 0 )
			{
				$score = 0;
			}
			else
			{
				$score = 0;#"1e-500";
				$zeroEVal = 1;
			}
		}
        	if ($topScore < $score) { $topScore = $score; }
	}
	
	if ( $include && $totalMagnitude )
	{
		$tree{'magnitude'}[$set] = $totalMagnitude;
	}
	
	if ( ! $combine )
	{
		$set++;
	}
}

if ( $zeroEVal )
{
	print "\nWARNING: Couldn't take log for e-values of 0.  Used 1e-413.\n\n";
}

my @attributeNames =
(
	'taxon',
	'rank',
	'score',
	'magnitude'
);

my $scoreName;

if ( $colorIdentity )
{
	$scoreName = 'Avg. % identity';
}
elsif ( $colorBitScore )
{
	$scoreName = 'Avg. bit score';
}
else
{
	$scoreName = 'Log avg. e-value';
}

my @attributeDisplayNames =
(
	'Taxon',
	'Rank',
	$scoreName,
	'Total'
);

print "Writing $outFile...\n";

writeTree
(
	\%tree,
	'Root',
	$outFile,
	$local,
	'magnitude',
	\@attributeNames,
	\@attributeDisplayNames,
	\@datasetNames,
	! ( $colorIdentity || $colorBitScore ),
	'score',
	$colorBitScore || $colorIdentity ? 0 : 120,
	$colorBitScore || $colorIdentity ? 120 : 0
);
