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

my $AMPHORA_TAX_LEVEL = "species";
my $totalMag;
my $outFile = 'report.krona.html';
my $include;
my $random;
my $combine;
my $local;
my $verbose;

GetOptions(
	'o=s' => \$outFile,
	'i'   => \$include,
	'r'   => \$random,
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
	my $magnitude = 0;
	
	while ( 1 )
	{
		my $line = <BLAST>;
		
		chomp $line;
                #print "$line";
		my
		(
			$taxID,
                        $taxLevel,
			$taxaName,
			$magnitude

		) = split /\s+/, $line; #split /\t/, $line;
                if ( ! defined $taxID )
                {
                        last; # EOF
                }
                if ( defined $taxID )
		{
                        if ($taxLevel eq $AMPHORA_TAX_LEVEL) {
			   # add the chosen hit
			   add($set, \%tree, $taxID, $magnitude);
                           $totalMagnitude += $magnitude;
			   $ties = 1;
                        }
		}
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
	'magnitude'
);

my @attributeDisplayNames =
(
	'Taxon',
	'Rank',
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
	\@datasetNames
);
