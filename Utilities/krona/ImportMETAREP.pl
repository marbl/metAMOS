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
use lib "$scriptPath/../lib";

use Getopt::Long;
use Krona;

my $outFile = 'metarep.krona.html';
my $random;
my $colorIdentity;
my $colorBitScore;
my $combine;
my $local;
my $verbose;

GetOptions(
	'o=s' => \$outFile,
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
   Infers taxonomic abundance from the best BLAST hits listed in the blast.tab
   files of METAREP data folders.  By default, separate datasets for each folder
   will be created and named after the folder (see -c).  The blast.tab files
   within the specified folders must be unzipped.

Usage:

ktImportMETAREP [options] <folder_1> <folder_2> ...

Input:

<folders>          METAREP data folders containing an unzipped blast.tab file.

Options:

   [-o <string>]  Output file name.  Default is metarep.krona.html.

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

my ($input) = @ARGV;

my $tree = newTree();

# taxonomy must be loaded for LCA

print "Loading taxonomy...\n";

loadTaxonomy();

my $lastReadID;
my $set = 0;
my @datasetNames;
my $zeroEVal;

foreach my $folder (@ARGV)
{
	$folder =~ /([^\/]+)/;
	
	if ( ! $combine )
	{
		push @datasetNames, $1;
	}
	
	print "Importing $folder...\n";
	
	open IN, "<$folder/blast.tab" or die
		"Couldn't open $folder/blast.tab.  Was it unzipped?";
	
	while ( my $line = <IN> )
	{
		chomp $line;
		my @values = split /\t/, $line;
		
		my $readID = $values[0];
		
		if ( $readID ne $lastReadID )
		{
			my $ties = 0;
			my $taxID;
			my $readLength = $values[2];
			
			while ( $values[15] =~ /taxon:(\d+)/g )
			{
				my $newTaxID = $1;
				
				if ( ! taxIDExists($newTaxID) )
				{
					print STDERR "$readID: Taxon $newTaxID does not exist; using root.\n";
					$newTaxID = 1;
				}
				
				if ( $random && int(rand(++$ties)) == 0 )
				{
					$taxID = $newTaxID;
				}
				elsif ( ! $random )
				{
					if ( $taxID )
					{
						$taxID = lowestCommonAncestor($taxID, $newTaxID);
					}
					else
					{
						$taxID = $newTaxID;
					}
				}
			}
			
			if ( $verbose )
			{
				print "$readID\ttaxID: $taxID\n";
			}
			
			my $score;
			
			if ( $colorIdentity )
			{
				$score = $values[11];
			}
			elsif ( $colorBitScore )
			{
				$score = $values[12];
			}
			else
			{
				if ( $values[19] > 0 )
				{
					$score = $values[19];
				}
				else
				{
					$score = "1e-500";
					$zeroEVal = 1;
				}
			}
			
			add($set, $tree, $taxID, $readLength, $score);
		}
		
		$lastReadID = $readID;
	}
	
	if ( ! $combine )
	{
		$set++;
	}
	
	close IN;
}

if ( $zeroEVal )
{
	print "\nWARNING: Couldn't take log for e-values of 0.  Used 1e-413.\n\n";
}

my $scoreName;

if ( $colorIdentity )
{
	$scoreName = 'Avg. % Identity';
}
elsif ( $colorBitScore )
{
	$scoreName = 'Avg. bit score';
}
else
{
	$scoreName = 'Log avg. e-value';
}

my @attributeNames =
(
	'magnitude',
	'taxon',
	'rank',
	'score',
);

my @attributeDisplayNames =
(
	'Total',
	'Taxon',
	'Rank',
	$scoreName
);

print "Writing $outFile...\n";

writeTree
(
	$tree,
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

