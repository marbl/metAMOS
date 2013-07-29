#! /usr/bin/perl

# Copyright © 2011, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.


use strict;

my $libPath = `ktGetLibPath`;
use lib `ktGetLibPath`;
use KronaTools;

use Getopt::Long;

my $includeConfidence;
my $totalMag;
my $outFile = 'report.krona.html';
my $include;
my $combine;
my $local;

GetOptions(
	'o=s' => \$outFile,
	'i'   => \$include,
	'c'   => \$combine,
	'p'   => \$includeConfidence,
	'l'   => \$local,
	);

setOption('out', $outFile);
setOption('local', $local);
setOption('combine', $combine);
setOption('name', "Root");

if
(
	@ARGV < 1
)
{
 	print '

Description:
   Infers taxonomic abundance from MetaPhyler results.  By default, each
   Metaphyler result file will create a separate dataset named after the file
   (see -c).

Usage:

ktImportMetaPhyler [options] \
   metaphyler_1:taxonomylevel \
   ...

Input:

   metaphyler_output  File containing MetaPhyler results in tabular format. If
                      running BLAST locally, subject IDs in the local database
                      must contain GI numbers in "gi|12345" format.
   
Options:

   [-o <string>]  Output file name.  Default is blast.krona.html.


   [-i]           Include queries with no hits in the total magnitude of the
                  chart.

   [-p]           Read and output confidence values. It is expected the file will be named the same as the input without the word epsilon. That is, if the main results are in proba.epsilon-nb_results.txt the confidence values should be in proba.nb_results.txt.

   [-c]           Combine input files into single dataset.

   [-l]           Create a local chart, which does not require an internet
                  connection to view (but will only work on this computer).

';
	exit;
}

my %tree;

# load taxonomy

print "Loading taxonomy...\n";
loadTaxonomy();


print "Loading name to taxonomy index...\n";
open INFO, "<$libPath/../taxonomy/taxonomy.tab" or die
	"Taxonomy not found.  Was updateTaxonomy.sh run?";
my %ids;
my @classes;
while ( my $line = <INFO> )
{
        chomp $line;
        my ($id, $depth, $parent, $rank, $name) = split /\t/, $line;
	$ids{lc($name)} = $id;
        $classes[$id] = $rank;
}

close INFO;

# parse MetaPhyler results

my $set = 0;
my @datasetNames;
my $zeroEVal;

foreach my $input (@ARGV)
{
	my $totalMagnitude;
	
	my ($fileName, $taxonomicLevel) = split /:/, $input;
	
	$fileName =~ /([^\/]+)\./;
	
	if ( ! $combine )
	{
		push @datasetNames, $1;
	}
 
        if ( ! defined $taxonomicLevel || length($taxonomicLevel) == 0) {
        	$taxonomicLevel = "class";
	}


	print "Importing $fileName...\n";
	 
	open metaPhyl, "<$fileName";

        my $taxName = undef;
	my $magnitude = 0;
        my $reading = 0;

	while ( 1 )
	{
		my $line = <metaPhyl>;
		
		chomp $line;

                if (!defined($line)) {
                   last;
                }

                if ($line =~ m/^>/) {
                   my $taxaLevel = substr $line, 1;
                   if ($taxaLevel eq lc($taxonomicLevel)) {
                      $reading = 1;
                      next;
                   } else {
                      $reading = 0;
                      next;
                   }
                }

                if ($reading == 1) {
		   my
		   (
                        $taxName,
			$magnitude

		    ) = split /\s/, $line;
                    if (!defined($taxName)) {
                       next;
                    }

                    # pick the lowest classified level to use
                    my $bestTaxa = $ids{lc($taxName)};
                    $totalMagnitude += $magnitude;
                    addByTaxID(\%tree, $set, $bestTaxa, 0, $magnitude, 1);
                    last;
                }
	}
        close metaPhyl;
 
        if ($totalMagnitude == 0) {
           addByTaxID(\%tree, $set, 1, undef, 0, 1, 0);
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
	'Confidence',
);

writeTree
(
	\%tree,
	\@attributeNames,
	\@attributeDisplayNames,
	\@datasetNames,
	0,
	120
);
