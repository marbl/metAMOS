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

my $includeConfidence;
my $totalMag;
my $outFile = 'report.krona.html';
my $include;
my $combine;
my $local;
my $verbose;

GetOptions(
	'o=s' => \$outFile,
	'i'   => \$include,
	'c'   => \$combine,
        'p'   => \$includeConfidence,
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
   Infers taxonomic abundance from FCP results.  By default, each
   FCP result file will create a separate dataset named after the file
   (see -c).

Usage:

ktImportFCP [options] \
   fcp_output_1:magnitute_file_1 \
   ...

Input:

   fcp_output     File containing FCP results in tabular format ("After epsilon filtering").
	            If running BLAST
                    locally, subject IDs in the local database must contain GI
                    numbers in "gi|12345" format.
   
Options:

   [-o <string>]  Output file name.  Default is blast.krona.html.


   [-i]           Include queries with no hits in the total magnitude of the
                  chart.

   [-p]           Read and output confidence values. It is expected the file will be named the same as the input without the word epsilon. That is, if the main results are in proba.epsilon-nb_results.txt the confidence values should be in proba.nb_results.txt.

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

print "Loading name to taxonomy index...\n";
# FCP outputs common names not GI numbers so we have to load the taxonomy and find the gi
open INFO, "<$scriptPath/taxonomy.tab" or die
	"Taxonomy not found.  Was updateTaxonomy.sh run?";
my %ids;
my @classes;
while ( my $line = <INFO> )
{
        chomp $line;
        my ($id, $depth, $parent, $rank, $name) = split /\t/, $line;
	$ids{$name} = $id;
        $classes[$id] = $rank;
}

close INFO;

# parse FCP results

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

        my $CONF;
        my %names;
        if ($includeConfidence) {
           print "Loading confidences...\n";
           my $confFile = $fileName;
           $confFile =~ s/epsilon-//g;
           open CONF, "<$confFile" or die $!;
        
           my $first = 1;
           while (my $line = <CONF>) {
              chomp $line;
              if ($first) {
                 my $count = 0;
                 foreach my $name (split/\t/, $line) {
                    $names{$name} = $count;
                    $count++;
                 }
                 $first = 0;
                 last;
              }
           }
        }

        my $annotFile = $fileName;
        $annotFile =~ s/epsilon-nb_results.txt/annots/;
        open ANNOTS, ">$annotFile" or die $!;	
	print "Importing $fileName...\n";
	 
	open FCP, "<$fileName";

        my $taxID = undef;
        my $taxonomy = undef;
        my $contigID = undef;
	my $magnitude = 0;
        my $taxa = 0;

        print ANNOTS "contigID\tclassID\n";
	while ( 1 )
	{
		my $line = <FCP>;
		
		chomp $line;
		my
		(
                        $contigID,
                        $magnitude,
			$taxonomy

		) = split /\s/, $line, 3;
                if (!defined($taxonomy)) {
                   # done parsing
                   last;
                } else { 
                   if (defined($magnitudes{$contigID})) {
                      $magnitude = $magnitudes{$contigID}
                    }
                    # pick the lowest classified level to use
                    my $bestTaxa = undef;
                    my $bestTaxaName = undef;
                    foreach my $taxa (split/;/, $taxonomy) {
                       if ($taxa eq "unclassified") {
                          last;
                       }
                       $bestTaxa = $ids{$taxa};
                       $bestTaxaName = $taxa;
                       $bestTaxaName =~ s/\s/_/g;

                       if ($classes[$bestTaxa] eq "class") {
                          print ANNOTS "$contigID\t$bestTaxa\n";
                       }
                    }

                    if (defined($bestTaxa)) {
                       my $confidence = 0;
                       if ($includeConfidence && defined($names{$bestTaxaName})) {
                          my $index = $names{$bestTaxaName};
                          while ( 1 ) {
                             my $confLine = <CONF>;
                             chomp $confLine;
                             my ($confID, $remainder) = split /\t/, $confLine, 2;
                             if ($confID eq $contigID) {
                                $confidence = ( split/\t/, $confLine )[$index]; 
                                last;
                             }
                          }
                       }
                       add($set, \%tree, $bestTaxa, $magnitude, $confidence);
                    }
                    $totalMagnitude += $magnitude;
                }
	}
        close ANNOTS;
        close FCP;
        if ($includeConfidence) {
           close CONF;
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

my @attributeDisplayNames =
(
	'Taxon',
	'Rank',
        'Confidence',
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
        0,
        'score',
        0,
	1
);
