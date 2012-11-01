#! /usr/bin/perl

# Copyright © 2011, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.


use strict;

my $libPath = `ktGetLibPath`;
use lib (`ktGetLibPath`);
use KronaTools;

use Getopt::Long;

my $includeConfidence;
my $totalMag;
my $outFile = 'report.krona.html';
my $include;
my $combine;
my $local;
my $url;
my $fastaClassified;

GetOptions(
	'o=s' => \$outFile,
	'i'   => \$include,
	'c'   => \$combine,
	'p'   => \$includeConfidence,
	'l'   => \$local,
	'u=s' => \$url,
        'f=s' => \$fastaClassified
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

';
	exit;
}

setOption('out', $outFile);
setOption('local', $local);
setOption('combine', $combine);
setOption('name', "Root");

if ( defined $url )
{
	setOption('url', $url);
}

my %tree;

# load taxonomy

print "Loading taxonomy...\n";
loadTaxonomy();

# load set of what we could have classified
my %classified;
if (defined($fastaClassified)) {
   my @files = split /:/, $fastaClassified;
   foreach my $file (@files) {
      open F, "<$file" or die $!;
      while (my $line = <F>) {
         chomp $line;
         if ($line =~ /^>(\S*)/) {
           $classified{$1} = 1;
         }
      }
      close(F)
   }
}

print "Loading name to taxonomy index...\n";
# FCP outputs common names not GI numbers so we have to load the taxonomy and find the gi
open INFO, "<$libPath/../taxonomy/taxonomy.tab" or die
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
	
	my ($fileName, $magFile, $taxonomicLevel) = split /:/, $input;
	
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
 
        if ( ! defined $taxonomicLevel || length($taxonomicLevel) == 0) {
        	$taxonomicLevel = "class";
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

        my %outputAnnot = { };
        my $annotFile = $fileName;
        $annotFile =~ s/epsilon-nb_results.txt/annots/;
        open ANNOTS, ">$annotFile" or die $!;	
	print "Importing $fileName...\n";
	 
	open FCP, "<$fileName";

        my $taxID = undef;
        my $taxonomy = undef;
        my $confidences = undef;
        my $contigID = undef;
	my $magnitude = 0;
        my $taxa = 0;

	while ( 1 )
	{
		my $line = <FCP>;
		
		chomp $line;
		my
		(
                        $contigID,
			$taxonomy,
                        $confidences

		) = split /\t/, $line;
                if (!defined($taxonomy)) {
                   # done parsing
                   last;
                } else { 

                   if (defined($magnitudes{$contigID})) {
                      $magnitude = $magnitudes{$contigID};
                    } else {
                      $magnitude = 1;
                    }
                    # pick the lowest classified level to use
                    my $bestIndex = -1;
                    my $bestTaxa = undef;
                    my $bestTaxaName = undef;
                    foreach my $taxa (split/;/, $taxonomy) {
                       if ($taxa eq "unclassified") {
                          last;
                       }
                       $bestTaxa = $ids{$taxa};
                       $bestTaxaName = $taxa;
                       $bestTaxaName =~ s/\s/_/g;
                       $bestIndex++;

                       if (!defined($outputAnnot{$contigID}) && $classes[$bestTaxa] eq $taxonomicLevel) {
                          print ANNOTS "$contigID\t$bestTaxa\n";
                          $outputAnnot{$contigID} = 1;
                       }
                    }

                    if (defined($bestTaxa)) {
                       my $confidence = 0;
                       my $conf2 = 0;
                       if (defined($confidences)) {
                          my @confs = split/;/, $confidences;
                          $confidence = $confs[$bestIndex];
                          $conf2 = log(-1*$confidence);
                       } elsif ($includeConfidence && defined($names{$bestTaxaName})) {
                          my $index = $names{$bestTaxaName};
                          while ( 1 ) {
                             my $confLine = <CONF>;
                             chomp $confLine;

                             my ($confID, $remainder) = split /\t/, $confLine, 2;
                             if ($confID eq $contigID) {
                                $confidence = ( split/\t/, $confLine )[$index]; 
                                $conf2 = log(-1*$confidence);
                                last;
                             }
                             if (eof(CONF)) {
                                die "Error confidence file is not in sync with output, could not find matching confidence value for $contigID\n";
                             }
                          }
                       }
                       addByTaxID(\%tree, $set, $bestTaxa, $contigID, $magnitude, $conf2);
                       $classified{$contigID} = 0;
                    }
                    else
                    {
                       addByTaxID(\%tree, $set, 1, $contigID, $magnitude, 0);
                       $classified{$contigID} = 0;
                    }
                    $totalMagnitude += $magnitude;
                }
	}
        close ANNOTS;
        close FCP;
        if ($includeConfidence) {
           close CONF;
        }
	
        foreach my $id (keys %classified) {
           if ($classified{$id} == 1) {
              my $magnitude = (defined($magnitudes{$id}) ? $magnitudes{$id} : 1);
              addByTaxID(\%tree, $set, 0, $id, $magnitude, 0);
              $totalMagnitude+=$magnitude;
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
	'magnitude',
	'count',
	'unassigned',
	'score',
	'taxon',
	'rank'
);

my @attributeDisplayNames =
(
	'Raw reads',
	'Contigs & singletons',
	'Unassigned',
	'Confidence',
	'Taxon',
	'Rank'
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
