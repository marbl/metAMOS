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
   Infers taxonomic abundance from results.  By default, each
   Result file will create a separate dataset named after the file
   (see -c).

Usage:

ktImportFCP [options] \
   fcp_output_1:magnitute_file_1 \
   ...

Input:

   output     File containing results in tabular format
   
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

# parse kraken results

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


        my %outputAnnot = { };
        my $annotFile = $fileName;
        $annotFile =~ s/hits/annots/;
        open ANNOTS, ">$annotFile" or die $!;	
	print "Importing $fileName...\n";
	 
	open GEN, "<$fileName";

        my $taxID = undef;
        my $contigID = undef;
	my $magnitude = 0;

	while ( 1 )
	{
		my $line = <GEN>;
		
		chomp $line;
		my
		(
                        $contigID
			$taxID,

		) = split /\t/, $line;
                if (!defined($taxID)) {
                   # done parsing
                   last;
                } else { 
                   if (defined($magnitudes{$contigID})) {
                      $magnitude = $magnitudes{$contigID};
                    } else {
                      $magnitude = 1;
                    }

                    if ($taxID != 0) {
                       my $parent = $taxID;
                       while (1) {
                          if ($parent == 1) {
                             last;
                          }
                          if (!defined($outputAnnot{$contigID}) && getTaxRank($parent) eq $taxonomicLevel) {
                             print ANNOTS "$contigID\t$parent\n";
                          }
                          $parent = getTaxParent($parent);
                       }
                       addByTaxID(\%tree, $set, $taxID, $contigID, $magnitude, 0);
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
        close GEN;
	
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
