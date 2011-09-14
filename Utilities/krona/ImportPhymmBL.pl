#! /usr/bin/perl

# Copyright © 2011, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.


use strict;

# get the path of this script; dependencies should be in the same directory
#
my $scriptPath;
BEGIN
{
	use Cwd 'abs_path';
	abs_path($0) =~ /(.*)\//;
	$scriptPath = $1;
}
use lib "$scriptPath/";

use Getopt::Long;
use Krona;

my $outFile;
my $name = 'all';
my $minConfidence = 0;
my $unclassified;
my $combine;
my $phymm;
my $local;

GetOptions(
	'o=s' => \$outFile,
	'n=s' => \$name,
	'm=f' => \$minConfidence,
	'c'   => \$combine,
	'p'   => \$phymm,
	'l'   => \$local
	);

if ( @ARGV < 1 )
{
	print '
ktImportPhymmBL [options] \
   <results_1>[:magnitude_file_1] \
   [<results_2>[:magnitude_file_2]]
   ...

Input:

   results           PhymmBL results files (results.03.*).  Results can also be
                     from Phymm alone (results.01.*), but -p must be specified.
                     By default, separate datasets will be created for each file
                     (see -c).

   [magnitude_file]  Optional file listing query IDs with magnitudes, separated
                     by tabs.  The can be used to account for read length or
                     contig depth to obtain a more accurate representation of
                     abundance.  By default, query sequences without specified
                     magnitudes will be assigned a magnitude of 1.  Magnitude
                     files for Newbler or Celera Assembler assemblies can be
                     created with getContigMagnitudesNewbler.pl or
                     getContigMagnitudesCA.pl

Options:

   [-o <string>]  Output name.  Default is phymmbl.krona.html or
                  phymm.krona.html.

   [-n <string>]  Name of the highest level.  Default is "all".

   [-m <number>]  Minimum confidence (0-1).  Each query sequence will only be
                  added to taxa that were predicted with a confidence score of
                  at least this value.

   [-c]           Combine input files into single dataset.

   [-p]           Input is phymm only (no confidence scores)

   [-l]           Create a local chart, which does not require an internet
                  connection to view (but will only work on this computer).

';
	
	exit;
}

if ( ! defined $outFile )
{
	if ( $phymm )
	{
		$outFile = 'phymm.krona.html';
	}
	else
	{
		$outFile = 'phymmbl.krona.html';
	}
}

my %all = ();

my %children = ();
my @magnitudes = ();

$all{'children'} = \%children;
$all{'magnitude'} = \@magnitudes;

my $set = 0;
my @datasetNames;

foreach my $input ( @ARGV )
{
	my ($fileName, $magFile) = split /:/, $input;
	
	$fileName =~ /([^\/]+)\./;
	
	if ( ! $combine )
	{
		push @datasetNames, $1;
	}
	
	my %magnitudes;
	
	if ( defined $magFile )
	{
		print "Loading magnitudes for $fileName...\n";
		
		open MAG, "<$magFile" or die $!;
		
		while ( my $line = <MAG> )
		{
			chomp $line;
			my ( $id, $magnitude ) = split /\t/, $line;
			$magnitudes{$id} = $magnitude;
			
			if ( $unclassified )
			{
				$all{'magnitude'}[$set] += $magnitude;
			}
		}
		
		close MAG;
	}
	
	print "Importing $fileName...\n";
	
	open INFILE, "<$fileName" or die $!;
	
	<INFILE>; # eat header
	
	while ( my $line = <INFILE> )
	{
		chomp $line;
		
		my $magnitude = 1;
		
		my
		(
			$readID,
			$bestMatch,
			$score,
			$genus,
			$genusConf,
			$family,
			$familyConf,
			$order,
			$orderConf,
			$class,
			$classConf,
			$phylum,
			$phylumConf
		);
		
		if ( $phymm )
		{
			(
				$readID,
				$bestMatch,
				$score,
				$genus,
				$family,
				$order,
				$class,
				$phylum
			) = split /\t/, $line;
			
			(
				$genusConf,
				$familyConf,
				$orderConf,
				$classConf,
				$phylumConf
			) = (1, 1, 1, 1, 1);
		}
		else
		{
			(
				$readID,
				$bestMatch,
				$score,
				$genus,
				$genusConf,
				$family,
				$familyConf,
				$order,
				$orderConf,
				$class,
				$classConf,
				$phylum,
				$phylumConf
			) = split /\t/, $line;
			
			if ( ! defined $phylumConf )
			{
				print STDERR
					"\nNot enough fields in $fileName.  Is it a PhymmBL result file?\n";
				exit;
			}
		}
		
		# return special characters in place of their Phymm placeholders
		#
		$bestMatch = decode($bestMatch);
		$genus = decode($genus);
		$family = decode($family);
		$order = decode($order);
		$class = decode($class);
		$phylum = decode($phylum);
		
		if ( defined %magnitudes )
		{
			$magnitude = $magnitudes{$readID};
			
			if ( ! defined $magnitude )
			{
				print STDERR "Warning: $readID doesn't exist in magnitude file\n";
			}
		}
		
		if
		(
			$phylumConf >= $minConfidence &&
			! ( defined $magFile && $unclassified )
		)
		{
			$all{'magnitude'}[$set] += $magnitude;
		}
		
		my $phylumHash = addPhymm($set, \%all, $phylum, $magnitude, $phylumConf, 'Phylum');
		my $classHash = addPhymm($set, $phylumHash, $class, $magnitude, $classConf, 'Class');
		my $orderHash = addPhymm($set, $classHash, $order, $magnitude, $orderConf, 'Order');
		my $familyHash = addPhymm($set, $orderHash, $family, $magnitude, $familyConf, 'Family');
		my $genusHash = addPhymm($set, $familyHash, $genus, $magnitude, $genusConf, 'Genus');
		addPhymm($set, $genusHash, $bestMatch, $magnitude, $genusConf, 'Species/Subspecies');#$score); # TODO: translate score to conf
	}
	
	close INFILE;
	
	if ( ! $combine )
	{
		$set++;
	}
}

# tree output

my @attributeNames =
(
	'magnitude',
	'rank'
);

my @attributeDisplayNames =
(
	'Total',
	'Rank'
);

my $hueName;

if ( ! $phymm )
{
	push @attributeNames, 'score';
	push @attributeDisplayNames, 'Avg. confidence';
	$hueName = 'score';
}

print "Writing $outFile...\n";

writeTree
(
	\%all,
	$name,
	$outFile,
	$local,
	'magnitude',
	\@attributeNames,
	\@attributeDisplayNames,
	\@datasetNames,
	0,
	$hueName,
	0,
	120
);

# subroutines

sub addPhymm
{
	my $result;
	my ($set, $node, $child, $magnitude, $confidence, $rank) = @_;
	
	if ( $child eq '' )
	{
		$child = 'unknown taxon';
	}
	
	if ( defined $node && $confidence >= $minConfidence )
	{
		if ( ! defined ${$node}{'children'}{$child} )
		{
			${$node}{'children'}{$child} = ();
			${$node}{'children'}{$child}{'children'} = ();
			
			my @rank = ($rank);
			$node->{'children'}{$child}{'rank'} = \@rank;
			$node->{'children'}{$child}{'magnitude'} = ();
			
			$node->{'children'}{$child}{'score'} = ();
			
			if ( $phymm )
			{
				push @{$node->{'children'}{$child}{'score'}}, 1;
			}
			else
			{
				$node->{'children'}{$child}{'scoreTotal'} = ();
				$node->{'children'}{$child}{'scoreCount'} = ();
			}
		}
		
		$result = ${$node}{'children'}{$child};
		${$result}{'magnitude'}[$set] += $magnitude;
		
		if ( ! $phymm )
		{
			${$result}{'scoreTotal'}[$set] += $confidence;
			${$result}{'scoreCount'}[$set]++;
		}
	}
	
	return $result;
}

sub decode
{
	my ($string) = @_;
	
	$string =~ s/_/ /g;
	$string =~ s/UNDERSCORE/_/g;
	$string =~ s/SLASH/\//g;
	$string =~ s/LPAREN/\(/g;
	$string =~ s/RPAREN/\)/g;
	$string =~ s/SINGLEQUOTE/'/g;
	$string =~ s/DOUBLEQUOTE/"/g;
	$string =~ s/COLONCOLON/:/g;
	$string =~ s/SEMICOLON/;/g;
	
	return $string;
}

