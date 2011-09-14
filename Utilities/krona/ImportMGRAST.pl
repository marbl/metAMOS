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
use lib "$scriptPath/../lib";

use Getopt::Long;
use Krona;

my $outFile = 'MG-RAST.krona.html';
my $name = 'all';
my $combine;
my $colorIdentity;
my $local;

GetOptions(
	'o=s' => \$outFile,
	'n=s' => \$name,
	'c'   => \$combine,
	'p'   => \$colorIdentity,
	'l'   => \$local
	);

if
(
	!defined $outFile ||
	@ARGV < 1
)
{
	print '

importMG-RAST.pl [options] table1.tsv [table2.tsv] ...

Creates a Krona chart from tables exported from MG-RAST sequence profiles.
The profiles can be metabolic or phylogenetic, but must be consistent between
files.  By default, separate datasets will be created for each file (see -c).

Options:

   [-o <string>]  Output file.  Default is MG-RAST.krona.html

   [-n <string>]  Name of the highest level.  Default is "all".

   [-c]           Combine input files into single dataset.

   [-p]           Color taxa by average percent identity instead of average log
                  e-value.

   [-l]           Create a local chart, which does not require an internet
                  connection to view (but will only work on this computer).

';
	
	exit;
}

my %all;
my @ranks;
my @datasetNames;
my $set = 0;

foreach my $fileName ( @ARGV )
{
	$fileName =~ /([^\/]+)\./;
	
	if ( ! $combine )
	{
		push @datasetNames, $1;
	}
	print "Importing $fileName...\n";
	
	open INFILE, "<$fileName" or die $!;
	
	my $header = <INFILE>;
	my $offset;
	my @fields = split /\t/, $header;
	
	# remove metagenome and source (if present)
	#
	shift @fields;
	#
	if ( $fields[0] eq 'source' )
	{
		shift @fields;
		$offset = 2;
	}
	else
	{
		$offset = 1;
	}
	
	my $i = 0;
	
	while ( $fields[0] ne 'abundance' )
	{
		$ranks[$i] = shift @fields;
		$i++;
	}
	
	while ( my $line = <INFILE> )
	{
		chomp $line;
		
		my @data = split /\t/, $line;
		my $magnitude = $data[$offset + @ranks];
		my $score;
		
		if ( $colorIdentity )
		{
			$score = $data[$offset + 2 + @ranks];
		}
		else
		{
			$score = $data[$offset + 1 + @ranks];
		}
		
		add($set, \%all, $magnitude, $score, 0, @data[$offset..($offset + @ranks - 1)]);
	}
	
	close INFILE;
	
	if ( ! $combine )
	{
		$set++;
	}
}

my @attributeNames =
(
	'magnitude',
	'score',
	'rank'
);

my @attributeDisplayNames =
(
	'Abundance',
	$colorIdentity ? 'Avg. % identity' : 'Avg. log e-value',
	'Rank'
);

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
	'score',
	$colorIdentity ? 0 : 120,
	$colorIdentity ? 120 : 0
);


sub add
{
	my $child;
	my ($dataset, $node, $magnitude, $score, $depth, @calls) = @_;
	
	if ( ! defined ${$node}{'children'} )
	{
#		${$node}{'children'} = ();
#		${$node}{'magnitude'} = ();
#		${$node}{'scoreTotal'} = ();
#		${$node}{'scoreCount'} = ();
#		print "\t${$node}{'magnitude'}\n";
		
		if ( $depth > 0 )
		{
#			${$node}{'rank'} = ();
		}
	}
	
	#print "${$node}{'magnitude'}\t$magnitude\t@calls\n";
	
	${$node}{'magnitude'}[$dataset] += $magnitude;
	${$node}{'scoreTotal'}[$dataset] += $score * $magnitude;
	${$node}{'scoreCount'}[$dataset] += $magnitude;
	
	if ( $depth > 0 )
	{
		${$node}{'rank'}[0] = $ranks[$depth - 1];
	}
	
	if ( @calls > 0 )
	{
		my $name;
		
		do
		{
			$name = shift @calls;
			$depth++;
		}
		while ( $name eq '' || $name eq '-' );
		
		if ( ! defined ${$node}{'children'}{$name} )
		{
			my %newHash = ();
			${$node}{'children'}{$name} = \%newHash;
		}
		
		$child = ${$node}{'children'}{$name};
		
		add($dataset, $child, $magnitude, $score, $depth, @calls);
	}
}
