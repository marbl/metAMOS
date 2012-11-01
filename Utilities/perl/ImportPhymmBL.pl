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
	'l'   => \$local,
	'u=s' => \$url,
        'f=s' => \$fastaClassified
	);

if ( @ARGV < 1 )
{
	setOption('out', 'phymm(bl).krona.html');
	
	printUsage
	(
		'Creates a Krona chart of Phymm or PhymmBL results.

Note: Since confidence scores are not given for species/subspecies
classifications, they inheret confidence scores from genus classifications.',
		'phymmbl_results',
		'PhymmBL results files (results.03.*). Results can also be from Phymm
alone (results.01.*), but ' . getOptionString('phymm') . ' must be specified.',
		1,
		1,
	);
	
	exit 0;
}

setOption('out', $outFile);
setOption('local', $local);
setOption('combine', $combine);
setOption('name', "Root");

if ( defined $url )
{
	setOption('url', $url);
}

my @ranks =
(
        'Phylum',
        'Class',
        'Order',
        'Family',
        'Genus',
        'Species/Subspecies'
);

my %tree;

# load taxonomy

print "Loading taxonomy...\n";
loadTaxonomy();

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

my $tree = newTree();
my $set = 0;
my @datasetNames;
my $useMag;

foreach my $input ( @ARGV )
{
        my ($fileName, $magFile, $taxonomicLevel, $name) = split /:/, $input;

        $fileName =~ /([^\/]+)\./;

	if ( ! getOption('combine') )
	{
		push @datasetNames, $name;
	}
	
	my %magnitudes;
	
	print "Importing $fileName...\n";
        my $annotFile = $fileName;
        $annotFile =~ s/phymm.out/annots/;
        open ANNOTS, ">$annotFile" or die $!;
	
	if ( defined $magFile )
	{
		print "   Loading magnitudes from $magFile...\n";
		loadMagnitudes($magFile, \%magnitudes);
		$useMag = 1;
	}
	
	print "   Reading classifications from $fileName...\n";
	
	open INFILE, "<$fileName" or die $!;
	
	<INFILE>; # eat header
	
	while ( my $line = <INFILE> )
	{
		chomp $line;
		
		my @values = split /\t/, $line;
		my @lineage;
		my $scores;
		
		my $readID = shift @values;
		
		if ( getOption('phymm') )
		{
			my ($species, $score);
			
			($species, $score, @lineage) = @values;
			
			$scores = $score;
			
			@lineage = reverse @lineage;
			
			push @lineage, $species;
		}
		else
		{
			if ( @values < 12 )
			{
				my $phymm = getOptionString('phymm');
				ktDie("Not enough fields in $fileName.  Is it a PhymmBL result file (see $phymm)?");
			}
			
			$scores = ();
			
			$values[1] = $values[3]; # use genus conf for species
			
			for ( my $i = 0; $i < @values; $i += 2 )
			{
				unshift @lineage, $values[$i];
				unshift @$scores, $values[$i + 1];
			}
		}
		

		my $printed = 0;
		for ( my $i = 0; $i < @lineage; $i++ )
		{
			$lineage[$i] = decode($lineage[$i]);

                        if (!$printed && defined($lineage[$i]) && $lineage[$i] ne "" && lc($ranks[$i]) eq $taxonomicLevel) {
                                print ANNOTS "$readID\t$ids{$lineage[$i]}\n";
				$printed = 1;
			}
		}
		
		map { if ( $_ eq '' ) { $_ = 'unknown' } } @lineage;
		addByLineage($tree, $set, \@lineage, $readID, $magnitudes{$readID}, $scores, \@ranks);
                $classified{$readID} = 0;
	}
	
	close INFILE;
        close ANNOTS;

        foreach my $id (keys %classified) {
           if ($classified{$id} == 1) {
              my $magnitude = (defined($magnitudes{$id}) ? $magnitudes{$id} : 1);
              addByTaxID(\%tree, $set, 0, $id, $magnitude, 0);
              $totalMag+=$magnitude;
           }
        }
	
	if ( ! getOption('combine') )
	{
		$set++;
	}
}

# tree output

my @attributeNames =
(
	'magnitude',
	'count',
	'unassigned',
	'rank',
	'score'
);

my @attributeDisplayNames =
(
	$useMag ? 'Magnitude' : undef,
	'Count',
	'Unassigned',
	'Rank',
	getOption('phymm') ? 'Avg. score' : 'Avg. confidence'
);

writeTree
(
	$tree,
	\@attributeNames,
	\@attributeDisplayNames,
	\@datasetNames,
	getOption('hueBad'),
	getOption('hueGood')
);


# subroutines

sub decode
{
	# return special characters in place of their Phymm placeholders
	
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

