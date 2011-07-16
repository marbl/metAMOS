#!/usr/bin/perl

use strict;

use Cwd;

$| = 1;

my $dataFile = shift;

$dataFile =~ s/\(/\\\(/g;

$dataFile =~ s/\)/\\\)/g;

my $outputPrefix = $dataFile;

$outputPrefix =~ s/\//\_/g;

$outputPrefix =~ s/\./\_/g;

die("Usage: $0 <file containing query reads>\n") if not $dataFile;

# Date variable for progress reports.

my $date;

my $blastFile;

my $blastScore;

my $blastMatch;

###########################################################################################
# 
# Save the query lengths for later confidence-score computation.
# 
###########################################################################################

print "Scanning query lengths...";

open IN, "<$dataFile" or die("Can't open $dataFile for reading.\n");

my $currentID = '';

my $currentLength = 0;

my $first = 1;

my $queryLengths = {};

my @readIDs = ();

while ( my $line = <IN> ) {
   
   if ( $line =~ /^>(\S+)/ and $first ) {

      $currentID = $1;

      push @readIDs, $currentID;

      $first = 0;

   } elsif ( $line =~ /^>(\S+)/ ) {

      my $tempID = $1;

      $queryLengths->{$currentID} = $currentLength;

      $currentLength = 0;

      $currentID = $tempID;
      
      push @readIDs, $currentID;

   } elsif ( $line =~ /^(\S+)$/ ) {

      my $data = $1;

      $currentLength += length($data);
   }
}

$queryLengths->{$currentID} = $currentLength;

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# Set / check BLAST environment variables.
# 
###########################################################################################

my $oldBlastDB = $ENV{'BLASTDB'};

my $oldBlastMat = $ENV{'BLASTMAT'};

die("BLAST error: environment variable 'BLASTMAT' not found; please check your local blast binary installation.\n") if not $oldBlastMat;

###########################################################################################
# 
# Grab the list of ICMs.
# 
###########################################################################################

print "Reading IMM list...";

# ICM list: 1: .genomeData.

opendir DOT, '.genomeData' or die("Can't open . for scanning.\n");

my @dirs = readdir DOT;

closedir DOT;

my @ICMs;

foreach my $dir ( @dirs ) {
   
   if ( -d ".genomeData/$dir" and $dir !~ /^\./ ) {
      
      &scanDir(".genomeData/$dir", \@ICMs);
   }
}

# ICM list: 2: .genomeData/.userAdded, if it exists.

my $userDir = '.genomeData/.userAdded';

if ( -e $userDir ) {
   
   opendir DOT, $userDir or die("Can't open $userDir for scanning.\n");

   my @dirs = readdir DOT;

   closedir DOT;

   foreach my $dir ( @dirs ) {
      
      if ( -d "$userDir/$dir" and $dir !~ /^\./ ) {
	 
	 &scanDir("$userDir/$dir", \@ICMs);
      }
   }
}

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# Create a reverse-complement copy of the query file.
# 
###########################################################################################

print "Creating reverse complement of $dataFile...";

system(".scripts/revCompFASTA.pl $dataFile");

my $revDataFile = $dataFile;

$revDataFile =~ /(\.[^\.]+)$/;

my $extension = $1;

$revDataFile =~ s/$extension$/\.revComp$extension/;

my $fileCheck = $revDataFile;

$fileCheck =~ s/\\//g;

if ( not -e $fileCheck ) {
   
   die("File renaming problem [revCompFASTA.pl]: could not detect revComp file \"$fileCheck\".\n");
}

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# Read organism names so we can generate results files comparable to what BLAST will output.
# 
###########################################################################################

print "Loading taxonomy...";

# Accession scan 1: RefSeq orgs.

my $accFile = '.taxonomyData/.0_accessionMap/accessionMap.txt';

open IN, "<$accFile" or die("Can't open $accFile for reading.\n");

my $speciesDirName = {};

while ( my $line = <IN> ) {
   
   chomp $line;

   (my $orgName, my $prefix, my $seqType, my $desc) = split(/\t/, $line);
   
   $speciesDirName->{$prefix} = $orgName;
}

close IN;

# Accession scan 2: user-added orgs, if they exist.

$accFile = '.taxonomyData/.0_accessionMap/accessionMap_userAdded.txt';

if ( -e $accFile ) {
   
   open IN, "<$accFile" or die("Can't open $accFile for reading.\n");

   while ( my $line = <IN> ) {
      
      chomp $line;

      (my $orgName, my $prefix) = split(/\t/, $line);

      $speciesDirName->{$prefix} = $orgName;
   }

   close IN;
}

###########################################################################################
# 
# Load full taxonomic metadata for all organisms in the database.
# 
###########################################################################################

# Taxonomy scan 1: RefSeq organisms.

my $taxFile = '.taxonomyData/.3_parsedTaxData/distributionOfTaxa.txt';

open IN, "<$taxFile" or die("Can't open $taxFile for reading.\n");

my $tax = {};

while ( my $line = <IN> ) {
   
   if ( $line =~ /^\S/ ) {
      
      chomp $line;

      (my $taxType, my $taxVal, my $prefixAndSpecies, my $dirName) = split(/\t/, $line);

      if ( $taxType eq 'phylum' or $taxType eq 'class' or $taxType eq 'order' or $taxType eq 'family' or $taxType eq 'genus' ) {
	 
	 $tax->{$dirName}->{$taxType} = $taxVal;
      }
   }
}

close IN;

# Taxonomy scan 2: user-added organisms, if they exist.

$taxFile = '.taxonomyData/.3_parsedTaxData/distributionOfTaxa_userAdded.txt';

if ( -e $taxFile ) {
   
   open IN, "<$taxFile" or die("Can't open $taxFile for reading.\n");

   while ( my $line = <IN> ) {
      
      if ( $line =~ /^\S/ ) {
	 
	 chomp $line;

	 (my $taxType, my $taxVal, my $prefixAndSpecies, my $dirName) = split(/\t/, $line);

	 if ( $taxType eq 'phylum' or $taxType eq 'class' or $taxType eq 'order' or $taxType eq 'family' or $taxType eq 'genus' ) {
	    
	    $tax->{$dirName}->{$taxType} = $taxVal;
	 }
      }
   }
}

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# Score the query data (using both forward and reverse directions) with Phymm.
# 
###########################################################################################

print "Scoring reads with Phymm...";

my $score = {};

my $topScoringICM = {};

my $rawPhymmFile = 'rawPhymmOutput_' . $outputPrefix . '.txt';

$rawPhymmFile =~ s/\\//g;

open OUT, ">$rawPhymmFile" or die("Can't open $rawPhymmFile for writing.\n");

print OUT "BEGIN_ICM_LIST\n";

foreach my $icm ( sort { $a cmp $b } @ICMs ) {
   
   print OUT "$icm\n";
}

print OUT "END_ICM_LIST\nBEGIN_READID_LIST\n";

foreach my $readID ( @readIDs ) {
   
   print OUT "$readID\n";
}

print OUT "END_READID_LIST\nBEGIN_DATA_MATRIX\n";

my $icmCount = @ICMs;

my $icmsFinished = 0;

my $logFile = $outputPrefix . '_progress.txt';

$logFile =~ s/\\//g;

open LOG, ">$logFile" or die("Can't open $logFile for writing.\n");

print LOG "Scoring reads with Phymm [$date]...\n\n";

foreach my $ICM ( sort { $a cmp $b } @ICMs ) {
   
   my $icmPrefix = $ICM;

   $icmPrefix =~ s/.+\/([^\/]+)$/$1/;

   $icmPrefix =~ s/\.icm$//;
   
   my $fullScore = {};
   
   my $command = '(.scripts/.icmCode/bin/simple-score -N ' . $ICM . ' < ' . $dataFile . ' > tempFwd_' . $outputPrefix . '.txt) >& errFile_' . $outputPrefix . '.txt';
   
   system($command);
   
   my $inFile = "tempFwd_${outputPrefix}.txt";

   $inFile =~ s/\\//g;

   open IN, "<$inFile" or die("Can't open $inFile for reading.\n");
   
   while ( my $line = <IN> ) {
      
      if ( $line =~ /(\S+)\s+(\S+)/ ) {
	 
	 my $queryID = $1;

	 my $queryScore = $2;
	 
	 $fullScore->{$queryID}->{$icmPrefix} = $queryScore;
	 
	 if ( not $score->{$queryID} or $queryScore > $score->{$queryID} ) {
	    
	    $score->{$queryID} = $queryScore;

	    $topScoringICM->{$queryID} = $icmPrefix;
	 }
      }
   }

   close IN;

   my $command = '(.scripts/.icmCode/bin/simple-score -N ' . $ICM . ' < ' . $revDataFile . ' > tempRev_' . $outputPrefix . '.txt) >& errFile_' . $outputPrefix . '.txt';
   
   system($command);
   
   $inFile = "tempRev_${outputPrefix}.txt";
   
   $inFile =~ s/\\//g;

   open IN, "<$inFile" or die("Can't open $inFile for reading.\n");
   
   while ( my $line = <IN> ) {
      
      if ( $line =~ /(\S+)\s+(\S+)/ ) {
	 
	 my $queryID = $1;

	 my $queryScore = $2;
	 
	 if ( $queryScore > $fullScore->{$queryID}->{$icmPrefix} ) {
	    
	    $fullScore->{$queryID}->{$icmPrefix} = $queryScore;
	 }
	 
	 if ( not $score->{$queryID} or $queryScore > $score->{$queryID} ) {
	    
	    $score->{$queryID} = $queryScore;

	    $topScoringICM->{$queryID} = $icmPrefix;
	 }
      }
   }

   close IN;
   
   my $firstDataPoint = 1;

   foreach my $queryID ( @readIDs ) {
      
      if ( $firstDataPoint ) {
	 
	 print OUT "$fullScore->{$queryID}->{$icmPrefix}";

	 $firstDataPoint = 0;

      } else {
	 
	 print OUT "\t$fullScore->{$queryID}->{$icmPrefix}";
      }
   }
   
   print OUT "\n";

   $icmsFinished++;

   if ( $icmsFinished % 50 == 0 ) {
      
      $date = `date`;

      chomp $date;

      print LOG "   ...finished scoring read set with $icmsFinished / $icmCount IMMs [$date]...\n";
   }
}

$date = `date`;

chomp $date;

print LOG "\n...done.  [$date]\n";

close LOG;

print OUT "END_DATA_MATRIX\n";

close OUT;

###########################################################################################
# 
# Write the Phymm output.
# 
###########################################################################################

my $phymmOut = 'results.01.phymm_' . $outputPrefix . '.txt';

$phymmOut =~ s/\\//g;

open OUT, ">$phymmOut" or die("Can't open $phymmOut for writing.\n");

print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";

foreach my $queryID ( sort { $a cmp $b } keys %$score ) {
   
   my $spName = $speciesDirName->{$topScoringICM->{$queryID}};
   
   print OUT "$queryID\t$spName\t$score->{$queryID}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
}

close OUT;

system("rm tempFwd_${outputPrefix}.txt");

system("rm tempRev_${outputPrefix}.txt");

system("rm errFile_${outputPrefix}.txt");

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# Score the query data using BLAST.
# 
###########################################################################################

print "Scoring reads with BLAST...";

my $newBlastDB = cwd;

$newBlastDB .= '/.blastData/';

$ENV{'BLASTDB'} = $newBlastDB;

my $blastCmd = "blastall -i $dataFile -o rawBlastOutput_${outputPrefix}.txt -d phymm_BLAST_DB -p blastn -a 2 -m 9";

system($blastCmd);

if ( -e "error.log" ) {
   
   system("rm error.log");
}

###########################################################################################
# 
# Parse the raw BLAST results file to generate a best-hit list.
# 
###########################################################################################

$blastFile = 'rawBlastOutput_' . $outputPrefix . '.txt';

$blastFile =~ s/\\//g;

$blastScore = {};

$blastMatch = {};

open IN, "<$blastFile" or die("Can't open $blastFile for reading.\n");

while ( my $line = <IN> ) {
   
   if ( $line !~ /^#/ ) {
      
      chomp $line;

      my @fields = split(/\t/, $line);

      my $queryID = $fields[0];

      my $matchName = $fields[1];

      $matchName =~ s/\.\d+$//;

      my $currentScore = $fields[10];

      if ( not $blastScore->{$queryID} or $blastScore->{$queryID} > $currentScore ) {
	 
	 $blastScore->{$queryID} = $currentScore;

	 $blastMatch->{$queryID} = $matchName;
      }
   }
}

close IN;

my $blastOut = 'results.02.blast_' . $outputPrefix . '.txt';

$blastOut =~ s/\\//g;

open OUT, ">$blastOut" or die("Can't open $blastOut for writing.\n");

print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";

foreach my $queryID ( sort { $a cmp $b } keys %$blastScore ) {
   
   my $spName = $blastMatch->{$queryID};
   
   print OUT "$queryID\t$spName\t$blastScore->{$queryID}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
}

close OUT;

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# Combine Phymm and BLAST scores to generate PhymmBL predictions.
# 
###########################################################################################

print "Combining scores...\n\n";

# Read the raw BLAST E-values.

$blastFile = 'rawBlastOutput_' . $outputPrefix . '.txt';

$blastFile =~ s/\\//g;

$blastScore = {};

open IN, "<$blastFile" or die("Can't open $blastFile for reading.\n");

$date = `date`;

chomp $date;

print "   Scanning BLAST output [$date]...";

while ( my $line = <IN> ) {
   
   if ( $line !~ /^#/ ) {
      
      chomp $line;

      my @fields = split(/\t/, $line);

      my $queryID = $fields[0];

      my $matchName = $fields[1];

      $matchName =~ s/\.\d+$//;

      my $currentScore = $fields[10];

      if ( not $blastScore->{$queryID}->{$matchName} or $blastScore->{$queryID}->{$matchName} > $currentScore ) {
	 
	 $blastScore->{$queryID}->{$matchName} = $currentScore;
      }
   }
}

close IN;

$date = `date`;

chomp $date;

print "done.  [$date]\n";

# Read the raw Phymm scores.

my $phymmFile = 'rawPhymmOutput_' . $outputPrefix . '.txt';

$phymmFile =~ s/\\//g;

open IN, "<$phymmFile" or die("Can't open $phymmFile for reading.\n");

$date = `date`;

chomp $date;

print "   Scanning Phymm output [$date]...";

# 1: create an ordered array of ICMs (and metadata for those ICMs)
# using the map which has been saved in the first section of the
# rawPhymmOutput file.

my @indexedICMs = ();

while ( my $line = <IN> ) {
   
   if ( $line =~ /END\_ICM/ ) {
      
      last;

   } elsif ( $line !~ /^BEGIN\_ICM/ ) {
      
      if ( $line =~ /^(\S+)$/ ) {
	 
	 my $lineText = $1;

	 $lineText =~ /([^\/]+)\/[^\/]+\.icm$/;

	 my $tempOrgName = $1;

	 push @indexedICMs, $tempOrgName;

      } else {
	 
	 print "\n   WARNING: Bad filename for ICM: $line\n";
      }
   }
}


# 2: save the ordered query IDs from the next section of the rawPhymmOutput file.

my @indexedReads = ();

while ( my $line = <IN> ) {
   
   if ( $line =~ /END\_READID/ ) {
      
      last;

   } elsif ( $line !~ /BEGIN\_READID/ ) {
      
      if ( $line =~ /^(\S+)$/ ) {
	 
	 my $readID = $1;

	 push @indexedReads, $readID;

      } else {
	 
	 print "\n   WARNING: Bad read ID: $line\n";
      }
   }
}

# 3: scan the stored ICM scores and save the best score assigned to each read
# by (any one of the [possibly multiple] ICMs assigned to) each organism.

# Format: $phymmScore->{$query_ID}->{$org_dirName} = $bestNumericScoreForThatOrg;

my $phymmScore = {};

my $currentIcmIndex = 0;

while ( my $line = <IN> ) {
   
   if ( $line =~ /END\_DATA\_MATRIX/ ) {
      
      last;

   } elsif ( $line !~ /BEGIN\_DATA\_MATRIX/ ) {
      
      chomp $line;

      my @currentICMScores = split(/\t/, $line);

      my $currentOrg = $indexedICMs[$currentIcmIndex];

      for ( my $currentReadIndex = 0; $currentReadIndex <= $#indexedReads; $currentReadIndex++ ) {
	 
	 my $currentScore = $currentICMScores[$currentReadIndex];

	 my $currentReadID = $indexedReads[$currentReadIndex];
	 
	 # If we haven't yet scored this read with an ICM from this organism, save the current score as the best seen so far.

	 if ( not $phymmScore->{$currentReadID}->{$currentOrg} ) {
	    
	    $phymmScore->{$currentReadID}->{$currentOrg} = $currentScore;

	 } elsif ( $phymmScore->{$currentReadID}->{$currentOrg} < $currentScore ) {
	    
	    # If we have scored this read with some other ICM from this organism, but this current score
	    # is better, than replace the existing score with the current one as the best seen so far.

	    $phymmScore->{$currentReadID}->{$currentOrg} = $currentScore;
	 }
      }

      $currentIcmIndex++;
   }
}

close IN;

$date = `date`;

chomp $date;

print "done.  [$date]\n";

# Iterate on each query ID.  Phymm will score everything;
# BLAST will (very occasionally) turn up no matches, so we
# use the $phymmScore keys as the reference index list.
# 
# Combine scores: the formula is [IMM score] + [1.2 * (4 - log(BLAST E-value))].

my $combinedScore = {};

my $combinedMatch = {};

$date = `date`;

chomp $date;

print "   Mixing scores [$date]...";

foreach my $queryID ( keys %$phymmScore ) {
   
   # Check each match.

   foreach my $matchName ( keys %{$phymmScore->{$queryID}} ) {
      
      # Make sure there's a corresponding BLAST score.

      if ( not $blastScore->{$queryID} or not $blastScore->{$queryID}->{$matchName} ) {
	 
	 # If there isn't, just find the best Phymm score.

	 if ( not $combinedScore->{$queryID} ) {
	    
	    $combinedScore->{$queryID} = $phymmScore->{$queryID}->{$matchName};

	    $combinedMatch->{$queryID} = $matchName;

	 } elsif ( $combinedScore->{$queryID} < $phymmScore->{$queryID}->{$matchName} ) {
	    
	    $combinedScore->{$queryID} = $phymmScore->{$queryID}->{$matchName};

	    $combinedMatch->{$queryID} = $matchName;
	 }

      } else {
	 
	 # There is a corresponding BLAST score: compute the combined score and compare.
	 
	 my $tempScore = $phymmScore->{$queryID}->{$matchName};
	 
	 if ( $blastScore->{$queryID}->{$matchName} == 0 ) {
	    
	    $tempScore += 1.2 * (4 - log(1e-180));
	    
	 } else {
	    
	    $tempScore += 1.2 * (4 - log($blastScore->{$queryID}->{$matchName}));
	 }
	 
	 # Update the best-hit record for this query if we've got
	 # something better than something we've seen so far.
	 
	 if ( not $combinedScore->{$queryID} or $tempScore > $combinedScore->{$queryID} ) {
	    
	    $combinedScore->{$queryID} = $tempScore;

	    $combinedMatch->{$queryID} = $matchName;
	 }

      } # end if ( there's a BLAST score for this query and this matching reference organism )

   } # end foreach ( reference organism to which this query was compared )

} # end foreach ( query key in the $phymmScore hash )

$date = `date`;

chomp $date;

print "done.  [$date]\n";

###########################################################################################
# 
# Output the combined (i.e., PhymmBL) results.
# 
###########################################################################################

my $combinedOut = 'results.03.phymmBL_' . $outputPrefix . '.txt';

$combinedOut =~ s/\\//g;

$date = `date`;

chomp $date;

print "   Writing PhymmBL results to $combinedOut [$date]...\n\n";

open OUT, ">$combinedOut" or die("Can't open $combinedOut for writing.\n");

print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tGENUS_CONF\tFAMILY\tFAMILY_CONF\tORDER\tORDER_CONF\tCLASS\tCLASS_CONF\tPHYLUM\tPHYLUM_CONF\n";

foreach my $queryID ( sort { $a cmp $b } keys %$combinedScore ) {
   
   my $matchName = $combinedMatch->{$queryID};
   
   my $genusConf = &estimateAccuracy('genus', $queryLengths->{$queryID}, $combinedScore->{$queryID});

   my $familyConf = &estimateAccuracy('family', $queryLengths->{$queryID}, $combinedScore->{$queryID});

   my $orderConf = &estimateAccuracy('order', $queryLengths->{$queryID}, $combinedScore->{$queryID});

   my $classConf = &estimateAccuracy('class', $queryLengths->{$queryID}, $combinedScore->{$queryID});

   my $phylumConf = &estimateAccuracy('phylum', $queryLengths->{$queryID}, $combinedScore->{$queryID});

   printf OUT ("$queryID\t$matchName\t$combinedScore->{$queryID}\t$tax->{$matchName}->{'genus'}\t%.3f\t$tax->{$matchName}->{'family'}\t%.3f\t$tax->{$matchName}->{'order'}\t%.3f\t$tax->{$matchName}->{'class'}\t%.3f\t$tax->{$matchName}->{'phylum'}\t%.3f\n", $genusConf, $familyConf, $orderConf, $classConf, $phylumConf);
}

close OUT;

$date = `date`;

chomp $date;

print "done.  [$date]\n\n";

###########################################################################################
# 
# SUBROUTINES
# 
###########################################################################################

# Scan a directory for all the .icm files in it, and save the local path for each into the array we received as an argument.

sub scanDir {
   
   my $dir = shift;

   my $arrayRef = shift;
   
   opendir DOT, $dir or die("Can't open $dir for scanning.\n");

   my @files = readdir DOT;

   closedir DOT;

   foreach my $file ( @files ) {
      
      if ( $file =~ /\.icm$/ ) {
	 
	 push @$arrayRef, "$dir/$file";
      }
   }
}

# Generate confidence estimates from the two-variable polynomials fitted to experimental accuracy data.

sub estimateAccuracy {
   
   my $clade = shift;

   my $readLength = shift;

   my $phymmblScore = shift;
   
   my $retVal;
   
   if ( $readLength < 1 ) {
      
      return 0.0;
   }

   if ( $clade eq 'genus' ) {
      
      $retVal = 0.86 - ((-$phymmblScore - (1.25*$readLength - 550))**(6 + 12/$readLength)) / 1.5e17 + 8/$readLength;
      
   } elsif ( $clade eq 'family' ) {
      
      $retVal = 0.97 - ((-$phymmblScore - (1.25*$readLength - 550))**(6 + 12/$readLength)) / 1.5e17 + 8/$readLength;

   } elsif ( $clade eq 'order' ) {
      
      $retVal = 0.98 - ((-$phymmblScore - (1.25*$readLength - 550))**(6 + 12/$readLength)) / 1.5e17 + 8/$readLength;

   } elsif ( $clade eq 'class' ) {
      
      $retVal = 0.99 - ((-$phymmblScore - (1.25*$readLength - 550))**(6 + 10/$readLength)) / 1.5e17 + 8/$readLength;

   } elsif ( $clade eq 'phylum' ) {
      
      $retVal = 0.99 - ((-$phymmblScore - (1.25*$readLength - 550))**(5.9 + 12/$readLength)) / 1.5e17 + 8/$readLength;
   }

   if ( $retVal eq 'nan' or $retVal > 1 ) {
	 
      $retVal = 0.99999;

   } elsif ( $retVal < 0 ) {
	 
      $retVal = 0.00001;
   }

   return $retVal;
}

