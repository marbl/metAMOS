#!/usr/bin/perl

### david.studholme@tsl.ac.uk

### Generates contigs (in FastA) and scaffolding information (in AGP) from Velvet 'contigs.fa' supercontigs file

### Use entirely at you own risk!! There may be bugs!

### modified by chienchi@lanl.gov  2010/09/13
### add flags -i -size -o
### change naming convention for HMP assembly purpose

### modified by mpop@umd.edu 2010/9/15
### added flag -name
### allowing to change names of sequences  

### modified by mpop@umd.edu 2010/9/27
### removed 10bp limit for inter-contig gaps

use strict;
use warnings ;
use Bio::SeqIO ;
use Getopt::Long;
use File::Basename;

my $file;
my $outDir;
my $scaf_size_cutoff=0;
my $name_prefix="PGA";

GetOptions( 
            'i=s' => \$file,     # scaf seq in fasta
            'size=i' => \$scaf_size_cutoff,
            'o=s'   => \$outDir,
            'name=s' => \$name_prefix);

die "Usage: $0 -i <sequence file> -o <out_dir> [-size <scaf_size_cutoff>] [-name <name>]\n" unless $file;
my ($file_name, $path, $suffix)=fileparse("$file", qr/\.[^.]*/);
my ($sample_name,$center)=split /_/,$file_name;

### Output file for contigs in Fasta format
### The various output file names can be tuned here.
my $contig_outfile = "$outDir/$name_prefix.contigs.fa";
my $scaffolds_outfile = "$outDir/$name_prefix.scaffolds.fa";
my $agp_outfile = "$outDir/$name_prefix.agp";


#open (FILE, ">$outdir/$contig_outfile") and
#   warn "Will write contigs to file '$contig_outfile' to $outDir Directory\n" or
#   die "Failed to write to file '$fasta_outfile'\n";

my $contig_out = Bio::SeqIO->new('-file' => ">$contig_outfile",
				 '-format' => 'Fasta') and
    warn "Will write contigs to file '$contig_outfile' to $outDir Directory\n" or
    die "Failed to write to file '$contig_outfile'\n"; 
   
open (AGP_FILE, ">$agp_outfile") and
   warn "Will write contigs to file $agp_outfile to $outDir Directory\n" or
   die "Failed to write to file $agp_outfile\n";

#open (SCAFF_FILE, ">$outDir/$scaffolds_outfile") and
my $scaff_out = Bio::SeqIO->new('-file' => ">$scaffolds_outfile",
				'-format' => 'Fasta') and
    warn "Will write scaffolds to file $scaffolds_outfile to $outDir Directory\n" or 
    die "Failed to write to file $scaffolds_outfile\n";

print AGP_FILE "# Generated from SOAPdenovo assembly file $file using script $0\n";

warn "Scaffold Sequence cutoff = $scaf_size_cutoff nt\n" if ($scaf_size_cutoff);


my $i = 0;# a counter, used for generating unique contig names

my $inseq = Bio::SeqIO->new('-file' => "<$file",
               '-format' => 'Fasta' ) ;


while (my $seq_obj = $inseq->next_seq ) {

   my $supercontig_id = $seq_obj->id ;
   my $supercontig_seq = $seq_obj->seq ;
   my $supercontig_desc = $seq_obj->description ;
   my $supercontig_length = length($supercontig_seq);

# only process the long supercontigs
   next if ($supercontig_length < $scaf_size_cutoff);


   ### NCBI do not allow coverage and length information in the FastA identifier
   ### e.g. NODE_1160_length_397673_cov_14.469489 is an illegal FastA ID
   ### So we will replace these with simple numbers
   if ($supercontig_id =~ m/NODE_(\d+)_length_\d+_cov_\d+/ or
       $supercontig_id =~ m/^(\d+)$/ or
       $supercontig_id =~ m/^scaffold(\d+)$/) {
       $supercontig_id = "${name_prefix}_scaffold_$1";
   }


#   print SCAFF_FILE ">$supercontig_id\n$supercontig_seq\n";
   my $scf = Bio::PrimarySeq->new(-seq => "$supercontig_seq",
				  -id => "$supercontig_id");

   $scaff_out->write_seq($scf);

   my $start_pos = 1; # keep track of whereabouts in this supercontig we are
     my %substring_sequences;
   foreach my $substring_sequence ( split /(N+)/i, $supercontig_seq ) {
   #warn "\n$substring_sequence\n" if $supercontig_id eq '1160'; for #debugging only

   ### Define the AGP column contents
   my $object1 = $supercontig_id;
   my $object_beg2 = $start_pos;
   my $object_end3 = $start_pos + length($substring_sequence) - 1;
   my $part_number4 = $i;
   my $component_type5;
   my $component_id6a;
   my $gap_length6b;
   my $component_beg7a;
   my $gap_type7b;
   my $component_end8a;
   my $linkage8b;
   my $orientation9a;
   my $filler9b;
     if (  $substring_sequence =~ m/^N+$/i ) {
       ### This is poly-N gap between contigs
       $component_type5 = 'N';
       $gap_length6b = length($substring_sequence);
       $gap_type7b = 'fragment';
       $linkage8b = 'yes';
       $filler9b = '';
         } elsif ( $substring_sequence =~ m/^[ACGT]+$/i ) {
       ### This is a contig
       $i++; # a counter, used for generating unique contig names
       $component_type5 = 'W';
       $component_id6a = "${name_prefix}_$i";
       $component_beg7a = 1;
       $component_end8a = length($substring_sequence);
       $orientation9a = '+';
             ### Print FastA formatted contig
#       print FILE ">$component_id6a\n$substring_sequence\n";
       my $ctg = Bio::PrimarySeq->new(-seq => "$substring_sequence",
				      -id => "$component_id6a");

       $contig_out->write_seq($ctg);


   } else {
       die "Illegal characters in sequence\n$substring_sequence\n";
   }
     $start_pos += length ($substring_sequence);
     if ($component_type5 eq 'N') {
       ### print AGP line for gap
       print AGP_FILE "$object1\t$object_beg2\t$object_end3\t$part_number4\t$component_type5\t$gap_length6b\t$gap_type7b\t$linkage8b\t$filler9b\n";
   } else {
       ### print AGP line for contig
       print AGP_FILE "$object1\t$object_beg2\t$object_end3\t$part_number4\t$component_type5\t$component_id6a\t$component_beg7a\t$component_end8a\t$orientation9a\n";
         }
   }
} 

$contig_out->close();
close AGP_FILE;
$scaff_out->close();
