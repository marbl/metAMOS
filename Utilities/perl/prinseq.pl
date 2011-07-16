#!/usr/bin/perl

#===============================================================================
#   Author: Robert SCHMIEDER, Computational Science Research Center @ SDSU, CA
#
#   File: prinseq-lite
#   Date: 2011-06-05
#   Version: 0.15 lite
#
#   Usage:
#      prinseq-lite [options]
#
#      Try 'prinseq-lite -h' for more information.
#
#    Purpose: PRINSEQ will help you to preprocess your genomic or metagenomic
#             sequence data in FASTA or FASTQ format. The lite version does not
#             require any non-core perl modules for processing.
#
#    Bugs: Please use http://sourceforge.net/tracker/?group_id=315449
#
#===============================================================================

use strict;
use warnings;

#use Data::Dumper; ###
use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile); #for output files
use Fcntl qw(:flock SEEK_END); #for log file
use Cwd;

$| = 1; # Do not buffer output

my $WINDOWSIZE = 64;
my $WINDOWSTEP = 32;
my $WORDSIZE = 3;
my $LOG62 = log(62);
my $LINE_WIDTH = 60;
my $TRIM_QUAL_WINDOW = 1;
my $TRIM_QUAL_STEP = 1;
my $TRIM_QUAL_TYPE = 'min';
my $TRIM_QUAL_RULE = 'lt';
my %MIDS = (ACGAGTGCGT => 0,
            ACGCTCGACA => 0,
            AGACGCACTC => 0,
            AGCACTGTAG => 0,
            ATCAGACACG => 0,
            ATATCGCGAG => 0,
            CGTGTCTCTA => 0,
            CTCGCGTGTC => 0,
            TAGTATCAGC => 0,
            TCTCTATGCG => 0,
            TGATACGTCT => 0,
            TACTGAGCTA => 0,
            CATAGTAGTG => 0,
            CGAGAGATAC => 0,
            ACACGACGACT => 0,
            ACACGTAGTAT => 0,
            ACACTACTCGT => 0,
            ACGACACGTAT => 0,
            ACGAGTAGACT => 0,
            ACGCGTCTAGT => 0,
            ACGTACACACT => 0,
            ACGTACTGTGT => 0,
            ACGTAGATCGT => 0,
            ACTACGTCTCT => 0,
            ACTATACGAGT => 0,
            ACTCGCGTCGT => 0);
my $MIDCHECKLENGTH = 15; #maximum MID length plus possible key length (by default 4 bp for 454)
my $VERSION = '0.15';
my $WHAT = 'lite';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print "PRINSEQ-$WHAT $VERSION\n"; exit; },
            'fastq=s',
            'fasta=s',
            'qual=s',
            'min_len=i',
            'max_len=i',
            'range_len=s',
            'min_gc=i',
            'max_gc=i',
            'range_gc=s',
            'min_qual_score=i',
            'max_qual_score=i',
            'min_qual_mean=i',
            'max_qual_mean=i',
            'ns_max_p=i',
            'ns_max_n=i',
            'noniupac',
            'seq_num=i',
            'derep=i',
            'lc_method=s',
            'lc_threshold=i',
            'trim_to_len=i',
            'trim_left=i',
            'trim_right=i',
            'trim_tail_left=i',
            'trim_tail_right=i',
            'trim_ns_left=i',
            'trim_ns_right=i',
            'trim_qual_left=i',
            'trim_qual_right=i',
            'trim_qual_type=s',
            'trim_qual_rule=s',
            'trim_qual_window=i',
            'trim_qual_step=i',
            'seq_case=s',
            'dna_rna=s',
            'line_width=i',
            'rm_header',
            'seq_id=s',
            'out_format=i',
            'out_good=s',
            'out_bad=s',
            'stats_len',
            'stats_dinuc',
            'stats_info',
            'stats_tag',
            'stats_dupl',
            'stats_ns',
            'stats_all',
            'log:s',
            'si13'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

PRINSEQ - PReprocessing and INformation of SEQuence data

=head1 VERSION

PRINSEQ-lite 0.15

=head1 SYNOPSIS

perl prinseq-lite.pl [-h] [-help] [-version] [-man] [-verbose] [-fastq input_fastq_file] [-fasta input_fasta_file] [-qual input_quality_file] [-min_len int_value] [-max_len int_value] [-range_len ranges] [-min_gc int_value] [-max_gc int_value] [-range_gc ranges] [-min_qual_score int_value] [-max_qual_score int_value] [-min_qual_mean int_value] [-max_qual_mean int_value] [-ns_max_p int_value] [-ns_max_n int_value] [-noniupac] [-seq_num int_value] [-derep int_value] [-lc_method method_name] [-lc_threshold int_value] [-trim_to_len int_value] [-trim_left int_value] [-trim_right int_value] [-trim_ns_left int_value] [-trim_ns_right int_value] [-trim_tail_left int_value] [-trim_tail_right int_value] [-trim_qual_left int_value] [-trim_qual_right int_value] [-trim_qual_type type] [-trim_qual_rule rule] [-trim_qual_window int_value] [-trim_qual_step int_value] [-seq_case case] [-dna_rna type] [-line_width int_value] [-rm_header] [-seq_id id_string] [-out_format int_value] [-out_good filename_prefix] [-out_bad filename_prefix] [si13] [stats_info] [stats_len] [stats_dinuc] [stats_tag] [stats_dupl] [stats_ns] [stats_all] [-log file]

=head1 DESCRIPTION

PRINSEQ will help you to preprocess your genomic or metagenomic sequence data in FASTA (and QUAL) or FASTQ format. The lite version does not require any non-core perl modules for processing.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-fastq> <file>

Input file in FASTQ format that contains the sequence and quality data. Use stdin instead of a file name to read from STDIN (-fasta stdin). This can be useful to process compressed files using Unix pipes.

=item B<-fasta> <file>

Input file in FASTA format that contains the sequence data. Use stdin instead of a file name to read from STDIN (-fastq stdin). This can be useful to process compressed files using Unix pipes.

=item B<-qual> <file>

Input file in QUAL format that contains the quality data.

=item B<-si13>

Quality data in FASTQ file is in Solexa/Illumina 1.3+ format and should be scaled to Phred quality scores ranging from 0 to 40. (Not required for Solexa/Illumina 1.5+, Sanger, Roche/454, Ion Torrent, PacBio data.)

=item B<***** OUTPUT OPTIONS *****>

=item B<-out_format> <integer>

To change the output format, use one of the following options. If not defined, the output format will be the same as the input format.

1 (FASTA only), 2 (FASTA and QUAL) or 3 (FASTQ)

=item B<-out_good> <string>

By default, the output files are created in the same directory as the input file containing the sequence data with an additional "_prinseq_good_XXXX" in their name (where XXXX is replaced by random characters to prevent overwriting previous files). To change the output filename and location, specify the filename using this option. The file extension will be added automatically. Use "-out_good null" to prevent the program from generating the output file(s) for data passing all filters. Use "-out_good stdout" to write data passing all filters to STDOUT (only for FASTA or FASTQ output files).

Example: use "file_passed" to generate the output file file_passed.fasta in the current directory

=item B<-out_bad> <string>

By default, the output files are created in the same directory as the input file containing the sequence data with an additional "_prinseq_bad_XXXX" in their name (where XXXX is replaced by random characters to prevent overwriting previous files). To change the output filename and location, specify the filename using this option. The file extension will be added automatically. Use "-out_bad null" to prevent the program from generating the output file(s) for data not passing any filter. Use "-out_bad stdout" to write data not passing any filter to STDOUT (only for FASTA or FASTQ output files).

Example: use "file_filtered" to generate the output file file_filtered.fasta in the current directory

Example: "-out_good stdout -out_bad null" will write data passing filters to STDOUT and data not passing any filter will be ignored

=item B<-log> <file>

Log file to keep track of parameters, errors, etc. The log file name is optional. If no file name is given, the log file name will be "inputname.log". If the log file already exists, new content will be added to the file.

=item B<***** FILTER OPTIONS *****>

=item B<-min_len> <integer>

Filter sequence shorter than min_len.

=item B<-max_len> <integer>

Filter sequence longer than max_len.

=item B<-range_len> <string>

Filter sequence by length range. Multiple range values should be separated by comma without spaces.

Example: -range_len 50-100,250-300

=item B<-min_gc> <integer>

Filter sequence with GC content below min_gc.

=item B<-max_gc> <integer>

Filter sequence with GC content above max_gc.

=item B<-range_gc> <string>

Filter sequence by GC content range. Multiple range values should be separated by comma without spaces.

Example: -range_gc 50-60,75-90

=item B<-min_qual_score> <integer>

Filter sequence with at least one quality score below min_qual_score.

=item B<-max_qual_score> <integer>

Filter sequence with at least one quality score above max_qual_score.

=item B<-min_qual_mean> <integer>

Filter sequence with quality score mean below min_qual_mean.

=item B<-max_qual_mean> <integer>

Filter sequence with quality score mean above max_qual_mean.

=item B<-ns_max_p> <integer>

Filter sequence with more than ns_max_p percentage of Ns.

=item B<-ns_max_n> <integer>

Filter sequence with more than ns_max_n Ns.

=item B<-noniupac>

Filter sequence with characters other than A, C, G, T or N.

=item B<-seq_num> <integer>

Only keep the first seq_num number of sequences (that pass all other filters).

=item B<-derep> <integer>

Type of duplicates to filter. Allowed values are 1, 2, 3, 4 and 5. Use integers for multiple selections (e.g. 124 to use type 1, 2 and 4). The order does not matter. Option 2 and 3 will set 1 and option 5 will set 4 as these are subsets of the other option.

1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4 (reverse complement exact duplicate), 5 (reverse complement 5'/3' duplicate)

=item B<-lc_method> <string>

Method to filter low complexity sequences. The current options are "dust" and "entropy". Use "-lc_method dust" to calculate the complexity using the dust method.

=item B<-lc_threshold> <integer>

The threshold value (between 0 and 100) used to filter sequences by sequence complexity. The dust method uses this as maximum allowed score and the entropy method as minimum allowed value.

=item B<***** TRIM OPTIONS *****>

=item B<-trim_to_len> <integer>

Trim all sequence from the 3'-end to result in sequence with this length.

=item B<-trim_left> <integer>

Trim sequence at the 5'-end by trim_left positions.

=item B<-trim_right> <integer>

Trim sequence at the 3'-end by trim_right positions.

=item B<-trim_tail_left> <integer>

Trim poly-A/T tail with a minimum length of trim_tail_left at the 5'-end.

=item B<-trim_tail_right> <integer>

Trim poly-A/T tail with a minimum length of trim_tail_right at the 3'-end.

=item B<-trim_ns_left> <integer>

Trim poly-N tail with a minimum length of trim_ns_left at the 5'-end.

=item B<-trim_ns_right> <integer>

Trim poly-N tail with a minimum length of trim_ns_right at the 3'-end.

=item B<-trim_qual_left> <integer>

Trim sequence by quality score from the 5'-end with this threshold score.

=item B<-trim_qual_right> <integer>

Trim sequence by quality score from the 3'-end with this threshold score.

=item B<-trim_qual_type> <string>

Type of quality score calculation to use. Allowed options are min, mean, max and sum. [default: min]

=item B<-trim_qual_rule> <string>

Rule to use to compare quality score to calculated value. Allowed options are lt (less than), gt (greater than) and et (equal to). [default: lt]

=item B<-trim_qual_window> <integer>

The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined, use a window size of 1. [default: 1]

=item B<-trim_qual_step> <integer>

Step size used to move the sliding window. To move the window over all quality scores without missing any, the step size should be less or equal to the window size. [default: 1]

=item B<***** REFORMAT OPTIONS *****>

=item B<-seq_case> <string>

Changes sequence character case to upper or lower case. Allowed options are "upper" and "lower". Use this option to remove soft-masking from your sequences.

=item B<-dna_rna> <string>

Convert sequence between DNA and RNA. Allowed options are "dna" (convert from RNA to DNA) and "rna" (convert from DNA to RNA).

=item B<-line_width> <integer>

Sequence characters per line. Use 0 if you want each sequence in a single line. Use 80 for line breaks every 80 characters. Note that this option only applies to FASTA output files, since FASTQ files store sequences without additional line breaks. [default: 60]

=item B<-rm_header>

Remove the sequence header. This includes everything after the sequence identifier (which is kept unchanged).

=item B<-seq_id> <string>

Rename the sequence identifier. A counter is added to each identifier to assure its uniqueness.

Example: "mySeq_10" will generate the IDs (in FASTA format) >mySeq_101, >mySeq_102, >mySeq_103, ...

=item B<***** SUMMARY STATISTIC OPTIONS *****>

The summary statistic values are written to STDOUT in the form: "parameter_name statistic_name value" (without the quotes). For example, "stats_info reads 10000" or "stats_len max 500". Only one statistic is written per line and values are separated by tabs.

If you specify any statistic option, no other ouput will be generated. To preprocess data, do not specify a statistics option.

=item B<-stats_info>

Outputs basic information such as number of reads (reads) and total bases (bases).

=item B<-stats_len>

Outputs minimum (min), maximum (max), range (range), mean (mean), standard deviation (stddev), mode (mode) and mode value (modeval), and median (median) for read length.

=item B<-stats_dinuc>

Outputs the dinucleotide odds ratio for AA/TT (aatt), AC/GT (acgt), AG/CT (agct), AT (at), CA/TG (catg), CC/GG (ccgg), CG (cg), GA/TC (gatc), GC (gc) and TA (ta).

=item B<-stats_tag>

Outputs the probability of a tag sequence at the 5'-end (prob5) and 3'-end (prob3) in percentage (0..100). Provides the number of predefined MIDs (midnum) and the MID sequences (midseq, separated by comma, only provided if midnum > 0) that occur in more than 34/100 (approx. 3%) of the reads.

=item B<-stats_dupl>

Outputs the number of exact duplicates (exact), 5' duplicates (5), 3' duplicates (3), exact duplicates with reverse complements (exactrevcom) and 5'/3' duplicates with reverse complements (revcomp), and total number of duplicates (total). The maximum number of duplicates is given under the value name with an additional "maxd" (e.g. exactmaxd or 5maxd).

=item B<-stats_ns>

Outputs the number of reads with ambiguous base N (seqswithn), the maximum number of Ns per read (maxn) and the maximum percentage of Ns per read (maxp). The maxn and maxp value are not necessary from the same sequence.

=item B<-stats_all>

Outputs all available summary statistics.


=back

=head1 AUTHOR

Robert SCHMIEDER, C<< <rschmieder_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <rschmieder_at_gmail_dot_com> >> or use http://sourceforge.net/tracker/?group_id=315449 so that I can make PRINSEQ better.

=head1 COPYRIGHT

Copyright (C) 2010-2011  Robert SCHMIEDER

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

#
################################################################################
## DATA AND PARAMETER CHECKING
################################################################################
#

my ($file1,$command,@dataread);

#Check if input file exists and check if file format is correct
if(exists $params{fasta} && exists $params{fastq}) {
    &printError('fasta and fastq cannot be used together');
} elsif(exists $params{fasta}) {
    $command .= ' -fasta '.$params{fasta};
    $file1 = $params{fasta};
    if($params{fasta} eq 'stdin') {
        if(exists $params{qual} && $params{qual} eq 'stdin') {
            &printError('input from STDIN is only allowed for either of the input files');
        } else {
            my $format = &checkInputFormat();
            unless($format eq 'fasta') {
                &printError('input data for -fasta is in '.uc($format).' format not in FASTA format');
            }
        }
    } elsif(-e $params{fasta}) {
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'fasta') {
            &printError('input file for -fasta is in '.uc($format).' format not in FASTA format');
        }
    } else {
        &printError("could not find input file \"".$params{fasta}."\"");
    }
} elsif(exists $params{fastq}) {
    $command .= ' -fastq '.$params{fastq};
    $file1 = $params{fastq};
    if($params{fastq} eq 'stdin') {
        my $format = &checkInputFormat();
        unless($format eq 'fastq') {
            &printError('input data for -fastq is in '.uc($format).' format not in FASTQ format');
        }
    } elsif(-e $params{fastq}) {
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'fastq') {
            &printError('input file for -fastq is in '.uc($format).' format not in FASTQ format');
        }
    } else {
        &printError("could not find input file \"".$params{fastq}."\"");
    }
} else {
    &printError("you did not specify an input file containing the query sequences");
}
if(exists $params{fastq} && exists $params{qual}) {
    &printError('fastq and qual cannot be used together');
} elsif(exists $params{qual}) {
    $command .= ' -qual '.$params{qual};
    if($params{qual} eq 'stdin') {
        &printError('QUAL data cannot be read from STDIN');
    } elsif(-e $params{qual}) {
        #check for file format
        my $format = &checkFileFormat($params{qual});
        unless($format eq 'qual') {
            &printError('input file for -qual is in '.uc($format).' format not in QUAL format');
        }
    } else {
        &printError("could not find input file \"".$params{qual}."\"");
    }
}

#check if stats_all
if(exists $params{stats_all}) {
    $params{stats_info} = 1;
    $params{stats_len} = 1;
    $params{stats_dinuc} = 1;
    $params{stats_tag} = 1;
    $params{stats_dupl} = 1;
    $params{stats_ns} = 1;
    delete($params{stats_all});
}

#check if anything todo
unless( exists $params{min_len} ||
        exists $params{max_len} ||
        exists $params{range_len} ||
        exists $params{min_gc} ||
        exists $params{max_gc} ||
        exists $params{range_gc} ||
        exists $params{min_qual_score} ||
        exists $params{max_qual_score} ||
        exists $params{min_qual_mean} ||
        exists $params{max_qual_mean} ||
        exists $params{ns_max_p} ||
        exists $params{ns_max_n} ||
        exists $params{noniupac} ||
        exists $params{seq_num} ||
        exists $params{derep} ||
        exists $params{lc_method} ||
        exists $params{lc_threshold} ||
        exists $params{trim_to_len} ||
        exists $params{trim_left} ||
        exists $params{trim_right} ||
        exists $params{trim_tail_left} ||
        exists $params{trim_tail_right} ||
        exists $params{trim_ns_left} ||
        exists $params{trim_ns_right} ||
        exists $params{trim_qual_left} ||
        exists $params{trim_qual_right} ||
        exists $params{trim_qual_type} ||
        exists $params{trim_qual_rule} ||
        exists $params{trim_qual_window} ||
        exists $params{trim_qual_step} ||
        exists $params{seq_case} ||
        exists $params{dna_rna} ||
        exists $params{line_width} ||
        exists $params{rm_header} ||
        exists $params{seq_id} ||
        exists $params{out_format} ||
        exists $params{stats_info} ||
        exists $params{stats_len} ||
        exists $params{stats_dinuc} ||
        exists $params{stats_tag} ||
        exists $params{stats_dupl} ||
        exists $params{stats_ns} ||
        exists $params{si13}
        ) {
    &printError('nothing to do with input data');
}
#prevent out of files for stats
if(exists $params{stats_info} || exists $params{stats_len} || exists $params{stats_dinuc} || exists $params{stats_tag} || exists $params{stats_dupl} || exists $params{stats_ns}) {
    $params{out_good} = 'null';
    $params{out_bad} = 'null';
    $params{stats} = 1;
} elsif(exists $params{out_good} && $params{out_good} eq 'null' && exists $params{out_bad} && $params{out_bad} eq 'null') {
    &printError('no output selected (both set to null)');
}

#check if FASTQ file is given for option si13
if(exists $params{si13}) {
    $command .= ' -si13';
    unless(exists $params{fastq}) {
        &printError('option -si13 can only be used for FASTQ input files');
    }
}

#check if output format is possible
if(exists $params{out_format}) {
    $command .= ' -out_format '.$params{out_format};
    if($params{out_format} =~ /\D/) {
        &printError('output format option has to be an integer value');
    } elsif($params{out_format} == 2 || $params{out_format} == 3) {
        unless(exists $params{fastq} || exists $params{qual}) {
            &printError('cannot use this output format option without providing quality data as input');
        }
    } elsif($params{out_format} != 1) {
        &printError('output format option not available');
    }
} else {
    if(exists $params{fastq}) {
        $params{out_format} = 3;
    } elsif(exists $params{fasta} && exists $params{qual}) {
        $params{out_format} = 2;
    } else {
        $params{out_format} = 1;
    }
}

#check if output names are different
if(exists $params{out_good} && exists $params{out_bad} && $params{out_good} eq $params{out_bad} && $params{out_good} ne 'null' && $params{out_good} ne 'stdout') {
    &printError('the output names for -out_good and -out_bad have to be different');
}
#check if output can be written to standard output
if($params{out_format} == 2 && ((exists $params{out_good} && $params{out_good} eq 'stdout') || (exists $params{out_bad} && $params{out_badd} ne 'stdout'))) {
    &printError('the output cannot be written to STDOUT for FASTA and QUAL file output. The option can only be used for FASTA only or FASTQ output');
}

#check dereplication option
#1 - exact dub, 2 - prefix, 3 - suffix, 4 - revcomp exact, 5 - revcomp prefix/suffix
my $derep = 0;
my %dereptypes;
if(exists $params{derep}) {
    $command .= ' -derep '.$params{derep};
    if($params{derep} < 0 || $params{derep} > 54321) {
        &printError('invalid option for dereplication');
    } else {
        my @tmp = split('',$params{derep});
        foreach(@tmp) {
            if($_ < 1 || $_ > 5) {
                &printError('invalid option '.$_.'for dereplication');
            } else {
                $derep = 1;
                $dereptypes{($_-1)} = 0;
            }
        }
    }
}
if(!exists $dereptypes{0} && (exists $dereptypes{1} || exists $dereptypes{2})) {
    $dereptypes{0} = 0;
}
if(!exists $dereptypes{3} && exists $dereptypes{4}) {
    $dereptypes{3} = 0;
}

#check for low complexity method
my $complval;
if(exists $params{lc_method}) {
    $command .= ' -lc_method '.$params{lc_method};
    unless($params{lc_method} eq 'dust' || $params{lc_method} eq 'entropy') {
        &printError('invalid low complexity method');
    }
    unless(exists $params{lc_threshold}) {
        &printError('the low complexity method requires a threshold value specified by -lc_threshold');
    }
    $command .= ' -lc_threshold '.$params{lc_threshold};
    $complval = $params{lc_threshold};
}
if(exists $params{lc_threshold} && !exists $params{lc_method}) {
    &printError('the low complexity threshold requires a method specified by -lc_method');
}

#check for quality trimming
my $trimscore;
if(exists $params{trim_qual_left} || exists $params{trim_qual_right}) {
    $command .= ' -trim_qual_right '.$params{trim_qual_right} if(exists $params{trim_qual_right});
    $command .= ' -trim_qual_left '.$params{trim_qual_left} if(exists $params{trim_qual_left});
    if(exists $params{trim_qual_type}) {
        unless($params{trim_qual_type} eq 'min' || $params{trim_qual_type} eq 'mean' || $params{trim_qual_type} eq 'max' || $params{trim_qual_type} eq 'sum') {
            &printError('invalid value for trim_qual_type');
        }
    } else {
        $params{trim_qual_type} = $TRIM_QUAL_TYPE;
    }
    $command .= ' -trim_qual_type '.$params{trim_qual_type};
    if(exists $params{trim_qual_rule}) {
        unless($params{trim_qual_rule} eq 'lt' || $params{trim_qual_rule} eq 'gt' || $params{trim_qual_rule} eq 'et') {
            &printError('invalid value for trim_qual_rule');
        }
    } else {
        $params{trim_qual_rule} = $TRIM_QUAL_RULE;
    }
    $command .= ' -trim_qual_rule '.$params{trim_qual_rule};
    unless(exists $params{trim_qual_window}) {
        $params{trim_qual_window} = $TRIM_QUAL_WINDOW;
    }
    $command .= ' -trim_qual_rule '.$params{trim_qual_rule};
    unless(exists $params{trim_qual_step}) {
        $params{trim_qual_step} = $TRIM_QUAL_STEP;
    }
    $command .= ' -trim_qual_window '.$params{trim_qual_window};
    $trimscore = 1;
}

#check sequence case
if(exists $params{seq_case}) {
    $command .= ' -seq_case '.$params{seq_case};
    unless($params{seq_case} eq 'upper' || $params{seq_case} eq 'lower') {
        &printError('invalid sequence case option');
    }
}

#check for dna/rna
if(exists $params{dna_rna}) {
    $command .= ' -dna_rna '.$params{dna_rna};
    unless($params{dna_rna} eq 'dna' || $params{dna_rna} eq 'rna') {
        &printError('invalid option for -dna_rna');
    }
}

#set remaining parameters
my $linelen;
if($params{out_format} == 3) {
    $linelen = 0;
} elsif(exists $params{line_width}) {
    $linelen = $params{line_width};
    $command .= ' -line_width '.$params{line_width};
} else {
    $linelen = $LINE_WIDTH;
}


if(exists $params{seq_id}) {
    $command .= ' -seq_id '.$params{seq_id};
    #remove spaces, ">" and quotes from sequence ids
    $params{seq_id} =~ s/\s//g;
    $params{seq_id} =~ s/\>//g;
    $params{seq_id} =~ s/[\"\'\`]//g;
}

my ($repAleft,$repTleft,$repAright,$repTright,$repNleft,$repNright);
if(exists $params{trim_tail_left}) {
    $command .= ' -trim_tail_left '.$params{trim_tail_left};
    $repAleft = 'A'x$params{trim_tail_left};
    $repTleft = 'T'x$params{trim_tail_left};
}
if(exists $params{trim_tail_right}) {
    $command .= ' -trim_tail_right '.$params{trim_tail_right};
    $repAright = 'A'x$params{trim_tail_right};
    $repTright = 'T'x$params{trim_tail_right};
}
if(exists $params{trim_ns_left}) {
    $command .= ' -trim_ns_left '.$params{trim_ns_left};
    $repNleft = 'N'x$params{trim_ns_left};
}
if(exists $params{trim_ns_right}) {
    $command .= ' -trim_ns_right '.$params{trim_ns_right};
    $repNright = 'N'x$params{trim_ns_right};
}

#add remaining to log command
if(exists $params{log}) {
    $command .= ' -log'.($params{log} ? ' '.$params{log} : '');
    if(exists $params{min_len}) {
        $command .= ' -min_len '.$params{min_len};
    }
    if(exists $params{max_len}) {
        $command .= ' -max_len '.$params{max_len};
    }
    if(exists $params{range_len}) {
        $command .= ' -range_len '.$params{range_len};
    }
    if(exists $params{min_gc}) {
        $command .= ' -min_gc '.$params{min_gc};
    }
    if(exists $params{max_gc}) {
        $command .= ' -max_gc '.$params{max_gc};
    }
    if(exists $params{range_gc}) {
        $command .= ' -range_gc '.$params{range_gc};
    }
    if(exists $params{min_qual_score}) {
        $command .= ' -min_qual_score '.$params{min_qual_score};
    }
    if(exists $params{max_qual_score}) {
        $command .= ' -max_qual_score '.$params{max_qual_score};
    }
    if(exists $params{min_qual_mean}) {
        $command .= ' -min_qual_mean '.$params{min_qual_mean};
    }
    if(exists $params{max_qual_mean}) {
        $command .= ' -max_qual_mean '.$params{max_qual_mean};
    }
    if(exists $params{ns_max_p}) {
        $command .= ' -ns_max_p '.$params{ns_max_p};
    }
    if(exists $params{ns_max_n}) {
        $command .= ' -ns_max_n '.$params{ns_max_n};
    }
    if(exists $params{noniupac}) {
        $command .= ' -noniupac';
    }
    if(exists $params{seq_num}) {
        $command .= ' -seq_num '.$params{seq_num};
    }
    if(exists $params{trim_to_len}) {
        $command .= ' -trim_to_len '.$params{trim_to_len};
    }
    if(exists $params{trim_left}) {
        $command .= ' -trim_left '.$params{trim_left};
    }
    if(exists $params{trim_right}) {
        $command .= ' -trim_right '.$params{trim_right};
    }
    if(exists $params{rm_header}) {
        $command .= ' -rm_header';
    }
    if(exists $params{stats_len}) {
        $command .= ' -stats_len';
    }
    if(exists $params{stats_dinuc}) {
        $command .= ' -stats_dinuc';
    }
    if(exists $params{stats_info}) {
        $command .= ' -stats_info';
    }
    if(exists $params{stats_tag}) {
        $command .= ' -stats_tag';
    }
    if(exists $params{stats_dupl}) {
        $command .= ' -stats_dupl';
    }
    if(exists $params{stats_ns}) {
        $command .= ' -stats_ns';
    }
    if(exists $params{verbose}) {
        $command .= ' -verbose';
    }
    if(exists $params{out_good}) {
        $command .= ' -out_good '.$params{out_good};
    }
    if(exists $params{out_bad}) {
        $command .= ' -out_bad '.$params{out_bad};
    }

    unless($params{log}) {
        $params{log} = join("__",$file1||'nonamegiven').'.log';
    }
    $params{log} = cwd().'/'.$params{log} unless($params{log} =~ /^\//);
    &printLog("Executing PRINSEQ with command: \"perl prinseq-".$WHAT.".pl".$command."\"");
}

#
################################################################################
## DATA PROCESSING
################################################################################
#

#order of processing
#seqnum
#trimleft
#trimright
#trimtails
#trimlen
#minlen
#maxlen
#rangelen
#mingc
#maxgc
#rangegc
#nsmaxp
#nsmaxn
#noniupac
#newparams
#compl
#derep
#seqid
#rmheaderl
#seqcase
#dnarna
#linelen


my $filename = $file1;
while($filename =~ /[\w\d]+\.[\w\d]+/) {
    $filename =~ s/\.[\w\d]+$//;
}

#create filehandles for the output data
my ($fhgood,$fhgood2,$fhbad,$fhbad2);
my ($filenamegood,$filenamegood2,$filenamebad,$filenamebad2);
my ($nogood,$nobad,$stdoutgood,$stdoutbad);
$nogood = $nobad = $stdoutgood = $stdoutbad = 0;
if(exists $params{out_good}) {
    if($params{out_good} eq 'null') {
        $nogood = 1;
    } elsif($params{out_good} eq 'stdout') {
        $stdoutgood = 1;
#        *$fhgood = *STDOUT; ###
    } else {
        open($fhgood,">".$params{out_good}.'.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
        $filenamegood = $params{out_good}.'.fast'.($params{out_format} == 3 ? 'q' : 'a');
    }
} else {
    $fhgood = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_XXXX',
                               SUFFIX => '.fast'.($params{out_format} == 3 ? 'q' : 'a'),
                               UNLINK => 0);
    $filenamegood = $fhgood->filename;
}
if(exists $params{out_bad}) {
    if($params{out_bad} eq 'null') {
        $nobad = 1;
    } elsif($params{out_bad} eq 'stdout') {
        $stdoutbad = 1;
    } else {
        open($fhbad,">".$params{out_bad}.'.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
        $filenamebad = $params{out_bad}.'.fast'.($params{out_format} == 3 ? 'q' : 'a');
    }
} else {
    $fhbad = File::Temp->new( TEMPLATE => $filename.'_prinseq_bad_XXXX',
                              SUFFIX => '.fast'.($params{out_format} == 3 ? 'q' : 'a'),
                              UNLINK => 0);
    $filenamebad = $fhbad->filename;
}
if($params{out_format} == 2) {
    if(exists $params{out_good}) {
        unless($nogood) {
            open($fhgood2,">".$params{out_good}.'.qual') or &printError('cannot open output file');
            $filenamegood2 = $params{out_good}.'.qual';
        }
    } else {
        $fhgood2 = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_XXXX',
                                    SUFFIX => '.qual',
                                    UNLINK => 0);
        $filenamegood2 = $fhgood2->filename;
    }
    if(exists $params{out_bad}) {
        unless($nobad) {
            open($fhbad2,">".$params{out_bad}.'.qual') or &printError('cannot open output file');
            $filenamebad2 = $params{out_bad}.'.qual';
        }
    } else {
        $fhbad2 = File::Temp->new( TEMPLATE => $filename.'_prinseq_bad_XXXX',
                                   SUFFIX => '.qual',
                                   UNLINK => 0);
        $filenamebad2 = $fhbad2->filename;
    }
}

my $numlines = 0;
my ($progress,$counter,$part);
$progress = $counter = $part = 1;
if(exists $params{verbose}) {
    print STDERR "Estimate size of input data for status report (this might take a while for large files)\n";
    $numlines = ($file1 eq 'stdin' ? 1 : &getLineNumber($file1));
    print STDERR "\tdone\n";

    #for progress bar
    $progress = 0;
    $counter = 1;
    $part = int($numlines/100);
}

#parse input data
print STDERR "Parse and process input data\n" if(exists $params{verbose});
print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
&printLog("Parsing and processing input data: \"".$file1."\"".(exists $params{qual} ? " and \"".$params{qual}."\"" : ""));
my $numseqs = 0;
my $goodcount = 0;
my $badcount = 0;
my ($tsvfile,$seqid,$header,$seq,$qual,$count,$numbases,$length);
$count = 0;
$qual = '';
$seq = '';
#stats data
my (%stats,%kmers,%odds,%counts);
#parse data
my $seqcount = 0;
my $seqbases = 0;
my $badbases = 0;
my (@seqs,@printtmp);

if($file1 eq 'stdin') {
    *FILE = *STDIN;
#    open(FILE, "<&=STDIN") or die "ERROR: Couldn't alias STDIN: $! \n";
} else {
    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file1 |") or die "ERROR: Could not open file $file1: $! \n";
}
if(exists $params{qual}) {
    if($params{qual} eq 'stdin') {
        &printError('QUAL data cannot be read from STDIN');
#        *FILE2 = *STDIN;
#        open(FILE2, "<&=STDIN") or die "ERROR: Couldn't alias STDIN: $! \n";
    } else {
        open(FILE2,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < ".$params{qual}." |") or die "ERROR: Could not open file ".$params{qual}.": $! \n";
    }
}

if(exists $params{fastq}) {
    foreach(@dataread) {
        chomp();
        if($count == 0 && /^\@(\S+)\s*(.*)$/) {
            $length = length($seq);
            if($seq && $length) {
                #remove space and dash from sequences
                $seq =~ s/[\s\-]//g;
                $numseqs++;
                $numbases += $length;
                &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"") unless($length == length($qual));
                if(exists $params{stats}) {
                    #calc summary stats
                    $seq = uc($seq);
                    &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
                    if(exists $params{stats_dupl}) {
                        push(@seqs,[$seq,$numseqs,$length]);
                    }
                } else {
                    #process data
                    &processData($seqid,$seq,$qual,$header);
                }
            }
            $seqid = $1;
            $header = $2 || '';
            $seq = '';
            $qual = '';
        } elsif($count == 1) {
            $seq = $_;
        } elsif($count == 3) {
            $qual = $_;
            $count = -1;
        }
        $count++;
    }
    while(<FILE>) {
        chomp();
        if($count == 0 && /^\@(\S+)\s*(.*)$/) {
            $length = length($seq);
            if($seq && $length) {
                #remove space and dash from sequences
                $seq =~ s/[\s\-]//g;
                $numseqs++;
                $numbases += $length;
                &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"") unless($length == length($qual));
                if(exists $params{stats}) {
                    #calc summary stats
                    $seq = uc($seq);
                    &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
                    if(exists $params{stats_dupl}) {
                        push(@seqs,[$seq,$numseqs,$length]);
                    }
                } else {
                    #process data
                    &processData($seqid,$seq,$qual,$header);
                }
            }
            $seqid = $1;
            $header = $2 || '';
            $seq = '';
            $qual = '';
        } elsif($count == 1) {
            $seq = $_;
        } elsif($count == 3) {
            $qual = $_;
            $count = -1;
        }
        $count++;
        #progress bar stuff
        $counter++;
        if($counter > $part) {
            $counter = 1;
            $progress++;
            $progress = 99 if($progress > 99);
            print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
        }
    }
    #add last one
    $length = length($seq);
    if($seq && $length) {
        #remove space and dash from sequences
        $seq =~ s/[\s\-]//g;
        $numseqs++;
        $numbases += $length;
        &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"") unless($length == length($qual));
        if(exists $params{stats}) {
            #calc summary stats
            $seq = uc($seq);
            &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
            if(exists $params{stats_dupl}) {
                push(@seqs,[$seq,$numseqs,$length]);
            }
        } else {
            #process data
            &processData($seqid,$seq,$qual,$header);
        }
    }
} elsif(exists $params{fasta}) {
    foreach(@dataread) {
        chomp();
        if(/^>(\S+)\s*(.*)$/) {
            $length = length($seq);
            if($seq && $length) {
                #get qual data if provided
                if(exists $params{qual}) {
                    while(<FILE2>) {
                        chomp();
                        last if(/^>/ && $qual);
                        next if(/^>/);
                        $qual .= $_.' ';
                    }
                    $qual = &convertQualNumsToAsciiString($qual);
                    &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"") unless($length == length($qual));
                }
                #remove space and dash from sequences
                $seq =~ s/[\s\-]//g;
                #print to tsv file
                $numseqs++;
                $numbases += $length;
                if(exists $params{stats}) {
                    #calc summary stats
                    $seq = uc($seq);
                    &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
                    if(exists $params{stats_dupl}) {
                        push(@seqs,[$seq,$numseqs,$length]);
                    }
                } else {
                    #process data
                    &processData($seqid,$seq,$qual,$header);
                }
                $qual = '';
            }
            $seqid = $1;
            $header = $2 || '';
            $seq = '';
        } else {
            $seq .= $_;
        }
    }
    while(<FILE>) {
        chomp();
        if(/^>(\S+)\s*(.*)$/) {
            $length = length($seq);
            if($seq && $length) {
                #get qual data if provided
                if(exists $params{qual}) {
                    while(<FILE2>) {
                        chomp();
                        last if(/^>/ && $qual);
                        next if(/^>/);
                        $qual .= $_.' ';
                    }
                    $qual = &convertQualNumsToAsciiString($qual);
                    &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"") unless($length == length($qual));
                }
                #remove space and dash from sequences
                $seq =~ s/[\s\-]//g;
                #print to tsv file
                $numseqs++;
                $numbases += $length;
                if(exists $params{stats}) {
                    #calc summary stats
                    $seq = uc($seq);
                    &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
                    if(exists $params{stats_dupl}) {
                        push(@seqs,[$seq,$numseqs,$length]);
                    }
                } else {
                    #process data
                    &processData($seqid,$seq,$qual,$header);
                }
                $qual = '';
            }
            $seqid = $1;
            $header = $2 || '';
            $seq = '';
        } else {
            $seq .= $_;
        }
        #progress bar stuff
        $counter++;
        if($counter > $part) {
            $counter = 1;
            $progress++;
            $progress = 99 if($progress > 99);
            print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
        }
    }
    #add last one
    $length = length($seq);
    if($seq && $length) {
        #get qual data if provided
        if(exists $params{qual}) {
            while(<FILE2>) {
                chomp();
                last if(/^>/ && $qual);
                next if(/^>/);
                $qual .= $_.' ';
            }
            $qual = &convertQualNumsToAsciiString($qual);
            &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"") unless($length == length($qual));
        }
        #remove space and dash from sequences
        $seq =~ s/[\s\-]//g;
        #print to tsv file
        $numseqs++;
        $numbases += $length;
        if(exists $params{stats}) {
            #calc summary stats
            $seq = uc($seq);
            &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
            if(exists $params{stats_dupl}) {
                push(@seqs,[$seq,$numseqs,$length]);
            }
        } else {
            #process data
            &processData($seqid,$seq,$qual,$header);
        }
        $qual = '';
    }
}

print STDERR "\r\tdone          \n" if(exists $params{verbose});
&printLog("Done parsing and processing input data");

if($derep && !exists $params{stats}) {
    &printLog("Remove duplicates");
    &derepSeqs();
    &printLog("Done removing duplicates");
}

close(FILE);
close(FILE2) if(exists $params{qual});

print STDERR "Clean up empty files\n" if(exists $params{verbose} && !exists $params{stats});
#close filehandles
close($fhgood) unless($nogood || $stdoutgood);
close($fhbad) unless($nobad || $stdoutbad);
if($params{out_format} == 2) {
    close($fhgood2) unless($nogood || $stdoutgood);
    close($fhbad2) unless($nobad || $stdoutbad);
}
#remove empty files
if($seqcount == 0 && !$nogood && !$stdoutgood) {
    unlink($filenamegood);
    if($params{out_format} == 2) {
        unlink($filenamegood2);
    }
}
if($badcount == 0 && !$nobad && !$stdoutbad) {
    unlink($filenamebad);
    if($params{out_format} == 2) {
        unlink($filenamebad2);
    }
}
print STDERR "\tdone\n" if(exists $params{verbose} && !exists $params{stats});
if(exists $params{verbose} && !exists $params{stats}) {
    print STDERR "Input and filter stats:\n";
    print STDERR "\tInput sequences: ".&addCommas($numseqs)."\n";
    print STDERR "\tInput bases: ".&addCommas($numbases)."\n";
    print STDERR "\tInput mean length: ".sprintf("%.2f",$numbases/$numseqs)."\n" if($numseqs);
    print STDERR "\tGood sequences: ".&addCommas($seqcount)." (".sprintf("%.2f",(100*$seqcount/$numseqs))."%)\n" if($numseqs);
    print STDERR "\tGood bases: ".&addCommas($seqbases)."\n" if($seqcount);
    print STDERR "\tGood mean length: ".sprintf("%.2f",$seqbases/$seqcount)."\n" if($seqcount);
    print STDERR "\tBad sequences: ".&addCommas($badcount)." (".sprintf("%.2f",(100*$badcount/$numseqs))."%)\n" if($numseqs);
    print STDERR "\tBad bases: ".&addCommas($badbases)."\n" if($badcount);
    print STDERR "\tBad mean length: ".sprintf("%.2f",$badbases/$badcount)."\n" if($badcount);
    &printLog("Input sequences: ".&addCommas($numseqs));
    &printLog("Input bases: ".&addCommas($numbases));
    &printLog("Input mean length: ".sprintf("%.2f",$numbases/$numseqs)) if($numseqs);
    &printLog("Good sequences: ".&addCommas($seqcount)." (".sprintf("%.2f",(100*$seqcount/$numseqs))."%)") if($numseqs);
    &printLog("Good bases: ".&addCommas($seqbases)) if($seqcount);
    &printLog("Good mean length: ".sprintf("%.2f",$seqbases/$seqcount)) if($seqcount);
    &printLog("Bad sequences: ".&addCommas($badcount)." (".sprintf("%.2f",(100*$badcount/$numseqs))."%)") if($numseqs);
    &printLog("Bad bases: ".&addCommas($badbases)) if($badcount);
    &printLog("Bad mean length: ".sprintf("%.2f",$badbases/$badcount)) if($badcount);
}

#print summary stats
if(exists $params{stats}) {
    if(exists $params{stats_info}) {
        $stats{stats_info}->{reads} = $numseqs;
        $stats{stats_info}->{bases} = $numbases;
    }
    if(exists $params{stats_len}) {
        $stats{stats_len} = &generateStats($counts{length});
    }
    if(exists $params{stats_dinuc}) {
        foreach my $i (keys %odds) {
	    $stats{stats_dinuc}->{lc($i)} = sprintf("%.9f",$odds{$i}/$numseqs);
	}
    }
    if(exists $params{stats_tag}) {
        #calculate frequency of 5-mers
        my $kmersum = &getTagFrequency(\%kmers);
        #check for frequency of MID tags
        my $midsum = 0;
        my $midcount = 0;
        my @midseqs;
        foreach my $mid (keys %MIDS) {
            $midsum += $MIDS{$mid};
            if($MIDS{$mid} > $numseqs/34) { #in more than 34/100 (approx. 3%) as this is estimated average error for MIDs
                $midcount++;
                push(@midseqs,$mid);
            }
        }
        $stats{stats_tag}->{midnum} = $midcount;
        if($midcount) {
            $stats{stats_tag}->{midseq} = join(",",@midseqs);
        }
        if($midsum > $kmersum->{5}) {
            $kmersum->{5} = $midsum;
        }
        foreach my $kmer (keys %$kmersum) {
            $stats{stats_tag}->{'prob'.$kmer} = sprintf("%d",(100/$numseqs*$kmersum->{$kmer}));
	}
    }
    if(exists $params{stats_dupl}) {
        #empty vars before n-plicate check
        %counts = %kmers = %odds = ();
        #0 - exact dub, 1 - prefix, 2 - suffix, 3 - revcomp exact, 4 - revcomp prefix/suffix
        my %types = (0 => 'exact', 1 => '5', 2 => '3', 3 => 'exactrevcomp', 4 => 'revcomp');
        my ($dupls,undef) = &checkForDupl(\@seqs,\%types,$numseqs);
        #set zero counts
        foreach my $s (keys %types) {
            $stats{stats_dupl}->{$types{$s}} = 0;
            $stats{stats_dupl}->{$types{$s}.'maxd'} = 0;
        }
        foreach my $n (keys %$dupls) {
	    foreach my $s (keys %{$dupls->{$n}}) {
		$stats{stats_dupl}->{$types{$s}} += $dupls->{$n}->{$s} * $n;
		$stats{stats_dupl}->{$types{$s}.'maxd'} = $n unless($stats{stats_dupl}->{$types{$s}.'maxd'} > $n);
		$stats{stats_dupl}->{total} += $dupls->{$n}->{$s} * $n;
	    }
	}
    }
    foreach my $type (sort keys %stats) {
        foreach my $value (sort keys %{$stats{$type}}) {
            print STDOUT join("\t",$type,$value,(defined $stats{$type}->{$value} ? $stats{$type}->{$value} : '-'))."\n";
        }
    }
}

##
#################################################################################
### MISC FUNCTIONS
#################################################################################
##

sub printError {
    my $msg = shift;
    print STDERR "\nERROR: ".$msg.".\n\nTry \'perl prinseq-".$WHAT.".pl -h\' for more information.\nExit program.\n";
    &printLog("ERROR: ".$msg.". Exit program.\n");
    exit(0);
}

sub printWarning {
    my $msg = shift;
    print STDERR "WARNING: ".$msg.".\n";
    &printLog("WARNING: ".$msg.".\n");
}

sub printLog {
    my $msg = shift;
    if(exists $params{log}) {
        my $time = sprintf("%02d/%02d/%04d %02d:%02d:%02d",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
        open(FH, ">>", $params{log}) or die "ERROR: Can't open file ".$params{log}.": $! \n";
        flock(FH, LOCK_EX) or die "ERROR: Cannot lock file ".$params{log}.": $! \n";
        print FH "[prinseq-".$WHAT."-$VERSION] [$time] $msg\n";
        flock(FH, LOCK_UN) or die "ERROR: cannot unlock ".$params{log}.": $! \n";
        close(FH);
    }
}

sub addCommas {
    my $num = shift;
    return unless(defined $num);
    return $num if($num < 1000);
    $num = scalar reverse $num;
    $num =~ s/(\d{3})/$1\,/g;
    $num =~ s/\,$//;
    $num = scalar reverse $num;
    return $num;
}

sub getLineNumber {
    my $file = shift;
    my $lines = 0;
    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    $lines += tr/\n/\n/ while sysread(FILE, $_, 2 ** 16);
    close(FILE);
    return $lines;
}


sub checkFileFormat {
    my $file = shift;

    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    while (<FILE>) {
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && /^[ACGTNacgtn]+/) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && /^[ACGTNacgtn]+/) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        }
    }
    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}

sub checkInputFormat {
    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    while (<STDIN>) {
        push(@dataread,$_);
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && /^[ACGTNacgtn]+/) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && /^[ACGTNacgtn]+/) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        }
    }

    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}

sub getArrayMean {
    my $array = shift;
    my $sum = 0;
    foreach(@$array) {
        $sum += $_;
    }
    return $sum/scalar(@$array);
}

sub getArrayMin {
    my $array = shift;
    my $min = $array->[0];
    foreach(@$array) {
        $min = $_ if($min > $_);
    }
    return $min;
}

sub getArrayMax {
    my $array = shift;
    my $max = $array->[0];
    foreach(@$array) {
        $max = $_ if($max < $_);
    }
    return $max;
}

sub getArraySum {
    my $array = shift;
    my $sum = 0;
    foreach(@$array) {
        $sum += $_;
    }
    return $sum;
}

sub convertQualNumsToAscii {
    my $qual = shift;
    my @ascii;
    $qual =~ s/^\s+//;
    $qual =~ s/\s+$//;
    my @nums = split(/\s+/,$qual);
    foreach(@nums) {
	push(@ascii,&convertQualNumToAscii($_));
    }
    return \@ascii;
}

sub convertQualNumsToAsciiString {
    my $qual = shift;
    my $ascii;
    $qual =~ s/^\s+//;
    $qual =~ s/\s+$//;
    my @nums = split(/\s+/,$qual);
    foreach(@nums) {
	$ascii .= &convertQualNumToAscii($_);
    }
    return $ascii;
}

sub convertQualAsciiToNums {
    my $qual = shift;
    my @nums;
    my @ascii = split('',$qual);
    foreach(@ascii) {
	push(@nums,&convertQualAsciiToNum($_));
    }
    return \@nums;
}

sub convertQualNumToAscii {
    my $qual = shift;
    return chr(($_ <= 93 ? $_ : 93) + 33);
}

sub convertQualAsciiToNum {
    my $qual = shift;
    return (ord($qual) - 33);
}

sub convertQualAsciiToNumSI13 {
    my $qual = shift;
    return ((10 * log(1 + 10 ** (ord($qual) - 64) / 10.0)) / log(10));
}

sub convertQualAsciiToNumsSI13 {
    my $qual = shift;
    my @nums;
    my @ascii = split('',$qual);
    foreach(@ascii) {
	push(@nums,&convertQualAsciiToNumSI13($_));
    }
    return \@nums;
}

sub convertQualArrayToString {
    my ($nums,$linelen) = @_;
    $linelen = 80 unless($linelen);
    my $str;
    my $count = 1;
    foreach my $n (@$nums) {
        $str .= ($n < 10 ? ' '.$n : $n).' ';
        if(++$count > $linelen) {
            $count = 1;
            $str =~ s/\s$//;
            $str .= "\n";
        }
    }
    $str =~ s/[\s\n]$//;
    return $str;
}

sub checkRange {
    my ($range,$val) = @_;
    my @ranges = split(/\,/,$range);
    foreach my $r (@ranges) {
	my @tmp = split(/\-/,$r);
	return 0 if($val < $tmp[0] || $val > $tmp[1]);
    }
    return 1;
}

#process sequence (and qual) data
sub processData {
    my ($sid,$seq,$qual,$header) = @_;
    #assume sequence is good ;-)
    my $good = 1;
    my $seqn = uc($seq);
    my $qualn = $qual;
    my $begin = 0;
    my $end = 0;
    my $length;
    my $bylength;

    #check for maximum number sequences requested
    $good = 0 if(exists $params{seq_num} && $params{seq_num} <= $seqcount);

    #trim sequence ends
    if($good && exists $params{trim_left}) {
        $begin += $params{trim_left};
        $seqn = substr($seqn,$begin);
        $qualn = substr($qualn,$begin) if(defined $qualn && length($qualn));
    }
    if($good && exists $params{trim_right}) {
        $end += $params{trim_right};
        $length = length($seqn);
        $seqn = substr($seqn,0,$length-$end);
        $qualn = substr($qualn,0,$length-$end) if(defined $qualn && length($qualn));
    }

    #check for quality scores
    if($good && defined $qualn && $trimscore) {
        my $quals;
        if(exists $params{si13}) { #scale data to Phred scale if necessary (e.g. $fileformat == 1)
            $quals = &convertQualAsciiToNumsSI13($qualn);
        } else {
            $quals = &convertQualAsciiToNums($qualn);
        }
        $length = scalar(@$quals);
        my $i = 0;
        my $begintmp = 0;
        my $endtmp = 0;
        my ($window,$val);
        #left
        if(exists $params{trim_qual_left}) {
            while($i < $length) {
                #calculate maximum window
                $window = ($i+$params{trim_qual_window} <= $length ? $params{trim_qual_window} : ($length-$i));
                #calculate value used to compare with given value
                if($window == 1) {
                    $val = $quals->[$i];
                } elsif($params{trim_qual_type} eq 'min') {
                    $val = &getArrayMin([@$quals[$i..($i+$window-1)]]);
                } elsif($params{trim_qual_type} eq 'max') {
                    $val = &getArrayMax([@$quals[$i..($i+$window-1)]]);
                } elsif($params{trim_qual_type} eq 'mean') {
                    $val = &getArrayMean([@$quals[$i..($i+$window-1)]]);
                } elsif($params{trim_qual_type} eq 'sum') {
                    last if($window < $params{trim_qual_window});
                    $val = &getArraySum([@$quals[$i..($i+$window-1)]]);
                } else {
                    last;
                }
                #compare values
                if(($params{trim_qual_rule} eq 'lt' && $val < $params{trim_qual_left}) || ($params{trim_qual_rule} eq 'gt' && $val > $params{trim_qual_left}) || ($params{trim_qual_rule} eq 'et' && $val == $params{trim_qual_left})) {
                    $begintmp += $params{trim_qual_step};
                    $i += $params{trim_qual_step};
                } else {
                    last;
                }
            }
            $seqn = substr($seqn,$begintmp);
            $qualn = substr($qualn,$begintmp);
            $begin += $begintmp;
        }
        #right
        if(exists $params{trim_qual_right}) {
            if($begintmp) {
                @$quals = @$quals[$begintmp..($length-1)];
                $length = scalar(@$quals);
            }
            @$quals = reverse(@$quals);
            $i = 0;
            while($i < $length) {
                #calculate maximum window
                $window = ($i+$params{trim_qual_window} <= $length ? $params{trim_qual_window} : ($length-$i));
                #calculate value used to compare with given value
                if($window == 1) {
                    $val = $quals->[$i];
                } elsif($params{trim_qual_type} eq 'min') {
                    $val = &getArrayMin([@$quals[$i..($i+$window-1)]]);
                } elsif($params{trim_qual_type} eq 'max') {
                    $val = &getArrayMax([@$quals[$i..($i+$window-1)]]);
                } elsif($params{trim_qual_type} eq 'mean') {
                    $val = &getArrayMean([@$quals[$i..($i+$window-1)]]);
                } elsif($params{trim_qual_type} eq 'sum') {
                    last if($window < $params{trim_qual_window});
                    $val = &getArraySum([@$quals[$i..($i+$window-1)]]);
                } else {
                    last;
                }
                #compare values
                if(($params{trim_qual_rule} eq 'lt' && $val < $params{trim_qual_right}) || ($params{trim_qual_rule} eq 'gt' && $val > $params{trim_qual_right}) || ($params{trim_qual_rule} eq 'et' && $val == $params{trim_qual_right})) {
                    $endtmp += $params{trim_qual_step};
                    $i += $params{trim_qual_step};
                } else {
                    last;
                }
            }
            $seqn = substr($seqn,0,$length-$endtmp);
            $qualn = substr($qualn,0,$length-$endtmp);
            $end += $endtmp;
        }
    }

    #check for tails with min trimtails char repeats
    if($good && exists $params{trim_tail_left}) {
        $length = length($seqn);
        my $begintmp = 0;
        if($seqn =~ m/^$repAleft/ || $seqn =~ m/^$repTleft/) {
            my @tmp = split(//,$seqn);
            my $tmpchar = $tmp[0]; #A or T
            $begintmp += $params{trim_tail_left};
            foreach ($params{trim_tail_left}..$length-1) {
                last unless($tmp[$_] eq $tmpchar || $tmp[$_] eq 'N');
                $begintmp++;
            }
            $seqn = substr($seqn,$begintmp);
            $qualn = substr($qualn,$begintmp) if(defined $qualn && length($qualn));
            $length = length($seqn);
            $begin += $begintmp;
        }
    }
    if($good && exists $params{trim_tail_right}) {
        $length = length($seqn);
        my $endtmp = 0;
        if($seqn =~ m/$repAright$/ || $seqn =~ m/$repTright$/) {
            my @tmp = split(//,$seqn);
            my $tmpchar = $tmp[$length-1]; #A or T
            $endtmp += $params{trim_tail_right};
            foreach (reverse 0..$length-$params{trim_tail_right}-1) {
                last unless($tmp[$_] eq $tmpchar || $tmp[$_] eq 'N');
                $endtmp++;
            }
            $seqn = substr($seqn,0,$length-$endtmp);
            $qualn = substr($qualn,0,$length-$endtmp) if(defined $qualn && length($qualn));
            $end += $endtmp;
        }
    }
    if($good && exists $params{trim_ns_left}) {
        $length = length($seqn);
        my $begintmp = 0;
        if($seqn =~ m/^$repNleft/) {
            my @tmp = split(//,$seqn);
            $begintmp += $params{trim_ns_left};
            foreach ($params{trim_ns_left}..$length-1) {
                last unless($tmp[$_] eq 'N');
                $begintmp++;
            }
            $seqn = substr($seqn,$begintmp);
            $qualn = substr($qualn,$begintmp) if(defined $qualn && length($qualn));
            $length = length($seqn);
            $begin += $begintmp;
        }
    }
    if($good && exists $params{trim_ns_right}) {
        $length = length($seqn);
        my $endtmp = 0;
        if($seqn =~ m/$repNright$/) {
            my @tmp = split(//,$seqn);
            $endtmp += $params{trim_ns_right};
            foreach (reverse 0..$length-$params{trim_ns_right}-1) {
                last unless($tmp[$_] eq 'N');
                $endtmp++;
            }
            $seqn = substr($seqn,0,$length-$endtmp);
            $qualn = substr($qualn,0,$length-$endtmp) if(defined $qualn && length($qualn));
            $end += $endtmp;
        }
    }

    #check if trim to certain length
    $length = length($seqn);
    if($good && exists $params{trim_to_len} && $length > $params{trim_to_len}) {
        $seqn = substr($seqn,0,$params{trim_to_len});
        $qualn = substr($qualn,0,$params{trim_to_len}) if(defined $qualn && length($qualn));
        $end += ($length-$params{trim_to_len});
    }

    #check for sequence length
    $length = length($seqn);
    $bylength = ($length ? 100/$length : 0);
    if($bylength == 0) {
        $good = 0;
    }
    if($good && exists $params{min_len} && $length < $params{min_len}) {
        $good = 0;
    }
    if($good && exists $params{max_len} && $length > $params{max_len}) {
        $good = 0;
    }
    if($good && exists $params{range_len} && !&checkRange($params{range_len},$length)) {
        $good = 0;
    }

    #check for quality scores
    if($good && defined $qualn && (exists $params{min_qual_score} || exists $params{max_qual_score} || exists $params{min_qual_mean} || exists $params{max_qual_mean})) {
        my $quals;
        if(exists $params{si13}) { #scale data to Phred scale if necessary (e.g. $fileformat == 1)
            $quals = &convertQualAsciiToNumsSI13($qualn);
        } else {
            $quals = &convertQualAsciiToNums($qualn);
        }
        if($good && exists $params{min_qual_score} && &getArrayMin($quals) < $params{min_qual_score}) {
            $good = 0;
        }
        if($good && exists $params{max_qual_score} && &getArrayMax($quals) < $params{max_qual_score}) {
            $good = 0;
        }
        if($good && exists $params{min_qual_mean} && &getArrayMean($quals) < $params{min_qual_mean}) {
            $good = 0;
        }
        if($good && exists $params{max_qual_mean} && &getArrayMean($quals) < $params{max_qual_mean}) {
            $good = 0;
        }
    }

    #check for GC content
    if($good && (exists $params{min_gc} || exists $params{max_gc} || exists $params{range_gc})) {
        my $gc = ($seqn =~ tr/G//);
        $gc += ($seqn =~ tr/C//);
        $gc = sprintf("%d",$gc*$bylength);
        $good = 0 if(exists $params{min_gc} && $gc < $params{min_gc});
        $good = 0 if(exists $params{max_gc} && $gc > $params{max_gc});
        $good = 0 if(exists $params{range_gc} && !&checkRange($params{range_gc},$gc));
    }

    #check for N's in sequence
    if($good && (exists $params{ns_max_p} || exists $params{ns_max_n})) {
        my $ns = ($seqn =~ tr/N//);
        $good = 0 if(exists $params{ns_max_p} && ($ns*$bylength) > $params{ns_max_p});
        $good = 0 if(exists $params{ns_max_n} && $ns > $params{ns_max_n});
    }

    #check for non IUPAC chars in sequence
    if($good && exists $params{noniupac}) {
        $good = 0 if($seqn =~ m/[^ACGTN]/);
    }

    #check for sequence complexity
    if($good && defined $complval) {
        my ($complCounts,$complNum,$rest,$steps,%vals);
        if($length <= $WINDOWSIZE) {
            $rest = $length;
            $steps = 0;
        } else {
            $steps = int(($length - $WINDOWSIZE) / $WINDOWSTEP) + 1;
            $rest = $length - $steps * $WINDOWSTEP;
            unless($rest > $WINDOWSTEP) {
                $rest += $WINDOWSTEP;
                $steps--;
            }
        }
        if($params{lc_method} eq 'dust') {
            foreach my $i (0..$steps) {
                ($complCounts,$complNum) = &count3mers(substr($seqn,($i * $WINDOWSTEP),($i == $steps ? $rest : $WINDOWSIZE)));
                $vals{&calcDust($complCounts,$complNum,$WINDOWSIZE)}++;
            }
            my $mean = 0;
            foreach my $val (keys %vals) {
                $mean += $val*$vals{$val};
            }
            $mean /= ($steps+1);
            $good = 0 if(int($mean * 100 / 31) > $complval);
        } elsif($params{lc_method} eq 'entropy') {
            foreach my $i (0..$steps) {
                ($complCounts,$complNum) = &count3mers(substr($seqn,($i * $WINDOWSTEP),($i == $steps ? $rest : $WINDOWSIZE)));
                $vals{&calcEntropy($complCounts,$complNum,$WINDOWSIZE)}++;
            }
            my $mean = 0;
            foreach my $val (keys %vals) {
                $mean += $val*$vals{$val};
            }
            $mean /= ($steps+1);
            $good = 0 if(int($mean * 100) < $complval);
        }
    }

    #check for read duplicates
    if($good && $derep) {
        push(@seqs,[$seqn,$goodcount++,$length]);

        #keep write data for possible duplicates
        if($params{out_format} == 1) { #FASTA
            push(@printtmp,[$sid,$header,$seq,$begin,$end,'']);
        } else { # FASTQ or FASTA+QUAL
            push(@printtmp,[$sid,$header,$seq,$begin,$end,$qual]);
        }

    } elsif($good && !$derep) { #passed filters
        $seqcount++;
        $seqbases += $length;
        return if($nogood);
        #check if change of sequence ID
        if(exists $params{seq_id}) {
            $sid = $params{seq_id}.$seqcount;
        }
        if(exists $params{rm_header}) {
            $header = undef;
        }
        #trim if necessary
        if($begin) {
            $seq = substr($seq,$begin);
            $qual = substr($qual,$begin) if(defined $qual && length($qual));
        }
        if($end) {
            $length = length($seq);
            $seq = substr($seq,0,$length-$end);
            $qual = substr($qual,0,$length-$end) if(defined $qual && length($qual));
        }
        #change case
        if(exists $params{seq_case}) {
            if($params{seq_case} eq 'lower') { #lower case
                $seq = lc($seq);
            } elsif($params{seq_case} eq 'upper') { #upper case
                $seq = uc($seq);
            }
        }
        #convert between DNA and RNA
        if(exists $params{dna_rna}) {
            if($params{dna_rna} eq 'dna') { #RNA to DNA
                $seq =~ tr/Uu/Tt/;
            } elsif($params{dna_rna} eq 'rna') { #DNA to RNA
                $seq =~ tr/Tt/Uu/;
            }
        }
        #set line length
        if($linelen && $params{out_format} != 3) {
            $seq =~ s/(.{$linelen})/$1\n/g;
            $seq =~ s/\n$//;
        }
        #write data
        if($params{out_format} == 1) { #FASTA
            if($stdoutgood) {
                print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
            } else {
                print $fhgood '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
            }
        } elsif($params{out_format} == 3) { # FASTQ
            &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
            if($stdoutgood) {
                print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '+'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $qual."\n";
            } else {
                print $fhgood '@'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
                print $fhgood '+'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood $qual."\n";
            }
        } elsif($params{out_format} == 2) { #FASTA+QUAL
            &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
            if($stdoutgood) {
                &printError('the output cannot be written to STDOUT for FASTA and QUAL file output. The option can only be used for FASTA only or FASTQ output');
            } else {
                print $fhgood '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
                print $fhgood2 '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
            }
        }
    } else { #filtered out
        $badcount++;
        $badbases += length($seq);
        return if($nobad);
        #set line length
        if($linelen && $params{out_format} != 3) {
            $seq =~ s/(.{$linelen})/$1\n/g;
            $seq =~ s/\n$//;
        }
        #write data
        if($params{out_format} == 1) { #FASTA
            if($stdoutbad) {
                print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
            } else {
                print $fhbad '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad $seq."\n";
            }
        } elsif($params{out_format} == 3) { # FASTQ
            if($stdoutbad) {
                &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '+'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $qual."\n";
            } else {
                &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                print $fhbad '@'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad $seq."\n";
                print $fhbad '+'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad $qual."\n";
            }

        } elsif($params{out_format} == 2) { #FASTA+QUAL
            if($stdoutbad) {
                &printError('the output cannot be written to STDOUT for FASTA and QUAL file output. The option can only be used for FASTA only or FASTQ output');
            } else {
                &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                print $fhbad '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad $seq."\n";
                print $fhbad2 '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
            }
        }
    }
}

#dereplicate sequences
sub derepSeqs {
    my $numseqs = scalar(@seqs);
    if($derep && $numseqs) {
        my ($sid,$seq,$qual,$header,$begin,$end);
        my ($dcounts,$dupls) = &checkForDupl(\@seqs,\%dereptypes,$numseqs);

        print STDERR "Write results to output file(s)\n" if(exists $params{verbose});
        #for progress bar
        my $progress = 0;
        my $counter = 1;
        my $part = int($numseqs/100);
        print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});

        foreach my $i (0..$numseqs-1) {
            #progress bar stuff
	    $counter++;
	    if($counter > $part) {
		$counter = 1;
		$progress++;
		$progress = 99 if($progress > 99);
		print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
	    }
            #get data
            $sid = $printtmp[$i]->[0];
            $seq = $printtmp[$i]->[2];
            $qual = $printtmp[$i]->[5];
            $header = $printtmp[$i]->[1];
            $begin = $printtmp[$i]->[3];
            $end = $printtmp[$i]->[4];
            #write data
            if(exists $dupls->{$i} || (exists $params{seq_num} && $params{seq_num} <= $seqcount)) { #bad
                $badcount++;
                $badbases += length($seq);
                next if($nobad);
                #set line length
                if($linelen && $params{out_format} != 3) {
                    $seq =~ s/(.{$linelen})/$1\n/g;
                    $seq =~ s/\n$//;
                }
                #write data
                if($params{out_format} == 1) { #FASTA
                    if($stdoutbad) {
                        print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                    } else {
                        print $fhbad '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad $seq."\n";
                    }
                } elsif($params{out_format} == 3) { # FASTQ
                    if($stdoutbad) {
                        &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                        print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                        print STDOUT '+'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $qual."\n";
                    } else {
                        &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                        print $fhbad '@'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad $seq."\n";
                        print $fhbad '+'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad $qual."\n";
                    }

                } elsif($params{out_format} == 2) { #FASTA+QUAL
                    if($stdoutbad) {
                        &printError('the output cannot be written to STDOUT for FASTA and QUAL file output. The option can only be used for FASTA only or FASTQ output');
                    } else {
                        &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                        print $fhbad '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad $seq."\n";
                        print $fhbad2 '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
                    }
                }
            } else { #good
                $seqcount++;
                $seqbases += (length($seq)-$begin-$end);
                next if($nogood);
                #check if change of sequence ID
                if(exists $params{seq_id}) {
                    $sid = $params{seq_id}.$seqcount;
                }
                if(exists $params{rm_header}) {
                    $header = undef;
                }
                #trim if necessary
                if($begin) {
                    $seq = substr($seq,$begin);
                    $qual = substr($qual,$begin) if(defined $qual && length($qual));
                }
                if($end) {
                    $length = length($seq);
                    $seq = substr($seq,0,$length-$end);
                    $qual = substr($qual,0,$length-$end) if(defined $qual && length($qual));
                }
                #change case
                if(exists $params{seq_case}) {
                    if($params{seq_case} eq 'lower') { #lower case
                        $seq = lc($seq);
                    } elsif($params{seq_case} eq 'upper') { #upper case
                        $seq = uc($seq);
                    }
                }
                #convert between DNA and RNA
                if(exists $params{dna_rna}) {
                    if($params{dna_rna} eq 'dna') { #RNA to DNA
                        $seq =~ tr/Uu/Tt/;
                    } elsif($params{dna_rna} eq 'rna') { #DNA to RNA
                        $seq =~ tr/Tt/Uu/;
                    }
                }
                #set line length
                if($linelen && $params{out_format} != 3) {
                    $seq =~ s/(.{$linelen})/$1\n/g;
                    $seq =~ s/\n$//;
                }
                #write data
                if($params{out_format} == 1) { #FASTA
                    if($stdoutgood) {
                        print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                    } else {
                        print $fhgood '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood $seq."\n";
                    }
                } elsif($params{out_format} == 3) { # FASTQ
                    &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                    if($stdoutgood) {
                        print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                        print STDOUT '+'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $qual."\n";
                    } else {
                        print $fhgood '@'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood $seq."\n";
                        print $fhgood '+'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood $qual."\n";
                    }
                } elsif($params{out_format} == 2) { #FASTA+QUAL
                    &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                    if($stdoutgood) {
                        &printError('the output cannot be written to STDOUT for FASTA and QUAL file output. The option can only be used for FASTA only or FASTQ output');
                    } else {
                        print $fhgood '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood $seq."\n";
                        print $fhgood2 '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
                    }
                }
            }
        }
        print STDERR "\r\tdone               \n" if(exists $params{verbose});
    }
}

#calculate summary statistics from sequences
sub calcSeqStats {
    my ($seq,$length,$stats,$kmers,$odds,$counts) = @_;

    #length related: min, max, range, mean, stddev, mode
    if(exists $params{stats_len}) {
        $counts->{length}->{$length}++;
    }

    #dinucleotide odds ratio related: aatt, acgt, agct, at, catg, ccgg, cg, gatc, gc, ta
    if(exists $params{stats_dinuc}) {
        &dinucOdds($seq,$length,$odds);
    }

    #tag related: probability of 5' and 3' tag sequence based on kmer counts
    if(exists $params{stats_tag}) {
        #get kmers
        if($length >= 5) {
            #get 5' and 3' ends
            my $str5 = substr($seq,0,5);
            my $str3 = substr($seq,$length-5);
            unless($str5 eq 'AAAAA' || $str5 eq 'TTTTT' || $str5 eq 'CCCCC' || $str5 eq 'GGGGG' || $str5 eq 'NNNNN') {
                $kmers->{5}->{$str5}++;
            }
            unless($str3 eq 'AAAAA' || $str3 eq 'TTTTT' || $str3 eq 'CCCCC' || $str3 eq 'GGGGG' || $str3 eq 'NNNNN') {
                $kmers->{3}->{$str3}++;
            }
        }
        #check for MID tags
        if($length >= $MIDCHECKLENGTH) {
            my $str5 = substr($seq,0,$MIDCHECKLENGTH);
            foreach my $mid (keys %MIDS) {
                if(index($str5,$mid) != -1) {
                    $MIDS{$mid}++;
                    last;
                }
            }
        }
    }

    #ambiguous base N related: seqswithn, maxp
    if(exists $params{stats_ns}) {
        my $bylength = 100/$length;
        my $ns = ($seq =~ tr/N//);
        $stats->{stats_ns}->{seqswithn}++ if($ns > 0);
        $stats->{stats_ns}->{maxn} = $ns if($ns > ($stats->{stats_ns}->{maxn}||0));
        $ns = ($ns > 0 && $ns*$bylength < 1 ? 1 : sprintf("%d",$ns*$bylength));
        $stats->{stats_ns}->{maxp} = $ns if($ns > ($stats->{stats_ns}->{maxp}||0));
    }
}

#dinucleotide odds ratio calculation
sub dinucOdds {
    my ($seq,$length,$odds) = @_;
    my ($mononum,$dinum,$i,$x,$y);
    $mononum = $dinum = 0;
    my %di = ('AA' => 0, 'AC' => 0, 'AG' => 0, 'AT' => 0, 'CA' => 0, 'CC' => 0, 'CG' => 0, 'CT' => 0, 'GA' => 0, 'GC' => 0, 'GG' => 0, 'GT' => 0, 'TA' => 0, 'TC' => 0, 'TG' => 0, 'TT' => 0);
    my %mono = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
    $i = 1;
    $x = substr($seq,0,1);
    while($i < $length) {
	$y = substr($seq,$i,1);
	if($y eq 'N') {
	    $i += 2;
	    if($i+1 < $length) {
		$x = substr($seq,$i+1,1);
	    }
	} elsif($x eq 'N') {
	    $i++;
	    $x = $y;
	} else {
	    $di{$x.$y}++;
	    $dinum++;
	    unless($i+1 < $length && substr($seq,$i+1,1) ne 'N') {
		$mono{$y}++;
		$mononum++;
	    }
	    $mono{$x}++;
	    $mononum++;
	    $i++;
	    $x = $y;
	}
    }
    if($dinum) {
	my $factor = 2 * $mononum * $mononum / $dinum;
	my $AT = $mono{'A'} + $mono{'T'};
	my $GC = $mono{'C'} + $mono{'G'};
	if($AT) {
	    my $AT2 = $factor / ($AT * $AT);
	    $odds->{'AATT'} += ($di{'AA'} + $di{'TT'}) * $AT2;
	    $odds->{'AT'}   +=          2 * $di{'AT'}  * $AT2;
	    $odds->{'TA'}   +=          2 * $di{'TA'}  * $AT2;
	    if($GC) {
		my $ATGC = $factor / ($AT * $GC);
		$odds->{'ACGT'} += ($di{'AC'} + $di{'GT'}) * $ATGC;
		$odds->{'AGCT'} += ($di{'AG'} + $di{'CT'}) * $ATGC;
		$odds->{'CATG'} += ($di{'CA'} + $di{'TG'}) * $ATGC;
		$odds->{'GATC'} += ($di{'GA'} + $di{'TC'}) * $ATGC;
		my $GC2 = $factor / ($GC * $GC);
		$odds->{'CCGG'} += ($di{'CC'} + $di{'GG'}) * $GC2;
		$odds->{'CG'}   +=          2 * $di{'CG'}  * $GC2;
		$odds->{'GC'}   +=          2 * $di{'GC'}  * $GC2;
	    }
	} elsif($GC) {
	    my $GC2 = $factor / ($GC * $GC);
	    $odds->{'CCGG'} += ($di{'CC'} + $di{'GG'}) * $GC2;
	    $odds->{'CG'}   +=          2 * $di{'CG'}  * $GC2;
	    $odds->{'GC'}   +=          2 * $di{'GC'}  * $GC2;
	}
    }
}

#calculate basic stats from an hash of number->count values
sub generateStats {
    my $counts = shift;
    my ($min,$max,$modeval,$mode,$mean,$count,$std,$x,$c,@vals,$num,$median,%stats);

    #min, max, mode and modeval
    $min = -1;
    $max = $modeval = $mean = $count = $std = $num = 0;
    while (($x, $c) = each(%$counts)) {
        if($min == -1) {
            $min = $x;
        } elsif($min > $x) {
            $min = $x;
        }
        if($max < $x) {
            $max = $x;
        }
        if($modeval < $c) {
            $modeval = $c;
            $mode = $x;
        }
        $mean += $x*$c;
        $count += $c;
        foreach(1..$c) {
            push(@vals,$x);
            $num++;
        }
    }

    #mean and stddev
    $mean /= $count;
    while (($x, $c) = each(%$counts)) {
        $std += $c*(($x-$mean)**2);
    }

    #median
    if($num == 1) {
        $median = $vals[0];
    } elsif($num == 2) {
        $median = ($vals[0]+$vals[1])/2;
    } else {
        @vals = sort {$a <=> $b} @vals;
        if($num % 2) {
            $median = $vals[($num-1)/2];
        } else {
            $median = ($vals[$num/2]+$vals[$num/2-1])/2;
        }
    }

    #save stats
    $stats{min} = $min;
    $stats{max} = $max;
    $stats{range} = $max-$min+1;
    $stats{modeval} = $modeval;
    $stats{mode} = $mode;
    $stats{mean} = sprintf("%.2f",$mean);
    $stats{stddev} = sprintf("%.2f",($std/$count)**(1/2));
    $stats{median} = $median;

    return \%stats;
}

sub count3mers {
    my $str = shift;
    my %counts = ();
    my $num = 0;
    foreach my $i (0..(length($str) - 3)) {
	$counts{substr($str,$i,3)}++;
	$num++;
    }
    return (\%counts,$num);
}

sub calcDust {
    my ($counts,$num,$w) = @_;
    if($num > 1) {
	my $score = 0;
	my $tmp;
	foreach(keys %$counts) {
	    $tmp = $counts->{$_};
	    $score += ($tmp * ($tmp - 1) / 2);
	}
        if($w > $num+2) {
	    return ($score / ($num-1)) * (($w - 2) / $num);
	} else {
	    return ($score / ($num-1));
	}
    } else {
	return 31; #to assign a maximum score based on the scaling factor 100/31
    }
}

sub calcEntropy {
    my ($counts,$num,$w) = @_;
    my $entropy = 0;
    my $tmp;
    foreach(keys %$counts) {
	$tmp = $counts->{$_};
	$entropy -= ($tmp / $num) * log($tmp / $num);
    }
    if($num <= 1) {
        return 0;
    } elsif($w > $num+2) {
	return ($entropy / log($num));
    } else {
	return ($entropy / $LOG62);
    }
}

#requires seqs array with [upper-case seq,array index, length] for each entry
sub checkForDupl {
    #requires seqs array with [upper-case seq, array index, length] for each entry
    my ($seqs,$types,$numseqs) = @_;
    my (@sort,$num,%dupls,$pretype,$precount,%counts);
    #precount = number duplicates for the same sequence

    print STDERR "Check for duplicates\n" if(exists $params{verbose});

    #for progress bar
    my $progress = 1;
    my $counter = 1;
    my $part = int($numseqs*4/100);
    print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});

    #exact duplicates and prefix duplicates
    if(exists $types->{0} || exists $types->{1} || exists $types->{2}) {
	$precount = 0;
	$pretype = -1;
        @sort = sort {$a->[0] cmp $b->[0]} @$seqs;
        foreach my $i (0..$numseqs-2) {
            if(exists $types->{0} && $sort[$i]->[2] == $sort[$i+1]->[2] && $sort[$i]->[0] eq $sort[$i+1]->[0]) {
                $dupls{$sort[$i]->[1]} = 0;
		if($pretype == 0) {
		    $precount++;
		} else {
		    if($pretype == 1 && $precount) {
			$counts{$precount}->{$pretype}++;
		    }
		    $pretype = 0;
		    $precount = 1;
		}
            } elsif(exists $types->{1} && $sort[$i]->[2] < $sort[$i+1]->[2] && $sort[$i]->[0] eq substr($sort[$i+1]->[0],0,$sort[$i]->[2])) {
                $dupls{$sort[$i]->[1]} = 1;
		if($pretype == 1) {
		    $precount++;
		} else {
		    if($pretype == 0 && $precount) {
			$counts{$precount}->{$pretype}++;
		    }
		    $pretype = 1;
		    $precount = 1;
		}
            } else {
		if($precount) {
		    $counts{$precount}->{$pretype}++;
		    $precount = 0;
		}
		$pretype = -1;
	    }
            $sort[$i] = undef;
	    #progress bar stuff
	    $counter++;
	    if($counter > $part) {
		$counter = 1;
		$progress++;
		$progress = 99 if($progress > 99);
		print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
	    }
        }
	if($precount) {
	    $counts{$precount}->{$pretype}++;
	}
    }
    #suffix duplicates
    if(exists $types->{2}) {
        $num = 0;
        @sort = ();
        foreach(@$seqs) {
            next if(exists $dupls{$_->[1]});
            push(@sort,[(scalar reverse $_->[0]),$_->[1],$_->[2]]);
            $num++;
        }
	if($num > 1) {
	    $precount = 0;
	    $pretype = -1;
	    @sort = sort {$a->[0] cmp $b->[0]} @sort;
	    foreach my $i (0..$num-2) {
		if($sort[$i]->[2] < $sort[$i+1]->[2] && $sort[$i]->[0] eq substr($sort[$i+1]->[0],0,$sort[$i]->[2])) {
		    $dupls{$sort[$i]->[1]} = 2;
		    if($pretype == 2) {
			$precount++;
		    } else {
			$pretype = 2;
			$precount = 1;
		    }
		} else {
		    if($precount) {
			$counts{$precount}->{$pretype}++;
			$precount = 0;
		    }
		    $pretype = -1;
		}
		$sort[$i] = undef;
		#progress bar stuff
		$counter++;
		if($counter > $part) {
		    $counter = 1;
		    $progress++;
		    $progress = 99 if($progress > 99);
		    print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
		}
	    }
	    if($precount) {
		$counts{$precount}->{$pretype}++;
	    }
	}
    }
    #reverse complement exact and prefix/suffix duplicates
    if(exists $types->{3} || exists $types->{4}) {
        $num = 0;
        @sort = ();
        foreach(@$seqs) {
            if(exists $dupls{$_->[1]}) {
		$counter++;
		next;
	    }
            push(@sort,[$_->[0],$_->[1],$_->[2],0]);
            push(@sort,[&revcompuc($_->[0]),$_->[1],$_->[2],1]);
            $num += 2;
        }
	if($num > 1) {
	    $precount = 0;
	    $pretype = -1;
	    @sort = sort {$a->[0] cmp $b->[0]} @sort;
	    foreach my $i (0..$num-2) {
		unless($sort[$i]->[3] == $sort[$i+1]->[3] || $sort[$i]->[1] eq $sort[$i+1]->[1] || exists $dupls{$sort[$i]->[1]}) { #don't check if both same (original or revcomp) or already counted as dubs
		    if(exists $types->{3} && $sort[$i]->[2] == $sort[$i+1]->[2] && $sort[$i]->[0] eq $sort[$i+1]->[0]) {
			$dupls{$sort[$i]->[1]} = 3;
			if($pretype == 3) {
			    $precount++;
			} else {
			    if($pretype == 4 && $precount) {
				$counts{$precount}->{$pretype}++;
			    }
			    $pretype = 3;
			    $precount = 1;
			}
		    } elsif(exists $types->{4} && $sort[$i]->[2] < $sort[$i+1]->[2] && $sort[$i]->[0] eq substr($sort[$i+1]->[0],0,$sort[$i]->[2])) {
			$dupls{$sort[$i]->[1]} = 4;
			if($pretype == 4) {
			    $precount++;
			} else {
			    if($pretype == 3 && $precount) {
				$counts{$precount}->{$pretype}++;
			    }
			    $pretype = 4;
			    $precount = 1;
			}
		    } else {
			if($precount) {
			    $counts{$precount}->{$pretype}++;
			    $precount = 0;
			}
			$pretype = -1;
		    }
		}
		$sort[$i] = undef;
		#progress bar stuff
		$counter++;
		if($counter > $part) {
		    $counter = 1;
		    $progress++;
		    $progress = 99 if($progress > 99);
		    print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
		}
	    }
	    if($precount) {
		$counts{$precount}->{$pretype}++;
	    }
	}
    }
    print STDERR "\r\tdone               \n" if(exists $params{verbose});
    return (\%counts,\%dupls);
}

#get the frequency of possible tags by shifting kmers by max 2 positions when aligned
sub getTagFrequency {
    my ($kmers) = @_;

    #find most abundant kmer counts
    my $percentone = $numseqs/100;
    my $percentten = $numseqs/10;
    my %most;
    foreach my $sp (keys %$kmers) {
	$most{$sp}->{max} = 0;
	foreach(keys %{$kmers->{$sp}}) {
#	    next if($_ eq 'A'x5 || $_ eq 'T'x5 || $_ eq 'C'x5 || $_ eq 'G'x5 || $_ eq 'N'x5);
	    if($kmers->{$sp}->{$_} >= $percentten) {
		$most{$sp}->{ten}++;
	    } elsif($kmers->{$sp}->{$_} >= $percentone) {
		$most{$sp}->{one}++;
	    }
	    #get max count
	    $most{$sp}->{max} = $kmers->{$sp}->{$_} if($most{$sp}->{max} < $kmers->{$sp}->{$_});
	}
    }

    #filter kmers by frequency - threshold of >10% occurrence -> max of 9 different kmers or more if there is non with >10% occurrence
    my $numseqssub = $numseqs/10;
    my $onecount = 2;
    foreach my $sp (keys %$kmers) {
	foreach(keys %{$kmers->{$sp}}) {
	    if(exists $most{$sp}->{ten} && $most{$sp}->{ten} > 0) {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} < $percentten);
	    } elsif(exists $most{$sp}->{one} && $most{$sp}->{one} > 0) {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} < $percentone);
	    } else {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} != $most{$sp}->{max});
	    }
	}
    }

    my (%kmersum,%kmershift);
    foreach my $sp (sort {$b <=> $a} keys %$kmers) { #5' before 3'
	#if more than one kmer in array, test if shifted by max 2 positions
	my $numkmer = scalar(keys %{$kmers->{$sp}});

	if($numkmer > 1) {
	    my @matrix;
	    my @kmersort = sort {$kmers->{$sp}->{$b} <=> $kmers->{$sp}->{$a}} keys %{$kmers->{$sp}};
	    foreach my $i (0..($numkmer-2)) {
		foreach my $j (($i+1)..($numkmer-1)) {
		    $matrix[$i]->[$j-($i+1)] = &align2seqs($kmersort[$j],$kmersort[$i]);
		}
	    }
	    my $countgood = 0;
	    foreach my $i (0..($numkmer-2)) {
		if(!defined @{$matrix[0]->[$i]}) { #not matching
		    my $count = 0;
		    foreach my $j (1..($numkmer-2)) {
			$count++;
			last if(defined $matrix[$j]->[$i-$j]); #found shift using other kmers
		    }
		    if($count < ($numkmer-1) && $i > 0) {
			my $sum = 0;
			my $sign;
			foreach my $j (0..$count) {
			    next unless(defined $matrix[$j] && defined $matrix[$j]->[$i-1]); #fix: 08/2010
			    if(defined $sign) {
				if(($sign < 0 && (defined $matrix[$j]->[$i-1]->[0] && $matrix[$j]->[$i-1]->[0] < 0)) || ($sign > 0 && (defined $matrix[$j]->[$i-1]->[0] && $matrix[$j]->[$i-1]->[0] > 0))) {
				    $sum += $matrix[$j]->[$i-1]->[0];
				} elsif(($sign < 0 && (defined $matrix[$j]->[$i-1]->[1] && $matrix[$j]->[$i-1]->[1] < 0)) || $sign > 0 && (defined $matrix[$j]->[$i-1]->[1] && $matrix[$j]->[$i-1]->[1] > 0)) {
				    $sum += $matrix[$j]->[$i-1]->[1];
				}
			    } elsif(defined $matrix[$j]->[$i-1]->[0]) {
				$sum += $matrix[$j]->[$i-1]->[0];
			    }
			    $sign = ((defined $matrix[$j]->[$i-1]->[0] && $matrix[$j]->[$i-1]->[0] < 0) ? -1 : 1);
			}
			$matrix[0]->[$i] = [$sum] if(defined $sign); #fix: 08/2010
		    }
		}
		if(!defined @{$matrix[0]->[$i]}) {
		    last;
		} else {
		    $countgood++;
		}
	    }
	    if($countgood) {
		my $min;
		if($sp == 3) { #3' prime end, 5 for 5' end
		    #find maximum shift to right (pos value)
		    $min = -100;
		    foreach my $i (0..($countgood-1)) {
			$min = ((defined $matrix[0]->[$i]->[0] && $min > $matrix[0]->[$i]->[0]) ? $min : $matrix[0]->[$i]->[0]);
		    }
		    if($min > 0) {
			$min = -$min;
		    } else {
			$min = 0;
		    }
		} else {
		    #find maximum shift to left (neg value)
		    $min = 100;
		    foreach my $i (0..($countgood-1)) {
			$min = ((defined $matrix[0]->[$i]->[0] && $min < $matrix[0]->[$i]->[0]) ? $min : $matrix[0]->[$i]->[0]);
		    }
		    if($min < 0) {
			$min = abs($min);
		    } else {
			$min = 0;
		    }
		}
#		$kmershift{$sp}->{$kmersort[0]} = $min;
		$kmersum{$sp} += $kmers->{$sp}->{$kmersort[0]};
		foreach my $i (0..($countgood-1)) {
#		    $kmershift{$sp}->{$kmersort[$i+1]} = $matrix[0]->[$i]->[0]+$min;
		    $kmersum{$sp} += $kmers->{$sp}->{$kmersort[$i+1]};
		}
	    } else {
		my $tmp = (sort {$kmers->{$sp}->{$b} <=> $kmers->{$sp}->{$a}} keys %{$kmers->{$sp}})[0];
#		$kmershift{$sp}->{$tmp} = 0;
		$kmersum{$sp} += $kmers->{$sp}->{$tmp};
	    }
	} elsif($numkmer == 1) {
	    my $tmp = (keys %{$kmers->{$sp}})[0];
#	    $kmershift{$sp}->{$tmp} = 0;
	    $kmersum{$sp} += $kmers->{$sp}->{$tmp};
	}
    }

   return \%kmersum;
}

sub align2seqs {
    my ($seq1,$seq2) = @_;
    my @shift;

    #get number of shifted positions
    if(substr($seq1,0,4) eq substr($seq2,1,4)) { #shift right by 1
	push(@shift,1);
    } elsif(substr($seq1,0,3) eq substr($seq2,2,3)) { #shift right by 2
	push(@shift,2);
    }
    if(substr($seq1,1,4) eq substr($seq2,0,4)) { #shift left by 1
	push(@shift,-1);
    } elsif(substr($seq1,2,3) eq substr($seq2,0,3)) { #shift left by 2
	push(@shift,-2);
    }

    return \@shift;
}

sub revcompuc {
    my $seq = shift;
    $seq = scalar reverse $seq;
    $seq =~ tr/GATC/CTAG/;
    return $seq;
}

sub compuc {
    my $seq = shift;
    $seq =~ tr/GATC/CTAG/;
    return $seq;
}
