#!/usr/bin/env perl
#
# script to download and package up a standalone version of amphora2
#
use strict;
use warnings;
`rm -rf Amphora-2`;
#`rm -rf bioperl-live`;
`git clone git://github.com/gjospin/Amphora-2.git`;
`rm -rf Amphora-2/.git`;
`git clone git://github.com/bioperl/bioperl-live.git`;
`mv bioperl-live/Bio* Amphora-2/lib`;
`wget http://search.cpan.org/CPAN/authors/id/G/GR/GROMMEL/Math-Random-0.71.tar.gz`;
`tar xvzf Math-Random-0.71.tar.gz`;
chdir("Math-Random-0.71");
`perl Makefile.PL`;
`make`;
`mv blib/arch/auto ../Amphora-2/lib/`;
`mv blib/lib/Math ../Amphora-2/lib/`;
chdir("..");

`wget http://search.cpan.org/CPAN/authors/id/R/RV/RVOSA/Bio-Phylo-0.45.tar.gz`;
`tar xvzf Bio-Phylo-0.45.tar.gz`;
chdir("Bio-Phylo-0.45");
`perl Makefile.PL`;
`make`;
`mv blib/lib/Bio/Phylo* ../Amphora-2/lib/Bio/`;
chdir("..");

`wget http://search.cpan.org/CPAN/authors/id/J/JE/JESSE/Locale-Maketext-Simple-0.21.tar.gz`;
`tar xvzf Locale-Maketext-Simple-0.21.tar.gz`;
chdir("Locale-Maketext-Simple-0.21");
`perl Makefile.PL`;
`make`;
`mv lib/Locale/ ../Amphora-2/lib/`;
chdir("..");

`wget http://search.cpan.org/CPAN/authors/id/J/JP/JPEACOCK/version-0.95.tar.gz`;
`tar xvzf version-0.95.tar.gz`;
chdir("version-0.95");
`perl Makefile.PL --perl_only`;
`make`;
`mv blib/lib/version* ../Amphora-2/lib/`;
chdir("..");

my @timerval = localtime();
my $datestr = (1900+$timerval[5]);
$datestr .= 0 if $timerval[4] < 9; 
$datestr .= ($timerval[4]+1);
$datestr .= 0 if $timerval[3] < 9; 
$datestr .= $timerval[3];

`rm -rf version-0.95*`;
`rm -rf Locale-Maketext-Simple-0.21*`;
`rm -rf Bio-Phylo-0.45*`;
`rm -rf Math-Random-0.71*`;
`rm -rf bioperl-live`;

`tar cvzf amphora2-$datestr.tar.gz Amphora-2/`
