#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. May 2017. Uppsala University ##################### #
#
# Reformats site freq file (Ding's & Simon's format) to format required by iqtree (-fs)
# It calculates e^^(-x) for each frequency and print out also the line number
#
##########################################################################################

my $usage = "sitefreq2iqtree.pl\n";
my $infile = $ARGV[0] or die "Can't open file: $!\n";

my $line_count = 0;

open (IN, "<", $infile) or die "Can't find file $infile!\n";

while ( my $line =<IN> ) {

    chomp $line;
    $line_count++;
    
    my @lines = split / / , $line;

    # transform to absolute frequencies (undo lnX) and print first line
    print $line_count;
    
    foreach my $freq ( @lines ) {

	my $inv = (-1) * $freq;
	my $exp_freq = exp($inv);
	# round to 6 decimals
	my $exp_freq_round = sprintf("%.6f", $exp_freq);
	
	print " $exp_freq_round";
    }
    print "\n";
}
	
