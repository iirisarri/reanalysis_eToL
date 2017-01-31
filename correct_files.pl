#!/usr/bin/perl

use strict;
use warnings;


# Small script to correct gene alignments from Katz and Grant 2015
# Some gene files included "n\t" after "Ba_pr_Bpen" and this scripts removes those

my $usage = "correct_files.pl infile.fa > STDOUT\n";
my $infile = $ARGV[0] or die $usage; 

open(IN, "<", $infile) or die;

while ( my $line =<IN> ) {

    chomp $line;
    if ( $line =~ /n\t(.+)/ ) {

	print "$1\n";
    }
    else {
	print "$line\n";
    }
}
