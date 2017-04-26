#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Apr 2017. Uppsala University ##################### #
#
# Splits Pfam database (file Pfam-A.hmm) into individual hmm profiles
#
##########################################################################################

my $usage = "split_pfam.pl\n";
my $file="Pfam-A.hmm";

my $outfile;

open (IN, "<", $file) or die "Can't find file $file!\n";

while ( my $line =<IN> ) {

    chomp $line;

    # skip first line
    next if ( $line =~ /^HMMER3.+/ );

    # get name
    if ( $line =~ /^NAME  (.+)/ ) {

#	my @lines = split "\t", $line;
	my $name = $1;
	$outfile = $name . ".hmm";

	# open output file
	open (OUT, ">", $outfile) or die "Can't print to file $outfile!\n";

	# print out first two lines
	print OUT "HMMER3/f [3.1b2 | February 2015]\n";
	print OUT "$line\n";
	next;
    }
    # close output file 
    elsif ( $line =~ /^\/\// ) {

	print OUT "$line\n";
	close(OUT);
	next;
    }
    # all other lines, print to file
    else {
	print OUT "$line\n";
    }
}

