#!/usr/bin/perl

use warnings;
use strict;

##########################################################################################
#
# #################### Iker Irisarri. Feb 2017. Uppsala University ##################### #
#
# modified from standarize_taxids_in_fasta.pl
#
# Simple script to swap two taxon names in fasta (hard-coded in %swaps)
#
#
##########################################################################################


my $usage = "standarize_taxids_in_fasta.pl infile.fa > STDOUT\n";
my $infile = $ARGV[0] or die $usage; 


# Sr_st_Bhom should be replaced by Op_fu_Beme & viceversa
my %swaps = (">Sr_st_Bhom" => ">Op_fu_Beme",
			">Op_fu_Beme" => ">Sr_st_Bhom");

# read-in fasta and replace 
open (IN, "<", $infile) or die "Cannot open file $infile\n";

while ( my $line =<IN> ) {
	chomp $line;

	if ( exists $swaps{$line} ) {
		
		print "$swaps{$line}\n";
	}
	# print sequence lines && rest of headers
	else {

		print "$line\n";
	} 
}
close(IN);

print STDERR "\ndone!\n\n";

