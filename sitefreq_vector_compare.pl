#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Math::Round;

##########################################################################################
#
# #################### Iker Irisarri. Jun 2017. Uppsala University ##################### #
#
# Compares vectors of AA frequencies (one per alignment colunn) 
# 	These can be obtained from sitefreq_from_fasta.pl or generated iqtree etc.
# 	It assumes iqtree format (sensu sitefreq_from_fasta.pl)
# 
##########################################################################################

my $usage = "sitefreq_vector_compare.pl input1 input2 > stdout\n";
my $input1 = $ARGV[0] or die $usage;
my $input2 = $ARGV[1] or die $usage;

# store info in hash
my $infile1_ref = open_and_read($input1);
my $infile2_ref = open_and_read($input2);

my %infile1 = %{ $infile1_ref };
my %infile2 = %{ $infile2_ref };

#print Dumper $infile1_ref;

# get number of keys (positions) and ensure they are the same in both cases
my $sites1 = scalar keys %{ $infile1_ref };
my $sites2 = scalar keys %{ $infile2_ref };

if ( $sites1 != $sites2 ) {
	die "Input files do not contain the same number of sites\n";
}

# compare hashes
for (my $i=1; $i < $sites1+1; $i++) { # $i refers to lines in input files (1-based)

	my $sum_diff_sq_sqrt = 0;

	# de-reference arrays for that position
	my @array1 = @{ $infile1{$i} };
	my @array2 = @{ $infile2{$i} };

	# compare elements of arrays
	for (my $j=0; $j < scalar @array1; $j++) { # $j refers to positions in arrays (0-based)
	
		# calculate sqrt of square differences
		my $diff = $array1[$j] - $array2[$j];
		my $diff_sq = $diff **2;					
		my $diff_sq_sqrt = sqrt $diff_sq;
		
		# sum up differences
		$sum_diff_sq_sqrt += $diff_sq_sqrt;
		# checking everything is calculated correctly
		#print "a:$array1[$j] b:$array2[$j] d:$diff sq:$diff_sq sq_sqrt:$diff_sq_sqrt sum:$sum_diff_sq_sqrt\n";
	}
	
	print "$i $sum_diff_sq_sqrt\n";
}


sub open_and_read {
	my $infile = shift;
	
	open (IN, "<", $infile) or die "Can't open file $infile!\n";

	my %hash;

	while ( my $line =<IN> ) {

		chop $line;
		my @lines = split (/ /, "$line");
		my $num = shift @lines;
		$hash{$num} = \@lines;
	}
	return \%hash;
}

		