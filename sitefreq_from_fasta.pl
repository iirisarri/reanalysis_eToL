#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Math::Round;

##########################################################################################
#
# #################### Iker Irisarri. Jun 2017. Uppsala University ##################### #
#
# Calculates AA frequencies per site given a fasta alignment
# Modif of calculate_aa_comp.pl to print out frequencies per site
# Gaps are ignored when calculating frequencies
# 
##########################################################################################

my $usage = "sitefreq_from_fasta.pl input.fa\n";
my $input = $ARGV[0] or die $usage;

my %fasta_hash;
my $aln_length = 0;
my $num_seqs = 0;

my %site_pattern = ();
my @order = qw (A R N D C Q E G H I L K M F P S T W Y V);

use Bio::SeqIO;

my $seqio_obj = Bio::SeqIO->new('-file' => "<$input",
								'-format' => "fasta");

# store sequences into hash
while (my $seq_obj = $seqio_obj->next_seq){
    my $seqname = $seq_obj->primary_id;
    my $seq =  $seq_obj->seq;
    $aln_length = length $seq;
    my @seq_array = split //, $seq;
    $fasta_hash{$seqname} = [@seq_array];
}

# loop through each alignment position

for (my $i=0; $i<$aln_length; $i++) {

	my $site = $i + 1;
	%site_pattern = (); # initialize with new site
	my $unambig_aa = 0;	# initialize with new site
	my $freq_sum = 0;

	# loop through the hash to get site pattern (count amino acids)
	foreach my $key ( sort keys %fasta_hash ) {

		my $aa = ${ $fasta_hash{$key} }[$i];
		
		# skip gaps
		next if $aa eq "-";

		# save pattern into %site_pattern
		if ( !exists $site_pattern{$aa} ) {
		
			$site_pattern{$aa} = 1;
			$unambig_aa++;
		}
		else {
		
			$site_pattern{$aa}++;
			$unambig_aa++;			
		}
	}
	#calculate frequencies and store in same hash
	foreach my $e ( keys %site_pattern ) {

		my $freq = $site_pattern{$e} / $unambig_aa;
		
		$site_pattern{$e} = $freq;
		$freq_sum += $freq;
	}

	# print out frequencies in "iqtree" format: 1 A R N D C Q E G H I L K M F P S T W Y V
	print "$site";
	
	foreach my $AA ( @order ) {
	
		if ( exists $site_pattern{$AA} ) {
			print " $site_pattern{$AA}";
		}
		else {
			print " 0";
		}
	}
	# checking that freqs within each line add up to 1 
	# print " SUM:$freq_sum\n";
}

print STDERR "\ndone!\n\n";

		