#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2017. Uppsala University ##################### #
# 
# modification of extract_homo_for_annotation.pl to extract instead the longest sequence
# 
# Script to extract reference taxa for annotation of gene alignments
# Typically be used into a for loop
# 	for f in *.fa; do perl extract_homo_for_annotation.pl $f ; done
#
# Taxon names are hard-coded into arrays. They contain the preferred taxon in position 1
# followed by other taxa that can replace 1 if not present for that particular gene
#
# Assumes gene files in fasta format
#
##########################################################################################

my $usage = "extract_longest_for_annotation.pl infile > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;

# define header for target taxa (human)
# and alternative headers in order of preference

# GET SEQUENCE

# read fasta file with SeqIO and store in %hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
	                	        '-format' => "fasta");

my %infile = ();
my $name_longest = ();
my $seq_longest = 0;
my $length = 0;
                		        
while (my $seq_obj = $seqio_obj->next_seq){

	# store sequences into %hash
    my $seqname = $seq_obj->primary_id;
	my $sequenc = $seq_obj->seq;
	# remove dashes if aligned
	$sequenc =~ s/-//g;
	my $new_length = length $sequenc;
	
	# if longer sequence is found, then reasign
	if ( $new_length > $length ) {
	
		$length = $new_length;
		$name_longest = $seqname;
		$seq_longest = $sequenc;
	}
}

# print out
if ( $length != 0 ) {

	    print ">$fasta\@$name_longest\n";
    	print "$seq_longest\n";
}
else {
	print STDERR "File $fasta: reference sequence not found!\n";
}

