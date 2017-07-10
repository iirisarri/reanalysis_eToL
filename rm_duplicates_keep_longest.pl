#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jul 2017. Uppsala University ##################### #
# 
# Script to remove sequences belonging to the same species (they have same header) but 
# 	keeping the longest sequence
#
# It will also remove gaps if present in the file (necessary to calculate lengths)
#
##########################################################################################

my $usage = "rm_duplicates_keep_longest.pl infile > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;

# read fasta file with SeqIO and store in %hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
	                	        '-format' => "fasta");

my %infile = ();

print STDERR "FILE: $fasta\n";
                		        
while (my $seq_obj = $seqio_obj->next_seq){

	# store sequences into %hash
    my $seqname = $seq_obj->primary_id;
	my $sequenc = $seq_obj->seq;
	
	# remove gaps
	$sequenc =~ s/-//g;

	# if sequence is not present, save it
	if ( !exists $infile{$seqname} ) {

			$infile{$seqname} = $sequenc;
	}
	# if present, check if length is longer than saved record
	else {
	
		my $length_current = length $sequenc;
		my $length_retrieved = length $infile{$seqname};
		print STDERR "duplicated record: $seqname\n";
		print STDERR "\tlengths current: $length_current previous: $length_retrieved\n";
		
		# if longer, overwrite with new record
		if ( $length_current > $length_retrieved ) {
		
			$infile{$seqname} = $sequenc;
			print STDERR "\trecord overwritten!\n";
		}
	}
}


# print out saved records ( containing longest sequence in case of duplicates)
foreach my $key ( sort keys  %infile ) {

	print ">$key\n";
    print $infile{$key}, "\n";
}
