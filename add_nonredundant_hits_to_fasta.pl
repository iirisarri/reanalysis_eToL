#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jul 2017. Uppsala University ##################### #
# 
# Script to add sequences to a reference fasta from a new fasta, only in cases where the
# 	species is not already present
#
# It allows adding multiple sequences for each species, only for the new added sequences
#	(e.g., multiple blast hits that can be then checked in gene trees)
#
# It will append sequence length to new sequences
#
# NOTE: this script also removed gaps from sequences
#
##########################################################################################

my $usage = "add_nonredundant_hits_to_fasta.pl reference.fa new.fa > STDOUT\n";
my $ref_fasta = $ARGV[0] or die $usage;
my $new_fasta = $ARGV[1] or die $usage;

# read fasta file with SeqIO and store in %infile
my $seqio_obj = Bio::SeqIO->new('-file' => "<$ref_fasta",
	                	        '-format' => "fasta");

my %infile = ();
my %taxa = ();

# save all sequences in reference fasta file 
while (my $seq_obj = $seqio_obj->next_seq){

	# store sequences into %hash
    my $seqname = $seq_obj->primary_id;
    my ( $taxon, $dataset) = split "\@", $seqname;
	my $seq = $seq_obj->seq;
	
	# remove gaps
	$seq =~ s/-//g;

	$infile{$seqname} = $seq;
	$taxa{$taxon} = 1;
}

# read in fasta with additional sequences and save only if taxa not present
my $seqio_obj_new = Bio::SeqIO->new('-file' => "<$new_fasta",
	                	        	'-format' => "fasta");

# save all sequences in reference fasta file 
while (my $seq_obj_new = $seqio_obj_new->next_seq){

	# store sequences into %hash
    my $seqname_new = $seq_obj_new->primary_id;
    my ( $taxon_new, $dataset_new) = split "\@", $seqname_new;
	my $seq_new = $seq_obj_new->seq;

	# remove gaps
	$seq_new =~ s/-//g;
	
	if ( !exists $taxa{$taxon_new} ) {
	
		my $length = length $seq_new;
		$seqname_new .= "LENGTH:$length";
		$infile{$seqname_new} = $seq_new;

	}
}

# print out saved records ( containing longest sequence in case of duplicates)
foreach my $key ( sort keys  %infile ) {

	print ">$key\n";
    print $infile{$key}, "\n";
}
