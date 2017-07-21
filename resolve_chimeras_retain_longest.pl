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
# Given an alignment and a list of sequences, it chooses the longest and renames it
#	as *_Genus_sp_N/A@
#
# This aims to resolve chimeras in Burki and Baurain
#
# It requires unaligned sequences
#
##########################################################################################


my $usage = "filter_fasta_published_taxa.pl infile.fa query_file > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;
my $query = $ARGV[1] or die $usage;

print STDERR "FILE: $fasta\n";

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
                	    		'-format' => "fasta");
                		        
my %queries;
my %saved;
my $chimera_epithet; # of final chimeric sequence
my $chimera;	# of final chimeric sequence

my $longest_taxa;
my $longest_seq;

open (IN , $query) or die "Can't open $query, $!\n";

while (my $line = <IN>){
	chomp $line;
	
	# Get epithet (in principle all taxa for chimera should have the same genus)
	# If not, then hard-code epithet above and comment out these lines
	my @lines = split "_", $line;
	# Genus and species should be fields 6 and 7
	$chimera_epithet = $lines[5] . "_sp"; # I will loose information after the @ sign for chimeras
	$chimera = $lines[0] . "_" . $lines[1] . "_" . $lines[2] . "_" . $lines[3] . "_" .
				 $lines[4] . "_" . $chimera_epithet . "_N/A";
	
	$queries{$line} = "1";	# this will save the sequence names to be chimerized
}



while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;
    my ( $species, $dadaset ) = split "\@", $seqname;
    my $length = length $seq_obj->seq;
    
    # save sequences not in query file
    if ( !exists ( $queries{$species} ) ) {
    
    	$saved{$seqname} = $seq_obj->seq;
    	next;

 	}
	# here dealing already only with species in query file ( chimeras )
	# use epithet to save chimeric sequence in %saved
	if ( !exists $saved{$chimera} ) {
	
		print STDERR "\tSaved: $length $seqname\n";
		$saved{$chimera} = $seq_obj->seq;
		next;	
    }
    # if sequence is saved, then replace if length is longer
    else {
    
    	# obtain length of saved record
    	my $length_saved = length $saved{$chimera};
    	
    	# and replace record if new sequence is longer
    	if ( $length > $length_saved ) {
    	
			print STDERR "\tReplaced by: $length $seqname\n";
    		$saved{$chimera} = $seq_obj->seq;
    	}
    	else {
			print STDERR "\tSkipped: $length $seqname\n";  		
    	
    	}
    }
}

# print out sequences
foreach my $key ( sort keys %saved ) {

	print ">$key@\n";
	print "$saved{$key}\n";
}
