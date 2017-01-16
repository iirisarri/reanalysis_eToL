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

my $usage = "extract_homo_for_annotation.pl infile > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;

# define header for target taxa (human)
# and alternative headers in order of preference

# Katz2015
my @reference = qw (Op_me_hsap Op_me_ptro Op_me_mmus Op_me_rnor Op_me_clup Op_me_ecab);
# Burki2016
#my @reference = qw (Homo_sapiens Danio_rerio Branchiostoma_floridae Lottia_gigantea Daphnia_pulex Nematostella_vectensis);
# Brown2013
#my @reference = qw (HomoSapi Nematostel Ampqueen);
#Baurain2010
#my @reference = qw (Hydra_magnipapillata Nematostella_vectensis Reniera_sp. Monosiga_brevicollis Monosiga_ovata);

# GET SEQUENCE

# read fasta file with SeqIO and store in %hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
	                	        '-format' => "fasta");

my %infile = ();
                		        
while (my $seq_obj = $seqio_obj->next_seq){

	# store sequences into %hash
    my $seqname = $seq_obj->primary_id;
	my $sequenc = $seq_obj->seq;
	# remove dashes if aligned
	$sequenc =~ s/-//g;
	
	$infile{$seqname} = $sequenc;
}

# search for references in order of preference

my $success = 0;

foreach my $tax ( @reference ) {

	if ( exists ( $infile{$tax} ) ) {
	    
	    print ">", $fasta, "@", $tax, "\n";
    	print $infile{$tax}, "\n";
    	$success++;
    	last;
    }
	# go to next item in @reference when array is not empty
    elsif ( scalar @reference !=0 && !exists ( $infile{$tax} ) ) {
    
    	next;
    }
}

# when none of the taxa in @reference are present
if ( $success == 0 ) {	
    	print STDERR "File $fasta: sequence for queried taxa do not exist\n";
	    print ">", $fasta, "_MISSING\n????????????????????????\n";
}

