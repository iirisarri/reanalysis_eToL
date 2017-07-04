#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jun 2017. Uppsala University ##################### #
# 
# Removes duplicated taxa from alignment according to a table of taxids provided as input
#
# This script is suited for this particular case and more general usage will need a a few
#	changes. This is because fasta headers are very irregular and I needed to treat
#	independently a few cases. Note that I only added fixes for the taxa that might be
#	affected by duplications in the current case! 
#
##########################################################################################

my $usage = "rm_redundant_taxa_table_fasta.pl infile_table fasta > STDOUT\n";
my $infile = $ARGV[0] or die $usage;
my $fasta = $ARGV[1] or die $usage;

# GET TAXA CORRESPONDENCE
my %taxa_correspondence = ();

open(IN, "<", $infile) or die "Can't open file $infile!\n";

while ( my $line =<IN> ) {

	chomp $line;
	my @lines = split "\t", $line;
	
	# save taxa in hash with first (preferred) element as key and the rest in an array
	my $first = shift @lines;	
	$taxa_correspondence{$first} = \@lines;
}

# GET SEQUENCE
# read fasta file with SeqIO and store in %hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
	                	        '-format' => "fasta");

my %fasta = ();	# contains native fasta
my %taxa = ();	# contains clean taxids for further comparison
                		        
while (my $seq_obj = $seqio_obj->next_seq){

	# store sequences into %fasta
    my $seqname = $seq_obj->primary_id;
	my $sequenc = $seq_obj->seq;
	$fasta{$seqname} = $sequenc;

	# get clean names and store in %taxa
	my @seqnames = split "\@", $seqname;
	my $clean = $seqnames[0];

	#print "\n$clean\n\n";
	
	# simplify Burki's headers in four cases (only for species with putative 
	# redundant sequences, because only these will be used in the next steps)
	if ( $clean =~ /^Dinophyta_Dinophyceae-(.+)/ || $clean =~ /^Cryptophyta_Cryptophyceae-(.+)/ ||
		 $clean =~ /^Ochrophyta_Xanthophyceae-(.+)/ || $clean =~ /^Rhodophyta_Rhodellophyceae-(.+)/ ) {
	
		$clean = $1;
		#print "\n$clean\n\n";
	}
	# exclude possible "gene info" e.g., >Homo_sapiens-tubb@ENSP00000341289	
	if ( $clean =~/([\w_]+)-.*/ ) {
	
		$clean = $1;
		#print "\n$clean\n\n";
	}

	#print "\n$clean\n\n";
	$taxa{$clean} = "1";

}

# loop through the %taxa_correspondence to find taxa that might be duplicated

print STDERR "FILE: $fasta\n"; # to help tracking results when run in a loop

foreach my $key ( keys %taxa_correspondence ) {

	#print "\n$key\n\n";
 
	# if this is found in the alignment...
	if ( exists $taxa{$key} ) {
	
		# ...get all "secondary" taxids...
		my @inner = @{ $taxa_correspondence{$key} };
		
		foreach my $e ( @inner ) {
		
			#print "\n$e\n\n";
			# ALL HEADERS THAT ARE NOT BURKI'S SHOULD HAVE THE 
			# EXACT SAME NAME AS STORED IN %taxa_correspondence
			if ( exists $fasta{$e} ) {
			
				delete $fasta{$e};
				print STDERR "\tremoved $e\n";
			}
		}
	}
}

# print out fasta after removing duplicates
foreach my $k ( keys %fasta ) {
	
	print ">$k\n";
	print "$fasta{$k}\n";
}


