#!/usr/bin/perl

use warnings;
use strict;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2017. Uppsala University ##################### #
#
# Very simple script to standarize taxonIDs in fasta files
#
# It requires a fasta input file and a file with tab-delimited fields, like:
# #New name		Old name
# Am_ac_acas	>Acanthamoeba_castellanii
# EE_ap_Amas	>Amastigomonas_sp
# Sr_st_Apla	>Aplanochytrium
#
##########################################################################################


my $usage = "standarize_taxids_in_fasta.pl infile.fa translation_table > STDOUT\n";
my $infile = $ARGV[0] or die $usage; 
my $table = $ARGV[1] or die $usage; 

# read-in and store taxid correspondence
open (IN1, "<", $table) or die "Cannot open file $table\n";

my %taxids = ();

while ( my $line =<IN1> ) {
	chomp $line;
	next if ( $line =~ /^#.*/);
	 
	my ( $new, $old ) = split ("\t", $line);

	$taxids{$old} = $new;
}
close(IN1);


# read-in fasta and replace 
open (IN2, "<", $infile) or die "Cannot open file $infile\n";

while ( my $line2 =<IN2> ) {
	chomp $line2;
	next if ( $line2 =~ /^#/);

	# get headers
	if ( $line2 =~ /^>.+/ ) {
	
		if ( exists $taxids{$line2} ) {
		
			print ">$taxids{$line2}\n";
		}
		else {
		
			print "$line2\n";
			print STDERR "This taxa has no correspodence: $line2\n";
		}
	}
	# print sequence lines
	else {

		print "$line2\n";
	} 
}
close(IN2);

print STDERR "\ndone!\n\n";

