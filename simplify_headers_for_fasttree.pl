#!/usr/bin/perl

use strict;
use warnings;

##########################################################################################
#
# #################### Iker Irisarri. Jul 2017. Uppsala University ##################### #
# 
# Some fasta headers are too long and give errors in FastTree. This script simplifies them
# 	by removing redundant TRINITY contig info
#	Example: >Rhizaria_N/A_N/A_N/A_Gromiidae_Gromia_sphaerica_SRR2003399@TRINITY_DN4487_c0_g1::TRINITY_DN4487_c0_g1_i1::g.7052::m.7052LENGTH_249
#
##########################################################################################

my $usage = "simplify_headers_for_fasttree.pl infile > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;

open (IN, "<", $fasta);

while (my $line =<IN> ){

	chomp $line;

	if ( $line =~ /^>.+/ ) {
	
		# work with transcript information
		my ( $taxon, $transcript ) = split "\@", $line;

        #print "\t$taxon\n\t$transcript\n";                                                                                                                       
		
        if ( $transcript =~ /TRINITY_DN\d+_c\d+_g\d+::TRINITY_DN\d+_c\d+_g\d+_i\d+.+/ ) {

			my @transcripts = split "::", $transcript;	
	
			my $header = $taxon . "\@" . $transcripts[1] . "_" . $transcripts[2] . "_" . $transcripts[3];
			print "$header\n";
			next;
		}
	# print header
	print "$line\n";
	}
	else {
		# print sequence
		print "$line\n";
	}
}
