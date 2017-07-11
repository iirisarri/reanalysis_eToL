#!/usr/bin/perl

use warnings;
use strict;

##########################################################################################
#
# #################### Iker Irisarri. Jul 2017. Uppsala University ##################### #
#
# It modifies the headers of the peptide output from transdecoder
#
# By specifying the full name of the species with Jurgen's system, the output will have
#	already the required format
#
##########################################################################################


my $usage = "rename_transdecoder_output.pl infile.fa Species_name > STDOUT\n";
my $infile = $ARGV[0] or die $usage; 
my $species = $ARGV[1] or die $usage; 


# read-in fasta and replace 
open (IN, "<", $infile) or die "Cannot open file $infile\n";

while ( my $line =<IN> ) {
    chomp $line;
    next if ( $line =~ /^#/);

    # get headers & reformat them
    if ( $line =~ /^>(.+)/ ) {
	
		# header format:
		# >Gene.16::CF259181.1::g.16::m.16 type:5prime_partial len:118 gc:universal CF259181.1:645-292(-)

		my @headers = split ":", $1;
		my $gene =$headers[2];
		
		print ">$species\@$gene\n";
	}
	else {
	    
	    # remove stop codons from sequences
	    $line =~ s/\*//g;
	    print "$line\n";
	}
}
close(IN);

print STDERR "\ndone!\n\n";
