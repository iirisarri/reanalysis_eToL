#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2017. Uppsala University ##################### #
#
# Very simple blast parser to extract Query and Hit names into a table
#
# It requires that only the best hit has been saved:
# -num_descriptions 1 -num_alignments 1
#
# By default, it parses gene name from query fasta header in the format "gene@taxa"
# 
# Output by default is filename ENSPXXX gene_symbol gene_name
#
##########################################################################################

my $usage = "parse_blastp_annotations.pl blast_report > STDOUT\n";
my $blast_report = $ARGV[0] or die $usage;

# read blast report
my $report = new Bio::SearchIO(-format => "blast",
			       -file   => "<$blast_report"
    );

while( my $result = $report->next_result ) {

    while( my $hit = $result->next_hit ) {

	my $query = $result->query_name;
	my $hit_name = $hit->name;
	my $hit_desc = $hit->description;
	    
	# get gene name
	my ($gene, $taxon) = split /\@/, $query;

	# get ENSP
	my @hit_names = split / /, $hit_name;
	my $ensp = $hit_names[0];

	# get gene info
	$hit_desc =~ /[\w\d\s]+?gene_symbol:(\w+)/;
	my $symbol = $1;
	$hit_desc =~ /[\w\d\s]+?description:([\w\s\-\,]+)/;
	my $description = $1;
	
	print "$gene\t$ensp\t$symbol\t$description\n";
    }
}


print STDERR "\nDone!\n\n";

