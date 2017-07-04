#!/usr/bin/perl

##USAGE: perl FigTree2seq_remove_red.pl nexus_file original_fasta_file
##This script removes the red sequences coloured in FigTree from a fasta file
##Last modified Aug 18, 2014. 
## Fabien Burki

use strict;
use warnings;
use Bio::SeqIO;

my $nexus = $ARGV[0];
chomp $nexus;
my $short_nexus = $nexus;
$short_nexus =~ s/\.\w+$//;

my $fasta = $ARGV[1];
chomp $fasta; #removes newline/whitespace
my $short_fasta = $fasta; # $short_fasta is a variable to store the name of the trimal file
#$short_fasta =~ s/\.\w+$//; # removed the last ext

open IN, "<$nexus";
my $in  = Bio::SeqIO->new(-file => "$fasta", -format => 'Fasta');
my $out = Bio::SeqIO->new(-file => ">$short_fasta".'.cleaned', -format => 'Fasta'); #modified outname

my $all_lines = do {local $/; <IN>};
my ($a, $b) = split ('taxlabels', $all_lines);
my ($c, $d) = split ("\n\;", $b);
$c =~ s/\n//g;
$c =~ s/\'//g;
my @ids = split ("\t", $c);

#For each sequence in the tre file, if the color id is red, it skips it
##Otherwise it adds the sequence to the newly created fasta file. 

my %hashids;
foreach my $id (@ids) {
	if ($id =~ /\S+\[\&\!color\=\#ff0000\]/) {
#		print "$id\n";
		next;
	}
	else {
		$hashids{$id} = 1;
#		print "$id\n";
	}
}

while (my $seq = $in->next_seq()) {
	my $id = $seq->id();
	foreach my $key (keys %hashids) {
		$key =~ s/\[\&.*$//; #remove the color tag so that the name matches the sequence id
#		print "$key\n";
		if ($key eq $id) {
			print "$key\n";
	   		$out->write_seq($seq); 
	  	}
	}
}


