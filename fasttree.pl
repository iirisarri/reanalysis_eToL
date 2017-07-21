#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::SeqIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jul 2017. Uppsala University ##################### #
# 
# Script run fasttree (fastest possible option) after trimming gappy positions (>50%) with trimal
#
##########################################################################################

my $usage = "fasttree.pl alignment > tree\n";
my $fasta = $ARGV[0] or die $usage;
my $trimmed_aln = $fasta . ".trim";

# run trimal (remove columns >50% missing                                                                                                           
print STDERR "Trimming columns with >50% gaps...\n";
system("/opt/Trimal_1.4/trimal -gt 0.5 -in $fasta -fasta -out $trimmed_aln");
print STDERR "done.\n";

# run fasttree
print STDERR "Running FastTree analysis... OPTIONS: -lg -fastest -nosupport\n";
system ("/opt/FastTree/version_2.1.10/FastTreeMP -lg -fastest -nosupport $trimmed_aln > $trimmed_aln.fasttree.tre");

system("rm $trimmed_aln");

print STDERR "\ndone!\n\n";
