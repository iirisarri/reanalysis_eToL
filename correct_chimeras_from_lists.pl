#!/usr/bin/perl

use warnings;
use strict;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2017. Uppsala University ##################### #
#
# Very simple script to replace taxon names to chimera names
#
# It requires a fasta input file and a file with tab-delimited fields, like:
# #New name		Old name
# Am_ac_acas	>Acanthamoeba_castellanii
# EE_ap_Amas	>Amastigomonas_sp
# Sr_st_Apla	>Aplanochytrium
#
##########################################################################################


my $usage = "correct_chimeras_from_lists.pl infile.txt > STDOUT\n";
my $infile = $ARGV[0] or die $usage; 

# Burki
#=pod
my %tax_to_chimera = (	"EE_br_Pbif" => "EE_br_Brev",
						"EE_br_bant" => "EE_br_Brev",
						"Ex_ja_Jlib" => "Ex_ja_Jako",
						"Ex_ja_Jbah" => "Ex_ja_Jako",
						"Ex_ma_Mjak" => "Ex_ma_Mala",
						"Ex_ma_Mcal" => "Ex_ma_Mala",
						"Sr_st_pram" => "Sr_st_Phyt",
						"Sr_st_Psoj" => "Sr_st_Phyt"
					);
#=cut

# Baurain
=pod
my %tax_to_chimera = (	"Sr_ap_Cpar" => "Sr_ap_Cryp",
						"Sr_ap_chom" => "Sr_ap_Cryp",
						"Sr_di_Atam" => "Sr_ch_Din1",
						"Sr_di_Kmic" => "Sr_ch_Din1",
						"Sr_ch_Acar" => "Sr_ch_Din1",
						"Sr_di_Lpol" => "Sr_ch_Din1",
						"Sr_di_Htri" => "Sr_ch_Din1",
						"Sr_ch_Kbre" => "Sr_ch_Din1",
						"Pl_gr_otau" => "Pl_gr_Ostr",
						"Pl_gr_Oluc" => "Pl_gr_Ostr"
					);
=cut

open (IN1, "<", $infile) or die "Cannot open file $infile\n";

my %taxa = ();
my %seen = ();

my $original = 0;
my $substitution = 0;
my $final = 0;

while ( my $line =<IN1> ) {

	chomp $line;
	$original++;
	
	if ( exists $tax_to_chimera{$line} ) {

		my $chimera = $tax_to_chimera{$line};
		$substitution++;
		
		$seen{$chimera} = 1;
		
		if ( !exists $seen{$chimera} ) {
		
			print "$chimera\n";
			$final++;
		}
	}
	else {
		
		print "$line\n";
		$final++;

		# store what has been printed
		$seen{$line} = 1;
	}
}

print STDERR "\nNumber of lines in infile: $original\n";
print STDERR "Number of lines in outfile: $final\n";
print STDERR "Number of substitutions: $substitution\n\n";


