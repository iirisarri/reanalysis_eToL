#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::TreeIO;

##########################################################################################
#
# #################### Iker Irisarri. Dec 2017. Uppsala University ##################### #
# 
# calc_delta_per_gene.pl
# 
# It will calculate likelihood differences between 3 trees:
# HA2 – unconstrained
# HA1 – constrained with PPE monophyly
# H0 – same as HA1 except it does not contain branch defining PPE monophyly 
# 
# likelihoods are calculated using best-fit models using IQTREE
#
# Requiremenets:	'iqtree-omp' must be on the path
#					BIOPERL
#
# Known problems: alignments with identical sequences will generate constrained iqtree to
#	abort & script will die. In such cases, repeat using the generated .uniqueseq.phy files
#
##########################################################################################

my $usage = "perl calc_delta_per_gene.pl *sdtaxids.fas\n";
my @infiles = @ARGV or die $usage; 
my $outfile = "_delta_results.tab";

open (OUT, ">", $outfile);

print OUT "#HA2: unconstrained\n#HA1: PPE monoph constrained\n#H0: same as HA1 without monophyly branch\n";
print OUT "#delta_A1 = lnL(X|H_A1) - lnL(X|H_0)\n#delta_A2 = lnL(X|H_A2) - lnL(X|H_0)\n";
print OUT "#Alignment\tLikelihood_H0\tLikelihood_HA1\tLikelihood_Ha2\tdelta_a1\tdelta_a2\n";

foreach my $fasta (@infiles) {

# 1. RUN UNCONSTRAINED TREE (HA2)
	system("iqtree-omp -nt 1 -s $fasta -m TEST");



# 2. RUN CONSTRAINED TREE (HA1)
	# make constraint file
	my $constr_file = $fasta . ".constr";
	my @PPE;
	my @nonPPE;

	my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta", 
					'-format' => "fasta"
					);

	while (my $seqio_obj = $seqio_obj->next_seq) {
	
		if ( $seqio_obj->primary_id =~ /^Pl_.+/ ) {
	
			push (@PPE, $seqio_obj->primary_id);
		}
		else {
	
			push (@nonPPE, $seqio_obj->primary_id);
		}
	}

	open(OUT2, ">", $constr_file);

	print OUT2 "((", join (",", @PPE), "),(", join (",", @nonPPE), "));\n";

	close(OUT2);

	# get best model
	my $logfile = $fasta . ".log";

	my $aicc_string = `grep "Corrected Akaike Information Criterion" $logfile`;

	chomp $aicc_string;

	my ($aicc, $model) = split ": ", $aicc_string;

	# run iqtree
	system("iqtree-omp -nt 1 -s $fasta -m $model -g $constr_file -pre $constr_file");



# 3. COMPUTE LIKELIHOOD OF H0

	# remove monophyly branch from tree HA1	
	# rm_branch_from_tree.pl
	my $input = $constr_file . ".treefile";
 
	my $in = new Bio::TreeIO(-file => $input,
							 -format => "newick");
						 
	my @node_obj;

	my $tree = $in->next_tree;

	# get node objects from PPE
	foreach my $i (@PPE) {
		my $node = $tree->find_node(-id => $i);
		push (@node_obj, $node); # $node is a Bio::Tree::Node object
	}

	# get LCA of PPE node objects
	my $lca = $tree->get_lca(-nodes => \@node_obj);

	# remove the LCA branch (defiles the monophyly)
	$tree->splice($lca);

	# write new tree to file
	my $constr_file2 = $input . "h0_constr2";

	my $h0_tree = new Bio::TreeIO(-file => ">$constr_file2",
								 -format => 'newick');

	# print out tree to file
	$h0_tree->write_tree($tree);


	# run iqtree
	system("iqtree-omp -nt 1 -s $fasta -m $model -t $constr_file2 -blscale -pre $constr_file2");

# 4. COMPUTE LIKELIHOOD DIFFERENCES AND PRINT
	
	# get likelihoods from output files
	# typical line: Log-likelihood of the tree: -4947.0596 (s.e. 312.5233)
	my $h0_string = `grep "Log-likelihood of the tree:" $constr_file2.iqtree`;
	chomp $h0_string;
	my $ha1_string = `grep "Log-likelihood of the tree:" $constr_file.iqtree`;
	chomp $ha1_string;
	my $ha2_string = `grep "Log-likelihood of the tree:" $fasta.iqtree`;
	chomp $ha2_string;

	$h0_string =~ /Log-likelihood of the tree: (-\d+\.\d+) \(s\.e\. \d+\.\d+\)/;
	my $lh_h0 = $1;

	$ha1_string =~ /Log-likelihood of the tree: (-\d+\.\d+) \(s\.e\. \d+\.\d+\)/;
	my $lh_ha1 = $1;

	$ha2_string =~ /Log-likelihood of the tree: (-\d+\.\d+) \(s\.e\. \d+\.\d+\)/;
	my $lh_ha2 = $1;
	
	# checking I am using the right numbers
	#print "$h0_string\n$ha1_string\n/$ha2_string\n";
	#print "lh H0: $lh_h0\nlh HA1: $lh_ha1\nlh HA2: $lh_ha2\n";
		
	# compute likelihood differences
	# \delta_A1 = lnL(X|H_A1) - lnL(X|H_0)
	# \delta_A2 = lnL(X|H_A2) - lnL(X|H_0)

	my $delta_a1 = $lh_ha1 - $lh_h0;
	my $delta_a2 = $lh_ha2 - $lh_h0;

	# print out stuff
	print OUT "$fasta\t$lh_h0\t$lh_ha1\t$lh_ha2\t$delta_a1\t$delta_a2\n";
}
	
