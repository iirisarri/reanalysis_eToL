#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;

# Iker Irisarri. Dec 2017. Uppsala University
# simple code to remove only the branch leading to monophyly of PPE (everything else stays as it was)
# to be used in delta_per_gene.pl

# THIS WORKS FINE TO REMOVE MRCA BRANCH OF PPE

my $input = shift;
 
my $in = new Bio::TreeIO(-file => $input,
           		         -format => "newick");
                         
my $tree = $in->next_tree;

my @nodes = qw (Pl_gl_Cpad Pl_gr_Cper Pl_gr_crei Pl_gr_micr Pl_gr_Mver Pl_gr_Ostr Pl_gr_ppat Pl_gr_Ptae Pl_rh_cmer Pl_rh_Gsul);
my @node_obj;

# get node objects from PPE
foreach my $i (@nodes) {
    my $node = $tree->find_node(-id => $i);
	push (@node_obj, $node); # $node is a Bio::Tree::Node object
}

# get LCA of PPE node objects
my $lca = $tree->get_lca(-nodes => \@node_obj);

# this seems to work
# set lca branch length to 0
$tree->splice($lca);

# remove that node from tree

my $h0_tree = new Bio::TreeIO(-file => ">$input.h0_constr.tre",
	                         -format => 'newick');

# print out tree to file
$h0_tree->write_tree($tree)


__END__

# Here I am trying to collapse all ancestral branches of each glauco, green and reds. It seems to work.
# However, these 3 lineages are not monophyletic in gene trees, so this would not always work

my $input = shift;
my $h0_tree_name = $input . "_h0";
 
my $in = new Bio::TreeIO(-file   => $input,
                            -format => "newick");

                         
my $tree = $in->next_tree;

my @glauco = qw (Pl_gl_Cpad);
my @green = qw (Pl_gr_Cper Pl_gr_crei Pl_gr_micr Pl_gr_Mver Pl_gr_Ostr Pl_gr_ppat Pl_gr_Ptae);
my @red = qw (Pl_rh_cmer Pl_rh_Gsul);

my @glauco_obj;
my @green_obj;
my @red_obj;

# get node objects from PPE

if ( scalar @glauco > 1) {
 
	foreach my $i (@glauco) {
		my $node_i = $tree->find_node(-id => $i);
		push (@glauco_obj, $node_i); # $node is a Bio::Tree::Node object
	}

	# get LCA of PPE node objects
	my $lca_glauco = $tree->get_lca(-nodes => \@glauco_obj);

	# get ancestor
	my $anc_glauco = $lca_glauco->ancestor;

	# remove that branch
	$tree->splice($anc_glauco);
}

# same for green
foreach my $j (@green) {
    my $node_j = $tree->find_node(-id => $j);
	push (@green_obj, $node_j); # $node is a Bio::Tree::Node object
}
my $lca_green = $tree->get_lca(-nodes => \@green_obj);

my $anc_green = $lca_green->ancestor;

$tree->splice($anc_green);

# and same for red
foreach my $k (@red) {
    my $node_k = $tree->find_node(-id => $k);
	push (@red_obj, $node_k); # $node is a Bio::Tree::Node object
}
my $lca_red = $tree->get_lca(-nodes => \@red_obj);

my $anc_red = $lca_red->ancestor;

$tree->splice($anc_red);

my $h0_tree = new Bio::TreeIO(-file => ">$h0_tree_name.tre",
                         -format => 'newick');

# print out tree to file
$h0_tree->write_tree($tree)

__END__


# This removes the branch and all its descendents. Not what I wanted...

my $input = shift;
my $h0_tree_name = $input . "_h0";
 
my $in = new Bio::TreeIO(-file   => $input,
                            -format => "newick");

                         
my $tree = $in->next_tree;

my @nodes = qw (Pl_gl_Cpad Pl_gr_Cper Pl_gr_crei Pl_gr_micr Pl_gr_Mver Pl_gr_Ostr Pl_gr_ppat Pl_gr_Ptae Pl_rh_cmer Pl_rh_Gsul);
my @node_obj;

# get node objects from PPE
foreach my $i (@nodes) {
    my $node = $tree->find_node(-id => $i);
	push (@node_obj, $node); # $node is a Bio::Tree::Node object
}

# get LCA of PPE node objects
my $lca = $tree->get_lca(-nodes => \@node_obj);

# remove that node from tree
$tree->remove_Node($lca);

my $h0_tree = new Bio::TreeIO(-file => ">$h0_tree_name",
                         -format => 'newick');

# print out tree to file
$h0_tree->write_tree($tree)