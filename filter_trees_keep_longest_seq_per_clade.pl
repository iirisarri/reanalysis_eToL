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
# Script to filter gene trees that contain too many sequences (added after blast) 
# 	It chooses longest sequence in every clade given that all child leaves belong to the 
#	same species. All redundant sequences are removed from the corresponding fasta file.
#
# Script will also exclude sequences that are not present in the tree but may be present
#	in the fasta file (but were removed by trimal when running fasttree.pl). To avoid this
#	behaviour, simply comment out line 144 in printing block
#
# This will still keep multiple sequences per species if these appear in different clades
#
# Script ignores the root of the tree (effect considered minor for very large trees)
#
##########################################################################################

my $usage = "filter_trees_keep_longest_seq_per_clade.pl infile.tre infile.fa > stdout 2> removed taxa\n";
my $tree = $ARGV[0] or die $usage;
my $fasta = $ARGV[1] or die $usage;

my $outgroup; # by default, script will re-root tree in first encountered prokaryote
my @taxa_to_rm_current; # to store shorter sequences of monophyletic groups, specific for each $node
my %taxa_to_remove; # to store shorter sequences of monophyletic groups for removal from fasta, global for $tree
my %leaves = ();

# read in trees
my $treeio = new Bio::TreeIO(-file   => "$tree",
						     -format => "newick");

my $node_count = 0; # will be used to give names to nodes

while( my $tree = $treeio->next_tree ) {

	# store all taxa in tree tree. Any sequence present in fasta but not in the tree 
	# will be excluded (e.g., excluded by trimal when using fasttree.pl)
	%leaves = %{ &get_all_leaves($tree->get_root_node) };

	NODE: for my $node ( $tree->get_nodes ) {
	
		$node_count++;

		# save @taxa_to_rm_current from previous into %taxa_to_remove
		# before initializing stuff
		# in case of nested clades all containing the same taxa, records will be overwritten
		# but this is not problematic since the longest sequence will be always kept
		foreach my $e ( @taxa_to_rm_current ) { # this will be skipped if array is empty
		
			$taxa_to_remove{$e} = 1;
		}
		my %descendents = ();
		my $previous_descendent = "0";
		my $previous_taxa = "0";
		my $previous_length = "0";
		@taxa_to_rm_current = ();

		for my $descendent ( $node->get_all_Descendents ) {
		
			# get all leaves for each node and save them
			if ( $descendent->is_Leaf ) {
			
				$descendents{$descendent->id} = "1";
			}
		}

		# just checks for debugging purpuses
		#print Dumper \%taxa_to_remove;
		#print Dumper \%descendents;

		# process descendents from that node
		foreach my $k ( keys %descendents ) {
		
			# get species and transcript info
			my ( $taxa, $transcript ) = split "\@", $k;
			
			# skip if any transcript does not contain length information
			if ( $transcript !~ /.+?LENGTH_(\d+)/ ) {
			
				@taxa_to_rm_current = (); # empty taxa to remove
				next NODE;
			}	
			
			$transcript =~ /.+?LENGTH_(\d+)/;
			my $length =  $1;
			
				# for first record in %descendent
				if ( $previous_descendent eq "0" ) {
				
					$previous_descendent = $k;
					$previous_taxa = $taxa;
					$previous_length = $length;
					next;
				}
				# for next records, check if taxa are different
				if ( $taxa ne $previous_taxa ) {
				
					@taxa_to_rm_current = (); # empty taxa to remove
					next NODE; # if leave is different, go to next node
				}
				# at this point $taxa eq $first_k
				# if shorter, keep current record for removal
				if ( $length < $previous_length ) {
				
					push (@taxa_to_rm_current, $k);
					next;
				}
				# if longer, remove previous record && update $previous_descendent and $previous_length to current record
				if ( $length > $previous_length ) {
				
					push (@taxa_to_rm_current, $previous_descendent);
					$previous_descendent = $k;  # update with current
					$previous_taxa = $taxa;		# update with current
					$previous_length = $length;	# update with current
				}					
			}
	}
}

# print to STDERR removed sequences
foreach my $to_remove ( sort keys %taxa_to_remove ) {

	print STDERR "$to_remove\n";
}

# modif from exclude_seq_from_fasta.pl

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
                	        	'-format' => "fasta",
                	        	'alphabet' => "protein");

while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;
    
    # exclude taxa is missing from the tree (e.g. removed by trimal)
    next if ( !exists $leaves{$seqname} );

    if ( !exists ( $taxa_to_remove{$seqname} ) ) {

	print ">", $seq_obj->primary_id, "\n";
	print $seq_obj->seq, "\n";
    }
}



## SUBROUTINES ##

sub get_all_leaves {
    #Esta rutina devuelve todas las "hojas" descendientes de un nodo
    my $initial_node = $_[0];
    my %nodes;
    if( $initial_node->is_Leaf ) {
	$nodes{ $initial_node->id } = 1;
	return \%nodes;
    }
    foreach my $node ( $initial_node->get_all_Descendents() ) {
		if( $node->is_Leaf ) { 
			# for example use below
			$nodes{$node->id} = 1;
		}
    }
    return \%nodes;
}

__END__

# NEED TO REROOT TREE?
# get all leaves
#    my %leaves = %{ &get_all_leaves($tree->get_root_node) };
    
    # root tree in one prokaryote sequence
#    REROOT: foreach my $leaf (keys %leaves) {
#
#		if ( $leaf =~ /^Prokaryota.+/ ) {
#		
#			$outgroup = $leaf;
#		    my $root = $tree->find_node( -id => $outgroup );
#		    $tree->reroot( $root );
#		    last REROOT; # scape foreach loop after rooting
