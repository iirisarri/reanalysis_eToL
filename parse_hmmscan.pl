#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Mar 2017. Uppsala University ##################### #
#
# Parses hmmscan output (domain information, in tabular format (--domtblout) and selects
#	non-overlapping domains. It also prints (STDERR) overlapping domains that should be 
#	excluded.
#
# Possibility to filter by e-value (should be i-evalue from domains). Default: 10. 
# 	Non-significant hits are also printed to STDERR. This filter is on the i-evalue.
#	To avoid having too many results to parse, run hmmscan with '--E 1e-3', the
#	"full sequence evalue", which is most meaningful
#
# Example of hmmscan run:
# hmmscan --cpu 8 -E 1e-3 -o output.pfam --domtblout output.tblout Pfam-A.hmm infile.fa
#
##########################################################################################


## WHAT HAPPENS WHEN WE DO NOT INCLUDE A FILE

my $usage = "parse_hmmscan.pl hmmscan.out > STDOUT 2> (overlapping domains)\n";
my $infile = $ARGV[0] or die $usage;
my $evalue_threshold = 10; # no evalue threshold

open (IN, "<", $infile) or die "Cann't open $infile!\n";

my %hash;
my %hash_no_overlap;

MAIN: while ( my $line =<IN> ) {

	chomp $line;
	
	# skip headers
	next if $line =~ /^#.*/;
	
	# convert to tab-delimited file
	$line =~ s/\s+/\t/g;
	
	my @lines = split ("\t", $line);
	
	my $hit_name = $lines[0];
	#my $hit_acc = $lines[1];
	#my $hit_len = $lines[2];     
	my $query_name = $lines[3];
	my $evalue = $lines[6]; # full sequence evalue
	my $bias = $lines[14];
	my $coord_hit_st = $lines[15];
	my $coord_hit_end = $lines[16];
	my $coord_query_st = $lines[17];
	my $coord_query_end = $lines[18];
	my $hsp_num = 1;

	# skip non-significant domains (too high evalue) 
	if ( $evalue > $evalue_threshold ) {
	
		print STDERR "QUERY=$query_name HIT=$hit_name ";
		print STDERR "FULL-SEQ EVALUE=$evalue non-significant hit\n";
		next MAIN;
	}

	# MAKE HITS UNIQUE (BY HSP)
	# append number to hit if it contains multiple hsp (high scoring pairs; 
	#	ie, >1 homologous stretches for present a query_hit pair)
	if ( exists $hash{$query_name}{$hit_name} ) {
	
		# get last hsp_num and add one
		my @hits = sort keys ${ $hash{$query_name} }{$hit_name};
		my $hsp_num_old = pop @hits;
		my $hsp_num = $hsp_num_old +1;
	}
	
	# SKIP OVERLAPPING (including single position overlaps)
	# skip domain if overlapping with a previously existing one (for same query)
	# this is a greedy algorithm: it favors domains lower evalue (more significant) 
	# because they appear in decreasing order in the infile
	if ( exists $hash{$query_name} ) {
	
		# get "old" coordinates saved in the %hash
		# loop through $hit_names's
		foreach my $k ( keys %{ $hash{$query_name} } ) {

			# loop through $hsp_num's
			foreach my $l ( keys %{ ${ $hash{$query_name} }{$k} } ) {

				# get "old" coordinates
				my $coord_query_st_OLD =  ${ ${ $hash{$query_name} }{$k} }{$l}[2];
				my $coord_query_end_OLD =  ${ ${ $hash{$query_name} }{$k} }{$l}[3];
	
				#print "$coord_query_st_OLD-$coord_query_end_OLD\n";
	
				# SKIP OVERLAPPING HITS (including those with single position overlaps)
				# nS < oS && nE > oS # new hit overlaps from left
				if ( $coord_query_st <= $coord_query_st_OLD && $coord_query_end >= $coord_query_st_OLD ) {
				
					print STDERR "QUERY=$query_name HIT=$hit_name ";
					print STDERR "COORDS OLD=$coord_query_st_OLD-$coord_query_end_OLD NEW=$coord_query_st-$coord_query_end left overlap\n";
					next MAIN;
				}
				# nS < oE && nE > oE # new hit overlaps from right
				if ( $coord_query_st <= $coord_query_end_OLD && $coord_query_end >= $coord_query_end_OLD ) {

					print STDERR "QUERY=$query_name HIT=$hit_name ";
					print STDERR "COORDS OLD=$coord_query_st_OLD-$coord_query_end_OLD NEW=$coord_query_st-$coord_query_end right overlap\n";
					next MAIN;
				}
				# nS < oE && nE < oE # new hit is included in old hit
				if ( $coord_query_st_OLD <= $coord_query_st && $coord_query_st <= $coord_query_end_OLD && $coord_query_end <= $coord_query_end_OLD ) {

					print STDERR "QUERY=$query_name HIT=$hit_name ";
					print STDERR "COORDS OLD=$coord_query_st_OLD-$coord_query_end_OLD NEW=$coord_query_st-$coord_query_end total overlap\n";
					next MAIN;
				}
			}
		}
	}

	# After skipping non-significant and overlapping domains and making hits unique (by hsp)...

	# SAVE RECORD INTO %hash
	$hash{$query_name}{$hit_name}{$hsp_num} = [$coord_hit_st, $coord_hit_end, $coord_query_st, $coord_query_end, $evalue, $bias];

	#print Dumper \%hash;
}	

#print Dumper \%hash;

# PRINT OUT SUMMARY

# <= CHECK HERE!!

print "\nQUERY\tHIT_DOMAIN\tHIT_START\tHIT_END\tQUERY_SART\tQUERY_END\tE-VALUE\tBIAS\n";
foreach my $m ( sort keys %hash ) {

	foreach my $n ( sort keys %{ $hash{$m} } ) {

		foreach my $o ( sort keys %{ ${ $hash{$m} }{$n} } ) {

			print "$m\t$n\t";
			my @array = @{ ${ ${ $hash{$m} }{$n} }{$o} };
			print join ("\t", @array ), "\n";
			last;
		}
	}
}

print STDERR "\ndone!\n\n";

__END__

Testing the overlap between domains
n=new; o=old; S=start, E=end

nS < oS && nE > oS # new hit overlaps from left
nS < oE && nE > oE # new hit overlaps from right
nS < oE && nE < oE # new hit is included in old hit

nE < oS # new hit on the left, non-overlapping
nS > oE # new hit on the right, non-overlapping


# Exclude manually from infile the lines that we want to ignore?
# Eg save the info to exclude to exclude from hmmscan output file?