#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2017. Uppsala University ##################### #
# 
# Script to extract reference taxa for annotation of gene alignments
# Typically be used into a for loop
# 	for f in *.fa; do perl extract_homo_for_annotation.pl $f ; done
#
# Taxon names are hard-coded into arrays. They contain the preferred taxon in position 1
# followed by other taxa that can replace 1 if not present for that particular gene
#
# Assumes gene files in fasta format
#
##########################################################################################

my $usage = "extract_homo_for_annotation.pl infile > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;

# define header for target taxa (human)
# and alternative headers in order of preference

# Katz2015
#my @reference = qw (Op_me_hsap Op_me_ptro Op_me_mmus Op_me_rnor Op_me_clup Op_me_ecab);
# Burki2016
#my @reference = qw (Homo_sapiens Danio_rerio Branchiostoma_floridae Lottia_gigantea Daphnia_pulex Nematostella_vectensis);
# Brown2013
#my @reference = qw (HomoSapi Nematostel Ampqueen);
#Baurain2010
#my @reference = qw (Hydra_magnipapillata Nematostella_vectensis Reniera_sp. Monosiga_brevicollis Monosiga_ovata);

# for jawed vertebrates
my @reference = qw (Homo_sapiens Homo_sapiens Mus_musculus Felis_catus Canis_familiaris Callorhinchus_milii Lepisosteus_oculatus Lepisosteus_platyrhinus Latimeria_chalumnae Macropus_eugenii Neoceratodus_forsteri Lepidosiren_paradoxa Protopterus_aethiopicus Protopterus_annectens Ornithorhynchus_anatinus Loxodonta_africana Xenopus_tropicalis Bombina_maxima Acipenser_baeri Alligator_sinensis Ambystoma_mexicanum Amia_calva Anas_platyrhynchos Andrias_davidianus Anolis_carolinensis Atelopus_zeteki Basiliscus_plumifrons Boa_constrictor Caiman_crocodilus Calotriton_asper Carcharodon_carcharias Caretta_caretta Carlia_rubrigularis Chamaeleo_chamaeleon Chelonoidis_nigra Chiloscyllium_griseum Chinemys_reevesii Crocodylus_niloticus Crotalus_adamanteus Cyclorana_alboguttata Cynops_pyrrhogaster Danio_rerio Dasypus_novemcinctus Discoglossus_pictus Dromaius_novaehollandiae Echis_coloratus Elgaria_multicarinata Emys_orbicularis Espadarana_prosoblepon Eublepharis_macularius Gallus_gallus Geotrypetes_seraphini Ginglymostoma_cirratum Hymenochirus_curticeps Hynobius_chinensis Iguana_iguana Lampropholis_coggeri Leucoraja_erinacea Megophrys_nasuta Meleagris_gallopavo Micrurus_fulvius Monodelphis_domestica Neotrygon_kuhlii Notophthalmus_viridescens Opheodrys_aestivus Ophiophagus_hannah Pantherophis_guttatus Pelodiscus_sinensis Pelophylax_lessonae Pelophylax_nigromaculatus Pelusios_castaneus Phrynops_hilarii Pipa_pipa Pleurodeles_waltl Podarcis_sp. Pogona_vitticeps Polypterus_senegalus Proteus_anguinus Python_regius Raja_clavata Rana_chensinensis Rhinatrema_bivittatum Salamandra_salamandra Saproscincus_basiliscus Sarcophilus_harrisii Sceloporus_undulatus Scincella_lateralis Scyliorhinus_canicula Siren_lacertina Sistrurus_miliarius Sphenodon_punctatus Sternotherus_odoratus Struthio_camelus Taeniopygia_guttata Takifugu_rubripes Tarentola_mauritanica Thamnophis_elegans Trachemys_scripta Tupinambis_teguixin Typhlonectes_compressicauda Typhlonectes_natans);

# GET SEQUENCE

# read fasta file with SeqIO and store in %hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
	                	        '-format' => "fasta");

my %infile = ();
                		        
while (my $seq_obj = $seqio_obj->next_seq){

	# store sequences into %hash
    my $seqname = $seq_obj->primary_id;
	my $sequenc = $seq_obj->seq;
	# remove dashes if aligned
	$sequenc =~ s/-//g;
	
	$infile{$seqname} = $sequenc;
}

# search for references in order of preference

my $success = 0;

foreach my $tax ( @reference ) {

	if ( exists ( $infile{$tax} ) ) {
	    
	    print ">", $fasta, "@", $tax, "\n";
    	print $infile{$tax}, "\n";
    	$success++;
    	last;
    }
	# go to next item in @reference when array is not empty
    elsif ( scalar @reference !=0 && !exists ( $infile{$tax} ) ) {
    
    	next;
    }
}

# when none of the taxa in @reference are present
if ( $success == 0 ) {	
    	print STDERR "File $fasta: sequence for queried taxa do not exist\n";
}

