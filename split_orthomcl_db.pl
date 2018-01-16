#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

##########################################################################################
#
# #################### Iker Irisarri. Jan 2018 Uppsala University ###################### #
#
# Splits orthoMCL database into individual fasta profiles for provided taxa & renames sequences
#
##########################################################################################

my $usage = "split_orthomcl_db.pl\n";
my $fasta="aa_seqs_OrthoMCL-5.fasta";


my %queries = (
    "aaeg" => "Opisthokonta-Metazoa_Arthropoda_Insecta_Diptera_Culicidae_Aedes_aegypti_OrthoMCL5",
    "agam" => "Opisthokonta-Metazoa_Arthropoda_Insecta_Diptera_Culicidae_Anopheles_gambiae_OrthoMCL5",
    "afum" => "Opisthokonta-Fungi_Ascomycota_Eurotiomycetes_Eurotiales_Aspergillaceae_Aspergillus_fumigatus_OrthoMCL5",
    "bmor" => "Opisthokonta-Metazoa_Arthropoda_Insecta_Lepidoptera_Bombycidae_Bombyx_mori_OrthoMCL5",
    "bmaa" => "Opisthokonta-Metazoa_Nematoda_Chromadorea_Spirurida_Onchocercidae_Brugia_malayi_OrthoMCL5",
    "cbri" => "Opisthokonta-Metazoa_Nematoda_Chromadorea_Rhabditida_Rhabditidae_Caenorhabditis_briggsae_OrthoMCL5",
    "clup" => "Opisthokonta-Metazoa_Chordata_Mammalia_Carnivora_Canidae_Canis_lupus-familiaris_OrthoMCL5",
    "cimm" => "Opisthokonta-Fungi_Ascomycota_Eurotiomycetes_Onygenales_N/A_Coccidioides_immitis_OrthoMCL5",
    "cpip" => "Opisthokonta-Metazoa_Arthropoda_Insecta_Diptera_Culicidae_Culex_pipiens_OrthoMCL5",
    "ecab" => "Opisthokonta-Metazoa_Chordata_Mammalia_Perissodactyla_Equidae_Equus_caballus_OrthoMCL5",
    "egos" => "Opisthokonta-Fungi_Ascomycota_Saccharomycetes_Saccharomycetales_Saccharomycetaceae_Eremothecium_gossypii_OrthoMCL5",
    "mmul" => "Opisthokonta-Metazoa_Chordata_Mammalia_Primates_Cercopithecidae_Macaca_mulatta_OrthoMCL5",
    "mgri" => "Opisthokonta-Fungi_Ascomycota_Sordariomycetes_Magnaporthales_Magnaporthaceae_Magnaporthe_grisea_OrthoMCL5",
    "mdom" => "Opisthokonta-Metazoa_Chordata_Mammalia_Didelphimorphia_Didelphidae_Monodelphis_domestica_OrthoMCL5",
    "oana" => "Opisthokonta-Metazoa_Chordata_Mammalia_Monotremata_Ornithorhynchidae_Ornithorhynchus_anatinus_OrthoMCL5",
    "ptro" => "Opisthokonta-Metazoa_Chordata_Mammalia_Primates_Hominidae_Pan_troglodytes_OrthoMCL5",
    "phum" => "Opisthokonta-Metazoa_Arthropoda_Insecta_Psocodea_Pediculidae_Pediculus_humanus_OrthoMCL5",
    "tnig" => "Opisthokonta-Metazoa_Chordata_Actinopteri_Tetraodontiformes_Tetraodontidae_Tetraodon_nigroviridis_OrthoMCL5",
    "anid" => "Opisthokonta-Fungi_Ascomycota_Eurotiomycetes_Eurotiales_Aspergillaceae_Emericella_nidulans_OrthoMCL5",
    "glae" => "Excavata-Metamonada_Metamonada_N/A_N/A_Hexamitidae_Giardia_lamblia_OrthoMCLv5",
    "gzea" => "Opisthokonta-Fungi_Ascomycota_Sordariomycetes_Hypocreales_Nectriaceae_Gibberella_zeae_OrthoMCL5",
    "lbra" => "Excavata-Discoba_Euglenozoa_N/A_Kinetoplastida_Trypanosomatidae_Leishmania_braziliensis_OrthoMCL5",
    "lmex" => "Excavata-Discoba_Euglenozoa_N/A_Kinetoplastida_Trypanosomatidae_Leishmania_mexicana_OrthoMCL5",
    "pchr" => "Opisthokonta-Fungi_Basidiomycota_Agaricomycetes_Polyporales_Phanerochaetaceae_Phanerochaete_chrysosporium_OrthoMCL5",
    "psti" => "Opisthokonta-Fungi_Ascomycota_Saccharomycetes_Saccharomycetales_Pichiaceae_Pichia_stipitis_OrthoMCL5",
    "tcon" => "Excavata-Discoba_Euglenozoa_N/A_Kinetoplastida_Trypanosomatidae_Trypanosoma_congolense_OrthoMCL5"
    );


my $outfile;


# read fasta file with SeqIO and store in %hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
				'-format' => "fasta");

while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;

    my @seqnames = split (/\|/, $seqname);

    my $code = $seqnames[0];
    #$seqnames[1] =~ /\s*([\S]+)\s*/; # remove gaps in this field
    my $seqid = $seqnames[1];

    # remove any

    if ( exists $queries{$code} ) {

	my $outfile = $code . "_OrthoMCL5.fa";

	open (OUT, ">>", $outfile); # append

	print OUT ">", $queries{$code}, "\@$seqid\n";
	print OUT $seq_obj->seq, "\n";

	close(OUT);
    }
}
