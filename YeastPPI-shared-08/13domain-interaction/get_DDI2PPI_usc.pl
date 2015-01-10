######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3

# Program to Extract Features of DDI derived PPI for protein pair list subset
# 
# Use two DDI2PPI data: 
##	ProteinAction_80SGDMIPS.txt        ( trained from the MIPS PPI ; Note: This feature overfits. So we would not use anymore )
#	ProteinAction_025_80SGDY2H.txt 		( trained from the Y2H PPI )
#
#  -------------------------------------------------------------------------------------------------
## ( Note: this trained from the MIPS PPI ; Note: This feature overfits. So we would not use anymore )
##  $ perl  get_DDI2PPI_usc.pl ./lists/sciencesubset.txt YeastGeneListOrfGeneName-106_pval_v9.0.txt ProteinAction_80SGDMIPS.txt ./lists/sciencesubset.MIPSddippi
##  Gene Name mapping protein - 6270;
##  MIPS  DDI2PPI Pair: 75468;
#
#  $ perl  get_DDI2PPI_usc.pl ./lists/sciencesubset.txt YeastGeneListOrfGeneName-106_pval_v9.0.txt ProteinAction_025_80SGDY2H.txt ./lists/sciencesubset.Y2Hddippi
#  Gene Name mapping protein - 6270;
#  Y2H DDI2PPI Pair: 125263;
#
#  ------------------------------------------------------------------------------------------------- 
#   The input pair list file format: “ORF1 ORF2 Flag” ( 0 rand or 1 postive)
#   The out put file format:  1 – real value feature, the last one is the class flag
#   This feature is between 0 and 1… 



use strict; 
die "Usage: command protein_pair_file yeast_name_mapping usc_ddi2ppi_file out_file_name\n" if scalar(@ARGV) < 4;

my ($int_file, $genename_file , $ddippi_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 



#--------------------- read in gene name mapping file  -------------------------------

my %gene_name = (); 

open( GEN, $genename_file) || die(" Can not open file(\"$genename_file\").\n"); 

while (<GEN>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\t', $_);
	
	$orf1 = $orf[0]; 
	$orf2 = substr($orf[1], 0, -1); 

	if ($orf2 eq "#REF!")  
		{ $gene_name{ "$orf1" } = $orf1 ;  }
	else 
		{ $gene_name{ "$orf1" } = $orf2 ;  }
}

close(GEN);



#--------------------- read in the ddi2ppi USC file  -------------------------------

open(DDI, $ddippi_file) || die(" Can not open file(\"$ddippi_file\").\n"); 

my %ddi_ppi =(); 

while (<DDI>)	
{
	chomp; 
	next if /^#/;			#ignore comments
	next if /^\s*$/; 			#ignore blank lines
	
	my @per_line = split(/\s+/, $_); 
	if ( $per_line[0] eq "")
		{ shift(@per_line); }; 
	$orf1 = $per_line[1]; 
	$orf2 = $per_line[3];
	
	my @temp = split(':', $orf2); 
	$orf2 = $temp[0]; 
	my $prob = $per_line[4]; 
	
	$ddi_ppi{"$orf1:$orf2"} = $prob ; 
}
close(DDI);



#--------------------- some statistic numbers  -------------------------------

my $proteinSize = scalar keys (%gene_name); 
my $pairsize = scalar keys (%ddi_ppi);

print "\nGene Name mapping protein - $proteinSize; \nDDI2PPI Pair: $pairsize; \n"; 


#--------------------- Begin to process the yeast_input_pair set and find if the pair in ddi2ppi set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");

my $score = 0; 
$count =  0; 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\s', $_);

	my $gene1 = $gene_name{ "$orf1" }; 
	my $gene2 = $gene_name{ "$orf2" }; 
	
	$pair_l = $gene1.":".$gene2; 
	$pair_r = $gene2.":".$gene1;
	
	if (defined $ddi_ppi{"$pair_l"}) 
		{$score = $ddi_ppi{"$pair_l"} ; }
	elsif (defined $ddi_ppi{"$pair_r"})
		{$score = $ddi_ppi{"$pair_r"} ; }
	else
		{$score = 0 ; }			
	
	print OUT "$score,$flag\n";
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);