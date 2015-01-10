######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# 
# This program is to make the pair list from the protein list 
#
# perl make_fullpair_4protein.pl protein_list.txt pos_list output_pos output_other 
# 
# 
# $ perl  make_fullpair_4protein.pl YeastGeneListOrfGeneName-106_pval_v9.0.txt Science03-pos_MIPS_complexes.txt mipsPosPair.txt mipsRandpair.txt
#==> There are 6270 unique proteins in original list.
#==>  There should have 19653315 pairs possibly generated totally !
# There are 8617 POS pairs originally .
# fullpairs has: 7390 POS pairs.
# fullpairs has: 19645925 RAND pairs.
#  ==>  There are 19653315 pairs generated !
# 


use strict; 
die "Usage: command gene_name_file pos_pairFile outPosPairFile outRandPairFile \n" if scalar(@ARGV) < 4; 

my ($gene_file, $pos_pairFile, $outPosFileName, $outRandFileName ) = @ARGV; 


#--------------------- read in the gene name file -------------------------------
my (@gene, $per_gene, $count); 
my %proteinLookup = ();                  # lookup table to gene list

open(GENE, $gene_file) || die(" Can not open file(\"$gene_file\").\n"); 
$count = 0;
@gene = (); 

while (<GENE>)	
{
	chomp $_;
	chop $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	my @cur_line = split('\t', $_); 
	$per_gene = $cur_line[0] ;
	$gene[$#gene+1] = $per_gene;
	
	if (defined $proteinLookup {"$per_gene"})	
	{
		print "-- $per_gene  Duplicate in protein list ! "
	}
	else {
		$proteinLookup{"$per_gene" } =  $count; 
		$count = $count + 1; 
	}
}
close(GENE);

print "==> There are $count unique proteins in original list. \n"; 
my $temptemp = $count * ( $count - 1) / 2; 
print "==>  There should have $temptemp pairs possibly generated totally ! \n"; 



#--------------------- read in the original pos protein pair list file -------------------------------

open(PR1, $pos_pairFile) || die(" Can not open file(\"$pos_pairFile\").\n");

my (@per_line, $orf1, $orf2, $flag, $index, %posPairLookup); 

$index = 0; 
%posPairLookup = (); 
while (<PR1>)
{
	chomp; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	@per_line = split('\t', $_); 
	
	$orf1 = $per_line[0]; 
	$orf2 = $per_line[1];

	my $temp = $orf1.":".$orf2; 
	$posPairLookup{"$temp"} = 1 ; 
	
	$index = $index +1; 
}
close(PR1); 
print "# There are $index POS pairs originally .\n";  



#--------------------- generate two files we need -------------------------------

open(OUTP, "> $outPosFileName") || die(" Can not open file(\"$outPosFileName\").\n"); 
open(OUTN, "> $outRandFileName") || die(" Can not open file(\"$outRandFileName\").\n"); 

my ($countp, $countn, $i, $j); 
$countp = 0; 
$countn = 0; 

for ($i = 0; $i < $count; $i ++)
{
	for ($j = $i+1; $j < $count; $j ++)
	{
		$orf1 = $gene[$i]; 
		$orf2 = $gene[$j];	
		
		my $temp1 = $orf1.":".$orf2; 		
		my $temp2 = $orf2.":".$orf1; 		
		if (defined $posPairLookup{"$temp1"}) 
		{
			$flag = $posPairLookup{"$temp1"}; 
			print OUTP "$orf1\t$orf2\t$flag\n";
			$countp = $countp +1; 
		}
		elsif ($posPairLookup{"$temp2"}) 
		{
			$flag = $posPairLookup{"$temp2"}; 
			print OUTP "$orf1\t$orf2\t$flag\n";
			$countp = $countp +1; 
		}
		else {
			print OUTN "$orf1\t$orf2\t0\n";
			$countn = $countn +1; 
		}
	}
}
print "# fullpairs has: $countp POS pairs.\n";  
print "# fullpairs has: $countn RAND pairs.\n";  

$temptemp = $countp + $countn; 
print "  ==>  There are $temptemp pairs generated ! \n"; 

close(OUTP);								
close(OUTN);

