######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to get protein essential feature for protein pair list 
# essential ORFs.txt : (http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt)
# 
# "	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
# "	( the co-essential feature is 3-value categorized feature: 0 means NE/EN, 1 means NN, 2 means EE)
# 
# command yeast_pr_pair_file_name Essential_ORFs_list.txt outputfile


use strict; 
die "Usage: command yeast_pr_pair_file_name Essential_ORFs_list.txt outputfile \n" if scalar(@ARGV) < 3 ;

my ($int_file, $e_file, $out_file) = @ARGV;

open(ES, $e_file) || die(" Can not open file(\"$e_file\").\n"); 
my %yeast_eorfs = (); 
my ($count, $orf, $remain); 

while (<ES>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf, $remain) = split('\s+', $_);

	if (defined $yeast_eorfs{"$orf"})
	{
		$count = $yeast_eorfs{"$orf"} + 1 ; 
		$yeast_eorfs{"$orf"} =  $count; 
	}
	else 
	{
		$count = 1; 	
		$yeast_eorfs{"$orf"} = $count; 
	}
}


#--------------------- Begin to add essential feature  -------------------------------


open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");

my ( $orf1, $orf2, $es1, $es2, $flag); 
while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	($orf1, $orf2, $flag) = split('\s', $_);
	$es1 = $yeast_eorfs{$orf1};
	$es2 = $yeast_eorfs{$orf2};
	
	my $score = 0; 
	if (((defined $es1))&&((defined $es2)))
	{
		$score = 2; 
	}
	elsif ((!(defined $es1))&&(!(defined $es2)))
	{
		$score = 1; 
	} 
	else {
		$score = 0; 
	}
	print OUT "$score,$flag\n"; 	
}
close(INT);					
close(OUT);
