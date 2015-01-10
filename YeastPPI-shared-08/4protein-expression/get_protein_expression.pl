######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to get protein expression (abundance) for protein pair list 
# Based on the genome-wide protein localization analysis published by 
# Huh et al. (2003) Nature 425:686-691 and protein abundance data published by 
# Ghaemmaghami et al. (2003) Nature 425:737-741
#
# perl get_protein_expression.pl temp.line.txt OSheaLocalization.WeissmanAbundance.tab temp.line.exp

use strict; 
die "Usage: command yeast_pr_pair_file_name nature03_expression_file out_file_name\n" if scalar(@ARGV) < 3;

my ($int_file, $nat_file, $out_file) = @ARGV;

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");


#--------------------- read in the nature03 expression & localization file -------------------------------

open(NAT, $nat_file) || die(" Can not open file(\"$nat_file\").\n"); 
my (@per_line, $orf, $expression, %nat_exp); 

%nat_exp = (); 
while (<NAT>)	
{
	chomp $_;
	chop $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	@per_line = split('\t', $_);
	$orf = $per_line[1];
	$expression = $per_line[6];	
	$nat_exp{$orf} = $expression; 
}
close(NAT);


#--------------------- Begin to add nature03 protein expression feature  -------------------------------

my $count = 0; 
my ( $orf1, $orf2, $exp1, $exp2, $flag); 
while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	($orf1, $orf2, $flag) = split('\s', $_);
	$exp1 = $nat_exp{$orf1};
	$exp2 = $nat_exp{$orf2};
	
	my $score = 0; 
	if ((!(defined $exp1))||(!(defined $exp2)))
	{
		$score = -100; 
	}
	elsif (( $exp1 =~ /[0-9]+/)&&( $exp2 =~ /[0-9]+/ ))
	{
		$score =  log(abs( $exp1 - $exp2) + 1 )/log(10) ;
	} 
	else {
		$score = -100; 
	}
	print OUT "$score,$flag\n"; 	
	$count = $count + 1;
}
print " $count pairs !"; 
close(INT);					
close(OUT);