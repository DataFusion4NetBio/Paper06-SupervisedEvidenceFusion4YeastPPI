######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to get gene expression (abundance) for protein pair list 
# Based on the gene expression data given by Ziv published by 
# Computational discovery of gene modules and regulatory networks. 
# Nature Biotechnology, 21(11) pp. 1337-42, 2003  
# 
# perl get_gene_expression.pl ./lists/Science03PPI.sublist.txt YeastGeneListOrfGeneName-106_pval_v9.0.txt all_expression_fixed_s4_csv.txt expressionYanjunSplit.txt 0.7 ./lists/Science03PPI.sublist.genexp
#
# Note: 
# - in the protein pair list file, 
# - for each protein pair, we have a flag representing the class label: "1" means postive pair, "0" means random pairs
#"	The input pair list file format: ORF1 ORF2 Flag
#			( flag is either 0 rand or 1 postive)
#"	The out put file format: 20 gene expression features separated by ",", the last one is the class flag


#use strict; 

die "Usage: command yeast_pr_pair_file_name yeast_name_list gene_expression_file expression_split_file percent_miss out_file_name\n" if scalar(@ARGV) < 5;

my ($pr_pair_file, $gene_list_file, $gene_expression_file, $expression_split_file, $percent_miss, $out_file_name) = @ARGV;



#--------------------- read in the gene name list  file -------------------------------

open(GENLIST, $gene_list_file) || die(" Can not open file(\"$gene_list_file\").\n"); 
my ( %gene_name ); 

%gene_name = (); 

$indexnum = 0; 
while (<GENLIST>)	
{
	chomp $_;
	chop $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	my @cur_per_line = (); 
	@cur_per_line = split('\t', $_);

	$indexnum = $indexnum + 1; 	
	$gene_name{"$cur_per_line[0]"} = $indexnum; 
}
close(GENLIST);



#--------------------- read in the gene expression file -------------------------------

open(GEN, $gene_expression_file) || die(" Can not open file(\"$gene_expression_file\").\n"); 
my ( %gene_exp, $indexnum ); 

%gene_exp = (); 
$indexnum = 0; 
while (<GEN>)	
{
	chomp $_;
	chop $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	my @cur_per_line = (); 
	@cur_per_line = split(',', $_);
	$indexnum = $indexnum + 1; 
	$gene_exp{$indexnum} = \@cur_per_line; 
}
close(GEN);



#--------------------- read in the gene expression list split file -------------------------------

my @splitlist = (); 
open(LIST, $expression_split_file) || die(" Can not open file(\"$expression_split_file\").\n"); 
while (<LIST>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
		
	my $curline = substr($_, 0, -1);
	push(@splitlist, $curline); ; 
}
close(LIST); 



#--------------------- Begin to generate ~20 gene expression feature  -------------------------------

open(INT, $pr_pair_file) || die(" Can not open file(\"$pr_pair_file\").\n"); 
open(OUT, "> $out_file_name") || die(" Can not open file(\"$out_file_name\").\n");

my $count = 0; 
my $feaCount = 0; 

my ( $orfs1, $orfs2, $orf1, $orf2, @exp1, @exp2, $flag); 
while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	($orfs1, $orfs2, $flag) = split('\s', $_);

	$orf1 = $gene_name{"$orfs1"} ; 
	$orf2 = $gene_name{"$orfs2"} ; 

	if (( ! (defined $orf1)  )||( ! (defined $orf2 )) )
	{
		for($feaCount = 0 ; $feaCount <= $#splitlist ; $feaCount++)
		{	print OUT "-100,";  }
		print OUT "$flag\n"; 
	}
	else {
		$exp1 = $gene_exp{$orf1};
		$exp2 = $gene_exp{$orf2};
	
		for($feaCount = 0 ; $feaCount < $#splitlist ; $feaCount++)
		{
			my $start = $splitlist[$feaCount] - 1; 
			my $endi = $splitlist[$feaCount + 1 ] - 1; 
			
			my @expf1 = (); 
			my @expf2 = ();  
			my $j=0; 
			for($j = $start ; $j < $endi ; $j++)
			{
				push(@expf1, $exp1->[$j]); 
				push(@expf2, $exp2->[$j]); 
			}
			
			my $temp =  &pearsoncc(\@expf1, \@expf2, $percent_miss); 
			print OUT "$temp,";  
		}
		my $start = $splitlist[$#splitlist] - 1; 
		my $endi = scalar @$exp2 - 1; 

		my @expf1 = (); 
		my @expf2 = ();  
		my $j=0; 
		for($j = $start ; $j < $endi ; $j++)
		{
			push(@expf1, $exp1->[$j]); 
			push(@expf2, $exp2->[$j]); 
		}

		my $temp =  &pearsoncc(\@expf1, \@expf2, $percent_miss); 
		print OUT "$temp,$flag\n";  
	}
	
	$count = $count + 1;
	
	if ( $count % 5000 == 0)
		{print "$count ";} 
}
print "\n\n ==> $count pairs !"; 

close(INT);					
close(OUT);


# ==================================================================================
# here we define the subroutine to get pearson correlation coefficients


sub pearsoncc 
{
   my(@a) = @{$_[0]};
   my(@b) = @{$_[1]};
   my($percent) = $_[2];

   my $i = 0;        
   my $revalue = 0; 
   
   my @af = (); 
   my @bf = (); 
   
   for($i = 0 ; $i <= $#a ; $i ++)
   {
	if (( $a[$i] < 90 ) && ( $b[$i] < 90))
	{
		push(@af, $a[$i]); 
		push(@bf, $b[$i]);
	}	
   }
   
   my $n = $#af + 1; 
   if (( $#af < $#a * ( 1- $percent )) || ($n <= 1))
   { 
   	$revalue =  -100;
   }    
   else 
   {
	   	my $sum_xy = 0; 
   		my $sum_x = 0; 
		my $sum_x2 = 0;    	
   		my $sum_y = 0; 
   		my $sum_y2 = 0; 
   	
   		for($i = 0 ; $i <= $#af ; $i ++)
   		{
   			$sum_xy = $sum_xy + $af[$i] * $bf[$i] ; 
   			$sum_x = $sum_x + $af[$i] ; 
			$sum_x2 = $sum_x2 + $af[$i] * $af[$i] ;
   			$sum_y = $sum_y + $bf[$i] ;
   			$sum_y2 = $sum_y2 + $bf[$i] * $bf[$i];  
   		}
		my $up = $sum_xy - $sum_x * $sum_y / $n ; 
        	my $low = ($sum_x2 - ($sum_x * $sum_x)/$n)*( $sum_y2 - ($sum_y * $sum_y)/$n);   	
   		$revalue = $up/sqrt($low + 0.00001); 
   } 
   return($revalue); 
}
