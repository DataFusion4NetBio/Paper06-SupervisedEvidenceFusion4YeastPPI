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
# Program to Extract summary Feature of TF group based binding features for protein pair list
# 
# The Binding information is based on the transcription factor features published by 
# Transcriptional Regulatory Networks in Saccharomyces cerevisiae Lee et al. Science 298:799-804 (2002)
# 
# Note: 
# 
# - for each protein pair, we have a flag representing the class label: "1" means postive pair, "0" means random pairs
#"	The input pair list file format: ORF1 ORF2 Flag ( 0 rand or 1 postive)
#"	The out put file format:  1 summary tf feature separated by "," with the last one is the class flag
# 
# - Note: summary tf feature is a discrete feature: , for a specific thresholding of the P-value
#	the discrete value means how many TFs overall binds to both the two genes for a certain threshodling P-value
# 
# perl get_tfGroupBinding_summary.pl ./lists/sciencesubset.txt pvalbygene_nature04.txt 0.05 ./lists/sciencesubset.tfsummary
# 


use strict; 
die "Usage: command yeast_pr_pair_file_name tf_file threshold out_file_name\n" if scalar(@ARGV) < 4;

my ($pr_pair_file, $tf_file, $threshold, $out_file_name) = @ARGV;

my @perline = (); 

#--------------------- read in the tf file -------------------------------

open(GEN, $tf_file) || die(" Can not open file(\"$tf_file\").\n"); 
my ( %gene_tf, $indexorf ); 
%gene_tf = (); 

my $numfea = 0; 
while (<GEN>)	
{
	chomp $_;
	chop $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	@perline = split('\t', $_);
	$indexorf = $perline[0]; 

	my @cur_per_line = (); 
	my $j = 0; 
	for ( $j = 3; $j <= $#perline; $j ++)
	{
		push(@cur_per_line, $perline[$j]); 
	}
		
	$numfea = $#cur_per_line + 1; 
	$gene_tf{$indexorf} = \@cur_per_line; 	
}
close(GEN);

print "- $numfea TF factors! \n"; 




#--------------------- Begin to generate tf feature  -------------------------------

open(INT, $pr_pair_file) || die(" Can not open file(\"$pr_pair_file\").\n"); 
open(OUT, "> $out_file_name") || die(" Can not open file(\"$out_file_name\").\n");

my $count = 0; 
my $feaCount = 0; 
my $grpCount = 0;

my ( $orf1, $orf2, @exp1, @exp2, $flag, $exp1, $exp2); 
while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	($orf1, $orf2, $flag) = split('\s', $_);

	$exp1 = $gene_tf{$orf1};
	$exp2 = $gene_tf{$orf2};

	
	if (!($flag =~ /[0-1]/))
		{die(" $count: $orf1, $orf2, $flag : this pair flag is wrong ! "); }
	
	if ((!(defined ($exp1))) || (!(defined ($exp2))))
	{
		print OUT "-100, $flag\n"; 
	}
	else {
		my $temp = 0; 
		for($feaCount = 0 ; $feaCount < $numfea ; $feaCount++)
		{
			my $p1 = $exp1->[ $feaCount ]  ; 
			my $p2 = $exp2->[ $feaCount ]  ; 
			
			if (($p1 eq 'NaN') || ($p2 eq 'NaN'))
				{$temp = $temp + 0; }
			elsif (($p1 <= $threshold ) && ($p2 <= $threshold ))
				{$temp = $temp + 1; }
			else
				{$temp = $temp + 0; }
		}
		print OUT "$temp,";  
		print OUT "$flag\n";  
	}
	
	$count = $count + 1;
}
print "- Input $count pairs !"; 

close(INT);					
close(OUT);