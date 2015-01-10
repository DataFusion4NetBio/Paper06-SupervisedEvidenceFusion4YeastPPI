######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to Extract Features of TF group based binding features for protein pair list subset
# 
# The TF grouping is based on the MIPS Class: 
# -	ftp://ftpmips.gsf.de/yeast/catalogues/protein_classes/
# -	Two files: classes.scheme and protein_classes 
# 
# The Binding information is based on the transcription factor features published by 
# Transcriptional Regulatory Networks in Saccharomyces cerevisiae Lee et al. Science 298:799-804 (2002)
# 
# Note: 
# 
# - for each protein pair, we have a flag representing the class label: "1" means postive pair, "0" means random pairs
#"	The input pair list file format: ORF1 ORF2 Flag ( 0 rand or 1 postive)
#"	The out put file format:  16 group tf features separated by ",", the last one is the class flag
# - Note: each tf feature is a discrete feature: , for a specific thresholding of the P-value
#	the discrete value means how many TFs in that group binds to both the two genes for a certain threshodling P-value
# 
# perl get_tfGroupBinding.pl ./lists/sciencesubset.txt pvalbygene_nature04.txt 204_pvalbygene_nature04_TFs.groupIndex 0.05 ./lists/sciencesubset.tfgroup
# 


use strict; 
die "Usage: command yeast_pr_pair_file_name tf_file tf_group_file threshold out_file_name\n" if scalar(@ARGV) < 5;

my ($pr_pair_file, $tf_file, $tf_group_file, $threshold, $out_file_name) = @ARGV;

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


#--------------------- read in the tf grouping file -------------------------------

open(GENGRP, $tf_group_file) || die(" Can not open file(\"$tf_group_file\").\n"); 

my @gene_tfgroup = (); 
my $numfeaGroup = 0; 

while (<GENGRP>)	
{
	chomp $_;
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	@perline = split('\t', $_);
	
	my $indexGroup = $perline[0]; 
	print "$indexGroup\t"; 

	shift( @perline ); 
	
	my @cur_per_line = @perline ; 
	$gene_tfgroup[ $numfeaGroup ] = \@cur_per_line; 
	
	$numfeaGroup = $numfeaGroup + 1; 
	
	my $size = $#cur_per_line + 1; 
	print "$size\n"; 	
}
close(GENGRP);
print "- $numfeaGroup TF groups !\n"; 



#--------------------- Begin to generate tf feature  -------------------------------

open(INT, $pr_pair_file) || die(" Can not open file(\"$pr_pair_file\").\n"); 
open(OUT, "> $out_file_name") || die(" Can not open file(\"$out_file_name\").\n");

my $count = 0; 
my $feaCount = 0; 
my $grpCount = 0;

my ( $orf1, $orf2, $flag, $exp1, $exp2); 
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
		for($grpCount = 0 ; $grpCount < $numfeaGroup ; $grpCount++)
		{	print OUT "-100,";  }
		print OUT "$flag\n"; 
	}
	else {
		for($grpCount = 0 ; $grpCount < $numfeaGroup ; $grpCount++)
		{
			my $p1 =0;
			my $p2 =0; 
			my $temp = 0; 
			my @curGrp = @{$gene_tfgroup[ $grpCount ]};

			for($feaCount = 0 ; $feaCount <= $#curGrp ; $feaCount++)
			{
				my $curIndex = $curGrp[ $feaCount ]; 
				$p1 = $exp1->[ $curIndex ]  ; 
				$p2 = $exp2->[ $curIndex ]  ; 
			
				if (($p1 eq 'NaN') || ($p2 eq 'NaN'))
					{	}
				elsif (($p1 <= $threshold ) && ($p2 <= $threshold ))
					{$temp = $temp + 1; }
				else
					{$temp = $temp + 0; }
			}
			print OUT "$temp,";  
		}
		print OUT "$flag\n";  
	}
	
	$count = $count + 1;
}
print "- Input $count pairs !"; 

close(INT);					
close(OUT);