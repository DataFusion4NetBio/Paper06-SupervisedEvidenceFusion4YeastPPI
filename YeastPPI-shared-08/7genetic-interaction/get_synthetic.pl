######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to Extract Features of synthetic lethal for protein pair list subset
# 
# Use two synthetic data: 
#"	Nature-compare-synthetic: 886 pairs
#"	Tong2004: 4624 pair
# together: Synthetic : protein - 1413; Pair: 5329;
# 
#  ==> in make the feature: 
#	all related proteins… pairs detected a 1
#	pairs not detected a 0
# 
#  -	$ perl get_synthetic.pl ./lists/sciencesubset.txt tong2004.tab nature-compare-synthetic-lethality.txt ./lists/sciencesubset.gensyn

use strict; 
die "Usage: command protein_pair_file Tong2004 Nature-compare-synthetic out_file_name\n" if scalar(@ARGV) < 4;

my ($int_file, $tong_file, $compare_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ($count_synthetic, $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 

my %syn_orf = (); 
my %syn_protein = (); 


#--------------------- read in the nature_analyze Y2H features -------------------------------

open(SYN, $tong_file) || die(" Can not open file(\"$tong_file\").\n"); 

while (<SYN>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\s+', $_);
	
	$orf1 = $orf[0]; 
	$orf2 = substr($orf[1], 7); 

	my $temp = $orf1.":".$orf2; 

	$syn_orf{$temp } = 1; 
	
	if (!(defined $syn_protein{$orf1}))
		{ $syn_protein{$orf1} = 1; }
	if (!(defined $syn_protein{$orf2}))
		{ $syn_protein{$orf2} = 1; }		
}

close(SYN);



#--------------------- read in the nature_compare_synthetic  hits features -------------------------------

open(NATSYN, $compare_file) || die(" Can not open file(\"$compare_file\").\n"); 

while (<NATSYN>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\s+', $_);
	
	$orf1 = $orf[0]; 
	$orf2 = $orf[1]; 

	my $temp = $orf1.":".$orf2; 

	$syn_orf{$temp } = 1; 
	
	if (!(defined $syn_protein{$orf1}))
		{ $syn_protein{$orf1} = 1; }
	if (!(defined $syn_protein{$orf2}))
		{ $syn_protein{$orf2} = 1; }			
}

close(NATSYN);

my $proteinSize = scalar keys (%syn_protein); 
my $synpairsize = scalar keys (%syn_orf);

print "Synthetic : protein - $proteinSize; Pair: $synpairsize; \n"; 
# ==> 

#--------------------- Begin to process the yeast_int set and find if the pair in nature_analzye feature set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");


$count =  0; 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\s', $_);

	$pair_l = $orf1.":".$orf2; 
	$pair_r = $orf2.":".$orf1;
	
	if ((! defined $syn_protein{$orf1}) || (! defined $syn_protein{$orf2}))
	{	$count_synthetic = -100;	}
	else 
	{
		if (defined $syn_orf{$pair_l}) 
			{$count_synthetic = 1 ; }
		elsif (defined $syn_orf{$pair_r})
			{$count_synthetic = 1 ; }
		else
			{$count_synthetic = 0 ; }			
	}
	
	print OUT "$count_synthetic,$flag\n";
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);