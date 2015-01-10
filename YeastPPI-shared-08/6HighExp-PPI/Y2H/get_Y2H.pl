######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to Extract Features of Y2H for protein pair list subset
# 
# Use two Y2H data: 
#"	- nature02_analyze y2h: 5614
#	- sgdlite-y2h_hits.txt: 2032 pairs
#
# Nov. 20. 2003 
# Changed for considering missing value and positive / negative case of the pairs within experiments result
#      ==>Y2H 
#	all related proteins… pairs detected a 1
#	pairs not detected a 0
#	Others as UNKNOWN ( -100 ) 
#
# perl get_Y2H.pl ./lists/sciencesubset.txt nature-analyze-y2h.txt sgdlite-y2h_hits.txt ./lists/sciencesubset.y2h

use strict; 
die "Usage: command protein_pair_file nature_analyze_y2h sgdlite-y2h_hits.txt out_file_name\n" if scalar(@ARGV) < 4;

my ($int_file, $y2h_file, $y2hsgd_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ($count_y2h, $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 

my %y2h_orf = (); 
my %y2h_protein = (); 


#--------------------- read in the nature_analyze Y2H features -------------------------------

open(Y2H, $y2h_file) || die(" Can not open file(\"$y2h_file\").\n"); 

while (<Y2H>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\t', $_);
	$orf1 = $orf[0]; 
	$orf2 = substr($orf[1], 0, -1); 

	my $temp = $orf1.":".$orf2; 

	$y2h_orf{$temp } = 1; 
	
	if (!(defined $y2h_protein{$orf1}))
		{ $y2h_protein{$orf1} = 1; }
	if (!(defined $y2h_protein{$orf2}))
		{ $y2h_protein{$orf2} = 1; }		
}

close(Y2H);



#--------------------- read in the sgd Y2H hits features -------------------------------

open(Y2HS, $y2hsgd_file) || die(" Can not open file(\"$y2hsgd_file\").\n"); 

while (<Y2HS>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\t', $_);
	$orf1 = $orf[0]; 
	$orf2 = $orf[1]; 

	my $temp = $orf1.":".$orf2; 

	$y2h_orf{$temp } = 1; 
	
	if (!(defined $y2h_protein{$orf1}))
		{ $y2h_protein{$orf1} = 1; }
	if (!(defined $y2h_protein{$orf2}))
		{ $y2h_protein{$orf2} = 1; }		
}

close(Y2HS);

my $proteinSize = scalar keys (%y2h_protein); 
my $y2hpairsize = scalar keys (%y2h_orf);

print "Y2H: protein - $proteinSize; Pair: $y2hpairsize; \n"; 
# ==> Y2H: protein - 4018; Pair: 7389;

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
	
	if ((! defined $y2h_protein{$orf1}) || (! defined $y2h_protein{$orf2}))
	{	$count_y2h = -100;	}
	else 
	{
		if (defined $y2h_orf{$pair_l}) 
			{$count_y2h = 1 ; }
		elsif (defined $y2h_orf{$pair_r})
			{$count_y2h = 1 ; }
		else
			{$count_y2h = 0 ; }			
	}
	
	print OUT "$count_y2h,$flag\n";
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);