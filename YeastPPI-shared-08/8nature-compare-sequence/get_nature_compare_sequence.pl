######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to Extract nature-compare-sequence based features for protein pair list subset
#
#-	What we use for this feature is the union of the above three methods derived pairs … totally 7446 interactions 
#-	      nature_compare_gene-fusion: 358 ( postive integer ) 
#-	      nature_compare_gene-neighborhood: 6387  ( postive integer ) 
#-	      nature_compare_gene-cooccurence: 997  ( 0/1 feature )
# 
# 
#  ==> in make the feature: 
#	all related proteins… pairs detected a 1
#	pairs not detected a 0
# 
# 
# -	$ perl get_nature_compare_sequence.pl ./lists/sciencesubset.txt nature-compare-gene-fusion.txt nature-compare-gene-neighborhood.txt nature-compare-gene-cooccurence.txt ./lists/sciencesubset.natCopSeq
# "	==> Nature-Compare-Sequence : protein - 1094; Pair: 7446;



use strict; 
die "Usage: command protein_pair_file gene_fusion gene_neighbor gene_cooccur out_file_name\n" if scalar(@ARGV) < 5;

my ($int_file, $fusion_file, $neighbor_file, $cooccur_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 

my %natCompSeque_orf = (); 
my %natCompSeque_protein = (); 


#--------------------- read in the nature_compare  features -------------------------------

open(NATF, $fusion_file) || die(" Can not open file(\"$fusion_file\").\n"); 

while (<NATF>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\s+', $_);
	
	$orf1 = $orf[0]; 
	$orf2 = $orf[1]; 

	my $temp = $orf1.":".$orf2; 

	$natCompSeque_orf{$temp } = 1; 
	$natCompSeque_protein{$orf1} = 1; 
	$natCompSeque_protein{$orf2} = 1; 		
}
close(NATF);


open(NATN, $neighbor_file) || die(" Can not open file(\"$neighbor_file\").\n"); 

while (<NATN>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\s+', $_);
	
	$orf1 = $orf[0]; 
	$orf2 = $orf[1]; 

	my $temp = $orf1.":".$orf2; 

	$natCompSeque_orf{$temp } = 1; 
	$natCompSeque_protein{$orf1} = 1; 
	$natCompSeque_protein{$orf2} = 1; 		
}
close(NATN);



open(NATC, $cooccur_file) || die(" Can not open file(\"$cooccur_file\").\n"); 

while (<NATC>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\s+', $_);
	
	$orf1 = $orf[0]; 
	$orf2 = $orf[1]; 

	my $temp = $orf1.":".$orf2; 

	$natCompSeque_orf{$temp } = 1; 
	$natCompSeque_protein{$orf1} = 1; 
	$natCompSeque_protein{$orf2} = 1; 		
}
close(NATC);




# ----------  print out some statistic number  ---------------------

my $proteinSize = scalar keys (%natCompSeque_protein); 
my $pairsize = scalar keys (%natCompSeque_orf);

print "Nature-Compare-Sequence : protein - $proteinSize; Pair: $pairsize; \n"; 


#--------------------- Begin to process the yeast_int set and find if the pair in nature_compare sequence feature set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");

my $count_natCompSeque = 0; 
$count =  0; 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\s', $_);

	$pair_l = $orf1.":".$orf2; 
	$pair_r = $orf2.":".$orf1;
	
	if ((! defined $natCompSeque_protein{$orf1}) || (! defined $natCompSeque_protein{$orf2}))
	{	$count_natCompSeque = 0;	}
	else 
	{
		if (defined $natCompSeque_orf{$pair_l}) 
			{$count_natCompSeque = 1 ; }
		elsif (defined $natCompSeque_orf{$pair_r})
			{$count_natCompSeque = 1 ; }
		else
			{$count_natCompSeque = 0 ; }			
	}
	
	print OUT "$count_natCompSeque,$flag\n";
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);