######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3

# Program to Extract SGD BestHits Similarity Features for protein pair list subset
# 
# Use :  ftp://genome-ftp.stanford.edu/pub/yeast/data_download/sequence_similarity/
#   - 	best_hits.tab ==>  	best_hits.sc.tab  ( ORF vs. ORF)   All BLASTP hits with E-value less than or equal to 0.01 are shown.
#
#
#"	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
#"	The out put file format:  1 - blastP E value feature, the last one is the class flag
#"	For all E value < 0.01, the feature is 0; For all others, the feature is the nature log of the E-value, which means the smaller the better
#
# perl get_scBestHits.pl ./lists/sciencesubset.txt best_hits.sc.tab ./lists/sciencesubset.besthits
 

use strict; 
die "Usage: command protein_pair_file sgd_bestHits_sc out_file_name\n" if scalar(@ARGV) < 3;

my ($int_file, $info_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 


#--------------------- read in the info features -------------------------------

my %info = (); 

open(INFO, $info_file) || die(" Can not open file(\"$info_file\").\n"); 

while (<INFO>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\t', $_);
	$orf1 = $orf[0]; 
	$orf2 = $orf[7]; 
	
	# we would use the negative natural log of aligned E-value 
	my $evalue = - $orf[6];  
	
	my $temp = $orf1.":".$orf2; 

	$info{$temp } = $evalue; 
}

close(INFO);

my $pairsize = scalar keys (%info);
print "Info Feature Pair: $pairsize; \n"; 




#--------------------- Begin to process feature set -------------------------------

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

	$pair_l = $orf1.":".$orf2; 
	$pair_r = $orf2.":".$orf1;
	
	if (defined $info{$pair_l}) 
		{$score =  $info{$pair_l}; }
	elsif (defined $info{$pair_r})
		{$score = $info{$pair_r} ; }
	else
		{$score = 0 ; }			
	
	print OUT "$score,$flag\n";
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);