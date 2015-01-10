######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to Extract Features of PPI mass specttrometry for protein pair list subset
# 
# Use MS data:    ==> Mass Spectrom
#
#	For HMS-PCI
#	 - nature02_analyze hms.matrix: 28252: protein - 1578
#	 - nature02_analyze hms.spoke: 3618:  protein - 1578
#	For TAP
#        - nature02_analyze tap.matrix: 18677: protein 1363
#        - nature02_analyze tap.spoke: 3225;  protein - 1363;

#
#    To get the MS features: 
#       -  All complexes pairng inputs 
#       -  All related Proteins 
#	all related proteins… pairs detected a 1
#	pairs not detected a 0
#	Others as UNKNOWN ( -100 ) 
#    But how to: 
#        Within a complex,  bait to others a 1
#        Within a complex,  others to others a 0 ( ?? )
#
# "	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
# "	The out put file format:  1-HMS , 1-TAP features, the last one is the class flag
#
#
#
#  $ perl get_PPIMS.pl ./lists/sciencesubset.txt nature-analyze-hmspci-spoke.txt nature-analyze-tap-spoke.txt ./lists/sciencesubset.msspoke
#  $ perl get_PPIMS.pl ./lists/sciencesubset.txt nature-analyze-hmspci-matrix.txt nature-analyze-tap-matrix.txt ./lists/sciencesubset.msmatrix


use strict; 
die "Usage: command protein_pair_file nature_analyze_hms nature_analyze_tap out_file_name\n" if scalar(@ARGV) < 4;

my ($int_file, $hms_file, $tap_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 

my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 

#--------------------- read in the nature_analyze MS features -------------------------------

my %hms_orf = (); 
my %hms_protein = (); 

open(HMS, $hms_file) || die(" Can not open file(\"$hms_file\").\n"); 

while (<HMS>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\t', $_);
	$orf1 = $orf[0]; 
	$orf2 = substr($orf[1], 0, -1); 

	my $temp = $orf1.":".$orf2; 

	$hms_orf{$temp } = 1; 
	
	if (!(defined $hms_protein{$orf1}))
		{ $hms_protein{$orf1} = 1; }
	if (!(defined $hms_protein{$orf2}))
		{ $hms_protein{$orf2} = 1; }		
}

close(HMS);


#--------------------- read in the another MS  features -------------------------------

my %tap_orf = (); 
my %tap_protein = (); 


open(TAP, $tap_file) || die(" Can not open file(\"$tap_file\").\n"); 

while (<TAP>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\t', $_);
	$orf1 = $orf[0]; 
	$orf2 = substr($orf[1], 0, -1);


	my $temp = $orf1.":".$orf2; 

	$tap_orf{$temp } = 1; 
	
	if (!(defined $tap_protein{$orf1}))
		{ $tap_protein{$orf1} = 1; }
	if (!(defined $tap_protein{$orf2}))
		{ $tap_protein{$orf2} = 1; }		
}

close(TAP);



my $hmsproteinSize = scalar keys (%hms_protein); 
my $hmspairsize = scalar keys (%hms_orf);
my $tapproteinSize = scalar keys (%tap_protein); 
my $tappairsize = scalar keys (%tap_orf);

print "HMS: protein - $hmsproteinSize; Pair: $hmspairsize; \n"; 
# ==> HMS: protein - 4018; Pair: 7389;
print "TAP: protein - $tapproteinSize; Pair: $tappairsize; \n"; 

#--------------------- Begin to process the yeast_int set and find if the pair in nature_analzye feature set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");

my $count_hms = 0; 
my  $count_tap = 0; 
$count =  0; 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\s', $_);

	$pair_l = $orf1.":".$orf2; 
	$pair_r = $orf2.":".$orf1;
	
	if ((! defined $hms_protein{$orf1}) || (! defined $hms_protein{$orf2}))
	{	$count_hms = -100;	}
	else 
	{
		if (defined $hms_orf{$pair_l}) 
			{$count_hms = 1 ; }
		elsif (defined $hms_orf{$pair_r})
			{$count_hms = 1 ; }
		else
			{$count_hms = 0 ; }			
	}

	if ((! defined $tap_protein{$orf1}) || (! defined $tap_protein{$orf2}))
	{	$count_tap = -100;	}
	else 
	{
		if (defined $tap_orf{$pair_l}) 
			{$count_tap = 1 ; }
		elsif (defined $tap_orf{$pair_r})
			{$count_tap = 1 ; }
		else
			{$count_tap = 0 ; }			
	}
	
	print OUT "$count_hms,$count_tap,$flag\n";
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);