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
# Program to Extract summary Features of MIPS ProteinClass and Phenotype for protein pair list subset
# Note: there is another similar program to extract the detailed features of MIPS prClass and Phenotetype
# 
# MIPS Class: 
# -	ftp://ftpmips.gsf.de/yeast/catalogues/protein_classes/
# -	Two files: classes.scheme and protein_classes 
# $ perl  get_MIPS_classPhe_summary.pl ./lists/temp.posFinal protein_classes ./lists/temp.posMIPSclassSum.csv 
# 
#  => Scheme First level: 25;   MIPS  : protein - 1017;
#"	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
#"	The out put file format:  1 summary MIPS protein Class features, the last one is the class flag
#
# ----------------------------------------------------------------------------------------------
#
# MIPS mutant phenotype:  
# -	ftp://ftpmips.gsf.de/yeast/catalogues/phenotype/
# -	From Two files: phencat_data_07102004 and phencat.scheme
#
# $ perl get_MIPS_classPhe_summary.pl ./lists/sciencesubset.txt phencat_data_07102004 ./lists/sciencesubset.phecatSummary
#"	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
#"	The out put file format:  1 summary MIPS phencat feature, the last one is the class flag
#
#  ==> in make the feature: 
#	if a pair's two protein has the same first levels in the MISP scheme  ...  pairs detected a 1
#	pairs not detected a 0
# 


use strict; 
die "Usage: command protein_pair_file MIPSlabel out_file_name\n" if scalar(@ARGV) < 3;

my ($int_file, $proteinLabel_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 

my %protein = (); 

#--------------------- read in the MIPS protein label features -------------------------------

open(PRO, $proteinLabel_file) || die(" Can not open file(\"$proteinLabel_file\").\n"); 

while (<PRO>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	@orf = split('\|', $_);
	
	$orf1 = uc($orf[0]); 
	$orf2 = $orf[1]; 
	
	if (defined $protein{ "$orf1" })
	{
		$ { $protein{"$orf1"} } {"$orf2"} = 1; 		
	}
	else {
		my %cur = (); 		
		$cur{"$orf2"} = 1;
		$protein{"$orf1" } = \%cur; 
	}
}

close(PRO);

my $proteinSize = scalar keys (%protein); 
print "MIPS  : protein labeled here:  $proteinSize;\n\n"; 


#--------------------- Begin to process the yeast_int set and find if the pair in feature set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");

$count =  0; 
my ( $temp1, $temp2 ); 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\s', $_);

	my $t1 = $protein{"$orf1"}; 
	my $t2 = $protein{"$orf2"}; 	
	my $i = 0; 
	
	if ((! defined $t1) || (! defined $t2 ))
	{	
		print OUT "-100,"; 
	}
	else { 
		my @p1Schemes = keys %{$t1};  
		my @p2Schemes = keys %{$t2};  
		my $levelScore = 0; 
		
		foreach $temp1 (@p1Schemes)
		{
			foreach $temp2 (@p2Schemes)
			{
				my @levels1 = split('\.', $temp1);
				my @levels2 = split('\.', $temp2);
		
				my $checkSize = scalar @levels1; 
				if ( $checkSize > scalar @levels2 )
					{ $checkSize = scalar @levels2; 	}
				my $t; 
		
				for ($t = 0; $t < $checkSize; $t ++ )
				{ 
					if ( $levels1[$t] == $levels2[$t] )
					{	$levelScore = $levelScore + 1; }
					else 
					{	$levelScore = $levelScore + 0;
						last; 
					}
				}
			}
		}
		print OUT "$levelScore,";
	}
	
	print OUT "$flag\n"; 	
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);