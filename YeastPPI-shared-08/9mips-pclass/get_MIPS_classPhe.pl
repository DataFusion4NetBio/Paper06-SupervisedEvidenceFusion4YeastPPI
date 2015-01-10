######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3

# Program to Extract Features of MIPS ProteinClass and Phenotype for protein pair list subset
# 
# MIPS Class: 
# -	ftp://ftpmips.gsf.de/yeast/catalogues/protein_classes/
# -	Two files: classes.scheme and protein_classes 
# $ perl get_MIPS_classPhe.pl ./lists/sciencesubset.txt classes.scheme protein_classes ./lists/sciencesubset.proClass > MIPSproteinClass.features.order
#  => Scheme First level: 25;   MIPS  : protein - 1017;
#"	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
#"	The out put file format:  25 -MIPS protein Class features, the last one is the class flag
#"	Each of the above features is a categorical feature with two possible features 0 or 1. 
#
# ----------------------------------------------------------------------------------------------
#
# MIPS mutant phenotype:  
# -	ftp://ftpmips.gsf.de/yeast/catalogues/phenotype/
# -	From Two files: phencat_data_07102004 and phencat.scheme
# $ perl get_MIPS_classPhe.pl ./lists/sciencesubset.txt phencat.scheme phencat_data_07102004 ./lists/sciencesubset.phecat > MIPS-Phencat.features.order
#"	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
#"	The out put file format:  11 -MIPS phencat features, the last one is the class flag
#"	Each of the above features is a categorical feature with two possible features 0 or 1. 
#
#  ==> in make the feature: 
#	if a pair's two protein has the same first levels in the MISP scheme  ...  pairs detected a 1
#	pairs not detected a 0
# 


use strict; 
die "Usage: command protein_pair_file MIPSscheme MIPTlabel out_file_name\n" if scalar(@ARGV) < 4;

my ($int_file, $scheme_file, $proteinLabel_file, $out_file) = @ARGV;

my ($orf1, $orf2, $flag); 
my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 

my %scheme = (); 


#--------------------- read in the MIPS scheme features -------------------------------

open(SYN, $scheme_file) || die(" Can not open file(\"$scheme_file\").\n"); 

while (<SYN>)	
{
	chomp $_; 
	next if /^\D/;				#ignore comments
	next if /^\s*$/; 			#ignore blank lines
	
	@orf = split('\s+', $_);
	
	if ($orf[0] =~ m/\./ )
		{}
	else {
		$orf1 = $orf[0]; 
		if (! defined $scheme{$orf1 } )
		{
			shift( @orf ); 
			$scheme{$orf1 } = join(" ",  @orf );
		}	
	}
}

close(SYN);


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

	if ($orf[1] =~ m/\./ )
	{
		my @tempp = split('\.', $orf[1] );
		$orf2 = $tempp[0]; 
	}
	else {
		$orf2 = $orf[1]; 
	}
	
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
my $synpairsize = scalar keys (%scheme);

print "MIPS  : protein - $proteinSize;   \n\n Scheme First level: $synpairsize; \n"; 
while ( my ($key, $value) = each(%scheme) ) 
	{    print "$key	$value\n"; 	}


#--------------------- Begin to process the yeast_int set and find if the pair in feature set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");


$count =  0; 

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
		for ( $i = 0; $i < $synpairsize; $i ++ )
		{
			print OUT "-100,"; 
		}	
		print OUT "$flag\n"; 
	}
	else { 
		for my $key ( keys %scheme ) 
		{
			if (( $t1 ->{"$key"} >0 ) && ( $t2 ->{"$key"} >0 ))
				{ print OUT "1,"; }
			else 
				{ print OUT "0,";}
		}
		print OUT "$flag\n";
	}
	$count = $count + 1; 
}

print " $count pairs ! ";

close(INT);					
close(OUT);