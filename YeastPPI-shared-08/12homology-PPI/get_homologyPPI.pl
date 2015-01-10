######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3

# Program to Extract Features of Homology PPI  for protein pair list subset
# 
# Two data used: 
# - 1. SGD PSI-BLAST sequence comparison result
# - 2. Other species' PPI features from DIP  
#
# -----------------------------------------------------------------------------------------------
#  ==> in make the feature: 
#"	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
#"	The out put file format:  1 -discrete value feature, the last one is the class flag
#"	This feature most time is 0 or 1… 
#
# -----------------------------------------------------------------------------------------------
#  Currently we are using the following four species to derive the homology PPI for the yeast protein pairs. 
#"	S. cerevisiae(Baker's yeast)
#"	C. elegans
#"	D. melanogaster (fruit fly)
#"	H. sapiens (Human)
# -----------------------------------------------------------------------------------------------
# $ perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.sc.tab ./dip-ppi/Scere20041003.tab ./lists/sciencesubset.scHmPPI
# $ perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.dmelanogaster.tab ./dip-ppi/Dmela20041003.tab ./lists/sciencesubset.DmeHmPPI
# $ perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.humans.tab ./dip-ppi/Hsapi20041003.tab ./lists/sciencesubset.HsapiHmPPI
# $ perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.celegan.tab ./dip-ppi/Celeg20041003.tab ./lists/sciencesubset.CeHmPPI


use strict; 
die "Usage: command protein_pair_file PBlastFile SpeciePPI out_file_name\n" if scalar(@ARGV) < 4;

my ($int_file, $homology_file, $hPPI_file, $out_file) = @ARGV;


my ($orf1, $orf2, $flag); 
my ( $pair_l, $pair_r, $count);
my ( @orf, $per_orf ); 



#--------------------- read in the homology PSI-BLAST SGD file  -------------------------------

my %homology = (); 

open(HGY, $homology_file) || die(" Can not open file(\"$homology_file\").\n"); 

while (<HGY>)	
{
	chomp $_; 
	next if /^\s*$/; 			#ignore blank lines
	
	@orf = split('\t', $_);
	
	$orf1 = $orf[0]; 
	my $eValue = $orf[6]; 
	my $temp = $orf[7];
	
	if( $temp =~ m/gi\|(\d+)/) 
	{	$orf2 = $1; 	}	
	
	if (! defined $homology{ "$orf1" } )
	{
		my @cur = (); 
		$cur[0] = $orf2 ; 
		$homology{ "$orf1" } = \@cur; 
	}	
	else {
		push( @ {$homology{ "$orf1" }}, $orf2); 	
	}
}

close(HGY);



#--------------------- read in the DIP some other species PPI file  -------------------------------

my %speciesppr = (); 

open(PPI, $hPPI_file) || die(" Can not open file(\"$hPPI_file\").\n"); 

while (<PPI>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	
	@orf = split('\t', $_);
	
	$orf1 = $orf[5]; 
	$orf2 = $orf[10]; 

	my $gi1 = ""; 
	my $gi2 = ""; 
	
	if( $orf1 =~ m/GI:(\d+)/) 
	{	$gi1 = $1; 	}
	if( $orf2 =~ m/GI:(\d+)/) 
	{	$gi2 = $1; 	}
	
	my $pair = $gi1.":".$gi2; 
	$speciesppr{ "$pair" } = 1; 
}

close(PPI);


my $proteinSize = scalar keys (%homology); 
my $pairsize = scalar keys (%speciesppr);

print "\n\n- Yeast SGD PSI-BLAST file: $homology_file related yeast proteins - $proteinSize;  \n\n- The other species $hPPI_file PPI pairs: $pairsize; \n"; 


#--------------------- Begin to process the yeast_int set and find if the pair in feature set -------------------------------

open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");


$count =  0; 
my $score = 0; 
my $homCatchNum = 0; 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\s', $_);

	$score = 0; 
	if ( (defined $homology{"$orf1"} ) &&  (defined $homology{"$orf2"} ))
	{
		my @proteins1 = @{$homology{"$orf1"}}; 
		my @proteins2 = @{$homology{"$orf2"}}; 	
		my ($curpro1, $curpro2); 
	
		$score = 0; 
		foreach $curpro1 (@proteins1) {
			foreach $curpro2 (@proteins2) 
			{
				my $pair_l = $curpro1.":".$curpro2; 
				my $pair_r = $curpro2.":".$curpro1; 
			
				if (defined $speciesppr{"$pair_l"} )
					{ $score = $score + 1;  
					  $homCatchNum = $homCatchNum + 1; }
				elsif (defined $speciesppr{"$pair_r"} )
					{ $score = $score + 1;  
					  $homCatchNum = $homCatchNum + 1; }
				else 
					{ $score = $score + 0;  }
			}
		}
	}
	print OUT "$score,$flag\n"; 
	$count = $count + 1; 
}

print "\n- Input Yeast protein pair file: $count pairs ! ; \n- $homCatchNum has homology PPI in $hPPI_file; ";

close(INT);					
close(OUT);