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
# This is a program to combine the features and labels files together into an integrated data set
# 
# This program is for the summary type of 10PPI feature set 
#
#!perl -w CombineFeaturesSummary.pl datasetPre outset

use strict; 
die "Usage: command InFilePre outfile \n" if scalar(@ARGV) < 2;
my ($inputPair, $out_file ) = @ARGV;

my $cmd1 = $inputPair.".genexpSumry" ; 
my $cmd2 = $inputPair.".tfSumry" ; 
my $cmd3 = $inputPair.".gofuncSumry " ; 
my $cmd3_2 = $inputPair.".goprocSumry " ; 
my $cmd3_3 = $inputPair.".gocompSumry " ;
my $cmd4 = $inputPair.".proexp" ; 
my $cmd5 = $inputPair.".ess" ; 
my $cmd6 = $inputPair.".msmatrix" ; 
my $cmd6_2 = $inputPair.".y2h" ; 
my $cmd7 = $inputPair.".gensyn" ; 
my $cmd8 = $inputPair.".natCopSeq" ; 
my $cmd9 = $inputPair.".proClassSumry" ; 
my $cmd10 = $inputPair.".phecatSumry" ; 
my $cmd11 = $inputPair.".besthits" ; 
my $cmd12 = $inputPair.".scHmPPI" ; 
my $cmd12_2 = $inputPair.".DmeHmPPI" ; 
my $cmd12_3 = $inputPair.".HsapiHmPPI" ; 
my $cmd12_4 = $inputPair.".CeHmPPI" ; 
my $cmd13 = $inputPair.".Y2Hddippi" ; 

open(F1, $cmd1) || die(" Can not open file(\"$cmd1\").\n"); 
open(F2, $cmd2) || die(" Can not open file(\"$cmd2\").\n"); 
open(F3, $cmd3) || die(" Can not open file(\"$cmd3\").\n"); 
open(F3_2, $cmd3_2) || die(" Can not open file(\"$cmd3_2\").\n"); 
open(F3_3, $cmd3_3) || die(" Can not open file(\"$cmd3_3\").\n"); 
open(F4, $cmd4) || die(" Can not open file(\"$cmd4\").\n"); 
open(F5, $cmd5) || die(" Can not open file(\"$cmd5\").\n"); 
open(F6, $cmd6) || die(" Can not open file(\"$cmd6\").\n"); 
open(F6_2, $cmd6_2) || die(" Can not open file(\"$cmd6_2\").\n"); 
open(F7, $cmd7) || die(" Can not open file(\"$cmd7\").\n"); 
open(F8, $cmd8) || die(" Can not open file(\"$cmd8\").\n"); 
open(F9, $cmd9) || die(" Can not open file(\"$cmd9\").\n"); 
open(F10, $cmd10) || die(" Can not open file(\"$cmd10\").\n"); 
open(F11, $cmd11) || die(" Can not open file(\"$cmd11\").\n"); 
open(F12, $cmd12) || die(" Can not open file(\"$cmd12\").\n"); 
open(F12_2, $cmd12_2) || die(" Can not open file(\"$cmd12_2\").\n"); 
open(F12_3, $cmd12_3) || die(" Can not open file(\"$cmd12_3\").\n"); 
open(F12_4, $cmd12_4) || die(" Can not open file(\"$cmd12_4\").\n"); 
open(F13, $cmd13) || die(" Can not open file(\"$cmd13\").\n"); 


open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n");

my ( $per_line, $line_num, @lines); 

$line_num = 0; 
my $fea_num = 0; 
my @temp = (); 

while (<F1>)
{
		$line_num = $line_num +1; 
		$fea_num = 0; 
		
		$per_line = $_; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		
		if ( $fea_num != 1 )
			{ print "Gene expression feature wrong: not $#temp features ! \n"; }; 
 
		$per_line = <F3>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 	
		if ( $fea_num != 2 )
			{ print "GO func feature wrong: not $#temp features ! \n"; }; 	
	
	
		$per_line = <F3_2>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 	
		if ( $fea_num != 3 )
			{ print "GO proc feature wrong: not $#temp features ! \n"; }; 	
	
	
		$per_line = <F3_3>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 		
		if ( $fea_num != 4 )
			{ print "GO comp feature wrong: not $#temp features ! \n"; }; 	


		$per_line = <F4>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 5 )
			{ print " Pro Exp feature wrong: not $#temp features ! \n"; }; 	


		$per_line = <F5>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 6 )
			{ print " Essential feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F6>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 8 )
			{ print " MS feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F6_2>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 9 )
			{ print " Y2H feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F7>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 10 )
			{ print " Synthetic feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F8>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 11 )
			{ print " Natu Comp sequence feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F11>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 12 )
			{ print " sequence besthits feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F12>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 13 )
			{ print " homology PPI SC feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F12_2>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 14 )
			{ print " homology PPI De feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F12_3>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 15 )
			{ print " homology PPI Hum feature wrong: not $#temp features ! \n"; }; 	
		

		$per_line = <F12_4>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 16 )
			{ print " homology PPI Ce feature wrong: not $#temp features ! \n"; }; 	

		

		$per_line = <F13>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 		
		if ( $fea_num != 17 )
			{ print " DDI Y2H feature wrong: not $#temp features ! \n"; }; 	
		
		$per_line = <F2>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 		 
		if ( $fea_num != 18 )
			{ print " TF Bind summary feature wrong: not $#temp features ! \n"; }; 	

		$per_line = <F9>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		pop(@temp); 
		print OUT join(',', @temp);  
		print OUT ",";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 19 )
			{ print " MIPS pro feature wrong: not $#temp features ! \n"; }; 	

		$per_line = <F10>; 
		chop($per_line); 
		@temp = split(/,/, $per_line); 
		print OUT join(',', @temp);  
		print OUT "\n";	
		$fea_num = $fea_num + $#temp + 1; 
		if ( $fea_num != 21 )
			{ print " MIPS phenotype wrong: not $#temp features ! \n"; }; 	
}
print "\n$line_num rows !   $fea_num features;"; 

close(F1);
close(F2);
close(F3);
close(F3_2);
close(F3_3);
close(F4);
close(F5);
close(F6);
close(F6_2);
close(F7);
close(F8);
close(F9);
close(F10);
close(F11);
close(F12);
close(F12_2);
close(F12_3);
close(F12_4);
close(F13);

close(OUT); 