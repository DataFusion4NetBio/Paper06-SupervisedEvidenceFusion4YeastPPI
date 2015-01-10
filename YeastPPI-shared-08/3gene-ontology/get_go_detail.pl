######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# Program to Extract gene ontology feature for protein pair list  
#
# perl get_go_detail.pl temp.line.txt go_slim_mapping04.tab go_slim_function04.txt go_slim_process04.txt go_slim_component04.txt temp.gofunc.txt temp.goproc.txt temp.gocomp.txt
#
# Nov. 2004 Version
# In this verion, we only have the second level categories for FUNC, PROC, COMP
# "	The input pair list file format: "ORF1 ORF2 Flag" ( 0 rand or 1 postive)
# "	The out put file format:  go features separated by ",", the last one is the class flag
#
# Found some component name not in the go_slim_component.txt in this new go_slim_mapping.tab (2003.11.5 version)
#1.	  we have 21 go function features
#2.	  we have 33 go process features
#3.	  we have 23 go component features


use strict; 
die "Usage: command protein_pair_list yeast_go_map go_function go_process go_component out_funcfile_name out_procfile_name out_compfile_name\n" if scalar(@ARGV) < 8;

my ($int_file, $map_file, $func_file, $process_file, $com_file, $out_funcfile, $out_procfile, $out_compfile) = @ARGV;


#--------------------- read in the go_slim_yeast_mapping  -------------------------------

open(MAP, $map_file) || die(" Can not open file(\"$map_file\").\n"); 

my (%yeast_map_function, %yeast_map_process, %yeast_map_component); 
my ($orf, $ontology, $ontology_cont); 

%yeast_map_function = (); 
%yeast_map_process = ();
%yeast_map_component = ();

while (<MAP>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	my @map_line = split('\t', $_);

	$orf = $map_line[0]; 
	$ontology =  $map_line[3];  
	$ontology_cont = $map_line[4]; 
	
	my %ontl_cont_ary = (); 
	
	if ($ontology eq 'F')	{
		if (defined $yeast_map_function{"$orf"})
		{
			${ $yeast_map_function{"$orf"} }{ $ontology_cont } =   1; 
		}
		else 
		{
			$ontl_cont_ary{$ontology_cont} = 1; 
			$yeast_map_function{"$orf"} = \%ontl_cont_ary; 
		}
	}
	elsif ($ontology eq 'P') {
		if (defined $yeast_map_process{"$orf"})
		{
			${ $yeast_map_process{"$orf"} }{ $ontology_cont } =   1; 			
		}
		else 
		{
			$ontl_cont_ary{$ontology_cont} = 1; 
			$yeast_map_process{"$orf"} = \%ontl_cont_ary; 
		}
	}
	elsif ($ontology eq 'C') {
		if (defined $yeast_map_component{"$orf"})
		{
			${ $yeast_map_component{"$orf"} }{ $ontology_cont } =   1; 
		}
		else 
		{
			$ontl_cont_ary{$ontology_cont} = 1; 
			$yeast_map_component{"$orf"} = \%ontl_cont_ary; 
		}
	}	
	else {
		die(" Wrong format go mapping file (\"$map_file\").\n");
	}
}

close(MAP);


#--------------------- read in the three go_slim_yeast ontologies  -------------------------------

my ( @line, $target, %func_2nd, %process_2nd, %comp_2nd, $func_num, $proc_num, $comp_num); 

@line = (); 
%func_2nd =  %process_2nd =  %comp_2nd = ();
$func_num = $proc_num = $comp_num = 0; 

open(FUNC, $func_file) || die(" Can not open file(\"$func_file\").\n");
while (<FUNC>)	
{
	chomp $_; 
	next if /^!/;			#ignore comments !
	next if /^$/; 			#ignore blank lines

	my @line = split(/;/, $_);
	chop($line[0]); 
	$target = $line[0]; 
	
	if (defined $func_2nd{"$target"})
		{ $func_2nd{"$target"} = $func_2nd{"$target"} + 1; }
	else 
		{ $func_2nd{"$target"} =  1; }
}
close(FUNC);

$func_num = scalar keys %func_2nd; 
print "\n -----   $func_num yeast slim functions !   ---- \n"; 
while ( my ($key, $value) = each(%func_2nd) ) {
        print "$key => $value\n";
}


open(PROCS, $process_file) || die(" Can not open file(\"$process_file\").\n");
while (<PROCS>)	
{
	chomp $_; 
	next if /^!/;			#ignore comments !
	next if /^$/; 			#ignore blank lines

	my @line = split(/;/, $_);
	chop($line[0]); 
	$target = $line[0]; 


	if (defined $process_2nd{"$target"})
		{ $process_2nd{"$target"} = $process_2nd{"$target"} + 1; }
	else 
		{ $process_2nd{"$target"} =  1; }
}
close(PROCS);
$proc_num = scalar keys %process_2nd; 

print "\n -----   $proc_num yeast slim processes !   ---- \n"; 
while ( my ($key, $value) = each(%process_2nd) ) {
        print "$key => $value\n";
}



open(COMP, $com_file) || die(" Can not open file(\"$com_file\").\n");
while (<COMP>)	
{
	chomp $_; 
	next if /^!/;			#ignore comments !
	next if /^$/; 			#ignore blank lines

	my @line = split(/;/, $_);
	chop($line[0]); 
	$target = $line[0]; 


	if (defined $comp_2nd{"$target"})
		{ $comp_2nd{"$target"} = $comp_2nd{"$target"} + 1; }
	else 
		{ $comp_2nd{"$target"} =  1; }
}
close(COMP);
$comp_num = scalar keys %comp_2nd; 

print "\n -----   $comp_num yeast slim components !   ---- \n"; 
while ( my ($key, $value) = each(%comp_2nd) ) {
        print "$key => $value\n";
}


#--------------------- Begin to process the yeast_int set and find if the pair in go_yeast_slim three ontologies -------------------------------
# In our 2003 version, we used the GO-slim hierarchy as up to three levels.  
# Here we only use yeast_slim all as the second levels 
# For each of them, we would have one feature. The feature would be a 2-value category variable... 
# 1 means the pair genes both have this func/comp/proc. 0 means otherwise
# 	- Special: if have the same component ontology term, but not in the same go_comp list. still give feature 1


open(INT, $int_file) || die(" Can not open file(\"$int_file\").\n"); 
open(OUTF, "> $out_funcfile") || die(" Can not open file(\"$out_funcfile\").\n");
open(OUTP, "> $out_procfile") || die(" Can not open file(\"$out_procfile\").\n");
open(OUTC, "> $out_compfile") || die(" Can not open file(\"$out_compfile\").\n");


my ( $orf1,  $orf2, $flag);  
my (%func1, %func2, %proc1, %proc2, %comp1, %comp2); 
my ($temp1, $temp2); 

while (<INT>)	
{
	chomp $_; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
	($orf1, $orf2, $flag) = split('\t', $_);
	my ( $i, $j, $p); 		

   	my $temp_score = 0; 
   	if ((defined $yeast_map_function{"$orf1"}) && (defined $yeast_map_function{"$orf2"})) 
   	{
			%func1 = %{ $yeast_map_function{"$orf1"} }; 
			%func2 = %{ $yeast_map_function{"$orf2"} }; 

			for my $key ( keys %func_2nd ) 
			{
				if ($key ne "molecular_function unknown")
				{
					$temp1 = $func1{$key}; 
					$temp2 = $func2{$key}; 
					if (($temp1 > 0)&&($temp2 > 0))
						{ print OUTF "1,"; }
					else 
						{ print OUTF "0,";}
				}
    			}
   	} 
     	else {
     		for ( $i = 0; $i < $func_num - 1 ; $i++ )
			{ print OUTF "-100,"; }
     	}
   	print OUTF "$flag\n"; 
   

   	if ((defined $yeast_map_process{"$orf1"}) && (defined $yeast_map_process{"$orf2"})) 
   	{	
			%proc1 = %{ $yeast_map_process{"$orf1"} }; 
			%proc2 = %{ $yeast_map_process{"$orf2"} }; 

			for my $key ( keys %process_2nd ) 
			{
				if ($key ne "biological_process unknown")
				{
					$temp1 = $proc1{$key}; 
					$temp2 = $proc2{$key}; 
					if (($temp1 > 0)&&($temp2 > 0))
						{ print OUTP "1,"; }
					else 
						{ print OUTP "0,";}
				}
    			}
   	} 
     	else {
     		for ( $i = 0; $i < $proc_num - 1 ; $i++ )
			{ print OUTP "-100,"; }
     	}
   	print OUTP "$flag\n"; 

   
   	if ((defined $yeast_map_component{"$orf1"}) && (defined $yeast_map_component{"$orf2"})) 
   	{	
			%comp1 = %{ $yeast_map_component{"$orf1"} }; 
			%comp2 = %{ $yeast_map_component{"$orf2"} }; 

			for my $key ( keys %comp_2nd ) 
			{
				if ($key ne "cellular_component unknown")
				{
					$temp1 = $comp1{$key}; 
					$temp2 = $comp2{$key}; 
					if (($temp1 > 0)&&($temp2 > 0))
						{ print OUTC "1,"; }
					else 
						{ print OUTC "0,";}
				}
    			}
   	} 
     	else {
     		for ( $i = 0; $i < $comp_num - 1 ; $i++ )
			{ print OUTC "-100,"; }
     	}
   	print OUTC "$flag\n";    
   
}

close(int); 
close(OUTF);
close(OUTP);
close(OUTC);
