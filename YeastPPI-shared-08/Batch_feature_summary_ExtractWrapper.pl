######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3

# This program is a yeast PPI summary-style feature extraction wrapper 
# 
# for the summary version  - Attention: we would not re-generate those features that are the same as the detailed version
#
# perl command inputPairlist

use strict; 
die "Usage: command inputPairFile \n" if scalar(@ARGV) < 1;
my ($inputPair ) = @ARGV;


print "\n-------------------------- 1gene-expression  summary   -----------------------------------------\n"; 

# -------------------   1gene-expression  summary ------------------------------
# # perl get_gene_expression_summary.pl ./lists/Science03PPI.sublist.txt YeastGeneListOrfGeneName-106_pval_v9.0.txt all_expression_fixed_s4_csv.txt 0.7 ./lists/Science03PPI.sublist.genexpsum

my $cmdPre = "perl ./1gene-expression/get_gene_expression_summary.pl  "; 
my $cmdPro = "./1gene-expression/YeastGeneListOrfGeneName-106_pval_v9.0.txt ./1gene-expression/all_expression_fixed_s4_csv.txt  0.6 "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".genexpSumry" ; 
print "$cmd\n"; 
system($cmd); 



print "\n------------------------------2tf-binding  summary -------------------------------------\n"; 

# -------------------   2tf-binding  summary ------------------------------

# perl get_tfGroupBinding_summary.pl ./lists/sciencesubset.txt pvalbygene_nature04.txt 0.05 ./lists/sciencesubset.tfsummary

my $cmdPre = "perl ./2tf-binding/get_tfGroupBinding_summary.pl  "; 
my $cmdPro = " ./2tf-binding/pvalbygene_nature04.txt  0.05 "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".tfSumry" ; 
print "$cmd\n"; 
system($cmd); 




print "\n------------------------------ 3gene-ontology summary   -------------------------------------\n"; 

# -------------------   3gene-ontology summary   ------------------------------

# perl get_go_summary.pl temp.line.txt go_slim_mapping04.tab go_slim_function04.txt go_slim_process04.txt go_slim_component04.txt temp.gofuncSum.txt temp.goprocSum.txt temp.gocompSum.txt

my $cmdPre = "perl ./3gene-ontology/get_go_summary.pl  "; 
my $cmdPro = " ./3gene-ontology/go_slim_mapping04.tab ./3gene-ontology/go_slim_function04.txt ./3gene-ontology/go_slim_process04.txt ./3gene-ontology/go_slim_component04.txt "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".gofuncSumry ".$inputPair.".goprocSumry ".$inputPair.".gocompSumry " ;

print "$cmd\n"; 
system($cmd); 




# -------------------   4protein-expression  ------------------------------


# -------------------   5essentiality  ------------------------------


# -------------------  6HighExp-PPI   ------------------- 


# -------------------   7genetic-interaction  ------------------------------


# -------------------   8nature-compare-sequence  ------------------------------


# -------------------   9mips-pclass summary ------------------------------

# $ perl  get_MIPS_classPhe_summary.pl ./lists/temp.posFinal protein_classes ./lists/temp.posMIPSclassSum.csv 

my $cmdPre = "perl ./9mips-pclass/get_MIPS_classPhe_summary.pl  "; 
my $cmdPro = "   ./9mips-pclass/protein_classes "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".proClassSumry" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   10mips-phenotype  ------------------------------
#
# perl get_MIPS_classPhe_summary.pl ./lists/sciencesubset.txt phencat_data_07102004 ./lists/sciencesubset.phecatSummary

my $cmdPre = "perl ./10mips-phenotype/get_MIPS_classPhe_summary.pl  "; 
my $cmdPro = "  ./10mips-phenotype/phencat_data_07102004 "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".phecatSumry" ; 
print "$cmd\n"; 
system($cmd); 



# -------------------   11sequence-similarity  ------------------------------


# -------------------   12homology-PPI ------------------------------


# -------------------   13domain-interaction  ------------------------------


# -------------------   combine features into one set  ------------------------------
# #!perl -w CombineFeaturesSummary.pl datasetPre outset
# perl  CombineFeaturesSummary.pl ./mipsRandpairSubsets/mipsRandpairSubset ./mipsRandpairSubsets/mipsRandpairSubset.summary.feature
# 

print "\n-------------------------------------------------------------------\n"; 


my $cmdPre = "perl ./train-set/CombineFeaturesSummary.pl  "; 

my $cmdPro = "   "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".sumry.fea "; 
print "$cmd\n"; 
system($cmd); 



