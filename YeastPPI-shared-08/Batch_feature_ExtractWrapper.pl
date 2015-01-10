######################################################################3
#
# copyright @ Yanjun Qi , qyj@cs.cmu.edu
# Please cite: 
# Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 
# Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, “Random Forest Similarity for Protein-Protein Interaction Prediction from Multiple source”, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 
# 
######################################################################3


# This program is a yeast PPI feature extraction wrapper 
# perl command inputPairlist


use strict; 
die "Usage: command inputPairFile \n" if scalar(@ARGV) < 1;
my ($inputPair ) = @ARGV;


print "\n--------------------------- 1gene-expression ----------------------------------------\n"; 

# -------------------   1gene-expression   ------------------------------

my $cmdPre = "perl ./1gene-expression/get_gene_expression.pl  "; 
my $cmdPro = "./1gene-expression/YeastGeneListOrfGeneName-106_pval_v9.0.txt ./1gene-expression/all_expression_fixed_s4_csv.txt  ./1gene-expression/expressionYanjunSplit.txt 0.6 "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".genexp" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   2tf-group-binding  ------------------------------

# perl -d get_tfGroupBinding.pl ./lists/sciencesubset.txt pvalbygene_nature04.txt 204_pvalbygene_nature04_TFs.groupIndex 0.05 ./lists/sciencesubset.tfgroup

my $cmdPre = "perl ./2tf-binding/get_tfGroupBinding.pl  "; 
my $cmdPro = " ./2tf-binding/pvalbygene_nature04.txt  ./2tf-binding/204_pvalbygene_nature04_TFs.groupIndex  0.05 "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".tfgroup" ; 
print "$cmd\n"; 
system($cmd); 





print "\n-------------------------------------------------------------------\n"; 

# -------------------   3gene-ontology  ------------------------------

# perl get_go_detail.pl temp.line.txt go_slim_mapping04.tab go_slim_function04.txt go_slim_process04.txt go_slim_component04.txt temp.gofunc.txt temp.goproc.txt temp.gocomp.txt

my $cmdPre = "perl ./3gene-ontology/get_go_detail.pl  "; 
my $cmdPro = " ./3gene-ontology/go_slim_mapping04.tab ./3gene-ontology/go_slim_function04.txt ./3gene-ontology/go_slim_process04.txt ./3gene-ontology/go_slim_component04.txt "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".gofunc ".$inputPair.".goproc ".$inputPair.".gocomp " ;
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   4protein-expression  ------------------------------
# 
# perl get_protein_expression.pl temp.line.txt OSheaLocalization.WeissmanAbundance.tab temp.line.exp

my $cmdPre = "perl ./4protein-expression/get_protein_expression.pl  "; 
my $cmdPro = " ./4protein-expression/OSheaLocalization.WeissmanAbundance.tab "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".proexp" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   5essentiality  ------------------------------

# perl get_Essential.pl Essential_ORFs_list.txt inputfile outfile 

my $cmdPre = "perl ./5essentiality/get_Essential.pl  "; 
my $cmdPro = " ./5essentiality/Essential_ORFs_list.txt "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".ess" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------  6HighExp-PPI   ------------------- 

# perl get_PPIMS.pl ./lists/sciencesubset.txt nature-analyze-hmspci-matrix.txt nature-analyze-tap-matrix.txt ./lists/sciencesubset.msmatrix

my $cmdPre = "perl ./6HighExp-PPI/MS/get_PPIMS.pl "; 
my $cmdPro = " ./6HighExp-PPI/MS/nature-analyze-hmspci-matrix.txt ./6HighExp-PPI/MS/nature-analyze-tap-matrix.txt "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".msmatrix" ; 
print "$cmd\n"; 
system($cmd); 


print "\n-------------------------------------------------------------------\n"; 

# perl get_Y2H.pl ./lists/sciencesubset.txt nature-analyze-y2h.txt sgdlite-y2h_hits.txt ./lists/sciencesubset.y2h

my $cmdPre = "perl ./6HighExp-PPI/Y2H/get_Y2H.pl "; 
my $cmdPro = " ./6HighExp-PPI/Y2H/nature-analyze-y2h.txt  ./6HighExp-PPI/Y2H/sgdlite-y2h_hits.txt  "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".y2h" ; 
print "$cmd\n"; 
system($cmd); 




print "\n-------------------------------------------------------------------\n"; 

# -------------------   7genetic-interaction  ------------------------------

# perl get_synthetic.pl ./lists/sciencesubset.txt tong2004.tab nature-compare-synthetic-lethality.txt ./lists/sciencesubset.gensyn

my $cmdPre = "perl ./7genetic-interaction/get_synthetic.pl  "; 
my $cmdPro = " ./7genetic-interaction/tong2004.tab  ./7genetic-interaction/nature-compare-synthetic-lethality.txt "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".gensyn" ; 
print "$cmd\n"; 
system($cmd); 




print "\n-------------------------------------------------------------------\n"; 

# -------------------   8nature-compare-sequence  ------------------------------

# perl get_nature_compare_sequence.pl ./lists/sciencesubset.txt nature-compare-gene-fusion.txt nature-compare-gene-neighborhood.txt nature-compare-gene-cooccurence.txt ./lists/sciencesubset.natCopSeq

my $cmdPre = "perl ./8nature-compare-sequence/get_nature_compare_sequence.pl  "; 
my $cmdPro = " ./8nature-compare-sequence/nature-compare-gene-fusion.txt ./8nature-compare-sequence/nature-compare-gene-neighborhood.txt ./8nature-compare-sequence/nature-compare-gene-cooccurence.txt "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".natCopSeq" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   9mips-pclass  ------------------------------

# perl get_MIPS_classPhe.pl ./lists/sciencesubset.txt classes.scheme protein_classes ./lists/sciencesubset.proClass 

my $cmdPre = "perl ./9mips-pclass/get_MIPS_classPhe.pl  "; 
my $cmdPro = " ./9mips-pclass/classes.scheme  ./9mips-pclass/protein_classes "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".proClass" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   10mips-phenotype  ------------------------------

# perl get_MIPS_classPhe.pl ./lists/sciencesubset.txt phencat.scheme phencat_data_07102004 ./lists/sciencesubset.phecat

my $cmdPre = "perl ./10mips-phenotype/get_MIPS_classPhe.pl  "; 
my $cmdPro = " ./10mips-phenotype/phencat.scheme ./10mips-phenotype/phencat_data_07102004 "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".phecat" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   11sequence-similarity  ------------------------------

# perl get_scBestHits.pl ./lists/sciencesubset.txt best_hits.sc.tab ./lists/sciencesubset.besthits

my $cmdPre = "perl ./11sequence-similarity/get_scBestHits.pl  "; 
my $cmdPro = " ./11sequence-similarity/best_hits.sc.tab  "; 

my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".besthits" ; 
print "$cmd\n"; 
system($cmd); 



print "\n-------------------------------------------------------------------\n"; 

# -------------------   12homology-PPI ------------------------------

# perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.sc.tab ./dip-ppi/Scere20041003.tab ./lists/sciencesubset.scHmPPI
# perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.dmelanogaster.tab ./dip-ppi/Dmela20041003.tab ./lists/sciencesubset.DmeHmPPI
# perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.humans.tab ./dip-ppi/Hsapi20041003.tab ./lists/sciencesubset.HsapiHmPPI
# perl   get_homologyPPI.pl ./lists/sciencesubset.txt ./psi_blast/psi_blast.celegan.tab ./dip-ppi/Celeg20041003.tab ./lists/sciencesubset.CeHmPPI


my $cmdPre = "perl   ./12homology-PPI/get_homologyPPI.pl  "; 

my $cmdPro = " ./12homology-PPI/psi_blast/psi_blast.sc.tab ./12homology-PPI/dip-ppi/Scere20041003.tab  "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".scHmPPI" ; 
print "$cmd\n"; 
system($cmd); 


print "\n-------------------------------------------------------------------\n"; 

my $cmdPro = " ./12homology-PPI/psi_blast/psi_blast.dmelanogaster.tab ./12homology-PPI/dip-ppi/Dmela20041003.tab  "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".DmeHmPPI" ; 
print "$cmd\n"; 
system($cmd); 


print "\n-------------------------------------------------------------------\n"; 

my $cmdPro = " ./12homology-PPI/psi_blast/psi_blast.humans.tab ./12homology-PPI/dip-ppi/Hsapi20041003.tab  "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".HsapiHmPPI" ; 
print "$cmd\n"; 
system($cmd); 


print "\n-------------------------------------------------------------------\n"; 

my $cmdPro = " ./12homology-PPI/psi_blast/psi_blast.celegan.tab ./12homology-PPI/dip-ppi/Celeg20041003.tab  "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".CeHmPPI" ; 
print "$cmd\n"; 
system($cmd); 


print "\n-------------------------------------------------------------------\n"; 

# -------------------   13domain-interaction  ------------------------------

# perl  get_DDI2PPI_usc.pl ./lists/sciencesubset.txt YeastGeneListOrfGeneName-106_pval_v9.0.txt ProteinAction_80SGDMIPS.txt ./lists/sciencesubset.MIPSddippi

my $cmdPre = "perl ./13domain-interaction/get_DDI2PPI_usc.pl  "; 

# perl  get_DDI2PPI_usc.pl ./lists/sciencesubset.txt YeastGeneListOrfGeneName-106_pval_v9.0.txt ProteinAction_025_80SGDY2H.txt ./lists/sciencesubset.Y2Hddippi
my $cmdPro = " ./13domain-interaction/YeastGeneListOrfGeneName-106_pval_v9.0.txt  ./13domain-interaction/ProteinAction_025_80SGDY2H.txt  "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".Y2Hddippi" ; 
print "$cmd\n"; 
system($cmd); 




print "\n-------------------------------------------------------------------\n"; 

# -------------------   combine features into one set  ------------------------------

# perl  CombineFeatures.pl ./mipsRandpairSubsets/mipsRandpairSubset ./mipsRandpairSubsets/mipsRandpairSubset.tfgroup.feature  G
my $cmdPre = "perl ./train-set/CombineFeatures.pl  "; 

my $cmdPro = "   "; 
my $cmd = $cmdPre." ".$inputPair." ".$cmdPro." ".$inputPair.".tfgroup.feature"."  G" ; 
print "\n$cmd\n"; 
system($cmd); 


