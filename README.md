EvidenceFusionPPIRanking
========================

A data-driven Information Integration Framework for Protein Interaction Prediction



-------------------------------------

If using codes in subdir: YeastPPI-shared-08, please cite : 

<BR>

==>  Major: 
Y. Qi, Z. Bar-Joseph, J. Klein-Seetharaman, "Evaluation of different biological data and computational classification methods for use in protein interaction prediction", PROTEINS: Structure, Function, and Bioinformatics. 63(3):490-500. 2006

<BR>

@article{qi2006evaluation,
  title={Evaluation of different biological data and computational classification methods for use in protein interaction prediction},
  author={Qi, Yanjun and Bar-Joseph, Ziv and Klein-Seetharaman, Judith},
  journal={Proteins: Structure, Function, and Bioinformatics},
  volume={63},
  number={3},
  pages={490--500},
  year={2006},
  publisher={Wiley Online Library}
}


<BR>

==> Two related papers: 
Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, "A mixture of feature experts approach for protein-protein interaction prediction", BMC Bioinformatics 8 (S10):S6, 2007 

Y. Qi, J. Klein-Seetharaman, Z. Bar-Joseph, ìRandom Forest Similarity for Protein-Protein Interaction Prediction from Multiple sourceî, Pacific Symposium on Biocomputing 10: (PSB 2005) Jan. 2005. 


<BR>

	Abstract
Protein-protein interactions play a key role in many biological systems. High-throughput methods can directly detect the set of interacting proteins in yeast, but the results are often incomplete and exhibit high false-positive and false-negative rates. Recently, many different research groups independently suggested using supervised learning methods to integrate direct and indirect biological data sources for the protein interaction prediction task. However, the data sources, approaches, and implementations varied. Furthermore, the protein interaction prediction task itself can be subdivided into prediction of (1) physical interaction, (2) co-complex relationship, and (3) pathway co-membership. To investigate systematically the utility of different data sources and the way the data is encoded as features for predicting each of these types of protein interactions, we assembled a large set of biological features and varied their encoding for use in each of the three prediction tasks. Six different classifiers were used to assess the accuracy in predicting interactions, Random Forest (RF), RF similarity-based k-Nearest-Neighbor, Naïve Bayes, Decision Tree, Logistic Regression, and Support Vector Machine. For all classifiers, the three prediction tasks had different success rates, and co-complex prediction appears to be an easier task than the other two. Independently of prediction task, however, the RF classifier consistently ranked as one of the top two classifiers for all combinations of feature sets. Therefore, we used this classifier to study the importance of different biological datasets. First, we used the splitting function of the RF tree structure, the Gini index, to estimate feature importance. Second, we determined classification accuracy when only the top-ranking features were used as an input in the classifier. We find that the importance of different features depends on the specific prediction task and the way they are encoded. Strikingly, gene expression is consistently the most important feature for all three prediction tasks, while the protein interactions identified using the yeast-2-hybrid system were not among the top-ranking features under any condition.



<BR>

	supplementary Web: http://www.cs.cmu.edu/%7Eqyj/papers_sulp/proteins05_PPI.html
