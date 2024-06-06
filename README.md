# CLADES (CLAssification based DElimitation of Species)
## CLADES is a python code designed for species delimitation problem.
### June 6, 2024: CLADES was written by Jingwen Pei, formal PhD student at University of Connecticut. Jingwen has graduated and left acadmia. If you have questions about CLADES, please contact Yufeng Wu (yufeng.wu@uconn.edu) and/or post your questions under "Issues".

==== Prerequisite ====

1. python 2

All python code is written in version python/2.7.6. To run python code:

  $ python <PYTHONCODE>.py argv1 argv2 ...
  
2. libsvm

Two trained CLADES model (All/AllM) in directory model/ are presented in libsvm format. To delimit species with trained model or train your own CLADES model, one should install libsvm. LIBSVM is a library for support vector machine designed by Chih-Chung Chang and Chih-Jen Lin. For more details and how to install libsvm, please see libsvm website: https://www.csie.ntu.edu.tw/~cjlin/libsvm/.

After installation of LIBSVM, make sure that the following commands are accessible to use; otherwise, add the path of following commands to $PATH.

  svm-predict  svm-scale    svm-train

==== Usage ====

1. model

In our paper, we conducted two groups of extensive simulation to build SVM model, one (All) for data without considering migration and one (AllM) for data with migration. The features of simulated data and trained models are listed under directory model/. The files with same prefix contain the information for the same model. For example, files with prefix 'All' are for simulation with no migration:

  -- All.sumstat :features collected from simulation data
  
  -- All.sumstat.scale :normalized features to range (0, 1)
  
  -- All.range :normalization parameter for each feature
  
  -- All.sumstat.scale.model :trained SVM model to delimit two populations
 
 The same for files with prefix 'AllM'.
 
Trained data is simulated with three various parameters:

  -- θ :population size parameter, range:(0.0005, 0.02)
  
  -- 𝜏 :species divergence time, range:(θ/10,θ)
  
  -- M :migration parameter, M=Nm(N is effective population size, m is migration rate per generation), range:(0,5)

Two models can be used for general analysis for species detetion. If there is no significant gene flow between two populations, we suggest to use model 'All'. Both models can be used for species delimitation. CLADES is not sensitive to three parameters, and it is able to identify two species, especially for hard cases (recent diverged or with a medium level of gene flow).

2. summary statistics

We used 5 summary statistics to represent the data, they are:

  -- folded SFS with 3-bin :folded site frequency spectrum
  
  -- private position
  
  -- pairwise difference ratio
  
  -- Fst
  
  -- shared tract percentage
  
For more details and how to compute these summary statistics, please refer to our paper.

3. species delimitation

CLADES is able to delimit multiple populations(n>=2). To do species delimitation:

  $ python CLADES.py argv1 argv2 argv3
  
  -- argv1 :prefix of input sequence data
  
  -- argv2 :path and prefix of model
  
  -- argv3 :path of libsvm commands. if path is contained in $PATH, one can leave this argv empty.
  
For example, to analyze sequence data in test_seq.txt:

  $ python CLADES.py test ../model/All <path to libsvm commands>
  
*4. train your own model

With prior knowledge of related populations, one can trained a more specific model as your own classifier. For example, if one has an estimate range of θ, 𝜏 and M, one can generate simulated sequence data for two-species model. For pairwise populations, compute summary statistics and normalize to range (0, 1). Then use 'svm-train' to train your own model.

==== Data ====

CLADES uses sequences data for multiple individuals from multiple populations. There is an example input sequence file test_seq.txt containing two locus information for 80 individuals from 8 populations.

Sequence data should have format:

10 20

a01^A1 GTTTGCAAGAA-ATCATGGA
a02^A1 GTTTGCAAGAA-ATCATGGA
a03^A1 GTTTGCAAGAA-ATCATGGA
a04^A1 GTTTGCAAGAA-ATCATGGA
a05^A1 GTTTGCAAGAA-ATCATGGA
b01^B1 GATTGCTAGAATATC-TGGA
b02^B1 GATTGCTAGAATATC-TGGA
b03^B1 GATTGCTAGAATATC-TGGA
b04^B1 GATTGCTAGAATATC-TGGA
b05^B1 GATTGCTAGAATATC-TGGA

Header contains the information of number of inidviduals, and length of such loci. Sequence data for each individual should be listed after the header line one by one, with format of '<individual id> <sequence data>'. <individual id> should contains a unique ID for such individual with letters and numbers, and delimited with '^', then combined with population ID. For example, 'a01^A1' stands for individual a01 belonging to population A1. Input data allows missing information with '-' to represent. One can use multiple locus sequence data, use an empty line to delimit two locus.

==== Output ====

CLADES outputs the best assignment for populations in sequence data. For example, to use either model (All or AllM) to delimit test_seq.txt, CLADES outputs:

The Best Assignment of species clusters are:
['C1,C2', 'B1,B2', 'D1,D2', 'A1,A2']

CLADES clusters A1 and A2 into the same cluster, same with (B1, B2), (C1, C2) and (D1, D2). That means, species delimitation for these 8 populations is 4 species, and each species has two populations.

CLADES computes summary statistics for pairwise populations to file <pop1>-<pop2>.sumstat. Then use model to determine whether two populations belong to the same species. The estimation information is contained in file <pop1>-<pop2>.out. File <pop1>-<pop2>.out looks like:

labels 1 -1
-1 0.000253429 0.999747
-1 0.00024749 0.999753

Label '1' means cluster 'different species', '-1' stands for 'same species' cluster. First line '-1 0.000253429 0.999747' indicates that, the first sample is classified to 'same species' with probability 0.9997 and can be classified to 'different species' with a very low probability 0.000253. If you have true labels for your dataset, you can compare the classification result with true labels to get estimation accuracy.


