# CLADES
# CLADES (CLAssification based DElimitation of Species) is a python code designed for species delimitation problem.

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
  
For example, to analyze sequence data in temp_seq.txt:

  $ python CLADES.py temp ../model/All <path to libsvm commands>
  
4*. train your own model

==== Data ====

1. Input data

CLADES uses sequences data for multiple individuals from multiple populations. There is an example input sequence file temp_seq.txt.

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

Header contains the information of number of inidviduals, and length of such loci. Sequence data for each individual should be listed after the header line one by one, with format of '<individual id> <sequence data>'. <individual id> should contains a unique ID for such individual with letters and numbers, and delimited with '^', then combined with population ID. For example, 'a01^A1' stands for individual a01 belonging to population A1. Input data allows missing information with '-' to represent. One can use multiple locus sequence data, use an empty line to delimit locus.



==== Output ====
