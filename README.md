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

4*. train your own model

==== Data ====

==== Output ====
