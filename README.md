# CLADES (CLAssification based DElimitation of Species)
CLADES is a python code designed for species delimitation problem.

### June 6, 2024: CLADES was written by Jingwen Pei, former PhD student at University of Connecticut. Jingwen has graduated and left academia. If you have questions about CLADES, please contact Yufeng Wu (yufeng.wu@uconn.edu) and/or post your questions under "Issues".

## Background
### Model
In our paper, we conducted two groups of extensive simulation to build SVM model, one (All) for data without considering migration and one (AllM) for data with migration. The features of simulated data and trained models are listed under directory `model/`. The files with the same prefix contain the information for the same model. For example, files with prefix 'All' are for simulation with no migration:
* All.sumstat :features collected from simulation data
* All.sumstat.scale :normalized features to range (0, 1)
* All.range :normalization parameter for each feature
* All.sumstat.scale.model :trained SVM model to delimit two populations
 
 The same for files with prefix 'AllM'.
 
Trained data is simulated with three various parameters:
* θ :population size parameter, range:(0.0005, 0.02)
* 𝜏 :species divergence time, range:(θ/10,θ)
* M :migration parameter, M=Nm(N is effective population size, m is migration rate per generation), range:(0,5)

Two models can be used for general analysis for species detetion. If there is no significant gene flow between two populations, we suggest to use model 'All'. Both models can be used for species delimitation. CLADES is not sensitive to three parameters, and it is able to identify two species, especially for hard cases (recent diverged or with a medium level of gene flow).

### Summary Statistics
We used 5 summary statistics to represent the data, they are:
* Folded SFS with 3-bin :folded site frequency spectrum
* Private position
* Pairwise difference ratio
* Fst
* Shared tract percentage
  
For more details and how to compute these summary statistics, please refer to our paper.

## Prerequisites
* [Python 3](https://www.python.org/downloads/)
* [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
  * Ensure that the following commands are accessible: `svm-predict`, `svm-scale`, `svm-train`

## Usage
```
usage: CLADES.py [-h] [-p MODEL_PATH] [-n MODEL_NAME] [-t] [-T] seq_data

CLADES (CLAssification based DElimitation of Species)

positional arguments:
  seq_data              The path to the sequence data file to be processed.

options:
  -h, --help            show this help message and exit
  -p, --model_path MODEL_PATH
                        Path to model. Defaults to `model/` being in the same directory as `CLADES.py`
  -n, --model_name MODEL_NAME
                        Name of the model. Defaults to `All`
  -t, --total_sumstat   Outputs a single file containing all of the summary statistics on top of pairwise sumstat.
  -T, --only_total_sumstat
                        Only outputs a single file containing all of the summary statistics. Does not output pairwise sumstat.
```

CLADES is able to delimit multiple populations (n >= 2). To do species delimitation:
```
python CLADES.py seq_data
```

For example, to analyze the sequence data in test_seq.txt:
```
python CLADES.py test_seq.txt
```
You may also specify the path to your model. In this case, it's assuming that the model name is 'All':
```
python CLADES.py -p path/to/model/ test_seq.txt
```


### Data
CLADES uses sequences data for multiple individuals from multiple populations.

Sequence data should have the following format:

* Header contains the information of number of inidviduals, and length of such loci
* Following that is the sequence data for each individual should be listed line by line, with format of `<individual id> <sequence data>`
  * `<individual id>` should contains a unique ID for such individual with letters and numbers, and delimited with '^', then combined with population ID. For example, `a01^A1` stands for individual `a01` belonging to population `A1`
  * `<sequence data>` is allowed to have missing information and can be represented with '-'
* One can use multiple locus sequence data, use an empty line to delimit two locus.

```
num_individuals loci_length

individual1^population1 sequence_data
individual2^population1 sequence_data
...
```

Here is an example input sequence from `test_seq.txt`, which contains two locus information for 80 individuals from 8 populations:

```
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
```

### Output
CLADES outputs the best assignment for populations in sequence data. For example, to use either model (All or AllM) to delimit test_seq.txt, CLADES outputs:

```
...
The Best Assignment of species clusters are:
['C1,C2', 'B1,B2', 'D1,D2', 'A1,A2']
```

CLADES clusters A1 and A2 into the same cluster, same with (B1, B2), (C1, C2) and (D1, D2). That means, species delimitation for these 8 populations is 4 species, and each species has two populations.

CLADES computes summary statistics for pairwise populations to file `<pop1>-<pop2>.sumstat`. Then uses the model to determine whether two populations belong to the same species. The estimation information is contained in file `<pop1>-<pop2>.out`. File `<pop1>-<pop2>.out` looks like:

```
labels 1 -1
-1 0.000253429 0.999747
-1 0.00024749 0.999753
```

Label '1' means cluster 'different species', '-1' stands for 'same species' cluster. First line '-1 0.000253429 0.999747' indicates that, the first sample is classified to 'same species' with probability 0.9997 and can be classified to 'different species' with a very low probability 0.000253. If you have true labels for your dataset, you can compare the classification result with true labels to get estimation accuracy.

## Train Your Own Model
With prior knowledge of related populations, one can trained a more specific model as your own classifier. For example, if one has an estimate range of θ, 𝜏 and M, one can generate simulated sequence data for two-species model. For pairwise populations, compute summary statistics and normalize to range (0, 1). Then use `svm-train` to train your own model.
