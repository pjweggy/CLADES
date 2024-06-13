__author__ = 'pjweggy'

import re
import numpy as np
from itertools import chain
import os
import subprocess
import sys
from pathlib import Path
import argparse

__all__ = ["compute_summary_statistics", "species_delimitation", "compute_best_assignment", "get_species_names", "get_class_labels", "write_sumstat", "write_total_sumstat"]

def compute_summary_statistics(seq_data, nbin = 3, delim = '\^'):
    """
    Compute summary statistics for sequence data.

    Args:
        seq_data: Sequence data input file.
        delim: A delimiter used for the header.
        nbin: Number of bins to collapse the folded-SFS.

    Returns:
        Dictionary containing summary statistics.
    """
    SS = dict()
    rawdata = dict()

    with seq_data.open() as f:
        for line in f:
            if line[0:2] != "\n":
                header = line.split()[0]
            else:
                header = ""

            indvName = RecogIndv(header, delim)
            if len(indvName) > 0:
                seq = GetSeq(line)
                rawdata[indvName] = seq
            else:
                SS = Process(rawdata, delim, nbin, SS)
                rawdata.clear()

        SS = Process(rawdata, delim, nbin, SS)
        rawdata.clear()

    return SS

def species_delimitation(SS, output_dir, model):
    """
    Perform species delimitation for pairwise populations.

    Args:
        SS: Summary statistics.
        output_dir: Directory where output files will be saved.
        model: Path to the model used for scaling and prediction.

    Returns:
        Dictionary with results of species delimitation.
    """
    Res = dict()
    for key in sorted(SS):
        sumstat_filepath = output_dir / f"{key}.sumstat"

        subprocess.run(
            f"svm-scale -r {model}.range {sumstat_filepath} > {sumstat_filepath}.scale",
            shell = True,
            stdout = subprocess.PIPE
        )

        key_out_filepath = output_dir / f"{key}.out"
        subprocess.run(
            f"svm-predict -b 1 -q {sumstat_filepath}.scale {model}.sumstat.scale.model {key_out_filepath}",
            shell = True,
            stdout = subprocess.PIPE
        )

        Out = np.loadtxt(key_out_filepath, comments = "labels")
        if len(Out.shape) == 1:
            Out = Out.reshape(1, Out.shape[0])

        res = np.mean(Out, axis = 0)[1:3]
        Res[key] = res

    seq_out_filepath = output_dir / f"{seq_data.name}.out"
    with seq_out_filepath.open("w") as f:
        f.write("labels +1 -1\n")
        for key in sorted(Res):
            f.write(f"{key} {Res[key][0]:.4f} {Res[key][1]:.4f}\n")

    return Res

def compute_best_assignment(Res):
    """
    Compute the probability of the best assignment of species clusters.

    Args:
        Res: Dictionary with results of species delimitation.

    Returns:
        List of species clusters, the best assignment of species clusters, and its probability.
    """
    clusters = []
    cluster = GetSPN(Res)
    prob = CompProb(cluster, Res)
    best_cluster = None
    coal = 1

    while coal:
        coalc = 0
        for i in range(len(cluster) - 1):
            for j in range(i + 1, len(cluster)):
                cur_cluster = CoalSPN(cluster, cluster[i], cluster[j])
                cur_prob = CompProb(cur_cluster, Res)
                clusters.append(f"{cur_cluster} {cur_prob}")

                if cur_prob > prob:
                    best_cluster = cur_cluster
                    prob = cur_prob
                    coalc += 1

        if coalc > 0:
            cluster = best_cluster
        else:
            coal = 0

    return clusters, cluster, prob

def get_species_names(SS):
    """
    Get the species names from summary statistics.

    Args:
        SS: Summary statistics.

    Returns:
        A list containing the species names.
    """
    species_names = set()

    for key in SS:
        specie1, specie2 = key.split("-")
        species_names.add(specie1)
        species_names.add(specie2)

    return sorted(list(species_names))

def get_class_labels(species_names):
    """
    Get the class labels from species names.

    Args:
        species_names: A list containing the species names.

    Returns:
        A list containing the class labels.
    """
    class_labels = list()

    for i in range(len(species_names) - 1):
        for j in range(i + 1, len(species_names)):
            if species_names[i][0] == species_names[j][0]:
                class_labels.append("-1")
            else:
                class_labels.append("+1")

    return class_labels

def write_sumstat(SS, output_dir):
    """
    Writes summary statistics to a file named {key}.sumstat in the specified output directory.
    
    Args:
        SS: Summary statistics.
        output_dir: Directory where the summary statistics file will be saved.
    """
    for key in sorted(SS):
        sumstat_filepath = output_dir / f"{key}.sumstat"
        with sumstat_filepath.open("w") as f:
            f.write(SS[key])

def write_total_sumstat(SS, output_dir):
    """
    Writes all summary statistics to a single file named total.sumstat in the specified output directory.
    
    Args:
        SS: Summary statistics.
        output_dir: Directory where the total summary statistics file will be saved.
    """
    species_names = get_species_names(SS)
    class_labels = get_class_labels(species_names)

    pairwise_sumstats = list()
    for key in sorted(SS):
        sumstats = [sumstat for sumstat in SS[key].split('\n') if sumstat]
        pairwise_sumstats.append(sumstats)

    total_sumstats = list(map(list, zip(*pairwise_sumstats)))

    total_sumstat_filepath = output_dir / "total.sumstat"
    with total_sumstat_filepath.open("w") as f:
        for group in total_sumstats:
            for i, sumstat in enumerate(group):
                f.write(f"{class_labels[i]} {sumstat}\n")

def RecogIndv(line, delim):
    pattern = fr"[A-Za-z]+[0-9]+{delim}[A-Za-z0-9]+"
    return line if re.match(pattern, line) else ""

def Process(rawdata, delim, nbin, SS):
    if len(rawdata) != 0:
        popsize, Allmat, spn = Dict2Mat(rawdata, delim)
        if isinstance(Allmat, list):
            print("Please have at least 2 individuals per locus.")
            sys.exit(1)
        nsites = Allmat.shape[1]
        Allmat, snplist = Seq2SNP(Allmat, rawdata)
        startpoint = GetSP(popsize)

        for i in range(len(popsize) - 1):
            for j in range(i + 1, len(popsize)):
                pos = list(chain(range(startpoint[i], startpoint[i] + popsize[i]), range(startpoint[j], startpoint[j] + popsize[j])))
                hap = Allmat[pos,:]
                sumstat = CollectSS2pop(hap, popsize[i], popsize[j], nbin, nsites, snplist)
                key = f"{spn[i]}-{spn[j]}"

                if key in SS:
                    SS[key] += AddSS(sumstat)
                else:
                    SS[key] = AddSS(sumstat)

    return SS

def AddSS(sumstat):
    s = ""

    for i in range(len(sumstat)):
        s += f"{i + 1}:{sumstat[i]:.6f} "
    s += "\n"

    return s

def CollectSS2pop(hap, popsize1, popsize2, nbin, nsites, snplist):
    hap1 = hap[0:popsize1,:]
    hap2 = hap[popsize1:,:]
    sfs = SFS_fold_bin(hap, nbin, nsites)
    #sfs=SFS_fold_bin(hap,nbin)
    #pdmean,pdstd=PairDiff(hap)
    pdmean1 = PairDiff(hap1)
    pdmean2 = PairDiff(hap2)
    pdmean3 = PairDiff2(hap1, hap2)
    pdmean4 = (pdmean1 + pdmean2) / 2
    n_prvpos = prvpos(hap, popsize1, popsize2, nsites)
    lsubseq = Lsubseq(hap, popsize1, popsize2, nsites, snplist)
    #n_prvpos=prvpos(hap,popsize1,popsize2)
    #lsubseq=Lsubseq(hap,popsize1,popsize2)

    if pdmean4 == 0:
        ktheta = 2
    else:
        ktheta =pdmean3 / pdmean4

    if pdmean3 == 0:
        Fst = 0
    else:
        Fst = (pdmean3 - pdmean4) / pdmean3
    ss = np.hstack((sfs, ktheta))
    ss = np.hstack((ss, Fst))
    #ss=np.hstack((ss,pdmean4))
    ss = np.hstack((ss, n_prvpos))
    ss = np.hstack((ss, lsubseq))

    return ss

def SFS_fold(hap, nsites):
    if isinstance(hap, list):
        hap_dim = list()
        hap_dim.append(1)
        hap_dim.append(len(hap))
    else:
        hap_dim = hap.shape

    if hap_dim[0] % 2 != 0:
        sfs = np.zeros((hap_dim[0] + 1) // 2)
    else:
        sfs = np.zeros(hap_dim[0] // 2 + 1)

    for i in range(hap_dim[1]):
        temp = (hap[:,i] == 1).sum()
        if temp <= (hap_dim[0] / 2):
            sfs[int(temp)] += 1
        else:
            sfs[hap_dim[0] - int(temp)] += 1

    return sfs / hap_dim[1]
    #return sfs/nsites

###Folded SFS with bin
def SFS_fold_bin(hap, nbin, nsites):
    sfs = SFS_fold(hap, nsites)
    nlen = sfs.shape[0]
    window = nlen // nbin
    sfs_bin = np.zeros(nbin)

    for i in range(nbin):
        if i != (nbin - 1):
            sfs_bin[i] = sum(sfs[window * i:window * (i + 1)])
        else:
            sfs_bin[i] = sum(sfs[window * i:])

    return sfs_bin

###Intro-population Pairwise difference
def PairDiff(hap):
    hap_dim = hap.shape

    pairdiff = list()
    for i in range(hap_dim[0] - 1):
        for j in range(i + 1, hap_dim[0]):
            pairdiff.append(pd(hap[i,:],hap[j,:]))
    #print pairdiff
    #pairdiff_reg=[x*(float(1.0)/hap_dim[1]) for x in pairdiff]  ##regularize pairwise difference
    if len(pairdiff) > 0:
        pdmean = np.mean(pairdiff)
    else:
        pdmean = 0.0

    #pdstd=std(pairdiff_reg)
    return pdmean #,pdstd

def pd(hap1, hap2):
    pd = 0
    for i in range(len(hap1)):
        if abs(hap1[i] - hap2[i]) > 0:
            pd += 1

    return pd

###Outra-population Pairwise difference
def PairDiff2(hap1, hap2):
    hap1_dim = hap1.shape[0]
    hap2_dim = hap2.shape[0]
    pairdiff = list()

    for i in range(hap1_dim):
        for j in range(hap2_dim):
            pairdiff.append(pd(hap1[i,:],hap2[j,:]))
    #print pairdiff
    #pairdiff_reg=[x*(float(1.0)/hap1.shape[1]) for x in pairdiff]
    pdmean = np.mean(pairdiff)
    #pdstd=std(pairdiff)
    return pdmean#,pdstd

#percentage of private positions overall SNPs
def prvpos(hap, popsize1, popsize2, nsites):
    hap_dim = hap.shape
    n_prvpos = 0
    for i in range(hap_dim[1]):
        n0hap1 = (hap[0:popsize1,i] == 0).sum()
        n1hap1 = (hap[0:popsize1,i] == 1).sum()
        n0hap2 = (hap[popsize1:,i] == 0).sum()
        n1hap2 = (hap[popsize1:,i] == 1).sum()

        if (n0hap1 == popsize1 or n1hap1 == popsize1) and (n0hap2 != popsize2 and n1hap2 != popsize2):
            #print i+1
            n_prvpos += 1
        elif (n0hap2 == popsize2 or n1hap2 == popsize2) and (n0hap1 != popsize1 and n1hap1 != popsize1):
            #print i+1
            n_prvpos +=1

    return float(n_prvpos) / hap_dim[1]
    #return float(n_prvpos)/nsites

def Lsubseq(hap, popsize1, popsize2, nsites, snplist):
    tract = list()
    lsubseq = list()
    nsnp = hap.shape[1]

    for i in range(popsize1):
        for j in range(popsize1, popsize1 + popsize2):
            nshare = 0
            ##compare two sequences
            for k in range(nsnp):
                if hap[i,k] == hap[j,k]:
                    if k == 0:
                        nshare += snplist[k] + 1
                    else:
                        nshare += snplist[k] - snplist[k - 1]
                    if k == nsnp - 1:
                        tract.append(nshare)
                        nshare = 0
                else:
                    if k==0:
                        nshare += snplist[0]
                    else:
                        nshare += snplist[k] - snplist[k - 1] - 1
                    tract.append(nshare)
                    if k == nsnp - 1:
                        tract.append(nsites - snplist[-1])
                    nshare = 0
            #print(tract)
            lsubseq.append(max(tract))
            tract = list()

    Lsubseq = max(lsubseq)
    #return float(Lsubseq)/hap.shape[1]
    return float(Lsubseq) / nsites


def GetSP(popsize):
    startpoint = [0]
    s = 0

    for i in range(len(popsize) - 1):
        s += popsize[i]
        startpoint.append(s)

    return startpoint

def Dict2Mat(rawdata, delim):
    speciesname = list()
    speciesinfo = dict()
    template = ""
    Allmat = []

    for key in sorted(rawdata):
        ###process key
        #r=re.compile('([a-d][1-2])([0-9]+)('+delim+')([A-D][1-2])')
        r = re.compile('([A-Za-z]+[0-9]+)('+delim+')([A-Za-z0-9]+)')
        m = r.match(key)
        spn = m.group(3)

        if spn in speciesinfo:
            speciesinfo[spn] += 1
        else:
            speciesname.append(spn)
            speciesinfo[spn] = 1
        ### process sequence
        if len(template) == 0:
            template = rawdata[key]
            Allmat = Gen01fromSeq(rawdata[key], template)
        else:
            #Update(template,rawdata[key])
            Allmat = np.vstack((Allmat, Gen01fromSeq(rawdata[key], template)))

    popsize = [speciesinfo[spn] for spn in speciesname]

    return popsize, Allmat, speciesname

def Seq2SNP(Allmat, rawdata):
    nr, nc = Allmat.shape
    snplist = [i for i in range(nc) if any(Allmat[:, i])]

    if len(snplist) == 0:
        print(rawdata)
        print("Error Message:Gene has no variant! Please delete this gene or concatenate multiple genes.")
        sys.exit(1)

    AllSNP = np.array(Allmat[:,snplist])
#    print AllSNP.shape
    return AllSNP, snplist

def Update(template, newtem):
    for i in range(len(template)):
        if (template[i] == '-' and newtem[i] != '-'):
            if i != (len(template) - 1):
                tmp = template[0:i] + newtem[i] + template[(i + 1):]
            else:
                tmp = template[0:i] + newtem[i]
            template = tmp

    return template

def Gen01fromSeq(seq, template):
    hap = list()

    for i in range(len(template)):
        if seq[i] == '-':
            hap.append(2)
        elif seq[i] == template[i]:
            hap.append(0)
        else:
            hap.append(1)

    return hap

def GetSeq(line):
    header = line.split()[0]
    line = line.replace(header, "")
    line = line.replace(" ", "")
    line = line.replace("\n", "")
    return line

def GetSPN(Res):
    spn = []

    for key in sorted(Res):
        spn1 = key.split('-')[0]
        spn2 = key.split('-')[1]
        if not (spn1 in spn):
            spn.append(spn1)
        if not (spn2 in spn):
            spn.append(spn2)

    return spn

def CoalSPN(CurC, spn1, spn2):
    NewC = [spn1 + ',' + spn2] + [spn for spn in CurC if spn != spn1 and spn != spn2]
    return NewC

def ProbSS(splist1, Res):
    prob = 1

    if len(splist1) != 1:
        for i in range(len(splist1) - 1):
            for j in range(i + 1, len(splist1)):
                pair_key1 = f"{splist1[i]}-{splist1[j]}"
                pair_key2 = f"{splist1[j]}-{splist1[i]}"

                if pair_key1 in Res:
                    prob *= Res[pair_key1][1]
                elif pair_key2 in Res:
                    prob *= Res[pair_key2][1]

    return prob

def ProbDS(splist1, splist2, Res):
    prob = 1

    for spn1 in splist1:
        for spn2 in splist2:
            pair_key1 = f"{spn1}-{spn2}"
            pair_key2 = f"{spn2}-{spn1}"

            if pair_key1 in Res:
                prob *= Res[pair_key1][0]
            elif pair_key2 in Res:
                prob *= Res[pair_key2][0]

    return prob

def CompProb(CurC, Res):
    prob = 1

    for i in range(len(CurC)-1):
        splist1 = CurC[i].split(",")
        prob *= ProbSS(splist1, Res)

        for j in range(i + 1, len(CurC)):
            splist2 = CurC[j].split(",")
            prob *= ProbDS(splist1, splist2, Res)

    if len(CurC) == 1:
        splist1 = CurC[0].split(",")
        prob += ProbSS(splist1, Res)

    return prob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "CLADES (CLAssification based DElimitation of Species)")
    
    # Positional arguments
    parser.add_argument("seq_data", type = Path, help = "The path to the sequence data file to be processed.")
    
    # Optional argumentS
    parser.add_argument("-p", "--model_path", type = Path, default = Path("model/"), help = "Path to model.")
    parser.add_argument("-n", "--model_name", type = str, default = "All", help = "Name of the model.")
    parser.add_argument("-t", "--total_sumstat", action = "store_true", help = "Outputs a single file containing all of the summary statistics on top of pairwise sumstat.")
    parser.add_argument("-T", "--only_total_sumstat", action = "store_true", help = "Only outputs a single file containing all of the summary statistics. Does not output pairwise sumstat.")

    args = parser.parse_args()

    if not args.seq_data.is_file():
        print(f"The sequence data '{args.seq_data}' does not exist.")
        sys.exit(1)

    if not args.model_path.is_dir():
        print(f"The model path '{args.model_path}' does not exist.")
        sys.exit(1)

    seq_data = args.seq_data
    model = args.model_path / args.model_name
    output_dir = seq_data.parent / "output/"
    output_dir.mkdir(parents = True, exist_ok = True)

    # Step 1
    SS = compute_summary_statistics(seq_data)

    # Step 2
    if args.total_sumstat or args.only_total_sumstat:
        write_total_sumstat(SS, output_dir)
        if args.only_total_sumstat:
            sys.exit(0)

    write_sumstat(SS, output_dir)

    # Step 3
    Res = species_delimitation(SS, output_dir, model)

    # Step 4
    clusters, cluster, prob = compute_best_assignment(Res)
    for c in clusters:
        print(c)

    print("\nThe Best Assignment of species clusters are:")
    print(f"{cluster} {prob}")
