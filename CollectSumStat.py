__author__ = 'pjweggy'

import re
import numpy as np
#from SumStat import *
#import math
import sys

def RecogIndv(line,delim):
    if re.match('[a-d][1-2][0-9A-Za-z]+'+delim+'[A-D][1-2]',line):
        return line
    else:
        return ''


def Process(fout,rawdata,delim,nbin,nsites):
    if len(rawdata) != 0:
        [popsize,Allmat,classlabel]=Dict2Mat(rawdata,delim)
        Allmat,snplist=Seq2SNP(Allmat)
        sumstat1=[]
        startpoint=GetSP(popsize)
        for i in range(len(popsize)-1):
            for j in range(i+1,len(popsize),1):
                pos=range(startpoint[i],(startpoint[i]+popsize[i]),1)+range(startpoint[j],(startpoint[j]+popsize[j]),1)
                hap=Allmat[pos,:]
                if np.size(sumstat1)==0:
                    sumstat1=CollectSS2pop(hap,popsize[i],popsize[j],nbin,nsites,snplist)
                else:
                    sumstat1=np.vstack((sumstat1,CollectSS2pop(hap,popsize[i],popsize[j],nbin,nsites,snplist)))
        WriteSS(fout,sumstat1,classlabel)


def WriteSS(fout,sumstat,classlabel):
    nsam,nft=sumstat.shape
    for i in range(nsam):
        s=classlabel[i]+' '
        for j in range(nft):
            s=s+str(j+1)+':'+str("%.6f" % sumstat[i,j])+' '
        s=s+'\n'
        fout.write(s)



def CollectSS2pop(hap,popsize1,popsize2,nbin,nsites,snplist):
    hap1=hap[0:popsize1,:]
    hap2=hap[popsize1:,:]
    sfs=SFS_fold_bin(hap,nbin,nsites)
    #sfs=SFS_fold_bin(hap,nbin)
    #pdmean,pdstd=PairDiff(hap)
    pdmean1=PairDiff(hap1)
    pdmean2=PairDiff(hap2)
    pdmean3=PairDiff2(hap1,hap2)
    pdmean4=(pdmean1+pdmean2)/2
    n_prvpos=prvpos(hap,popsize1,popsize2,nsites)
    lsubseq=Lsubseq(hap,popsize1,popsize2,nsites,snplist)
    #n_prvpos=prvpos(hap,popsize1,popsize2)
    #lsubseq=Lsubseq(hap,popsize1,popsize2)

    if pdmean4==0:
        ktheta=2
    else:
        ktheta=pdmean3/pdmean4
    if pdmean3==0:
        Fst=0
    else:
        Fst=(pdmean3-pdmean4)/pdmean3
    ss=np.hstack((sfs,ktheta))
    ss=np.hstack((ss,Fst))
    #ss=np.hstack((ss,pdmean4))
    ss=np.hstack((ss,n_prvpos))
    ss=np.hstack((ss,lsubseq))
    return ss


def SFS_fold(hap,nsites):
    if isinstance(hap,list):
        hap_dim=list()
        hap_dim.append(1)
        hap_dim.append(len(hap))
    else:
        hap_dim=hap.shape
    if hap_dim[0]%2!=0:
        sfs=np.zeros((hap_dim[0]+1)/2)
    else:
        sfs=np.zeros(hap_dim[0]/2+1)
    for i in range(hap_dim[1]):
        temp=(hap[:,i] == 1).sum()
        if temp <= (hap_dim[0]/2):
            sfs[int(temp)] += 1
        else:
            sfs[hap_dim[0]-int(temp)] += 1
    return sfs/hap_dim[1]
    #return sfs/nsites

###Folded SFS with bin
def SFS_fold_bin(hap,nbin,nsites):
    sfs=SFS_fold(hap,nsites)
    nlen=sfs.shape[0]
    window=nlen/nbin
    sfs_bin=np.zeros(nbin)
    for i in range(nbin):
        if i != (nbin-1):
            sfs_bin[i]=sum(sfs[window*i:window*(i+1)])
        else:
            sfs_bin[i]=sum(sfs[window*i:])
    return sfs_bin

###Intro-population Pairwise difference
def PairDiff(hap):
    hap_dim=hap.shape

    pairdiff=list()
    for i in range(hap_dim[0]-1):
        for j in range(i+1,hap_dim[0],1):
            pairdiff.append(pd(hap[i,:],hap[j,:]))
    #print pairdiff
    #pairdiff_reg=[x*(float(1.0)/hap_dim[1]) for x in pairdiff]  ##regularize pairwise difference
    pdmean=np.mean(pairdiff)
    if len(pairdiff)==0:
        pdmean=0.0
    #pdstd=std(pairdiff_reg)
    return pdmean #,pdstd

def pd(hap1,hap2):
    pd=0
    for i in range(len(hap1)):
        if abs(hap1[i]-hap2[i])>0:
            pd+=1
    return pd

###Outra-population Pairwise difference
def PairDiff2(hap1,hap2):
    hap1_dim=hap1.shape[0]
    hap2_dim=hap2.shape[0]
    pairdiff=list()
    for i in range(hap1_dim):
        for j in range(hap2_dim):
            pairdiff.append(pd(hap1[i,:],hap2[j,:]))
    #print pairdiff
    #pairdiff_reg=[x*(float(1.0)/hap1.shape[1]) for x in pairdiff]
    pdmean=np.mean(pairdiff)
    #pdstd=std(pairdiff)
    return pdmean#,pdstd


#percentage of private positions overall SNPs
def prvpos(hap,popsize1,popsize2,nsites):
    hap_dim=hap.shape
    n_prvpos=0
    for i in range(hap_dim[1]):
        n0hap1=(hap[0:popsize1,i]==0).sum()
        n1hap1=(hap[0:popsize1,i]==1).sum()
        n0hap2=(hap[popsize1:,i]==0).sum()
        n1hap2=(hap[popsize1:,i]==1).sum()
        if (n0hap1==popsize1 or n1hap1==popsize1) and (n0hap2!=popsize2 and n1hap2 != popsize2):
            #print i+1
            n_prvpos += 1
        elif (n0hap2==popsize2 or n1hap2==popsize2) and (n0hap1 != popsize1 and n1hap1 !=popsize1):
            #print i+1
            n_prvpos +=1
    return float(n_prvpos)/hap_dim[1]
    #return float(n_prvpos)/nsites

def Lsubseq(hap,popsize1,popsize2,nsites,snplist):
    tract=list()
    lsubseq=list()
    nsnp=hap.shape[1]
    for i in range(popsize1):
        for j in range(popsize1,(popsize1+popsize2),1):
            nshare=0
            ##compare two sequences
            for k in range(nsnp):
                if hap[i,k] == hap[j,k] :
                    if k==0:
                        nshare += (snplist[k]+1)
                    else:
                        nshare += snplist[k]-snplist[k-1]
                    if k == nsnp-1:
                        tract.append(nshare)
                        nshare=0
                else:
                    if nshare != 0:
                        tract.append(nshare)
                        nshare = 0
            lsubseq.append(max(tract))
            tract=list()
    Lsubseq=max(lsubseq)
    #return float(Lsubseq)/hap.shape[1]
    return float(Lsubseq)/nsites


def GetSP(popsize):
    startpoint=list()
    startpoint.append(0)
    s=0
    for i in range(len(popsize)-1):
        s+=popsize[i]
        startpoint.append(s)
    return startpoint

def Dict2Mat(rawdata,delim):
    speciesname=list()
    speciesinfo=dict()
    template=''
    Allmat=[]
    for key in sorted(rawdata.iterkeys()):
        ###process key
        r=re.compile('([a-b][1-2])([0-9]+)('+delim+')([A-B][1-2])')
        m=r.match(key)
        spn=m.group(1)
        if spn in speciesinfo.keys():
            speciesinfo[spn]+=1
        else:
            speciesname.append(spn)
            speciesinfo[spn]=1
        ### process sequence
        if len(template)==0:
            template=rawdata[key]
            Allmat=Gen01fromSeq(rawdata[key],template)
        else:
            #Update(template,rawdata[key])
            Allmat=np.vstack((Allmat,Gen01fromSeq(rawdata[key],template)))

    popsize=list()
    for spn in speciesname:
        popsize.append(speciesinfo[spn])

    classlabel=list()
    for i in range(len(speciesname)-1):
        for j in range(i+1,len(speciesname),1):
            if speciesname[i][0]==speciesname[j][0]:
                classlabel.append('-1')
            else:
                classlabel.append('+1')
    return popsize,Allmat,classlabel

def Seq2SNP(Allmat):
    nr,nc=Allmat.shape
    snplist=list()
    for i in range(nc):
        if sum(Allmat[:,i]) != 0:
            snplist.append(i)
    AllSNP=np.array(Allmat[:,snplist])
#    print AllSNP.shape
    return AllSNP,snplist

def Update(template,newtem):
    for i in range(len(template)):
        if (template[i]=='-' and newtem[i]!='-'):
            if i!=(len(template)-1):
                tmp=template[0:i]+newtem[i]+template[(i+1):]
            else:
                tmp=template[0:i]+newtem[i]
            template=tmp
    return template

def Gen01fromSeq(seq,template):
    hap=list()
    for i in range(len(template)):
        if seq[i]=='-':
            hap.append(2)
        elif seq[i]==template[i]:
            hap.append(0)
        else:
            hap.append(1)
    return hap

def GetSeq(line):
    header=line.split()[0]
    line=line.replace(header,'')
    line=line.replace(' ','')
    line=line.replace('\n','')
    return line



prefix=sys.argv[1]
#prefix='temp'
'''
theta=float(sys.argv[2])
tau=float(sys.argv[3])
if tau<0.0001:
    seqprefix=prefix+'_'+str(theta)+'_'+str(("%.5f" % tau))
else:
    seqprefix=prefix+'_'+str(theta)+'_'+str(tau)
'''
delim='\^'
rawdata=dict()
nbin=3
nsites=10000

#f1=open(seqprefix+'_seq.txt','rb')
#fout=open(seqprefix+'.sumstat','ab')
f1=open(prefix+'_seq.txt','rb')
fout=open(prefix+'.sumstat','ab')
for line in f1:
    if line[0:2]!='\n':
        header=line.split()[0]
    else:
        header=''
    indvName=RecogIndv(header,delim)
    if len(indvName)>0:
        seq=GetSeq(line)
        rawdata[indvName]=seq
    else:
        Process(fout,rawdata,delim,nbin,nsites)
        rawdata.clear()

Process(fout,rawdata,delim,nbin,nsites)
rawdata.clear()
#print rawdata

f1.close()
fout.close()


