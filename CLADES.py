__author__ = 'pjweggy'

import re
import numpy as np
from itertools import chain
#from SumStat import *
#import math
import os
from subprocess import *
import sys

def RecogIndv(line,delim):
    if re.match('[A-Za-z]+[0-9]+'+delim+'[A-Za-z0-9]+',line):
        return line
    else:
        return ''


def Process(rawdata,delim,nbin,SS):
    if len(rawdata) != 0:
        [popsize,Allmat,spn]=Dict2Mat(rawdata,delim)
        nsites=Allmat.shape[1]
        Allmat,snplist=Seq2SNP(Allmat)
        startpoint=GetSP(popsize)
        for i in range(len(popsize)-1):
            for j in range(i + 1, len(popsize)):
                pos = list(chain(range(startpoint[i], startpoint[i]+popsize[i]), range(startpoint[j],startpoint[j]+popsize[j])))
                hap=Allmat[pos,:]
                sumstat=CollectSS2pop(hap,popsize[i],popsize[j],nbin,nsites,snplist)
                key=spn[i]+'-'+spn[j]
                if key in SS:
                    SS[key]=SS[key]+AddSS(sumstat)
                else:
                    SS[key]=AddSS(sumstat)

    return SS
        #WriteSS(fout,sumstat1,classlabel)

def AddSS(sumstat):
    nsam=len(sumstat)
    s=''
    for i in range(nsam):
        s=s+str(i+1)+':'+str("%.6f" % sumstat[i])+' '
    s=s+'\n'
    return s

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
        sfs=np.zeros((hap_dim[0]+1) // 2)
    else:
        sfs=np.zeros(hap_dim[0] // 2+1)
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
    window=nlen // nbin
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
                    if k==0:
                        nshare+=snplist[0]
                    else:
                        nshare+=snplist[k]-snplist[k-1]-1
                    tract.append(nshare)
                    if k==nsnp-1:
                        tract.append(nsites-snplist[-1])
                    nshare = 0
            #print(tract)
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
    for key in sorted(rawdata.keys()):
        ###process key
        #r=re.compile('([a-d][1-2])([0-9]+)('+delim+')([A-D][1-2])')
        r=re.compile('([A-Za-z]+[0-9]+)('+delim+')([A-Za-z0-9]+)')
        m=r.match(key)
        spn=m.group(3)
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

#    classlabel=list()
#    for i in range(len(speciesname)-1):
#        for j in range(i+1,len(speciesname),1):
#            if speciesname[i][0]==speciesname[j][0]:
#                classlabel.append('-1')
#            else:
#                classlabel.append('+1')

    return popsize,Allmat,speciesname

def Seq2SNP(Allmat):
    nr,nc=Allmat.shape
    snplist=list()
    for i in range(nc):
        if sum(Allmat[:,i]) != 0:
            snplist.append(i)
    if len(snplist)==0:
        print('Error Message:Gene has no variant! Please delete this gene or concatenate multiple genes.')
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

def GetSPN(Res):
    spn=[]
    for key in sorted(Res.keys()):
        spn1=key.split('-')[0]
        spn2=key.split('-')[1]
        if not (spn1 in spn):
            spn.append(spn1)
        if not (spn2 in spn):
            spn.append(spn2)
    return spn

def CoalSPN(CurC,spn1,spn2):
    NewC=list()
    NewC.append(spn1+','+spn2)
    for spn in CurC:
        if spn !=spn1 and spn !=spn2:
            NewC.append(spn)
    return NewC

def ProbSS(splist1,Res):
    prob=1
    if len(splist1)==1:
        prob=1
    else:
        for i in range(len(splist1)-1):
            for j in range(i+1,len(splist1),1):
                if (splist1[i]+'-'+splist1[j]) in Res:
                    prob=prob*Res[splist1[i]+'-'+splist1[j]][1]
                elif (splist1[j]+'-'+splist1[i]) in Res:
                    prob=prob*Res[splist1[j]+'-'+splist1[i]][1]
    return prob

def ProbDS(splist1,splist2,Res):
    prob=1
    for spn1 in splist1:
        for spn2 in splist2:
            if (spn1+'-'+spn2) in Res:
                prob=prob*Res[spn1+'-'+spn2][0]
            elif (spn2+'-'+spn1) in Res:
                prob=prob*Res[spn2+'-'+spn1][0]
    return prob


def CompProb(CurC,Res):
    prob=1
    for i in range(len(CurC)-1):
        splist1=CurC[i].split(',')
        prob=prob*ProbSS(splist1,Res)
        for j in range(i+1,len(CurC),1):
            splist2=CurC[j].split(',')
            prob=prob*ProbDS(splist1,splist2,Res)
    if len(CurC)==1:
        splist1=CurC[0].split(',')
        prob=prob*ProbSS(splist1,Res)
    return prob

prefix=sys.argv[1]
#prefix='test'

delim='\^'
rawdata=dict()
nbin=3
SS=dict()

##1. Compute Summary Statistics for data
f1=open(prefix+'_seq.txt','r')

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
        SS=Process(rawdata,delim,nbin,SS)
        rawdata.clear()

SS=Process(rawdata,delim,nbin,SS)
rawdata.clear()
f1.close()

##2.Species delimitation for pairwise pops
#model='/Users/pjweggy/Documents/academy/SpecD/code/model2/All'
#path1='/Users/pjweggy/Downloads/apps/libsvm-3.22/'
model=sys.argv[2]

if len(sys.argv)>3:
    path1=sys.argv[3]
else:
    path1=''

#path1=''
Res=dict()
for key in sorted(SS.keys()):
    try:
        os.remove(key+'.sumstat')
        os.remove(key+'.sumstat.scale')
    except OSError:
        pass
    fss=open(key+'.sumstat','ab')
    fss.write(SS[key].encode('utf-8'))
    fss.close()
    cmd1='{2}svm-scale -r {0}.range {1} > {1}.scale'.format(model,key+'.sumstat',path1)
    cmd2='{3}svm-predict -b 1 -q {1}.scale {0}.sumstat.scale.model {2}.out'.format(model,key+'.sumstat',key,path1)
    #print cmd1
    #print cmd2
    Popen(cmd1, shell = True, stdout = PIPE).communicate()
    while not os.path.exists(key+'.sumstat.scale'):
        continue
    Popen(cmd2, shell = True, stdout = PIPE).communicate()
    while not os.path.exists(key+'.out'):
        continue
    Out=np.loadtxt(key+'.out',comments='labels')
    res=np.mean(Out,axis=0)[1:3]
    Res[key]=res

try:
    os.remove(prefix+'.out')
except OSError:
    pass
fout=open(prefix+'.out','w')
fout.write('labels +1 -1\n')
for key in sorted(Res.keys()):
    fout.write(key+' '+str(("%.4f" % Res[key][0]))+' '+str(("%.4f" % Res[key][1]))+'\n')
fout.close()



##3.Compute probability of best assignment
'''
f=open('test.out','rb')
Res=dict()
for line in f:
    if line[0:5]=='label':
        continue
    else:
        key=line.split()[0]
        spn=[float(line.split()[1]),float(line.split()[2])]
        Res[key]=spn
f.close()
'''

spn=GetSPN(Res)
CurC=spn
curprob=CompProb(CurC,Res)
coal=1

#print CurC,curprob
while coal:
    coalc=0
    for i in range(len(CurC)-1):
        for j in range(i+1,len(CurC),1):
            NewC=CoalSPN(CurC,CurC[i],CurC[j])
            prob=CompProb(NewC,Res)
            print('{0} {1}'.format(NewC,prob))
            if prob>curprob:
                CurC_temp=NewC
                curprob=prob
                coalc += 1
    if coalc>0:
        CurC=CurC_temp
    else:
        coal=0

print('The Best Assignment of species clusters are:')
print('{0} {1}'.format(CurC,curprob))
