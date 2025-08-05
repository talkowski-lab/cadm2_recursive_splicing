#!/apps/lab/miket/anaconda/4.0.5/envs/dg520_py2/bin/python

from __future__ import division
import shlex, subprocess
import commands
import os, sys
import collections


def make(a):
    fa=a+"/cadm2_iN_all_transcripts.fa"
    #ap=a+"/alltrans_appearance.txt"
    #dic={}
    #d=open(ap,"r").readlines()
    #for x in d:
    #    t=x[:-1].split("\t")
    #    trans=t[0]
    #    tis=t[1].split(",")
    #    for y in tis:
    #        if (y[7:] in dic) == False:
    #            dic[y[7:]]=[trans]
    #        else:
    #            dic[y[7:]].append(trans)
    
    #print len(dic)
    d1=open(fa,"r").readlines()
    cmd="mkdir "+a+"/rsem_index_Uniform/"
    stat,out=commands.getstatusoutput(cmd)
    #for x in dic:
    o1=a+"/cadm2_iN_all_transcripts.fa"
    o2=a+"/rsem_index_Uniform/gene_trans.map"
    f1=open(o1,"w")
    f2=open(o2,"w")
    #for y in dic[x]:
    i=0
    while i<len(d1):
        #if d1[i][1:-1]==y:
        f1.write(d1[i])
        f1.write(d1[i+1])
        f2.write("CADM2\t"+d1[i][1:-1]+"\n")
        i+=2
    f1.close()
    f2.close()
    cmd="rsem-prepare-reference --bowtie2 --transcript-to-gene-map "+o2+" "+o1+" "+o1
    stat,out=commands.getstatusoutput(cmd)
    
    print(cmd)
    return 1
    


make("/data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_index")
