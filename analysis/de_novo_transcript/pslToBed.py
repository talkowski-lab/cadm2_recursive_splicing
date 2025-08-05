def construct(b,c):
    #d1=open(a,"r").readlines()
    #dic={}
    #for x in d1:
    #    t=x[:-1].split("\t")
    #    dic[t[0]]=t[2]
    d=open(b,"r").readlines()
    q=[]
    f=open(c,"w")
    #f.write("X\t70660351\t70662978\tSVA\t0\t-\t70660351\t70662978\t0\t1\t2627,\t0,\n")
    for x in d:
        t=x.split("\t")
        n=t[18]
        st=t[9]
        size=t[-3]
        cand=t[-1][:-2].split(",")
        cand1=map(int,cand)
        #cand0=t[-2]
        name=t[0]
        start=int(t[16])
        end=int(t[17])
        cand2=[y-start for y in cand1] 
        start=str(start)     #0-based [a,b)
        end=str(end)         #0-based [a,b)
        cand0=",".join(map(str,cand2))+","
        line="3\t"+start+"\t"+end+"\t"+name+"\t0\t"+st+"\t"+start+"\t"+end+"\t0\t"+n+"\t"+size+"\t"+cand0+"\n"
        f.write(line)
                    
    f.close()
    
    return 1

#construct('new_transcripts_IGV.psl','final_transcripts_merged.bed')
#construct('iN_0.8_merged.psl','test.fasta.bed')
#construct("iN_all_final.psl","iN_all_final.bed")
#construct("iPSC_all_final.psl","iPSC_all_final.bed")
construct("test.psl","test.bed")
