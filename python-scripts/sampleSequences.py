# April 15, 2019
# Tanmay Tanna

# Ensure that the input fasta file has no newlines within the features. Use the following command if there are any:
# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < file.fa > out.fa 

from __future__ import division
from collections import Counter
import sys, os, argparse, operator, fuzzysearch, random


parser = argparse.ArgumentParser()

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-p', '--inFile', help='path to .fasta file.', dest='file_inFile')
required.add_argument('-l1', '--minlen', help='minimum length of sequence to be sampled', dest='min_len')
required.add_argument('-l2', '--maxlen', help='maximum length of sequence to be sampled', dest='max_len')
required.add_argument('-o', '--outPath', help='path to output directory.', dest='path_outPath')


## user inputs optional
optional.add_argument('-n', '--num', help='maximum total number of sequences to be sampled', dest='seq_num', default=0)
optional.add_argument('-n1', '--minseq', help='minimum number of sequences to be sampled from each feature in reference', dest='minseq', default=0)
optional.add_argument('-n2', '--maxseq', help='maximum number of sequences to be sampled from each feature in reference', dest='maxseq', default=100)
optional.add_argument('-gc1', '--minGC', help='minimum gc content in sequences', dest='mingc', default=30)
optional.add_argument('-gc2', '--maxGC', help='maximum gc content in sequences', dest='maxgc', default=70)
optional.add_argument('-ngc', '--nonGC', help='percentage probability of sequence to be accepted if it is outside of gc range. If this is >0, +-10 percent GC content of the limits of the GC range will be selected for twice this probability, so that the distribution is gaussian-like', dest='nonGC', default=0)
optional.add_argument('-z', '--numZero', help='percentage probability of features in reference to not have any sampling', dest='noseq', default=50)
optional.add_argument('-q', '--outName', help='name of output file', dest='outName', default='NULL')

args = parser.parse_args()

# assign arguments to variables

inFile = str(args.file_inFile)
outPath = str(args.path_outPath)+'/'
l1=int(args.min_len)
l2=int(args.max_len)
n1=int(args.minseq)
n2=int(args.maxseq)
n=int(args.seq_num)
if n==0:
    n = float("inf")
maxgc=int(args.maxgc)
mingc=int(args.mingc)
nongc=int(args.nonGC)
noseq=int(args.noseq)
outName = str(args.outName)
if 'NULL' in outName:
    outName = inFile.split("/")[-1]
    if outName.endswith('.fasta'):
        outName = outName[:-6]

n_sect=0
## De-comment if a 5' or 3' bias is to be introduced in sequence sampling
# n_sect=10 #change if the sequence is to be split into sub-sequences and a bias for sub-sequences from specific positions is to be introduced
# skew=5 # can be 5 for 5' or 0 for no skew

if ('.fasta' in inFile or '.fa' in inFile):
    
    # open inFile for reading/writing and report file being processed
    F = open(inFile,mode='rU')
    G = open(outPath+outName+'.randomseq.info.txt',mode='w')
    H = open(outPath+outName+'.randomseq.fasta',mode='w') 
    I =  open(outPath+outName+'.counts.tsv', mode = 'w')
    G.write("Sequence"+'\t'+"Sequence_Length"+'\t'+"GC_content"+'\n')
    I.write("Label"+'\t'+"Counts"+'\n')
    i=0 # number of sampled sequences
    geneCounts={}
    for L in F:
        L=L.strip()
    	if i>n:
    		break
    	if '>' in L:
            seqname=L.strip()
            continue

        ## This section splits the sequence into 'n_sect' sub-sequences to introduce 5' or 3' bias
        # if n_sect is not 0:
        #     if int(len(L)/n_sect)<l2:
        #         continue
        #     L_sect=[] # array of sub-sequences
        #     bias = [] # array of bias towards corresponding subsequence
        #     for sect in range(0,n_sect):
        #         L_sect.append(L[(sect*int(len(L)/n_sect)):((sect+1)*int(len(L)/n_sect))])
        #         if skew is 0:
        #             bias.append(100/n_sect)
        #         if skew is 5:
        #             if sect<((n_sect/2)+1):
        #                 bias.append(((100*(n_sect-1)/n_sect)+1) - (100*sect/n_sect))
        #             else:
        #                 sectdiff=sect-(n_sect/2)
        #                 bias.append(((100*(n_sect-1)/n_sect)+1) - 50 + (50*sectdiff/n_sect))

        ## end of section

        if n_sect is 0 and len(L)<l2:
            continue
            
        if noseq>0:
            k=random.randint(1,100)
            if k<noseq:
                continue

        nseq=random.randint(n1,n2) # randomly choosing a number of sequences to be picked from current input, this is the max final "count" mapping to this input
        for num in range(1, nseq):
            l=random.randint(l1,l2)
            if n_sect is 0:
                L_skew=L
            # else:
            #     for sect in range(0,n_sect):
            #         prob = random.randint(1,100)
            #         if prob<bias[sect]:
            #             L_skew=L_sect[sect]
            #             continue
            #         if not 'L_skew' in locals():
            #             L_skew=L_sect[0]

            start=random.randint(0,len(L_skew)-l)
            seq = L_skew[start:start+l]
            count = Counter(seq)
            gc = (count['G']+count['C'])*100/len(seq)
            if gc>maxgc or gc<mingc:
                if nongc is 0:
                    continue
                if gc-10>maxgc or gc+10<mingc: 
                    prob=random.randint(1,100)
                    if prob>nongc:
                        continue
                    if seqname in geneCounts:
                        geneCounts[seqname]+=1
                    else:
                        geneCounts[seqname]=1

                    G.write(str(seq)+'\t'+str(len(seq))+'\t'+str(gc)+'\n')
                    H.write(str(seqname)+'_sampledSequence_'+str(i)+'\n'+str(seq)+'\n')
                    i+=1
                else:
                    prob=random.randint(1,100)
                    if prob>nongc*2:
                        continue
                    if seqname in geneCounts:
                        geneCounts[seqname]+=1
                    else:
                        geneCounts[seqname]=1

                    G.write(str(i)+'\t'+str(len(seq))+'\t'+str(gc)+'\n')
                    H.write(str(seqname)+'_spacer_'+str(i)+'\n'+str(seq)+'\n')

            if seqname in geneCounts:
                geneCounts[seqname]+=1
            else:
                geneCounts[seqname]=1

            G.write(str(seq)+'\t'+str(len(seq))+'\t'+str(gc)+'\n')
            H.write(str(seqname)+'_spacer_'+str(i)+'\n'+str(seq)+'\n')
            i+=1

    for seqname in geneCounts:
        I.write(str(seqname)+'\t'+str(geneCounts[seqname])+'\n')

    F.close()
    G.close()
    H.close()
    I.close()
  



