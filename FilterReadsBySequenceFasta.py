# October 5, 2018
# Tanmay Tanna



from __future__ import division
import sys, os, argparse, operator, fuzzysearch


parser = argparse.ArgumentParser()


required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-p', '--inFile', help='path to .fasta file.', dest='file_inFile')
required.add_argument('-s', '--seq', help='sequence to search for and remove', dest='sequence')
required.add_argument('-o', '--outPath', help='path to output directory.', dest='path_outPath')
required.add_argument('-f', '--filter', help='filter "in" or "out"', dest='filter')
required.add_argument('-m', '--misMatch', help='number of allowed mismatches', dest='mismatch')
required.add_argument('-n', '--readNumber', help='maximum number of reads in destination file', dest='readNumber')
args = parser.parse_args()

# assign arguments to variables

inFile = str(args.file_inFile)
outPath = str(args.path_outPath)+'/'
sequence=str(args.sequence)
outName = inFile.split("/")[-1]
readNumber = args.readNumber
mismatch = args.mismatch
filter = str(args.filter)
if outName.endswith('.fasta'):
    outName = outName[:-6]


if ('.fasta' in inFile or '.fa' in inFile):
    
    # open inFile for reading/writing 
    F = open(inFile,mode='rU')
    G = open(outPath+outName+'.filtered.fasta',mode='w')

    rawReads=0
    filteredReads=0
    ReadsWithSequence=0
    ReadsWithoutSequence=0
    VeryShortReads=0
    i=0

    os.system(str("echo '##################################################'"))
    os.system(str('echo '+"'"+inFile+' accepted for processing'+"'"))


    for L in F: # loop through reads in file
        if filteredReads>readNumber:
            break
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            rawReads+=1
            continue
        L=L.strip()

        if len(L) < len(sequence):
            VeryShortReads+=1
            continue

        if sequence in L:
            ReadsWithSequence+=1
            if filter=="in":
                filteredReads+=1
                G.write(readName+'\n'+L+'\n')
                
        if mismatch>0:
            d = fuzzysearch.find_near_matches(sequence1, L, max_l_dist = mismatch)

        if d:
            ReadsWithsequence1+=1
            continue
            if filter=="in":
                filteredReads+=1
                G.write(readName+'\n'+L+'\n')
        else:
            ReadsWithoutSequence+=1
            if filter=="out":
                filteredReads+=1
                G.write(readName+'\n'+L+'\n')



    F.close()
    G.close()

    # os.system(str('echo '+str('rawReads = '+str(rawReads)+' filteredReads = '+str(filteredReads)+' ReadsWithsequence1 = '+str(ReadsWithsequence1)+' VeryShortReads = '+str(VeryShortReads)))
