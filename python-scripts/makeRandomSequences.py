import sys, os, argparse, random


parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-l', '--len', help='length of random sequences', dest='length')
required.add_argument('-n', '--num', help='number of random sequences', dest='number')
args = parser.parse_args()
k=int(args.number)
l=int(args.length)
dna = ["A","G","C","T"]
database=[]
while len(database)<k:
	randseq=""
	for i in range(0,l):
		randseq+=random.choice(dna)
	if randseq not in database:
		database.append(randseq)

print(database)




