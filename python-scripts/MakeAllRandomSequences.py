import random
k=4 # replace 4 with the length of random sequence desired
dna = ["A","G","C","T"]
database=[]
while len(database)<4**k: # replace 4**k with any number if you want a particular number of sequences
	randseq=""
	for i in range(0,k):
		randseq+=random.choice(dna)
	if randseq not in database:
		database.append(randseq)

print(database)




