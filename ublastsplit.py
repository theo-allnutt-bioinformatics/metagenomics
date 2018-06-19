from Bio import SeqIO
import sys
import os
from StringIO import StringIO # Python 2


	#count reads

inputfile = sys.argv[1]
chunksize = int(str(sys.argv[2]))


count = SeqIO.index(inputfile, "fasta")

c= len(count)

#print c
a= int(c)
numchunks = (c/chunksize)

#print"Number of chunks=",numchunks

#sys.stdout.write(str(numchunks))

t=0
f=1
print"writing files"
outfile = "/OSM/HOME-MEL/all29c/scripts/tmp/ubin"+str(f)+".fna"
g=open(outfile,'w')
for i in SeqIO.parse(inputfile,"fasta"):
	
	t=t+1
	if t>=chunksize:
		f=f+1
		outfile = "/OSM/HOME-MEL/all29c/scripts/tmp/ubin"+str(f)+".fna"
		g=open(outfile,'w')
		t=0
	print"file",f,"record",t
	print i.seq
	raw_input()
	SeqIO.write(i,g,"fasta")
	
print
print "python files written"