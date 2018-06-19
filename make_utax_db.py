#!/usr/bin/python

import sys

from Bio import SeqIO

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[3],'w')
h = open(sys.argv[2],'r')
level=int(sys.argv[4])

print "loading taxonomy"

tax={}
labels1=["d:","p:","c:","o:","f:","g:","s:"]

for i in h:
	
	i=i.replace(",","_")
	i=i.replace(" ","_")
	acc = i.split("\t")[0].rstrip(" ")
	
	lineage = i.split("\t")[1].rstrip("\n").split(";")
	
	#if "Viridiplantae" not in lineage:
	
	lin=[";tax="]	
	s=-1
	for j in lineage:
		s=s+1
		
		if "NA" not in j and "uncultured" not in j and "sp." not in j and "incertae" not in j and "incerti" not in j and "unknown" not in j and "unclassified" not in j and "chloroplast" not in j:

			lin.append(labels1[s]+j+",")
			
	r=1
	t=0
	for p in labels1[:level]:
		t=t+1
		try:
			if p not in lin[t]:
				r=0
		except:
			r=0
		
	if r==1 and len(lin)>=level:
		
		tax[acc]="".join(str(v) for v in lin[:level+1])[:-1]+";"
	

c=0
n=0
for x in SeqIO.parse(infile,'fasta'):
	id=str(x.id).split(" ")[0]
	
	seq=str(x.seq)
	
	try:
		outfile.write(">"+id+tax[id]+"\n"+seq+"\n")
		n=n+1
	except:
		c=c+1
outfile.close()
print n,"taxonomies found", c, "not found"

'''
print "sorting by fullest taxonomies"

tosort=open(sys.argv[3],'r')

seqs=SeqIO.to_dict(SeqIO.parse(tosort,'fasta'))

t=[]
for i in seqs.keys():
	t.append((i,i.count(",")))
t.sort(key=lambda p: -p[1])

outfile2=open(sys.argv[3],'w')
for i in t:
	
	outfile2.write(">"+str(seqs[i[0]].description)+"\n"+str(seqs[i[0]].seq)+"\n")
'''	
print "done"

