#!/usr/bin/python

import sys
import os
import re
import glob

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

f = open(sys.argv[1],'r') #normalised otu table (floats)
				  
folder = sys.argv[2]

g = open(sys.argv[3],'w') 

filelist=glob.glob(folder+"/*")

filelist.sort(key=tokenize)

print filelist

c=0
t=0
#make collated habitat dict
data={}
habitats=[]

for i in filelist:
	f1=open(i,'r')
	
	otu = f1.readline().split(": ")[1].rstrip("\n")
	f1.readline()
	
	k = f1.readline().split("\t")
	pchit = float(k[0].split("%")[0])/100
	data[otu]={}
	
	for j in k[1:]:
		
		hab = j.split(":")[0]
		pchab=float(j.split(" ")[1].split("%")[0])/100
		
		data[otu][hab]=pchit*pchab
		if hab not in habitats:
			habitats.append(hab)

habitats.sort()	
		

samples=[]
sampledata={}
for i in f:
	k = i.split("\t")
	k[-1]=k[-1].rstrip("\n")
		
	if i[0]=="#": #get sample names
	
		for x in k[1:]:
			samples.append(x)
			sampledata[x]={}

	else:
		otu=k[0]
		c=-1
		for x in k[1:]:
			c=c+1
			sampledata[samples[c]][otu]=float(x)


title="#habitat\t"+"\t".join(str(p) for p in samples)+"\n"

g.write(title)

for x in habitats:

	g.write(x)
	
	for i in samples:
		
		v=0
		
		for j in data.keys(): #loop otus
			
			if x in data[j].keys():
				
				v = v + (data[j][x] * sampledata[i][j])
			
		g.write("\t"+str(v))
		

	g.write("\n")

for i in data.keys():

	print i, data[i].keys()







		
		
		
		
		
		
		
		
		
		
		
		