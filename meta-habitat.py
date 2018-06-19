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
				  
folder = sys.argv[1]

g = open(sys.argv[2],'w') 

filelist=glob.glob(folder+"/*")

filelist.sort(key=tokenize)

print filelist

c=0
t=0
for i in filelist:
	f=open(i,'r')
	

	
	otu = f.readline().split(": ")[1].rstrip("\n")
	f.readline()
	tophit=f.readline()
	p1 = tophit.split("\t")[0]
	p2=	tophit.split("\t")[1]
	a=float(p1.split("%")[0])/100
	b=float(p2.split(": ")[1].split("%")[0])/100
	#gut: 100.00%
	h1=p2.split(": ")[0]
	pc=a*b*100
	
	print otu+" "+h1+" "+str(pc)
	g.write(otu+"\t"+h1+" "+str(pc)+"\n")