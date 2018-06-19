#!/usr/bin/python
import sys
import os
import re
import glob
import subprocess

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))
				  
folder = sys.argv[1] #working folder

g=sys.argv[2] #output folder

filelist=glob.glob(folder+"/*" )

filelist.sort(key=tokenize)

print filelist

for i in filelist:
	print i
	name =i.split("/")[-1]
	
	outname=g+"/"+name+".txt"
	
	p1=subprocess.Popen("~/bin/metaphlan1/metaphlan.py %s --blastdb ~/bin/metaphlan1/blastdb/mpa -o %s --input_type multifasta --nproc 24 --tax_lev a --tmp_dir ./ --no_map" %(i,outname),shell=True).wait()
	

			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			