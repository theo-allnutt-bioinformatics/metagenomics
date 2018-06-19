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
				  
#makes an otu table from insput species lists

folder = sys.argv[1] #working folder
outfile=sys.argv[2] #output file
fmt=sys.argv[3] #file extension minus . dot

g=open(outfile,'w')

filelist=glob.glob(folder+"/%s" %fmt)

filelist.sort(key=tokenize)


print filelist

data={}
allspecies=[]
filenames=[]
for i in filelist:

	file1 = open(i,'r')
	filename = i.split("/")[-1].split(".")[0]
	filenames.append(filename)
	data[filename]={}
	
	for line in file1:
		if line[0]<>"#":
			sp1 = line.split("\t")[0]
			
			freq= line.split("\t")[1].rstrip("\n").rstrip("\r")

			if sp1 not in allspecies:
				allspecies.append(sp1)
				
			if sp1 not in data[filename]:
					
				data[filename][sp1]=""
			else:
			
				data[filename][sp1]=freq
		
	file1.close()
			
allspecies.sort()

print "\n".join(str(x) for x in allspecies)	

g.write("\t"+"\t".join(str(x) for x in filenames)+"\n")


for i in allspecies:
	out = i+"\t"
	for name in filenames:
		
		if i in data[name].keys():
			out = out + str(data[name][i])+"\t"
			
		else:
			out = out + "0" +"\t"
	
	g.write(out+"\n")

