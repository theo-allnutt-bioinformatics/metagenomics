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

level=sys.argv[4]
delim=sys.argv[5] #taxonomy output delimiter
file_delim=sys.argv[6] #delimter in file name to obtain proper sample name
g=open(outfile,'w')

filelist=glob.glob(folder+"/*.%s" %fmt)

filelist.sort(key=tokenize)


print filelist

data={}
allspecies=[]

for i in filelist:

	file1 = open(i,'r')
	filename = i.split("/")[-1].split("file_delim")[0]
	
	data[filename]={}
	
	for line in file1:
		if line[0]<>"#":
			sp1 = line.split("\t")[0]
			sp=sp1.split(delim)[-1]
			
			lev = sp.split("__")[0]
			
			freq= float(line.split("\t")[1].rstrip("\n").rstrip("\r"))
			
			if sp1=="unclassified":
				if sp1 not in allspecies:
					allspecies.append(sp1)
				
				if sp1 not in data[filename]:
					
					data[filename][sp1]=""
					
				data[filename][sp1]=freq	
			
			if lev==level:
			
				if sp1 not in allspecies:
					allspecies.append(sp1)
				
				if sp1 not in data[filename]:
					
					data[filename][sp1]=""
					
				data[filename][sp1]=freq
		
	file1.close()
			
allspecies.sort()

print "\n".join(str(x.split(delim)[-1]) for x in allspecies)	

g.write("\t"+"\t".join(str(x) for x in allspecies)+"\n")

sorted_data=[]
for i in data.keys():
	sorted_data.append(i)
	
sorted_data.sort(key=tokenize)
print sorted_data
sums={}
sum_names=[]
for name in sorted_data:

	out = name+"\t"
	sum_name=name.split(".")[0] #get first part of file name of bins to make a sum over all bins for that sample
	if sum_name not in sums.keys():
		sums[sum_name]={}
	if sum_name not in sum_names:
		sum_names.append(sum_name)
		for i in allspecies: #make sums dict
			sums[sum_name][i]=0
		
	for i in allspecies:
		
		if i in data[name].keys():
			out = out + str(data[name][i])+"\t"
			sums[sum_name][i]=sums[sum_name][i]+data[name][i]
			
		else:
			out = out + "0" +"\t"
	
	g.write(out+"\n")
#record the summary information
sum_names.sort(key=tokenize)

for sample in sum_names:
	g.write(sample)
	for i in allspecies:
	
		g.write("\t"+str(sums[sample][i]))
		
	g.write("\n")
	
g.close()
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			