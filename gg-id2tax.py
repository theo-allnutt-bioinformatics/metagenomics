

import sys



f=open(sys.argv[1],'r')

g=open(sys.argv[2],'r')

h=open(sys.argv[3],'w')

tax_map={}
for x in g:
	

	tax_map[x.split("\t")[0]]=x.split("\t")[1].rstrip("\n").rstrip("\t")
		
	

for j in f:
	line = j.split("\t")
	if j[0]=="#":
		h.write("\t".join(str(a) for a in line))
	else:	
		id=j.split("\t")[0].split(".")[0]
	
		try:
			h.write(tax_map[id]+"\t"+"\t".join(str(a) for a in line[1:]))
		except:
			print "not found", i
			h.write(tax_map[id]+"\t".join(str(a) for a in line))
			
			
			