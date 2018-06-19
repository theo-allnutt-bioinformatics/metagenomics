from Bio import SeqIO
import sys



f = open(sys.argv[1],'r') #locus tag list

h = open(sys.argv[2],'r') #genbank file

g= open(sys.argv[3],'w') #annotated tag list


record = SeqIO.to_dict(SeqIO.parse(h, 'genbank'))

a=record.keys()
a.sort()

tag=[]
for line in f:

	tag.append(line.rstrip("\n"))
		
	
for i in record.keys()
	data={}
	for feature in record[i].features:
		data[feature.qualifiers["locus_tag"][0]]=feature.qualifiers["exact_location"][0]
	
for t in tag:
	if t in data.keys():

		print tag, feature.qualifiers["exact_location"][0]
		g.write(tag+"\t"+feature.qualifiers["exact_location"][0]+"\n")
	else:
		print tag, "no hit"
		g.write(tag+"\t"+"no hit"+"\n")
				
			
	