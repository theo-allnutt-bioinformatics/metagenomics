import sys
import os
import subprocess
import re
from random import randrange

#taxonfetch.py theo allnutt 2014. Gets the taxonomy ID number from a local copy of the NCBI GI/taxonomy list file #and retrieves the scientific lineage field from nodes and names file using tax_trace.pl. ID and lineage are then appended to the end of the submitted #blast output (tab format, 6) file.

#this software has no guarantees of working or doing anything and should not be used for anything.


inputfile = sys.argv[1] #tab format blast output (6) GIs in second column
outputfile = sys.argv[2] #taxa are added to end columns - name of appended file
taxdb = sys.argv[3] #specify "p" or "n" for protein or nucleotide taxon list


f = open(inputfile,'r')
g = open(outputfile,'w')
subprocess.Popen("rm -rf ~/tmp.txt",shell=True).wait()


c=0
dataout = ""
t=0
 #get the GIs
gi=[]
dataout=[]
ids=[]

r2=str(randrange(100000))
split_place=1
for line2 in f:
	c=c+1
	
	gi.append(line2.split("\t")[1].split('|')[int(sys.argv[4])])
print str(c)+"   GIs"
#print gi
c=0
gitax=[]

names =[]
t=0

#CHANGE THESE PATHS TO THOSE FOR YOUR GI - TAX_ID FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if taxdb == "p":
	print "reading protein taxonomy \n"
	taxaid = "/OSM/MEL/AAHL_MicBio/data/theo/data/db/gi_taxid_prot.dmp"
if taxdb == "n":
	print "reading nucleotide taxonomy"
	taxaid = "/OSM/MEL/AAHL_MicBio/data/theo/data/db/gi_taxid_nucl.dmp"

print "searching.."
size = long(os.path.getsize(taxaid))

print "TaxID File Size: ", size
t = open(taxaid,'r')

for gis in gi: #binary search #########################################
	c=c+1
	#print gis
	found=False
	offset=0
	chunk=size
	pos=chunk/2
	while found == False and chunk>0:
		#print "posn: ", pos
		chunk = chunk/2
		#print "chunk:", chunk
		#print "offset: ", offset
		t.seek(pos)
		t.readline()
		entry = t.readline().split("\t")
		#print entry[0]
		if entry[0]:
			filegi = int(entry[0])
			filetax = int(entry[1].rstrip("\n"))
		
		#print gis, filegi, filetax
		#raw_input()
		
		if filegi == int(gis):
			answer = filetax
			#print filegi, "  FOUND"
			found = True
			#print 'chunk:',chunk,'   offset:',offset,'   posn:',pos
			print c,":",filegi, answer
		elif filegi > int(gis):
			pos = offset +(chunk/2)
		
		elif filegi < int(gis):
			offset = offset+chunk
			pos = pos + (chunk/2)
	
	if found == False:
		answer = "no taxonomy" #"32644"
		print c,":",filegi, answer
	
	
	gitax.append(str(answer))
	
#print gitax	#list of Gi's taxids
cc1=-1
v1=-1
p=""

#save gitax to tempfile

tmp1 = open('/OSM/HOME-MEL/all29c/scripts/tmp2/%stmp1.txt' %r2,'w')#change this to your preferred tempfile location

cnt1=-1

for i in gitax:
	cnt1=cnt1+1
	tmp1.write(str(gi[cnt1])+"\t"+str(i)+"\n")
	
tmp1.close()

#change this line to your location for tax_trace.pl and nodes and names and temp files
p1=subprocess.Popen("perl ~/scripts/tax_trace.pl ~/d/db/nodes.dmp  ~/d/db/names.dmp ~/scripts/tmp2/%stmp1.txt ~/scripts/tmp2/%stmp2.txt " %(r2,r2),shell=True).wait()
	

f.seek(0)
#change line to temp file location
p3=open("/OSM/HOME-MEL/all29c/scripts/tmp2/%stmp2.txt" %r2,'r')
p4=[]



for i in p3:
	names.append(i.rstrip("\n")) #list of lineages
	
#print "\n".join(str(x) for x in names)
#print len(names)
n=-1
for line4 in f:
	
	n=n+1
	
	
	if len(names[n])==1:
		names[n][0]="unclassified"
	
	lineage=line4.rstrip('\n')+"\t"+gitax[n].rstrip("\n")+"\t"+names[n]+"\n"
	#print lineage
	g.write(lineage)	


#change line to temp file location
#subprocess.Popen("rm -f ~/scripts/tmp2/%stmp1.txt" %r2,shell=True).wait()
#subprocess.Popen("rm -f ~/scripts/tmp2/%stmp2.txt" %r2,shell=True).wait()
g.close()
f.close()
