from Bio import SeqIO
import sys
import os
from StringIO import StringIO # Python 2
import subprocess


	#count reads

inputfile = sys.argv[1]
outputfile=sys.argv[2]
chunksize = int(str(sys.argv[3])) #number of reads in each division of the input file
threads = int(sys.argv[4])
fmt = sys.argv[5] #input format


subprocess.Popen("rm -rf ~/scripts/tmp2", shell=True).wait()
subprocess.Popen("mkdir ~/scripts/tmp2", shell=True).wait()

out1=open('/OSM/HOME-MEL/all29c/scripts/tmp2/finaloutput.ublast','w')

count = SeqIO.index(inputfile, fmt)

c= len(count)

#print c
a= int(c)
print "num reads=",a
if chunksize>a:
	chunksize=a
numchunks = (a/chunksize)

#print"Number of chunks=",numchunks

#sys.stdout.write(str(numchunks))

t=0
f=1
print"Writing split files..."
outfile = "/OSM/HOME-MEL/all29c/scripts/tmp2/ubin"+str(f)+".fna"
g=open(outfile,'w')
for i in SeqIO.parse(inputfile,fmt):
	
	t=t+1
	if t>chunksize:
		f=f+1
		outfile = "/OSM/HOME-MEL/all29c/scripts/tmp2/ubin"+str(f)+".fna"
		g=open(outfile,'w')
		t=1
	#print"file",f,"record",t
	#print i.seq
	#raw_input()
	SeqIO.write(i,g,"fasta")
	
print "Number of chunks=",numchunks
print "Files written"
print "Multithreading ublast"

#numt = number of times to run multithread 
p1=[];p2=[];p3=[];p4=[];p5=[]

q=0

numreps = int(numchunks/threads)
if numreps<1:
	numreps=1
	threads=numchunks

for setp in range(1,threads+2):
	p1.append("")
	p2.append("")
	p3.append("")
	p4.append("")
	p5.append("")

	
	print len(p1),"lenp1"
	
for rep in range(1,numreps+1):
	
	for x in range(1,threads+1):
		
		q=q+1 #file number
		
		p1[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/greengenes/gg1.udb -evalue 1e-10 -strand both -id 0.97 -accel 0.5 -maxaccepts 20 -maxrejects 200 -maxhits 1 -threads 1 -userout ~/scripts/tmp2/1ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		
		p2[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/greengenes/gg2.udb -evalue 1e-10 -strand both -id 0.97 -accel 0.5 -maxaccepts 20 -maxrejects 200 -maxhits 1 -threads 1 -userout ~/scripts/tmp2/2ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		
		p3[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/greengenes/gg3.udb -evalue 1e-10 -strand both -id 0.97 -accel 0.5 -maxaccepts 20 -maxrejects 200 -maxhits 1 -threads 1 -userout ~/scripts/tmp2/3ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		
		p4[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/greengenes/gg4.udb -evalue 1e-10 -strand both -id 0.97 -accel 0.5 -maxaccepts 20 -maxrejects 200 -maxhits 1 -threads 1 -userout ~/scripts/tmp2/4ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		
		p5[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/greengenes/gg5.udb -evalue 1e-10 -strand both -id 0.97 -accel 0.5 -maxaccepts 20 -maxrejects 200 -maxhits 1 -threads 1 -userout ~/scripts/tmp2/5ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		
	
	#wait for n threads to finish
	for x1 in range(1,threads+1):
		print "waiting for usearch to finish........."
		p1[x1].wait(); p2[x1].wait();p3[x1].wait(); p4[x1].wait();p5[x1].wait()
	
	for x2 in range(1,threads+1):		
		
		print "concatenating.....", x2
			
		subprocess.Popen("cat ~/scripts/tmp2/1ubin%dout.txt ~/scripts/tmp2/2ubin%dout.txt ~/scripts/tmp2/3ubin%dout.txt ~/scripts/tmp2/4ubin%dout.txt ~/scripts/tmp2/5ubin%dout.txt >> ~/scripts/tmp2/finaloutput.ublast" %(x2,x2,x2,x2,x2), shell=True).wait()
	
#subprocess.Popen("rm -rf ~/scripts/tmp2/ubin* ~/scripts/tmp2/1ubin* ~/scripts/tmp2/2ubin* ~/scripts/tmp2/3ubin* ~/scripts/tmp2/4ubin* ~/scripts/tmp2/5ubin* ", shell=True).wait()

print"adding taxonomy.." #taxon id is given by green genes instad of GI
subprocess.Popen("python ~/scripts/taxonfetchgreen.py ~/scripts/tmp2/finaloutput.ublast %s" %(outputfile), shell=True).wait()

#############otufreqgreen.py

f2 = open(outputfile,'r')
g1=open(outputfile+".family",'w')
g2=open(outputfile+".species",'w')
g3=open(outputfile+".genera",'w')

species=[]
genus=[]
family=[]
k=[]
ksp=""
kgn=""
kfm=""
c=0
print "Getting taxonomy.."
for line in f2:
	c=c+1
	print c
	#k=line.rstrip("\n").split("\t")[13].split(" ")
	k=line.rstrip("\n").split("\t")
	
	#get spp if present
	kfm = str(k[17])
	if kfm<>"":
		if kfm[0]=="[":
			kfm=kfm[1:]
		if kfm[:-1]=="]":
			kfm=kfm[:-1]
		family.append(kfm)
	else:
		family.append("unknown")
	
	kgn = str(k[18])
	if kgn<>"":
		if kgn[0]=="[":
			kgn=kgn[1:]
		if kgn[:-1]=="]":
			kgn=kgn[:-2]
		genus.append(kgn)
	else:
		genus.append("unknown")
	
	ksp = str(k[19])
	if ksp<>"":
		if ksp[0]=="[":
			ksp=ksp[1:]
		if ksp[:-1]=="]":
			ksp=ksp[:-1]
		species.append(kgn+" "+ksp)
	else:
		species.append("unknown")
	
	
	print str(k[0]), kfm, kgn, ksp
	
species_set=set(species)
family_set=set(family)
genus_set=set(genus)
print "Counting taxa.."
for i in species_set:
	g2.write(i+"\t"+str(species.count(i))+"\n")
for i in genus_set:
	g3.write(i+"\t"+str(genus.count(i))+"\n")
for i in family_set:
	g1.write(i+"\t"+str(family.count(i))+"\n")

f2.close()
g1.close()
g2.close()
g3.close()