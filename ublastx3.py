from Bio import SeqIO
import sys
import os
from StringIO import StringIO # Python 2
import subprocess


	#count reads

inputfile = sys.argv[1]
chunksize = int(str(sys.argv[2]))
threads = int(sys.argv[3])

subprocess.Popen("rm -rf ~/scripts/tmp2", shell=True).wait()
subprocess.Popen("mkdir ~/scripts/tmp2", shell=True).wait()

out1=open('/OSM/HOME-MEL/all29c/scripts/tmp2/finaloutput.ublast','w')

count = SeqIO.index(inputfile, "fasta")

c= len(count)

#print c
a= int(c)
numchunks = (c/chunksize)

#print"Number of chunks=",numchunks

#sys.stdout.write(str(numchunks))

t=0
f=1
print"writing files"
outfile = "/OSM/HOME-MEL/all29c/scripts/tmp2/ubin"+str(f)+".fna"
g=open(outfile,'w')
for i in SeqIO.parse(inputfile,"fasta"):
	
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
	
print "number of chunks=",numchunks
print "files written"
print "multithreading ublast"

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
		
		p1[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted1.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/1ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		print "p 1chunk=",q
		print "usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted1.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/1ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits"
		p2[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted2.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/2ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		print "p 1chunk=",q
		print "usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted2.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/2ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits"
		p3[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted3.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/3ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		print "p 1chunk=",q
		print "usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted3.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/3ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits"
		p4[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted4.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/4ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		print "p 1chunk=",q
		print "usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted4.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/4ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits"
		p5[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted5.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/5ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits" %(q,q), shell=True)
		print "p 1chunk=",q
		print "usearch -ublast ~/scripts/tmp2/ubin%d.fna -db ~/data/db/genpept/sorted5.udb -evalue 1e-10 -accel 0.5 -maxhits 1 -maxaccepts 100 -maxrejects 1024 -threads 1 -userout ~/scripts/tmp2/5ubin%dout.txt -userfields query+target+id+qlo+qhi+tlo+thi+ql+tl+alnlen+evalue+bits"

	#wait for n threads to finish
	for x1 in range(1,threads+1):
		print "waiting for usearch to finish........."
		p1[x1].wait(); p2[x1].wait();p3[x1].wait(); p4[x1].wait();p5[x1].wait()
	
	for x2 in range(1,threads+1):		
		
		print "concatenating.....", x2
			
		subprocess.Popen("cat ~/scripts/tmp2/1ubin%dout.txt ~/scripts/tmp2/2ubin%dout.txt ~/scripts/tmp2/3ubin%dout.txt ~/scripts/tmp2/4ubin%dout.txt ~/scripts/tmp2/5ubin%dout.txt >> ~/scripts/tmp2/finaloutput.ublast" %(x2,x2,x2,x2,x2), shell=True).wait()
	
subprocess.Popen("rm -rf ~/scripts/tmp2/ubin* ~/scripts/tmp2/1ubin* ~/scripts/tmp2/2ubin* ~/scripts/tmp2/3ubin* ~/scripts/tmp2/4ubin* ~/scripts/tmp2/5ubin* ", shell=True).wait()

print"adding taxonomy.."
subprocess.Popen("python ~/scripts/taxonfetch.py ~/scripts/tmp2/finaloutput.ublast ~/scripts/tmp2/finaloutput.taxa.ublast p", shell=True).wait()
