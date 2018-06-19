from Bio import SeqIO
import sys
import os
from StringIO import StringIO # Python 2
import subprocess
from random import randrange



if os.path.isfile('/db2/genpept.udb')== False:
	print "loading db to ram"

	subprocess.Popen("sudo mkdir -p /db2 && sudo mount -t tmpfs tmpfs /db2 && sudo chmod 777 /db2 && cp -n ~/data/db/genpept/genpept.udb /db2/genpept.udb", shell=True).wait()
else:
	print "db already exists"
	
def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


inputfile = sys.argv[1]
outputfile=sys.argv[2]
chunksize = int(str(sys.argv[3])) #number of reads in each division of the input file
threads = int(sys.argv[4])
fmt = sys.argv[5] #input format
fetchsize = sys.argv[6]
maxhits = sys.argv[7]
meta=sys.argv[8]

r2=str(randrange(100000))

subprocess.Popen("rm -rf ~/scripts/tmp2", shell=True).wait()
subprocess.Popen("mkdir ~/scripts/tmp2", shell=True).wait()

out1=open('/OSM/HOME-MEL/all29c/scripts/tmp2/finaloutput%s.ublast' %(r2),'w')

f2=open('/OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa' %(r2),'w')
#split large reads/contigs if >50kb
print "splitting large seq records if >50kb..."
for i in SeqIO.parse(inputfile,fmt):
	
	if len(i.seq)>50000:
		seqlist=list(chunkstring(i.seq,50000))
		c=-1
		for t in seqlist:
			c=c+1
			print i.id, c+1
			f2.write(">"+i.id+"-part-"+str(c+1)+"\n"+str(t)+"\n")
	else:
		f2.write(">"+i.id+"\n"+str(i.seq)+"\n")

	
f2.close()
inputfile2=open('/OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa' %(r2),'r')
count = SeqIO.index('/OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa' %(r2), fmt)

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
outfile = "/OSM/HOME-MEL/all29c/scripts/tmp2/%subin" %(r2)+str(f)+".fna"
g=open(outfile,'w')
for i in SeqIO.parse(inputfile2,fmt):
	
	t=t+1
	if t>chunksize:
		f=f+1
		outfile = "/OSM/HOME-MEL/all29c/scripts/tmp2/%subin" %(r2)+str(f)+".fna"
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
	
	print len(p1),"lenp1"
	
for rep in range(1,numreps+1):
	
	for x in range(1,threads+1):
		
		q=q+1 #file number
		
		p1[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/%subin%d.fna -db /db2/genpept.udb -evalue 1e-10 -accel 0.5 -maxhits %s -maxaccepts 100 -maxrejects 1024 -userout ~/scripts/tmp2/%subin%dout.txt -userfields query+target+id+ql+qlo+qhi+tlo+thi+tl+alnlen+evalue+bits" %(r2,q,maxhits,r2,q), shell=True)



	#wait for n threads to finish
	for x1 in range(1,threads+1):
		print "waiting for usearch to finish........."
		p1[x1].wait()
	
	for x2 in range(1,threads+1):		
		
		print "concatenating.....", x2
			
		if meta =="m":
			subprocess.Popen("cat ~/scripts/tmp2/%subin%dout.txt >> ~/scripts/tmp2/finaloutput%s.ublast" %(r2,x2,r2), shell=True).wait()
		else:
			subprocess.Popen("cat ~/scripts/tmp2/%subin%dout.txt >> %s" %(r2,x2,outputfile), shell=True).wait()
	
	
#subprocess.Popen("rm -rf ~/scripts/tmp2/ubin* ~/scripts/tmp2/1ubin* ~/scripts/tmp2/2ubin* ~/scripts/tmp2/3ubin* ~/scripts/tmp2/4ubin* ~/scripts/tmp2/5ubin* ", shell=True).wait()

if meta =="m":

	print"adding taxonomy.."
	subprocess.Popen("python ~/scripts/taxonfetch.py ~/scripts/tmp2/finaloutput%s.ublast %s p %s" %(r2,outputfile,fetchsize), shell=True).wait()

	#############otufreq.py


	f2 = open(outputfile,'r')
	g2=open(outputfile+".species",'w')
	g3=open(outputfile+".genera",'w')
	# read output file species into list
	species=[]
	genus=[]
	k=[]
	print "Getting species.."
	for line in f2:
		
		k=line.rstrip("\n").split("\t")[13].split(" ")
		print k
		if len(k)>1:
			
			species.append(" ".join(str(x5) for x5 in k[0:2])) #nb makes two element list into a single string
		if len(k)==1:
			species.append(k[0])
			
		genus.append(line.split("\t")[13].split(" ")[0])

	species_set=set(species)
	genus_set=set(genus)
	print "Counting species.."
	for i in species_set:
		g2.write(i+"\t"+str(species.count(i))+"\n")
	for i in genus_set:
		g3.write(i+"\t"+str(genus.count(i))+"\n")

subprocess.Popen("rm /OSM/HOME-MEL/all29c/scripts/tmp2/finaloutput%s.ublast" %(r2), shell=True).wait()
subprocess.Popen("rm /OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa" %(r2), shell=True).wait()
subprocess.Popen("rm -r /OSM/HOME-MEL/all29c/scripts/tmp2/%subin*" %(r2), shell=True).wait()


f2.close()
g2.close()
g3.close()