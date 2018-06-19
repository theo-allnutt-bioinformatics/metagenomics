from Bio import SeqIO
import sys
import os
from StringIO import StringIO # Python 2
import subprocess
from random import randrange
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-i", "--infile", dest="filename", help="input sequence file", metavar="FILE")
parser.add_option("-o", "--outfile", dest="outfile", default="ublast.out", help="output file name")
parser.add_option("-t", "--threads", dest="threads", default="8", help="output file name")
parser.add_option("-f", "--format", dest="infmt", default="fasta", help="input file format, e.g. fasta, fastq, genbank")
parser.add_option("-m", "--maxhits", dest="maxhits", default="10", help="max number of hits to return for each query")
parser.add_option("-d", "--identity", dest="id", default="0.9", help="percent identity threshold")
parser.add_option("-b", "--database", dest="db", default="/db/gbnt.udb", help="database, default = /db/gbnt.udb")
parser.add_option("-e", "--evalue", dest="evalue", default="1.0e-10", help="e value threshold")
parser.add_option("-c", "--cons", dest="cons", default="0.15", help="specify consensus taxon scoring threshold 0 to 1, default = 0.75")
parser.add_option("-s", "--stype", dest="stype", default="p", help="specify the sequence type. Default = protein, 'p'")

(options, args) = parser.parse_args()

r2=str(randrange(100000))
db=options.db
dbname = db.split("/")[-1]


filesize=int(os.path.getsize(db))

if filesize < 10000000000:

	if os.path.isfile('/db/%s' %dbname)== False:
		print "loading db to ram"
		subprocess.Popen("sudo mkdir -p /db && sudo mount -t tmpfs tmpfs /db && sudo chmod 777 /db && cp -n %s /db/%s" %(db,dbname), shell=True).wait()
		db = "/db/%s" %dbname
	else:
		db = "/db/%s" %dbname
		print "db already exists"
else:
	print "db too large for memory.. usearch to use directly, will be a little slower", "size = ", str(filesize/1000000000),"Gb"

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

	#count reads

inputfile = options.filename #sys.argv[1]
outputfile=options.outfile #sys.argv[2]
threads = int(options.threads)
fmt = options.infmt
maxhits= options.maxhits
id= options.id
evalue=options.evalue
cons = float(options.cons)
stype=options.stype

maxacc = 5
'''
subprocess.Popen("rm -rf ~/scripts/tmp2", shell=True).wait()
subprocess.Popen("mkdir ~/scripts/tmp2", shell=True).wait()
'''

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


c1= len(count)

if threads >c1:
	chunksize = c1
	threads = c1
else:
	chunksize = int(c1/threads)+1
seqs=[]

print "threads=",threads
print "input file=", c1
print "Chunksize =", chunksize

for i in SeqIO.parse(inputfile2,fmt):
	seqs.append(i)

print len(seqs)

if len(seqs)>chunksize:
	seqlist=list(chunkstring(seqs,chunksize))
	c=-1
	for t in seqlist:
		c=c+1
		
		SeqIO.write(t,'/OSM/HOME-MEL/all29c/scripts/tmp2/%subin%d.fna' %(r2,c),fmt)
else:
	SeqIO.write(seqs,'/OSM/HOME-MEL/all29c/scripts/tmp2/%subin0.fna' %(r2),fmt)
				

print "Files written"
print "Multithreading ublast"

#numt = number of times to run multithread 
p1=[]

q=-1

for i in range(0,threads):
	p1.append("")
	
	
for x in range(0,threads):
		
	q=q+1 #file number
	
	print""
	p1[x] = subprocess.Popen("usearch -ublast ~/scripts/tmp2/%subin%d.fna -db %s -lopen 3 -lext 1 -strand both  -maxaccepts %s -maxrejects 255 -maxhits %s -userout ~/scripts/tmp2/%subin%dout.txt -userfields  query+target+id+ql+qlo+qhi+tl+tlo+thi+alnlen+evalue+bits -id %s -evalue %s" %(r2,q,db,maxacc,maxhits,r2,q,id,evalue), shell=True)

	
	#wait for n threads to finish
for x1 in range(threads):
	print "waiting for usearch to finish........."
	p1[x1].wait()
	
for x2 in range(0,threads):		
		
	print "concatenating.....", x2
			
	subprocess.Popen("cat ~/scripts/tmp2/%subin%dout.txt  >> ~/scripts/tmp2/finaloutput%s.ublast" %(r2,x2,r2), shell=True).wait()
	
#subprocess.Popen("rm -rf ~/scripts/tmp2/ubin* ~/scripts/tmp2/1ubin* ~/scripts/tmp2/2ubin* ~/scripts/tmp2/3ubin* ~/scripts/tmp2/4ubin* ~/scripts/tmp2/5ubin* ", shell=True).wait()




print"adding lineage.."
subprocess.Popen("python ~/scripts/linfetch.py ~/scripts/tmp2/finaloutput%s.ublast %s %s" %(r2,outputfile,stype), shell=True).wait()

print "finding lca.."

subprocess.Popen("python ~/scripts/lca.py %s %s %s " %(outputfile,outputfile+".taxa",cons), shell=True).wait()

subprocess.Popen("rm /OSM/HOME-MEL/all29c/scripts/tmp2/finaloutput%s.ublast" %(r2), shell=True).wait()
subprocess.Popen("rm /OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa" %(r2), shell=True).wait()
subprocess.Popen("rm -r /OSM/HOME-MEL/all29c/scripts/tmp2/%subin*" %(r2), shell=True).wait()

