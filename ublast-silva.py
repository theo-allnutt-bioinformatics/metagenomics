from __future__ import division
from Bio import SeqIO
import sys
import os
from StringIO import StringIO # Python 2
import subprocess
from random import randrange
from optparse import OptionParser
import re
import glob

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

parser = OptionParser()
parser.add_option("-i", "--infile", dest="filename", help="input sequence file", metavar="FILE")
parser.add_option("-o", "--outfile", dest="outfile", default="ublast.out", help="output file name")
parser.add_option("-T", "--tax_lev", dest="tax_lev", default="a", help="report taxonomy level - default is all ranks")
parser.add_option("-f", "--format", dest="infmt", default="fasta", help="input file format, e.g. fasta, fastq, genbank")
parser.add_option("-m", "--maxhits", dest="maxhits", default="100", help="max number of hits to return for each query")
parser.add_option("-d", "--identity", dest="id", default="0.9", help="percent identity threshold")
parser.add_option("-b", "--database", dest="db", default="~/d/db/genpept/genpept.udb", help="database, default = ~/d/db/genpept/genpept.udb")
parser.add_option("-e", "--evalue", dest="evalue", default="1.0e-10", help="e value threshold")
parser.add_option("-c", "--cons", dest="cons", default="0.15", help="specify LCA scoring radius 0 to 1, default = 0.15")
parser.add_option("-s", "--stype", dest="stype", default="p", help="specify the sequence type. Default = protein, 'p'")
parser.add_option("-O", "--offset", dest="offset1", default="1", help="specify the split offset, default = 1, '1'")
parser.add_option("-a", "--accel", dest="accel", default="0.5", help="ublast accel parameter, default = 0.5, '0.5'")

(options, args) = parser.parse_args()

r2=str(randrange(100000))
db=options.db
folder = options.filename #sys.argv[1]
#outputfile=options.outfile #sys.argv[2]

acc = options.accel
fmt = options.infmt
maxhits= options.maxhits
id= options.id
evalue=options.evalue
cons = float(options.cons)
stype=options.stype
tax_lev=options.tax_lev
offset1=options.offset1

#017

'''
subprocess.Popen("rm -rf ~/scripts/tmp2", shell=True).wait()
subprocess.Popen("mkdir ~/scripts/tmp2", shell=True).wait()
'''		  
#folder = sys.argv[1]
				  			  
filelist=glob.glob(folder+"/*")
filelist.sort(key=tokenize)


print filelist	

for inputfile in filelist:
	output_name= inputfile.split("/")[-1].split(".")[0]
	inputname='/OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa' %(r2)
	f1=open(inputfile,'r')
	f2=open(inputname,'w')
	outputfile = output_name+".ublast"
	#split large reads/contigs if >50kb
	print "splitting large seq records if >50kb..."
	for i in SeqIO.parse(f1,fmt):
		
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

	p1= subprocess.Popen("usearch -ublast %s -db %s -strand both -maxhits %s -maxaccepts %s -userout %s -userfields  query+target+id+ql+qlo+qhi+tl+tlo+thi+alnlen+evalue+bits -id %s -evalue %s -accel %s" %(inputname,db,maxhits,maxhits,outputfile,id,evalue,acc), shell=True).wait() #-ka_ungapped_lambda 0.318 -ka_gapped_lambda 0.291 -ka_gapped_k 0.075 -ka_ungapped_k 0.134

	outputfile2=output_name+".taxa"

	print"adding lineage.."
	subprocess.Popen("python ~/scripts/linfetch-silva.py %s %s %s %s" %(outputfile,outputfile2, stype, offset1), shell=True).wait()


	subprocess.Popen("python ~/scripts/lca2.py %s %s %s %s" %(outputfile2,output_name+".lca",cons,tax_lev), shell=True).wait()

	subprocess.Popen("rm /OSM/HOME-MEL/all29c/scripts/tmp2/inp%s.fa" %(r2), shell=True).wait()


