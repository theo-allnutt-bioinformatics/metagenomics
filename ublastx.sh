#!/bin/bash

#len1=$(wc -l < ~/data/014-MiSeq/14.7/14.5d1.100.fasta)
#echo
#




#split the input file $1 into $2 chunks send output to ub.temp
#$3 = number of threads
#get number of chunks files: $1==in fasta

len1=$(grep -c ">" $1)
echo "Fasta file size="$len1
echo "Chunk size="$2
echo "number of threads="$3

threads=$3
rm -rf ~/scripts/tmp
mkdir ~/scripts/tmp


numchunks=`python ~/scripts/ublastsplit.py $1 $2` #later use start and end not $1 and $2
echo
echo "bash after files"


#run loop of ublast on splt python files threads times