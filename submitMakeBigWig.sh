#!/usr/bin/sh
#$ -cwd


#sort -k1,1 -k2,2n $FILE > $FILE.sorted
module load python

/magnuson-lab/jraab/scripts/bedGraphToBigWig $FILE /magnuson-lab/jraab/annotations/genome_hg19.txt $FILE.bw
