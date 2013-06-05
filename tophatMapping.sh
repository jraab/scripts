#!/bin/bash

# the name of the job
#$ -N jr_rnaMAP

# the name of the file that gets STDOUT

#$ -o STDOUT.txt
# the name of the file that gets STDERR
#$ -e STDERR.txt
#$ -m bae
#$ -M jesse.r.raab@gmail.com

# use the current working directory.
#$ -cwd

## This file should not be changed
# it called from a python script that sets up the mapping parameters
# the options to tophat can be changed if necessary, but these basic ones work

 
module load base
module load samtools/0.1.18 
module load bowtie2	
module load tophat

BASE=
FILE=
PE=
SE=

if [$PE eq TRUE]; 
	then 
	tophat -p 8 --library-type fr-firststrand --no-coverage-search -r 150 -o ${BASE}/MAPPING1/${FILE} /magnuson-lab/jraab/indexes/hg19 ${FILE}.fastq.gz  

fi 



if [$SE eq TRUE]; 
	then 
	tophat -p 8 --library-type fr-firststrand --no-coverage-search -r 150 -o ${BASE}/MAPPING1/${FILE} /magnuson-lab/jraab/indexes/hg19 ${FILE}.R1.fastq.gz ${FILE}.R2.fastq.gz

fi 









