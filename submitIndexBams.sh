#!/usr/bin/sh


#$ -o indexed.out
#$ -e index.err

module load samtools

samtools index $FILE
