#!/bin/sh
#$ -cwd 
#$ -e qc.err
#$ -o qc.out
module load fastqc


OUTDIR="/magnuson-lab/jraab/labwork/analysis/"#alter this to where you want file to land
fastqc -o ${OUTDIR} -t 2 -q ${FILE}
