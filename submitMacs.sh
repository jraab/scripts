#!/bin/bash
#$ -cwd
#$ -o macs.OUT
#$ -e macs.ERR

module load macs

macs -t ~/seqdata/ENCODE/datafiles/HepG2/combined/wgEncodeBroadHistoneHepg2Ezh239875Aln.sorted.merged.bam -c  ~/seqdata/ENCODE/datafiles/HepG2/combined/wgEncodeBroadHistoneHepg2ControlStdAln.sorted.merged.bam -f BAM -g hs -n ezh2.hepg2.combined. 



