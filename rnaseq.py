#! /usr/bin/env python

#purpose of this script is to take fastq.gz files and run through the 
#standard RNAseq qc, mapping, sorting, and counting of reads
#it should be called once per fastq.gz file and run from codon because it runs
# a bunch of codon jobs

# the basic idea is to run a few pipelines
#  fastqc -> read quality stats
#  tophat | samtools sort | htseq-count | -> counts for genes

#this script simply calls shell scripts to run jobs so a separate script for 
# each of the above portions of the pipeline is necessary

#it should be moved to the ~/analysis/RNAseq/<EXPT>/scripts folder and
#the pointer to the files should be changed in that copy

##############################################################
import argparse 
import subprocess
import os
###############################################################
 
def findFiles(searchPath): 
   os.chdir(searchPath)
   outFiles=[]
   for f in os.listdir("."):
      if f.endswith(".fastq.gz"):
         fileName=searchPath+f
         outFiles.appendi([fileName])

   return(outFiles)


def qc(fa,resultsDir): 
   FILE='FILE='+fa
   OUT='OUT='+resultsDir+fa
   subprocess.call(['qsub','-v',FILE,OUT,'/magnuson-lab/jraab/scripts/qc.sh']) 


################################################################
parser = argparse.ArgumentParser(description='main program for initial rnaseq processing') 
parser.add_argument('file', help='input file(s)')
parser.add_argument('-q', help='run fastqc on codon')
parser.add_argument('-t', help='run tophat on codon')
parser.add_argument('-s', help='run samtools sort on codon')
parser.add_argument('-c', help='run htseq-count on codon')

args = parser.parse_args()

################################################################
files=args.file
resultsDir=''

if (args.q):
   for f in files:
      fname=
      qc(f, resultsDir+fname)

if(args.t):
   for f in files:
      tophat(f)

if(args.s): 
   for f in files:
      samtoolsSort(f)

if(args.c): 
   for f in files:
      htCount(f)



      


