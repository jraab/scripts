#! /usr/bin/python

import sys
import os
import pysam
import glob

#read a list of file names from a directory and combine all replicates of a bam

def readList(f): 
#may need to change file to an ls command
   allfiles = []
   for l in f: 
      allfiles.append(l.split(".")[0])   
   roots = set()
   for name in allfiles:  
      roots.add(name.split("Rep")[0])
   return(roots) 
names = readList(open(sys.argv[1]))


def mergeFiles(name): 
   numeps = 2 
   f = glob.glob(name+'*.bam')
   print f
   samfile = pysam.Samfile(name+'Rep1.bam')
   outfile = pysam.Samfile('combined/'+name+'.merged.bam','wb', template=samfile)
   for i in f: 
       samfile = pysam.Samfile(i)  
       for read in samfile.fetch(): 
             outfile.write(read)
       samfile.close()
   outfile.close()

def dothework(filelist):
   for n in names: 
      mergeFiles(n)
      pysam.sort('combined/'+n+'.merged.bam', 'combined/'+n+'.sorted.merged')
      pysam.index('combined/'+n+'.sorted.merged.bam')

names=readList(open(sys.argv[1]))
dothework(names)
