#!/usr/bin/python 
import HTSeq
import itertools
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import argparse
import os
import sys
#########################################################################################################################
def read_feature(gtffile):
    tsspos = set ()
    for feature in ( gtffile ):
        line = feature.split()
        i = 0
        if len(line) < 6:
            line.extend('.')
	if len(line) < 5:
            line.append('.')
        if len(line) < 4:
            line.append ( '.')
        if len(line) < 3:
            sys.exit( "this file isn't even close to a valid bed\n")
        chr,start,end,name,score,strand = line    
	#if chr == "#":
        #    continue
	if strand is '+' or strand is '-':
	    strand = strand
   else :
      strand = "."
   start = int(start)
   end = int(end) 
	midpt = start+end/2
   gi  = HTSeq.GenomicInterval(chr,midpt,midpt+1,strand)
   tsspos.add( gi.start_d_as_pos )
    return(tsspos)	
def coverage_plot(tsspos,bamfile,halfwinwidth,fragmentsize):
    profile = numpy.zeros (2*halfwinwidth,dtype = "float")
    nump = 1
    for p in tsspos: 
        window = HTSeq.GenomicInterval ( p.chrom, p.pos - halfwinwidth - fragmentsize, p.pos + halfwinwidth + fragmentsize, "." )	 
        if window.start < 0:
            continue
        for almnt in bamfile [window]:
            if p.strand == "+":
                start_in_window = almnt.iv.start - p.pos + halfwinwidth
                end_in_window   = almnt.iv.end - p.pos + halfwinwidth 
            else :
                start_in_window = p.pos + halfwinwidth - almnt.iv.end
                end_in_window = p.pos + halfwinwidth - almnt.iv.start
            start_in_window = max (start_in_window, 0 )
	    end_in_window = min ( end_in_window, 2*halfwinwidth)
            if start_in_window >= 2*halfwinwidth or end_in_window < 0:
                continue
            profile [ start_in_window : end_in_window ] += 1
        nump += 1 
    return(profile,nump)

def aligned_reads (bamfile):
    count = 0
    for almnt in bamfile:
        if almnt.aligned:
            count+=1    
    return(count)
###########################################################################################################################
parser = argparse.ArgumentParser(description='Make average coverage plots over a region of interest')
parser.add_argument('-t', help='bed file of regions of interest')
parser.add_argument('-b', help='bam files of mapped reads', nargs='+')
parser.add_argument('-n', help='name for output files')
args = parser.parse_args()
roiFile  =  open(args.t)
name = args.n
halfwinwidth = 5000
fragmentsize = 200
features = read_feature(roiFile)
num_features = len(features)
bamfile = args.b 
for file in bamfile:
    bam = file.split('.')[0].split('/')[-1]
    print 'processing file - '+bam
    bf_reader_1 = HTSeq.BAM_Reader(file)
    num_reads_1 = aligned_reads(bf_reader_1)
    # num_reads_1=30000000
    print 'file:'+bam[0]+':'+str(num_reads_1) 

    bamfile1_coverage,nump = coverage_plot(features,bf_reader_1,halfwinwidth,fragmentsize)     
    b_coverage1 = (bamfile1_coverage/num_features)/num_reads_1 * 1000000

    pyplot.plot(numpy.arange(-halfwinwidth, halfwinwidth), b_coverage1 )
    ax = [-halfwinwidth, halfwinwidth, 0, 0.5 ] 
    pyplot.axis(ax)

    figName = name+"_"+bam+"coveragePlot.png"   
    pyplot.savefig(figName)
    pyplot.clf()


