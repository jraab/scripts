#!/usr/bin/env python
#########################
#Given a bed file of interval info and a bamfile of reads 
#calculate the enrichment over each interval - if named, intervals will 
#retain name, else get a random name
#####################################
import HTSeq
import pandas as pd
import numpy as np
import argparse 
import sys
import os
import gzip 
import pysam
import tables 
##############################
#arguments
parser= argparse.ArgumentParser(description='calculate enrichment at a set of intervals') 
parser.add_argument('-p', help='peakfile')
parser.add_argument('-b', help='bamfile')
parser.add_argument('-i', nargs='?', help='input file to normalize enrichment to')
parser.add_argument('-bamname', help='bamfile name - for naming output')
parser.add_argument('-pkwin', default=1000, help='size of enrichment window')
parser.add_argument('-o', nargs='?', default='.',  help='output directory')
parser.add_argument('-w', default = 5000, help='bp to output in coverage file')
parser.add_argument('-bg', default = 50000, help='size of the background window to calculate enrichment')

args = parser.parse_args()
#############################################
#function defs
def readGenes (fname): 
   """
   Take a file name that is bed-like (needs at least 3 columns) - will conserve strand and name if they exist
   Name will be the dictionary key for the output coverages
   """
   genes = {}
   i = 1
   with open(fname) as fn: 
      for line in fn: 
         if line.startswith('#'): 
            continue
         try: 
            chrom,start,end,name,score,strand = line.strip().split('\t')[:6] 
         except: 
            try: 
               a = line.strip().split('\t')
               length = len(a) 
               ext = 6-length
               a.extend(['.'] * ext)
               chrom, start, end, name, score, strand = a
            except: 
               sys.error(fname, 'Not a valid bed')
         if name == '.': 
            name = 'peak'+str(i)
         iv = HTSeq.GenomicInterval(chrom, int(start), int(end), strand)  
         genes[str(name)] = iv
         i += 1
   return(genes) 

def bamCoverage (pos_dic, bamfile, halfwinwidth, fragmentsize = 200 ):  
   frame = {}
   bf = pysam.Samfile(bamfile, 'rb')
   for key in pos_dic:
      p = pos_dic[key] 
      profile = np.zeros(2*halfwinwidth, dtype='float')
      window = HTSeq.GenomicInterval(p.chrom, p.start_d-halfwinwidth-fragmentsize, p.start_d+halfwinwidth+fragmentsize, '.') 
      if window.start < 0: 
         continue
      for a in bf.fetch(p.chrom, window.start, window.end): 
         if a.is_reverse: 
            strand = '-'
         else: 
            strand = '+'
         almnt = HTSeq.GenomicInterval(window.chrom, a.pos, a.pos+fragmentsize, strand) 
         if p.strand == '+' or p.strand == '.': 
            start_in_window = almnt.start - p.start_d + halfwinwidth
            end_in_window = almnt.end - p.start_d + halfwinwidth
         else: 
            start_in_window = p.start_d + halfwinwidth - almnt.end
            end_in_window = p.start_d + halfwinwidth - almnt.start
         start_in_window = max(start_in_window, 0)
         end_in_window = min(end_in_window, 2*halfwinwidth)
         if start_in_window >= 2*halfwinwidth or end_in_window < 0 : 
            continue
         profile[start_in_window:end_in_window] += 1
      frame[key] = profile
   return(frame)
    
   

def excludeRegions(genedict): 
   """
   take a dictionary of gene-intervals and check them against known blacklisted regions - returning only those 
   that do not overlap the restricted areas
   """
   e1 = gzip.open('/magnuson-lab/jraab/annotations/wgEncodeDukeMapabilityRegionsExcludable.bed.gz')
   e2 = gzip.open('/magnuson-lab/jraab/annotations/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
     
   excludedRegions = []
   cleanGenes ={} 
   for l1 in e1: 
      chr,start,end,name,score,strand=l1.split('\t')
      interval = HTSeq.GenomicInterval(chr, int(start), int(end), '.')
      excludedRegions.append(interval)
   for l2 in e2: 
      chr,start,end,name,score,strand=l2.split()
      interval = HTSeq.GenomicInterval(chr, int(start), int(end), '.')
      excludedRegions.append(interval)
   for key in genedict: 
      gene = genedict[key]
      hit = 0
      for interval in excludedRegions: 
         if interval.overlaps(gene): 
            hit=1
      if hit==0: 
         cleanGenes[key] = gene
   return(cleanGenes)      


def readDepth(bamfile): 
   counter = 0 
   bf = pysam.Samfile(bamfile, 'rb') 
   for almnt in bf.fetch(): 
      counter += 1
   return(counter) 

def rpm(val, reads): 
   return( val * 1e6/reads)

def enrichmentScore(numpy_vec, background_window, peak_window):
   """
   Given a numpy array calculate how much enrichment over the background window is
   found in the peak_window
   Simplistic socring of log2(median(small_window)+1/median(background)+1)
   """
   mid = int(len(numpy_vec)/2)
   subset = numpy_vec[mid-peak_window: mid+peak_window] 
   num = np.median(subset) + 1
   den = np.median(numpy_vec) + 1
   enrichment = np.log2(num/den)
   return(enrichment)

def saveHDF(from_dict, fname, group): 
   length_vals = len(from_dict.values()[0])
   class Coverage(tables.IsDescription): 
      name = tables.StringCol(24)
      vals = tables.Float64Col(shape=(1,length_vals) )
   f = tables.open_file(fname, 'w')
   t = f.create_table(f.root, 'coverages', Coverage)
   cov = t.row
   for k,v in from_dict.iteritems(): 
      cov['name'] = k
      cov['vals'] = v
      cov.append()
   t.flush()
   f.close()

##############################################################
def main(): 
   genes = readGenes(args.p)
   genes = excludeRegions(genes)
   peak_depth  = readDepth(args.b)
   peak_coverage = bamCoverage(genes, args.b, int(args.bg)/2)
   peak_coverage = {k:rpm(v, peak_depth) for k,v in peak_coverage.iteritems() } 



   #get info for the input file if included
   if args.i: 
      input_depth = readDepth(args.i)
      input_coverage = bamCoverage(genes, args.i, int(args.bg)/2)
      input_coverage = {k:rpm(v, input_depth) for k,v in input_coverage.iteritems() } 
   

   #normalize to input (as IP+1/Input+1)   
   output_coverage = {}
   for k,v in peak_coverage.iteritems(): 
      try: 
         input = input_coverage[k]
         out_val = (v+1)/(input+1)
      except: 
         out_val = v
      output_coverage[k] = out_val


   #calculate enrichment score over TSS 
   enrichments = {k:enrichmentScore(v, args.bg, args.pkwin) for k,v in output_coverage.iteritems() } 

   #convert to dataframe for saving
   #df = pd.DataFrame.from_dict(output_coverage, orient='index')
   enrich_df = pd.DataFrame.from_dict(enrichments, orient='index')
   cov_dir = args.o+'/coverages/'
   enrich_dir = args.o+'/enrichments/'

   for a in cov_dir,enrich_dir: 
      if not os.path.exists(a): 
         os.mkdir(a)

   out_fn = cov_dir+ args.bamname+'_coverage.h5'
   enrich_fn = enrich_dir+args.bamname+'_enrichment.csv'
   saveHDF(output_coverage, out_fn, 'group')
   enrich_df.to_csv(enrich_fn, sep=',', header=False)

if __name__ == '__main__': 
   main()



   
   
  

    
     















