#!/usr/bin/env/python 
#Given a data frame containing interval information about peaks, and bamfile of some histone mod or transcription factors signal
#calculate how enriched that signal is over each interval. 
###################################
import HTSeq
import pandas as pd
import numpy as np
#import prettyplotlib as ppl
#from matplotlib import pyplot as plt
import argparse 
import sys
import os
import gzip
import pysam
import tables
from matplotlib import rcParams
import brewer2mpl
#######################################
#arguments
parser = argparse.ArgumentParser(description='calculate enrichment at a set of intervals')
parser.add_argument('-p', help='peakfile')
parser.add_argument('-b', help='bamfile', nargs='?')
parser.add_argument('-bamname', help='bamfile name - for naming plots and files')
parser.add_argument('--plot', action='store_true', default=False, help='output plot of enrichment and scores')
parser.add_argument('-w', default=500, help='size of the enrichment window')
parser.add_argument('-bg', default=50000, help='size of the background window')
parser.add_argument('--storesize', default=10000, help='How much coverage info should be saved to hdf5')
# I took out the argument to normalize and will always do this to RPM, before calculating enrichments, should help scaling
#parser.add_argument('-o', nargs='?', help='directory to send files to')
#parser.add_argumnet('-i', nargs='?', help='bamfile of input') # add this at some point
args = parser.parse_args()
#########################################
#
print args
#define any constant things, these may later take variable form. 
out_fn = 'test_hdf5.hdf5'
if args.bamname:
   ab = args.bamname
else: 
   ab = 'test' #or change me
outdir = '/magnuson-lab/jraab/analysis/swi_snf/swi_snf_chipseq/expt2/processed_data/summaries/'
peakname = os.path.basename(args.p).split('-')[1] #this is specific for my naming scheme which is Lane-Antibody-replicate-bunchofstuff
#######################
#better matplotlib defaults
dark2_colors = brewer2mpl.get_map('Set2', 'Qualitative', 7).mpl_colors
rcParams['figure.figsize'] = (10,6)
rcParams['figure.dpi'] = 150
rcParams['axes.color_cycle'] = dark2_colors
rcParams['axes.facecolor'] = 'white'
rcParams['grid.linewidth'] = 0 
rcParams['font.size'] = 14
rcParams['patch.edgecolor'] = 'white'
rcParams['patch.facecolor'] = dark2_colors[0]
rcParams['font.family'] = 'StixGeneral'

def remove_border(axes = None, top=False, right=False, left=True, bottom=True):
   """
   Minimize chartjunk by stripping out unnecessary plot borders and axis ticks
   The top/right/left/bottom keywords toggle whether corresponding plot border is drawn
   """
   ax = axes or plt.gca()
   ax.spines['top'].set_visible(top)
   ax.spines['right'].set_visible(right)
   ax.spines['left'].set_visible(left)
   ax.spines['bottom'].set_visible(bottom)

   #turn off all ticks
   ax.yaxis.set_ticks_position('none')
   ax.xaxis.set_ticks_position('none')
   #turn off grid
   ax.grid=False

   #now re-enable
   if top:
      ax.xaxis.tick_top()
   if bottom:
      ax.xaxis.tick_bottom()
   if left:
      ax.yaxis.tick_left()
   if right:
      ax.yaxis.tick_right()
      
#########################
#function defs
def coverage (pos_dic, bamfile, halfwinwidth, fragmentsize=200):
   """
   Function calculates the coverage over a dict of HTSeq.GenomicIntervals,keyed on peakname, when given a bamreader object
   returns a dict where the peakname is the key
   """
   frame = {}
   bf = pysam.Samfile(bamfile, 'rb')
   for key in pos_dic: 
      p = pos_dic[key]
      profile = np.zeros(2*halfwinwidth, dtype='float')
      window = HTSeq.GenomicInterval(p.chrom, p.start_d-halfwinwidth-fragmentsize, p.start_d+halfwinwidth+fragmentsize, '.')
      if window.start <0: 
         continue
      for a in bf.fetch(p.chrom, window.start, window.end): 
         if a.is_reverse:
            strand = '-'
         else: 
            strand = '+'
         almnt = HTSeq.GenomicInterval(window.chrom, a.pos, a.pos+fragmentsize, strand)  
         if p.strand == '+' or p.strand =='.': 
            start_in_window = almnt.start - p.start_d + halfwinwidth
            end_in_window = almnt.end - p.start_d + halfwinwidth
         else: 
            start_in_window = p.start_d+halfwinwidth-almnt.end
            end_in_window = p.start_d+halfwinwidth-almnt.start
         start_in_window = max(start_in_window, 0)
         end_in_window = min(end_in_window, 2*halfwinwidth)
         if start_in_window >= 2*halfwinwidth or end_in_window < 0:
            continue
         profile[start_in_window:end_in_window] += 1
      frame[key] = profile
   return(frame)

def read_peaks(filename,stranded=False):
   """
   Take a file name that is approximately a bed file (needs at least 4 columns) - and will truncate any other info
   name column will be conserved as a dictionary key throughout this program
   """
   genes = {}
   with open(filename) as fn:  
      for line in fn:
         if line.startswith('#'):
            continue
         try: 
            chrom, start, end, name = line.strip().split('\t')[:4]
         except: 
            sys.exit(fn, 'not a valid bedfile')
         if stranded == True:
            strand = line.strip().split('\t')[5] # for this to work must be valid 6 column bed
         else:
            strand = '.'
         iv = HTSeq.GenomicInterval(chrom, int(start), int(end), strand)
         peaknum = os.path.basename(name).split('_')[-1] #this works for summit.bed files produced by macs, because they make a big long name 
         #based on original file name
         name = peakname+'_summit_'+peaknum
         genes[name] = iv
         j +=1
      return(genes)

def exclude_regions(genedict): 
   """
   take a dictionary of gene-intervals and check them against known blacklisted regions - returning only those
   that do not overlap the restricted areas
   """
   e1 = gzip.open('/magnuson-lab/jraab/annotations/wgEncodeDukeMapabilityRegionsExcludable.bed.gz')
   e2 = gzip.open('/magnuson-lab/jraab/annotations/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
   #e1 = gzip.open('/media/data/annotations/wgEncodeDukeMapabilityRegionsExcludable.bed.gz')
   #e2 = gzip.open('/media/data/annotations/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
   excludedRegions = []
   cleanGenes = {}
   for l1 in e1: 
      chrom, start, end, name, score, strand = l1.split('\t')
      interval = HTSeq.GenomicInterval(chrom, int(start), int(end), '.')
      excludedRegions.append(interval)
   for l2 in e2: 
      chrom, start, end, name, score, strand = l2.split('\t')
      interval = HTSeq.GenomicInterval(chrom, int(start), int(end), '.')
      excludedRegions.append(interval)
   for key in genedict: 
      gene = genedict[key] 
      hit = 0 
      for interval in excludedRegions: 
         if interval.overlaps(gene):
            hit = 1
      if hit ==0: 
         cleanGenes[key] = gene
   return(cleanGenes)

def enrichment_score(numpy_vec, background_window, peak_window): 
   """
   Given a numpy array calculate how much enrichment over the background window is found int he peak_window
   Simplistic scorign of log2(median(small_window)+1/median(background) +1 )
   """
   mid = int(len(numpy_vec)/2)
   subset = numpy_vec[mid-peak_window: mid+peak_window] 
   num = np.median(subset) +1
   den = np.median(numpy_vec) + 1
   enrichment = np.log2(num/den)
   return(enrichment)

def read_depth(bamfile): 
   counter = 0
   bf = pysam.Samfile(bamfile, 'rb')
   for almnt in bf.fetch(): 
      counter +=1
   return(counter)

def rpm(val, reads): 
   return(val*1e6/reads)

def save_npdicts_tohdf5(dicts, fn): 
   """
   given a dictionary of numpy arrays, save the result to hdf5
   """

def make_plots(adf, enrichments,title, label, outdir): 
   """
   make a 1x2 plot of enrichment around a peak, and overall enrichment scores
   stuck using matplotlib here for now, can't seem to update it on cocdon
   There seems to be issues with this on codon, checked out a new branch to run through all enrichment and coverages without plotting
   """
   win = adf.shape[1]/2  
   adf_collapse = adf.mean(axis=0)
   plt.figure()
   plt.suptitle(title, fontsize=20)
   plt.subplot(121)
   plt.plot(np.arange(-win,win), adf_collapse, label=label)
   plt.subplot(122)
   plt.hist(enrichments.values(),label=label)
   plt.legend (fontsize=18)
   remove_border()
   plt.savefig(outdir+'plots/'+label+'_'+ab+'_enrichment.png') 
   

   
def main():
   genes = read_peaks(args.p)
   genes = exclude_regions(genes) 
   total_reads = read_depth(args.b) 
   enrichwin = int(args.w)
   background = int(args.bg)
   all_coverage = coverage(genes, args.b, background/2)
   all_coverage = {k:rpm(v, total_reads) for k,v in all_coverage.iteritems() }
   enrichments = {k:enrichment_score(v,args.bg, args.w) for k,v in all_coverage.iteritems() } 
   subsets = {}
   for k,v in all_coverage.iteritems(): 
      midpt = int(len(v)/2)
      sub_win  = [midpt-args.storesize, midpt+args.storesize]
      subsets[k] = v[sub_win[0]:sub_win[1]] 
   sub_df = pd.DataFrame.from_dict(subsets, orient='index')
   enrichments_df = pd.DataFrame.from_dict(enrichments, orient='index')
   print len(sub_df.keys())
   #need to sort out how to move from a big dictionary of numpy arrays into a structure I can load into R or python 
   #save things I want
   if args.plot: 
      title = ab+' at all ' +label +' peaks'
      make_plots(sub_df, enrichments, title, label, outdir) 
   #save frame
   label=peakname
   e_fn = outdir+'enrichments/'+label+'_'+ab+'_enrichments.csv'
   s_fn = outdir+'coverages/'+label+'_'+ab+'_coverage.h5'
   #sub_df.to_csv(s_fn, sep=',')
   enrichments_df.to_csv(e_fn, sep=',')
   store = pd.HDFStore(s_fn)
   store['df'] = sub_df
   store.close()
   
if __name__ == '__main__': 
   main()
