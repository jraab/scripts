#!/usr/bin/python
################
#general script to take a bam file and make a bigwig file
# do not store intermediates - since these take pu stupid amounts of space

#program requires the ucsc utility bedgraphToBigWig to be somewhere in the $PATH
import argparse
import HTSeq
import numpy
import pybedtools

parser = argparse.ArgumentParser(description = 'Convert BAM files to wiggle tracks' ) 
parser.add_argument('-b', help='input bam' )
parser.add_argument('--norm', help = 'Normalize by read depth' , action='store_true')
parser.add_argument('-s', help = 'separate files for + and - strand', action='store_true', default=False)
parser.add_argument('-g', help = 'genome name e.g hg19')
args = parser.parse_args()

################################
def calc_norm_factor(bamfile):	
	bam_reader = HTSeq.BAM_Reader(args.b)
	for read in bam_reader: 
		if read.aligned: 
			read_counter +=1
	return( 1e6 /read_counter) 



###############################
#this file must exist
chrom_sizes = '/magnuson-lab/jraab/annotations/genome_hg19.txt'
#chrom_sizes = '/magnuson-lab/shared/jraab/annotations/genome_mm9.txt'
chrom_lens = {} 

name = args.b.split('/')[-1]
name = name.split('.')[0]

if args.norm:
	normfactor = calc_norm_factor(args.b)

bt = pybedtools.BedTool(args.b)
if args.s:	
        output_plus = name+'_plus.bw'
        output_minus = name+'_minux.bw'
        if args.norm:
		bedgraph_plus = bt.genome_coverage(genome=args.g, bg=True, strand="+", scale=normfactor)
		bedgraph_minus =bt.genome_coverage(genome=args.g, bg=True, strand='-', scale=normfactor)
	else:
		bedgraph_plus = bt.genome_coverage(genome=args.g, bg=True, strand="+")
		bedgraph_minus = bt.genom_coverage(genome=args.g, bg=True, strand="-")
        cmds = ['bedGraphToBigWig', bedgraph_plus.fn, args.g, output_plus] 
        os.system(' '.join(cmds))
        cmds = ['bedGraphToBigWig', bedgraph_minus.fn, args.g, output_minus] 
        os.system(' '.join(cmds))
else:
	output = name+'.bw'
	if args.norm: 
		bedgraph = bt.genome_coverage(genome=args.g, bg=True, strand=True, scale=normfactor)
	else: 
		bedgraph = bt.genome_coverage(genome=args.g, bg=True, strand=True)
	cmds = ['bedGraphToBigWig', bedgraph.fn, args.g, output]
        os.system(' '.join(cmds)







 
#with open(chrom_sizes) as f: 
#   for l in f:
#      chr,length = l.split('\t')
#      chrom_lens[chr] = int(length)
#
#if args.s == True: 
#   ga = HTSeq.GenomicArray(chrom_lens, stranded = True, typecode='d')
#else :
#   ga = HTSeq.GenomicArray(chrom_lens, stranded = False, typecode = 'd')
#
#
#bam_reader = HTSeq.BAM_Reader(args.b)
#
#read_counter = 0 
#for read in bam_reader: 
#   if read.aligned :
#      read_counter += 1 
#      ga[read.iv] += 1
#
#
#
#if args.norm:
#   for iv,value in ga.steps(): 
#      ga[iv].apply(lambda x: x*1e6/read_counter) 
#   if args.s == True:
#   	ga.write_bedgraph_file(name+'_plus.bg', strand='+')
#        ga.write_bedgraph_file(name+'_minus.bg', strand='-')
#   else: 
#        ga.write_bedgraph_file(name+'.bg', strand='.')

#else: 
#   for iv,value in ga.step():
#      ga[iv].apply(lambda x: x*1e6/read_counter)
#      if args.s==True:
#         ga.write_bedgraph_file(name+'_plus.bg', strand='+')
#         ga.write_bedgraph_file(name+'_minus.bg', strand='-')
#      else:
#         ga.write_bedgraph_file(name+'.bg', strand='.') 



