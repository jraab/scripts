#!/usr/bin/env python
import os
import sys 
############################
#Usage: python getCoords.py geneNames annotation > output.bed

##########################
def lookup (genesfile,annotations):
  genenames = set()
  output={}
  with open(genesfile) as f:
    for line in f: 
      line=line.strip('\n')
      gene_id  = line.split('\t')[0]
      gene_id  = gene_id.replace('\"', '')
      genenames.add(gene_id)
    for gene in genenames: 
      try:  
        fields = annotations[gene] 
        output[gene]=fields
      except:
        continue
    return(output)

def annotation_reader(annotationfile): 
#annotations should be in bed file
  anno = {}
  with open(annotationfile) as a: 
    for anno_line in a: 
      anno_line = anno_line.strip('\n')
      anno_line = anno_line.replace('\"', '')
      try:  
        chr, start, end, name, score, strand = anno_line.split('\t')
        anno[name] = [chr, start, end, name,score, strand]
      except: 
        print 'Annotations are not a valid bed file'
  return(anno)
genesfile = sys.argv[1]
annotationfile = sys.argv[2]
annotations = annotation_reader(annotationfile)
output = lookup(genesfile, annotations) 
print 'chr\tstart\tend\tname\tscore\tstrand'
for line in output.values(): 
  print '\t'.join(map(str,line))
  

