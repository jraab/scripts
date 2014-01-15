#!/usr/bin/python 
#script takes a target sequence and identifies guideRNAs
#usable with the Cas9 genome editing system
#gRNAs are in the form GN20GG
import string
import argparse 
import re

parser = argparse.ArgumentParser(description = 'search a sequence or file for gRNA')
parser.add_argument('file', help='input file(s)')
#parser.add_argument('-c', help='sequence on comand line')
args = parser.parse_args()
###############################################
def rc(dna):
    complements = string.maketrans('acgtACGT',
'tgcaTGCA')
    rcseq = dna.translate(complements)[::-1]
    return rcseq
##############################################
fin = open(args.file, 'r')
data= []
for line in fin:
    if line.startswith('>'):
        continue
    data.append(line)
data = ''.join(data)
sequence = data.replace('\n','')
reversecomp = rc(sequence) 
pattern = 'G[A,G,C,T]{20}GG' 
fMatch = re.findall(pattern,sequence.upper())
rMatch = re.findall(pattern,reversecomp.upper())

for i in fMatch:
    print i
for i in rMatch: 
    print i


