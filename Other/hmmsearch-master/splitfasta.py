#!/usr/bin/python
# SplitFASTA
# Scripts for the examination of neighborhood genes of bacterial genome
# Require locally installed Hmmer 3.0, CD-HIT
#
import re
import sys

capture = False

f = open(sys.argv[1])
gi = re.compile('>gi\|(\S+)\|ref\|(\S+)\| (.*) \[(.*)\]')
# gid = ''
desc=''
organism = ''
id = ''
refseq = ''
buffer=[]
while True:
    line = f.readline()
    if not line:
        break
    match = gi.match(line)
    if match:
        if capture:
            print refseq + '.fasta'
            fw = open(refseq + '.fasta', 'w')
            fw.writelines(buffer)
            fw.close
            capture = False
            refseq = ''
            buffer = []
        else:
            id = match.group(2).strip()
            capture = True
            refseq = id
        
    if capture:
        buffer.append(line)

if capture:
    
    fw = open(refseq + '.fasta', 'w')
    fw.writelines(buffer)
    fw.close
    capture = False
    refseq = ''
         