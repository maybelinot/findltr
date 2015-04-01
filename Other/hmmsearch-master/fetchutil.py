#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import re
import urllib2


def readFasta(filename):

    #
    # Simple FASTA reader. return sequence name with sequence as a string
    #

    sequence = ''
    sequenceName = ''
    if os.path.exists(filename):
        f = open(filename)
        while True:
            line = f.readline().strip()
            if not line:
                break
            match = re.match('^>(.*)', line)
            if match:
                sequenceName = match.group(1)
            else:
                sequence = sequence + line
        f.close()
        return (sequenceName, sequence)
    else:
        return None


def readAccession(filename):

    #
    # From FASTA header, guess accession number and source
    #

    refseqRegex = re.compile('>gi\|(\S+)\|ref\|(\S+)\|')
    uniprotRegex = re.compile('>(sp|tr)\|(\S+)\|')
    db = ''
    accession = ''

    if os.path.exists(filename):
        f = open(filename)
        line = f.readline().strip()
        if not line:
            return None
        refSeqMatch = refseqRegex.match(line)
        if refSeqMatch:
            db = 'refseq'
            accession = refSeqMatch.group(2)
        else:
            uniprotMatch = uniprotRegex.match(line)
            if uniprotMatch:
                db = 'uniprot'
                accession = uniprotMatch.group(2)
        f.close()
    return (db, accession)


def fetchfromUniprot(uniprotId):

    #
    # fetch fasta file using Rest interface of uniprot
    #  http://www.uniprot.org/faq/28
    #

    uniprotURL = \
        'http://www.uniprot.org/uniprot/{0}.fasta'.format(uniprotId)
    fileName = uniprotId + '.fasta'
    try:
        req = urllib2.Request(uniprotURL)
        u = urllib2.urlopen(req)
        sequence = u.read()
        fileName = uniprotId + '.fasta'
        f = open(fileName, 'w')
        f.write(sequence)
        f.close()
        return fileName
    except urllib2.URLError:

        print '{0} is not valid uniprot id'.format(uniprotId)
    return None

    return fileName


