#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Hmmer3 based Protein Domain Searcher/Visualizer
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#

from __future__ import division
import urllib
import urllib2
import os
import xml.etree.ElementTree as ET
import argparse
import re
import sys
import subprocess
from hmmer import HmmerScanRunner


def extractRefSeq(data):
    gi = re.compile('gi\|(\S+)\|ref\|(\S+)\|')
    strip = re.compile('(\S+).\S+')
    refseq = ''
    mat = gi.match(data)
    if mat:
        mat2 = strip.match(mat.group(2))
        if mat2:
            refseq = mat2.group(1)
        else:

                 # gid = mat.group(1)

            refseq = mat.group(2)
    else:
        refseq = data
    return refseq


def extractMultiFasta(multifasta, searchId):

    capture = False
    fileNames = {}
    buffer = []
    f = open(multifasta)
    no = 1
    header = re.compile('^>(.*)')

    # gid = ''

    refseq = ''
    while True:
        line = f.readline()
        if not line:
            break
        match = header.match(line)
        if match:
            if capture:
                print refseq + '.fasta'
                fw = open(refseq + '.fasta', 'w')
                fw.writelines(buffer)
                fw.close
                fileNames[id] = refseq + '.fasta'
                capture = False
                refseq = ''
                buffer = []

            id = match.group(1).strip()
            if id in searchId:
                capture = True
                refseq = extractRefSeq(id)
                if id == refseq:
                    refseq = str(no)
                    no += 1
        if capture:
            buffer.append(line)

    if capture:
        print refseq + '.fasta'
        fw = open(refseq + '.fasta', 'w')
        fw.writelines(buffer)
        fw.close
        capture = False
        refseq = ''
        buffer = []
        fileNames[id] = refseq + '.fasta'

    return fileNames


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


def readHmm(filename):
    hmm = ''
    if os.path.exists(filename):
        f = open(filename)
        while True:
            line = f.readline()
            if not line:
                break
            hmm = hmm + line
        f.close()
        return hmm
    else:
        return None


def fetchFasta(proteinId, db):

    #
    # Fetch list (proteinId) of FASTA from the selected db.
    #

    fileNames = {}

    # if db can be downloadable in uniprot, download it.

    current = 1
    no = len(proteinId)
    if db in [
        'uniprotkb',
        'swissprot',
        'uniprotrefprot',
        'rp15',
        'env_nr',
        'pdb',
        'rp35',
        'rp55',
        'rp75',
        'nr',
        'env_nr',
        'pdb',
        ]:
        for singleId in proteinId:
            if db in [
                'uniprotkb',
                'swissprot',
                'uniprotrefprot',
                'rp15',
                'env_nr',
                'pdb',
                'rp35',
                'rp55',
                'rp75',
                ]:
                dbformat = 'http://www.uniprot.org/uniprot/{0}.fasta'
            if db in ['nr', 'env_nr']:
                dbformat = \
                    'http://www.ncbi.nlm.nih.gov/protein/{0}?report=fasta&log$=seqview&format=text'

            fileName = fetchfromDB(singleId, dbformat)
            print '{0}/{1} downloading {2}'.format(current, no,
                    fileName)
            if fileName:
                if not singleId in fileNames:
                    fileNames[singleId] = fileName

            if db in ['pdb']:
                fileName = fetchfromPDB(singleId)
                if not fileName:
                    if not proteinId in fileNames:
                        fileNames[singleId] = fileName
            current += 1

    return fileNames


def fetchfromPDB(pdbId):

    #
    # Not Yet Implemented
    #

    return None


def fetchfromDB(dbId, dbformat):

    #
    # fetch fasta file using Rest interface of uniprot
    #  http://www.uniprot.org/faq/28
    #

    uniprotURL = dbformat.format(dbId)
    fileName = dbId + '.fasta'
    if os.path.exists(fileName):
        print '{0} is exist. Skip downloading..'
        return fileName
    else:
        try:
            req = urllib2.Request(uniprotURL)
            u = urllib2.urlopen(req)
            sequence = u.read()
            fileName = dbId + '.fasta'
            f = open(fileName, 'w')
            f.write(sequence)
            f.close()
            return fileName
        except urllib2.URLError:

            print '{0} is not valid uniprot id'.format(dbId)
        return None

    return fileName


class SmartRedirectHandler(urllib2.HTTPRedirectHandler):

    def http_error_302(
        self,
        req,
        fp,
        code,
        msg,
        headers,
        ):
        return headers


class HmmFetch(object):

    #
    # Hmmscan wrapper class
    #

    def __init__(self, hmmdb, domainName):

        self.hmmdb = hmmdb
        self.domainName = domainName

    def run(self):
        hmmFileName = self.domainName + '.hmm'
        p = subprocess.Popen(['hmmfetch', '-o', hmmFileName,
                             self.hmmdb, self.domainName],
                             stdout=subprocess.PIPE)
        p_stdout = p.stdout.read()
        regex = re.compile('Retrieved HMM')
        match = regex.search(p_stdout)
        if match:
            print '{0} is generated'.format(hmmFileName)
            return hmmFileName
        else:
            print 'Problems. {0} is not generated'.format(hmmFileName)
            os.remove(hmmFileName)
            return None


class HmmerDomainHit(object):

    #
    # Class for the Hmmer domain hits
    #

    def __init__(self, **kwargs):
        self.name = kwargs.get('name')
        self.acc = kwargs.get('acc')
        self.desc = kwargs.get('desc')
        self.evalue = kwargs.get('evalue')
        self.ievalue = kwargs.get('ievalue')
        self.cevalue = kwargs.get('cevalue')
        self.start = int(kwargs.get('start'))
        self.end = int(kwargs.get('end'))
        self.bitscore = kwargs.get('bitscore')
        self.exclude = False
        self.color = None
        self.label = True
        self.border = True
        self.startshow = True
        self.endshow = True
        self.gradient = False


class HmmerSearchHit(object):

    def __init__(self, **kwargs):

        self.target = kwargs.get('target')
        self.query = kwargs.get('query')
        self.desc = kwargs.get('desc')
        self.acc = kwargs.get('acc')
        self.acc2 = kwargs.get('acc2')
        self.evalue = kwargs.get('evalue')
        self.bitscore = kwargs.get('bitscore')
        self.domainNo = kwargs.get('domainNo')
        self.species = kwargs.get('species')
        self.domains = []


class HmmerSearch(object):

    #
    # Hmmsearch wrapper class
    #

    def __init__(self, **kwargs):

        self.file = kwargs.get('file')
        self.db = kwargs.get('db')
        self.cutoff = kwargs.get('evalue')
        self.localHmmDB = kwargs.get('localHmmDB')
        self.threshold = kwargs.get('threshold')
        self.species = kwargs.get('species')
        self.hits = []

    def runLocal(self):

        #
        # Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store
        #

        if os.path.exists(self.db):
            proteinOutputFile = self.file + '.tb'
            domainOutputFile = self.file + '.dtb'
            hmmeroutFile = self.file + '.hmmer'
            print 'Processing {0}'.format(self.file)

            #
            # hmmscan --domtbout outputFile --cut_ga, DBfile, queryFile
            #

            if self.threshold == 'cut_ga':
                p = subprocess.Popen([
                    'hmmsearch',
                    '--domtblout',
                    domainOutputFile,
                    '--tblout',
                    proteinOutputFile,
                    '--cut_ga',
                    self.file,
                    self.db,
                    ], stdout=subprocess.PIPE)
            else:
                p = subprocess.Popen([
                    'hmmsearch',
                    '--domtblout',
                    domainOutputFile,
                    '--tblout',
                    proteinOutputFile,
                    self.file,
                    self.db,
                    ], stdout=subprocess.PIPE)

            p_stdout = p.stdout.read()
            fw = open(hmmeroutFile, 'w')
            fw.write(p_stdout)
            fw.close()

            f = open(proteinOutputFile, 'r')
            data = f.readlines()
            f.close()

            #
            # Parse table formatted domain datafile generated by '--domtblout'.
            #

            for line in data:
                if line[0:8] == '# target':
                    desclocation = line.find('description')
                if line[0:1] != '#':
                    splited = line.split()
                    if float(splited[4]) < self.cutoff:
                        target = splited[0]
                        query = splited[2]
                        acc = splited[3]
                        evalue = splited[4]
                        bitscore = splited[5]
                        domainNo = splited[17]
                        desc = line[desclocation:]
                        hit = HmmerSearchHit(
                            target=target,
                            query=query,
                            desc=desc,
                            acc=acc,
                            bitscore=bitscore,
                            evalue=evalue,
                            domainNo=domainNo,
                            )
                        self.hits.append(hit)
            return True
        else:
            print '{0} file is not exist'.format(self.db)
            sys.exit()
            return False

    def runRemote(self):

        #
        # Using Hmmsearch in Hmmer3 web service, find locations of domains in the Fasta sequence and store into class
        #

        hmm = readHmm(self.file)
        if not hmm:
            print 'Fatal error:'
            print '{0} does not exist.'.format(self.file)
            sys.exit()

        # print hmm

        opener = urllib2.build_opener(SmartRedirectHandler())
        urllib2.install_opener(opener)
        print 'search {0} in remote database'.format(self.file)
        if not self.db in [
            'uniprotkb',
            'swissprot',
            'nr',
            'uniprotrefprot',
            'rp15',
            'env_nr',
            'pdb',
            'rp35',
            'rp55',
            'rp75',
            ]:
            print 'invalid db. It should be uniprotkb,swissprot,nr, uniprotrefprot, rp15, env_nr, pdb, rp35, rp55 or rp75'
            sys.exit()

        parameters = {'seqdb': self.db, 'threshold': self.threshold,
                      'seq': hmm}
        enc_params = urllib.urlencode(parameters)

        # post the seqrch request to the server

        request = \
            urllib2.Request('http://hmmer.janelia.org/search/hmmsearch'
                            , enc_params)

        # get the url where the results can be fetched from

        results_url = urllib2.urlopen(request).getheader('location')

        # modify the range, format and presence of alignments in your results here

        res_params = {'output': 'xml'}

        # add the parameters to your request for the results

        enc_res_params = urllib.urlencode(res_params)
        modified_res_url = results_url + '?' + enc_res_params

        # send a GET request to the server

        results_request = urllib2.Request(modified_res_url)
        data = urllib2.urlopen(results_request)

        # print out the results

        result = data.read()
        if result:

            #
            # Parse using ElementTree Modules (http://docs.python.org/2/library/xml.etree.elementtree.html)
            #

            f = open(self.file + '.xml', 'w')
            f.write(result)
            f.close()
            root = ET.fromstring(result)
            for child in root.iter('hits'):

                name = child.get('name')
                desc = child.get('desc')
                acc = child.get('acc')
                acc2 = child.get('acc2')
                evalue = child.get('evalue')
                domainNo = child.get('nreported')
                bitscore = child.get('score')
                species = child.get('species')
                query = ''
                if float(evalue) < float(self.cutoff):
                    if self.species:
                        speciesFilter = self.species.split('+')
                        if species in speciesFilter:
                            hit = HmmerSearchHit(
                                target=name,
                                query=query,
                                desc=desc,
                                acc=acc,
                                acc2=acc2,
                                bitscore=bitscore,
                                evalue=evalue,
                                domainNo=domainNo,
                                species=species,
                                )
                            self.hits.append(hit)
                    else:
                        hit = HmmerSearchHit(
                            target=name,
                            query=query,
                            desc=desc,
                            acc=acc,
                            acc2=acc2,
                            bitscore=bitscore,
                            evalue=evalue,
                            domainNo=domainNo,
                            species=species,
                            )
                        self.hits.append(hit)

            return True
        else:
            print 'Failed to retrieve results'
            return False


def generateInputFile(
    proteinNames,
    fastaFileNames,
    inputFileName,
    argument,
    ):

    f = open(inputFileName, 'w')
    for (id, name) in proteinNames.items():
        line = 'FILE\t{0}\t{1}\n'.format(name, fastaFileNames[id])
        f.write(line)

    if argument.localscan:
        f.write('SETUP\tLOCAL\t{0}\n'.format(argument.hmmdb))
    else:
        f.write('SETUP\tREMOTE\t{0}\n'.format('pfam'))

    if argument.domains:
        fileNameBase = ''.join(argument.domains)
        f.write('SETUP\tOUTPUTHTML\t{0}.html\n'.format(fileNameBase))
        f.write('SETUP\tOUTPUTSVG\t{0}.svg\n'.format(fileNameBase))

    f.close()

    return inputFileName


def main(argument):

    hmmFiles = []
    if len(argument.domains) > 0:
        domainNames = argument.domains
        hmmdb = argument.hmmdb
        for domainName in domainNames:
            singleHmm = HmmFetch(hmmdb, domainName)
            hmmFile = singleHmm.run()
            if hmmFile:
                hmmFiles.append(hmmFile)
    if len(argument.hmmFiles) > 0:
        for hmmFile in argument.hmmFiles:
            if os.path.exists(hmmFile):
                hmmFiles.append(hmmFile)
    if argument.threshold:
        threshold = 'cut_ga'
    else:
        threshold = 'No'

    ids = []
    proteinNames = {}
    if len(hmmFiles) > 0:
        for hmmFile in hmmFiles:
            hmmSearch = HmmerSearch(file=hmmFile,
                                    db=argument.proteindb,
                                    evalue=argument.evalue,
                                    threshold=threshold,
                                    species=argument.species)
            if argument.local:
                hmmSearch.runLocal()
                for hit in hmmSearch.hits:
                    id = hit.target + ' ' + hit.desc.strip()
                    ids.append(id)
                    if not id in proteinNames:
                        proteinNames[id] = hit.desc.strip()
            else:
                hmmSearch.runRemote()
                for hit in hmmSearch.hits:
                    if hit.acc2:
                        id = hit.acc2
                    else:
                        id = hit.acc
                    ids.append(id)
                    if not id in proteinNames:
                        proteinNames[id] = hit.desc.strip() + '[' \
                            + hit.species + ']'

        # # Bug fix needed,

        print '# of hits :{0}'.format(len(hmmSearch.hits))
        inputFileName = ''.join(argument.domains) + '.INP'
        if argument.local:
            fastaFileNames = extractMultiFasta(argument.proteindb, ids)
            inputFileName = generateInputFile(proteinNames,
                    fastaFileNames, inputFileName, argument)
        else:
            fastaFileNames = fetchFasta(ids, argument.proteindb)
            inputFileName = generateInputFile(proteinNames,
                    fastaFileNames, inputFileName, argument)
        hmmerscan = HmmerScanRunner(inputFileName=inputFileName)
        hmmerscan.run()
    else:
        print 'Nothing to do!'
        sys.exit()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',
        '--domains',
        nargs='+',
        dest='domains',
        default=[],
        help='Domains to search',
        )
    parser.add_argument('-hd', '--hmm', dest='hmmFiles', default=[],
                        help='HMM query File')
    parser.add_argument('-p', '--proteindb', dest='proteindb',
                        default='swissprot', help='protein db')
    parser.add_argument('-b', '--hmmdb', dest='hmmdb',
                        default='Pfam-A.hmm', help='protein db')
    parser.add_argument('-e', '--evalue_cutoff', dest='evalue',
                        default=1e-5, help='E-value cutoff')
    parser.add_argument(
        '-l',
        '--local',
        action='store_true',
        dest='local',
        default=False,
        help='run local Hmmer',
        )
    parser.add_argument(
        '-ls',
        '--local_scan',
        action='store_true',
        dest='localscan',
        default=False,
        help='run local Hmmerscan',
        )
    parser.add_argument(
        '-t',
        '--no_threshold',
        dest='threshold',
        action='store_false',
        default=True,
        help='Turn of Pfam gathering threshold. Enable to look up more weak(unreliable) domains'
            ,
        )
    parser.add_argument('-s', '--species', dest='species',
                        help='Confine search species. ex: Arabidopsis thaliana+Homo sapiens'
                        )
    results = parser.parse_args()

    main(results)

